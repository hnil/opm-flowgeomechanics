#!/usr/bin/env python3
"""Run or inspect the small geomechanics physics baselines (standard library only)."""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
import xml.etree.ElementTree as ET

from frac_area import area_of

REPO = Path(__file__).resolve().parents[1]
PROFILES = REPO / "tests/physics_baselines/cases.json"
REPOSITORIES = ("opm-common", "opm-grid", "opm-simulators", "opm-upscaling",
                "opm-flowgeomechanics")
WARNING_PATTERNS = {
    "missing_connection": r"Warning: Well connection with global index .* not found",
    "mechanical_bc": r"(?:Mechanical BC configuration:|Error in boundary condition specification)",
    "unsupported_keywords": r"Warning: Unsupported keywords or keyword items:",
    "refined_summary": r"Warning: Summary block vector .* refined away",
    "lgr_well_averages": r"Warning: Well .* completed in an LGR; its block average pressures",
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


HARNESS_SHA256 = sha256(Path(__file__))
GEOMETRY_READER_SHA256 = sha256(REPO / "scripts/frac_area.py")


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def capture(command, cwd=None):
    result = subprocess.run(command, cwd=cwd, capture_output=True, text=True,
                            timeout=60, check=True)
    return result.stdout


def read_cache(build):
    cache = {}
    for line in (build / "CMakeCache.txt").read_text().splitlines():
        match = re.match(r"([^/#][^:]*):[^=]+=(.*)", line)
        if match:
            cache[match[1]] = match[2]
    return cache


def provenance(build, profile_path, inputs, command):
    cache = read_cache(build)
    revisions = {}
    for name in REPOSITORIES:
        repo = REPO.parent / name
        diff = capture(["git", "diff", "HEAD", "--binary"], cwd=repo)
        untracked = capture(["git", "ls-files", "--others", "--exclude-standard",
                             "--", "opm", "cmake", "scripts", "tests", "docs", "data",
                             "CMakeLists.txt"], cwd=repo).splitlines()
        revisions[name] = {
            "commit": capture(["git", "rev-parse", "HEAD"], cwd=repo).strip(),
            "branch": capture(["git", "branch", "--show-current"], cwd=repo).strip(),
            "tracked_diff_sha256": hashlib.sha256(diff.encode()).hexdigest(),
            "untracked_sources": {p: sha256(repo / p) for p in untracked},
        }
    keys = {k: v for k, v in cache.items() if (
        k in {"CMAKE_BUILD_TYPE", "CMAKE_CXX_COMPILER", "CMAKE_CXX_FLAGS",
              "CMAKE_CXX_FLAGS_RELEASE", "CMAKE_PREFIX_PATH", "MPIEXEC_EXECUTABLE",
              "MPIEXEC_NUMPROC_FLAG", "MPIEXEC_PREFLAGS", "MPIEXEC_POSTFLAGS",
              "HYPRE_LIBRARY", "HYPRE_INCLUDE_DIRS"}
        or k.startswith("USE_") or k.endswith("_DIR") and
        k.startswith(("dune-", "CGAL", "Boost", "fmt")))}
    versions = {}
    compiler = cache.get("CMAKE_CXX_COMPILER")
    for label, tool in (("compiler", compiler), ("cmake", "cmake"),
                        ("mpi_launcher", cache.get("MPIEXEC_EXECUTABLE"))):
        if tool and not tool.endswith("NOTFOUND"):
            try:
                versions[label] = capture([tool, "--version"]).strip()
            except (OSError, subprocess.SubprocessError) as error:
                versions[label] = {"unavailable": str(error)}
    for key, directory in keys.items():
        if key.startswith("dune-") and key.endswith("_DIR"):
            version_file = Path(directory) / (key[:-4] + "-config-version.cmake")
            if version_file.is_file():
                match = re.search(r'set\(PACKAGE_VERSION "([^"]+)"', version_file.read_text())
                if match:
                    versions[key[:-4]] = match[1]
    for label, path, pattern in (
        ("CGAL", Path(cache.get("CGAL_DIR", "")) / "CGALConfigVersion.cmake",
         r'set\(CGAL_VERSION_PUBLIC_RELEASE_VERSION "([^"]+)"'),
        ("fmt", Path(cache.get("fmt_DIR", "")) / "fmt-config-version.cmake",
         r'set\(PACKAGE_VERSION "([^"]+)"'),
        ("HYPRE", Path(cache.get("HYPRE_INCLUDE_DIRS", "")) / "HYPRE_config.h",
         r'#define HYPRE_RELEASE_VERSION "([^"]+)"'),
    ):
        if path.is_file():
            match = re.search(pattern, path.read_text())
            if match:
                versions[label] = match[1]
    binaries = [build / "opm-flowgeomechanics/bin/flow_energy_geomech",
                build / "opm-common/bin/summary"]
    return {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "provenance_kind": "fresh_run",
        "repositories": revisions,
        "inputs": {str(p): sha256(p) for p in inputs},
        "profiles_sha256": sha256(profile_path),
        "harness_sha256": HARNESS_SHA256,
        "geometry_reader_sha256": GEOMETRY_READER_SHA256,
        "binaries": {str(p): sha256(p) for p in binaries},
        "cmake_cache": keys,
        "cmake_cache_sha256": sha256(build / "CMakeCache.txt"),
        "tool_versions": versions,
        "platform": {"system": platform.platform(), "machine": platform.machine(),
                     "cpu_count": os.cpu_count(), "python": platform.python_version()},
        "environment": {key: os.environ[key] for key in
                        ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
                        if key in os.environ},
        "command": command,
        "cwd": str(REPO.parent),
    }


def parse_summary(output):
    rows = []
    for line in output.splitlines():
        fields = line.split()
        if len(fields) != 5:
            if fields and re.match(r"^[+-]?(?:\d|\.\d)", fields[0]):
                raise ValueError("Malformed numeric summary row")
            continue
        try:
            values = [float(value) for value in fields]
        except ValueError:
            continue
        if not all(math.isfinite(value) for value in values):
            raise ValueError("Non-finite summary value")
        day, bhp, rate, fracture_rate, total = values
        if rows and day <= rows[-1]["day"]:
            raise ValueError("Summary times are not strictly increasing")
        rows.append(dict(day=day, bhp_bar=bhp, injection_m3_per_day=rate,
                         fracture_injection_m3_per_day=fracture_rate,
                         total_injected_m3=total,
                         fracture_flow_fraction=fracture_rate / rate if rate > 0 else None))
    if not rows:
        raise ValueError("No numeric summary rows: required vectors may be missing")
    return rows


def geometry_history(run_dir):
    collections = sorted(run_dir.glob("*Fracture_on_*_nr.pvd"))
    if not collections:
        raise ValueError("No fracture PVD collection; enable-write-all-solutions is required")
    if len(collections) != 1:
        raise ValueError("These one-fracture baselines require exactly one fracture collection")
    rows = []
    for entry in ET.parse(collections[0]).getroot().iter("DataSet"):
        day = float(entry.attrib["timestep"]) / 86400.0
        path = (run_dir / entry.attrib["file"]).resolve()
        if not path.is_relative_to(run_dir.resolve()):
            raise ValueError("Fracture PVD references a file outside the run directory")
        piece = ET.parse(path).getroot().find(".//Piece")
        widths = piece.find(".//DataArray[@Name='FractureWidth']") if piece is not None else None
        if widths is None or not widths.text or not all(
                math.isfinite(float(w)) and float(w) >= 0 for w in widths.text.split()):
            raise ValueError(f"Missing, negative, or non-finite fracture width in {path.name}")
        area, opened, volume = area_of(path)
        if not all(math.isfinite(v) and v >= 0 for v in (day, area, opened, volume)):
            raise ValueError("Invalid fracture geometry metric")
        if rows and day < rows[-1]["day"]:
            raise ValueError("Fracture collection times decrease")
        rows.append(dict(day=day, fracture_area_m2=area, fracture_open_area_m2=opened,
                         fracture_volume_m3=volume, file=path.name))
    if not rows:
        raise ValueError("Empty fracture collection")
    return rows


def log_metrics(log):
    clean = re.sub(r"\x1b\[[0-9;]*m", "", log)
    warnings = {name: len(re.findall(pattern, clean, re.I))
                for name, pattern in WARNING_PATTERNS.items()}
    other = [line for line in clean.splitlines() if re.search(r"\bwarning\s*:", line, re.I)
             and not any(re.search(pattern, line, re.I) for pattern in WARNING_PATTERNS.values())]
    warnings["other"] = len(other)
    totals = {}
    for name in ("nonlinear_iterations", "linear_solves", "linear_iterations",
                 "linear_solve_failures", "ladder_rescues", "closed_cell_toggles"):
        found = re.findall(r"\b" + name + r"=\S+ \(total (\d+)\)", clean)
        if found:
            totals[name] = max(map(int, found))
    return {
        "end_of_simulation": "End of simulation" in clean,
        "nonconvergence_messages": len(re.findall(
            r"did not converge|failed to converge|converged=false", clean, re.I)),
        "warning_occurrences": warnings,
        "other_warning_samples": other[:10],
        "reported_max_cumulative_solver_counts": totals,
        "solver_count_scope": "maximum emitted cumulative count, not a sum across MPI ranks",
    }


def assess(summary_rows, geometry_rows, log, profile, strict_warnings=False):
    final = summary_rows[-1]
    geometry = geometry_rows[-1]
    metrics = {**final, **{k: v for k, v in geometry.items() if k not in ("day", "file")},
               "geometry_day": geometry["day"],
               "geometry_lag_days": final["day"] - geometry["day"],
               "area_growth_ratio": geometry["fracture_area_m2"] /
               max(geometry_rows[0]["fracture_area_m2"], 1e-30)}
    checks = {}
    checks["final_day"] = abs(final["day"] - profile["expected_final_day"]) <= 1e-6
    checks["geometry_timestamp"] = -1e-6 <= metrics["geometry_lag_days"] <= profile["max_geometry_lag_days"]
    checks["fracture_growth"] = metrics["area_growth_ratio"] >= profile["minimum_area_growth_ratio"]
    for name, (lower, upper) in profile["bands"].items():
        value = metrics.get(name)
        checks[name] = value is not None and lower <= value <= upper
    diagnostics = log_metrics(log)
    checks["end_of_simulation"] = diagnostics["end_of_simulation"]
    checks["no_reported_nonconvergence"] = diagnostics["nonconvergence_messages"] == 0
    totals = diagnostics["reported_max_cumulative_solver_counts"]
    checks["solver_diagnostics_present"] = bool(totals)
    checks["no_linear_solve_failures"] = totals.get("linear_solve_failures", -1) == 0
    checks["no_solver_rescues"] = totals.get("ladder_rescues", -1) == 0
    if strict_warnings:
        checks["no_warnings"] = not any(diagnostics["warning_occurrences"].values())
    warnings_present = any(diagnostics["warning_occurrences"].values())
    unexpected = [name for name, count in diagnostics["warning_occurrences"].items()
                  if count and name not in profile["known_warnings"]]
    checks["no_unexpected_warning_types"] = not unexpected
    status = "fail" if not all(checks.values()) else "pass_with_warnings" if warnings_present else "pass"
    return dict(schema_version=1, status=status, metrics=metrics, checks=checks,
                diagnostics=diagnostics, known_warning_types=profile["known_warnings"],
                unexpected_warning_types=unexpected,
                acceptance_bands=profile["bands"])


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def analyze(run_dir, summary_tool, log_path, profile, strict_warnings=False, report_dir=None,
            profile_path=PROFILES):
    report_dir = report_dir or run_dir
    stem = Path(profile["deck"]).stem
    spec = run_dir / (stem + ".SMSPEC")
    well = profile["well"]
    output = capture([str(summary_tool), str(spec), "TIME", f"WBHP:{well}",
                      f"WWIR:{well}", f"WWIRFRAC:{well}", f"WWIT:{well}"])
    summary_rows = parse_summary(output)
    geometry_rows = geometry_history(run_dir)
    report = assess(summary_rows, geometry_rows, log_path.read_text(errors="replace"),
                    profile, strict_warnings)
    files = [spec, spec.with_suffix(".UNSMRY"), log_path]
    files += list(run_dir.glob("*Fracture_on_*_nr.pvd"))
    files += [run_dir / row["file"] for row in geometry_rows]
    report["analyzed_inputs"] = {str(p.resolve()): sha256(p) for p in files}
    report["profiles_sha256"] = sha256(profile_path)
    report["harness_sha256"] = HARNESS_SHA256
    report["geometry_reader_sha256"] = GEOMETRY_READER_SHA256
    write_csv(report_dir / "physics-summary.csv", summary_rows)
    write_csv(report_dir / "physics-geometry.csv", geometry_rows)
    write_json(report_dir / "physics-report.json", report)
    return report


def run_case(name, profile, args, profile_path):
    run_dir = (args.output_root / name).resolve()
    run_dir.mkdir(parents=True, exist_ok=False)
    build = args.build_dir.resolve()
    cache = read_cache(build)
    binary = build / "opm-flowgeomechanics/bin/flow_energy_geomech"
    summary = build / "opm-common/bin/summary"
    inputs = [REPO / profile[key] for key in ("deck", "fracture_parameters", "mechanics_solver")]
    command = []
    if cache.get("USE_MPI", "OFF").upper() in ("ON", "TRUE", "1"):
        launcher = args.mpi_launcher or cache.get("MPIEXEC_EXECUTABLE")
        if not launcher or launcher.endswith("NOTFOUND"):
            raise ValueError("MPI build needs --mpi-launcher or MPIEXEC_EXECUTABLE")
        command = [launcher, cache.get("MPIEXEC_NUMPROC_FLAG", "-n"), str(profile["mpi_ranks"])]
        command += [v for v in cache.get("MPIEXEC_PREFLAGS", "").split(";") if v]
        command += args.launcher_arg
    elif profile["mpi_ranks"] != 1:
        raise ValueError("The WELLREF baseline requires an MPI-enabled build")
    command += [str(binary)]
    if command[0] != str(binary):
        command += [v for v in cache.get("MPIEXEC_POSTFLAGS", "").split(";") if v]
    command += [str(inputs[0]), f"--fracture-param-file={inputs[1]}",
                f"--linear-solver-mech={inputs[2]}", "--enable-write-all-solutions=true",
                "--threads-per-process=1", f"--output-dir={run_dir}"]
    if name == "wellref":
        command += ["--edge-conformal=true", "--parsing-strictness=low"]
    manifest = provenance(build, profile_path, inputs, command)
    manifest["case"] = name
    manifest["mpi_ranks"] = profile["mpi_ranks"]
    manifest["threads_per_rank"] = 1
    manifest["acceptance_profile"] = profile
    write_json(run_dir / "physics-manifest.json", manifest)
    print(f"Running {name}: {profile['mpi_ranks']} rank(s), outputs in {run_dir}", flush=True)
    start = time.monotonic()
    with (run_dir / "simulator.log").open("w") as log:
        process = subprocess.Popen(command, cwd=REPO.parent, stdout=log,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = process.wait(timeout=args.timeout)
        except (subprocess.TimeoutExpired, KeyboardInterrupt) as error:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            manifest.update(timed_out=isinstance(error, subprocess.TimeoutExpired),
                            interrupted=isinstance(error, KeyboardInterrupt),
                            returncode=process.returncode,
                            elapsed_wall_seconds=time.monotonic() - start)
            write_json(run_dir / "physics-manifest.json", manifest)
            raise RuntimeError(f"{name} interrupted or exceeded {args.timeout} seconds")
    manifest.update(returncode=code, elapsed_wall_seconds=time.monotonic() - start,
                    log_sha256=sha256(run_dir / "simulator.log"))
    write_json(run_dir / "physics-manifest.json", manifest)
    if code:
        raise RuntimeError(f"{name} simulator exited with {code}; see {run_dir / 'simulator.log'}")
    return analyze(run_dir, summary, run_dir / "simulator.log", profile, args.strict_warnings,
                   profile_path=profile_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles", type=Path, default=PROFILES)
    sub = parser.add_subparsers(dest="mode", required=True)
    run = sub.add_parser("run", help="Run pinned SIMPLE/WELLREF profiles into new directories")
    run.add_argument("--build-dir", type=Path, required=True)
    run.add_argument("--output-root", type=Path, required=True)
    run.add_argument("--cases", nargs="+", choices=("simple", "wellref"), default=["simple", "wellref"])
    run.add_argument("--mpi-launcher")
    run.add_argument("--launcher-arg", action="append", default=[])
    run.add_argument("--timeout", type=float, default=900)
    inspect = sub.add_parser("analyze", help="Inspect existing outputs without claiming their run provenance")
    inspect.add_argument("--case", choices=("simple", "wellref"), required=True)
    inspect.add_argument("--run-dir", type=Path, required=True)
    inspect.add_argument("--summary-tool", type=Path, required=True)
    inspect.add_argument("--log", type=Path, required=True)
    inspect.add_argument("--report-dir", type=Path, required=True,
                         help="New directory for analysis reports; historical outputs are unchanged")
    for command in (run, inspect):
        command.add_argument("--strict-warnings", action="store_true")
    args = parser.parse_args()
    profiles = json.loads(args.profiles.read_text())
    if profiles.get("schema_version") != 1:
        parser.error("Unsupported profile schema version")
    if args.mode == "run":
        if not math.isfinite(args.timeout) or args.timeout <= 0:
            parser.error("--timeout must be finite and positive")
        if len(set(args.cases)) != len(args.cases):
            parser.error("--cases must not contain duplicates")
        try:
            args.output_root.mkdir(parents=True, exist_ok=False)
        except FileExistsError:
            parser.error("--output-root must be a new directory; existing evidence is never overwritten")
    else:
        try:
            args.report_dir.mkdir(parents=True, exist_ok=False)
        except FileExistsError:
            parser.error("--report-dir must be a new directory")
    results = {}
    names = args.cases if args.mode == "run" else [args.case]
    for name in names:
        profile = {**profiles["common"], **profiles["cases"][name]}
        try:
            if args.mode == "run":
                report = run_case(name, profile, args, args.profiles)
            else:
                report = analyze(args.run_dir.resolve(), args.summary_tool.resolve(),
                                 args.log.resolve(), profile, args.strict_warnings,
                                 args.report_dir.resolve(), args.profiles)
                report["provenance_kind"] = "analysis_only_original_run_provenance_not_verified"
                write_json(args.report_dir / "physics-report.json", report)
            results[name] = report
            print(f"{name}: {report['status']}; failed checks: "
                  f"{[key for key, passed in report['checks'].items() if not passed]}", flush=True)
        except (OSError, ValueError, RuntimeError, subprocess.SubprocessError, ET.ParseError) as error:
            results[name] = {"status": "fail", "error": str(error)}
            print(f"{name}: fail: {error}", file=sys.stderr, flush=True)
    if args.mode == "run":
        args.output_root.mkdir(parents=True, exist_ok=True)
        write_json(args.output_root / "physics-suite.json", results)
    return int(any(result["status"] == "fail" for result in results.values()))


if __name__ == "__main__":
    sys.exit(main())
