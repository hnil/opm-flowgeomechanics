"""Fast harness checks; these do not require the simulator or MPI."""

import copy
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location("physics_baselines", SCRIPTS / "physics_baselines.py")
harness = importlib.util.module_from_spec(spec)
spec.loader.exec_module(harness)


class BaselineChecks(unittest.TestCase):
    def setUp(self):
        profiles = json.loads(harness.PROFILES.read_text())
        self.profile = {**profiles["common"], **profiles["cases"]["simple"]}
        self.summary = harness.parse_summary(
            "TIME WBHP WWIR WWIRFRAC WWIT\n"
            "0 0 0 0 0\n106 220.98 20000 19276 2120000\n")
        self.geometry = [
            dict(day=0, fracture_area_m2=76, fracture_open_area_m2=76,
                 fracture_volume_m3=0.01, file="seed.vtu"),
            dict(day=101, fracture_area_m2=11695, fracture_open_area_m2=11473,
                 fracture_volume_m3=293.23, file="grown.vtu"),
        ]
        self.log = ("End of simulation\nFracture solve stats: nonlinear_iterations=4 (total 8), "
                    "linear_solve_failures=0 (total 0), ladder_rescues=0 (total 0), converged=true\n")

    def assess(self, summary=None, geometry=None, log=None, strict=False):
        return harness.assess(summary or self.summary, geometry or self.geometry,
                              self.log if log is None else log, self.profile, strict)

    def test_completed_fracturing_case_passes_with_distinct_timestamps(self):
        result = self.assess()
        self.assertEqual(result["status"], "pass")
        self.assertEqual(result["metrics"]["day"], 106)
        self.assertEqual(result["metrics"]["geometry_day"], 101)
        self.assertIsNone(self.summary[0]["fracture_flow_fraction"])

    def test_finished_case_without_fracturing_fails(self):
        geometry = copy.deepcopy(self.geometry)
        geometry[-1].update(fracture_area_m2=76, fracture_open_area_m2=76,
                            fracture_volume_m3=0.01)
        result = self.assess(geometry=geometry)
        self.assertFalse(result["checks"]["fracture_growth"])
        self.assertEqual(result["status"], "fail")

    def test_low_fracture_flow_fails_even_with_large_geometry(self):
        summary = copy.deepcopy(self.summary)
        summary[-1]["fracture_flow_fraction"] = 0.1
        self.assertFalse(self.assess(summary=summary)["checks"]["fracture_flow_fraction"])

    def test_truncated_run_or_stale_geometry_fails(self):
        summary = copy.deepcopy(self.summary)
        summary[-1]["day"] = 101
        self.assertFalse(self.assess(summary=summary)["checks"]["final_day"])
        geometry = copy.deepcopy(self.geometry)
        geometry[-1]["day"] = 96
        self.assertFalse(self.assess(geometry=geometry)["checks"]["geometry_timestamp"])

    def test_solver_failure_is_not_hidden_by_completion(self):
        result = self.assess(log=self.log + "linear_solve_failures=1 (total 1)\n")
        self.assertFalse(result["checks"]["no_linear_solve_failures"])
        self.assertFalse(self.assess(log=self.log + "did not converge\n")
                         ["checks"]["no_reported_nonconvergence"])
        self.assertFalse(self.assess(log="End of simulation\n")
                         ["checks"]["solver_diagnostics_present"])

    def test_known_warnings_are_visible_and_strict_mode_rejects_them(self):
        log = self.log + "Warning: Unsupported keywords or keyword items:\n"
        self.assertEqual(self.assess(log=log)["status"], "pass_with_warnings")
        self.assertEqual(self.assess(log=log, strict=True)["status"], "fail")

    def test_missing_lgr_completion_is_no_longer_allowlisted(self):
        profiles = json.loads(harness.PROFILES.read_text())
        for case in profiles["cases"].values():
            self.profile["known_warnings"] = case["known_warnings"]
            log = self.log + ("Warning: Well connection with global index 6562 "
                              "in grid 1 for well B-3H not found in schedule step 0\n")
            result = self.assess(log=log)
            self.assertEqual(result["status"], "fail")
            self.assertEqual(result["unexpected_warning_types"], ["missing_connection"])

    def test_mechanical_constraint_warnings_are_no_longer_allowlisted(self):
        profiles = json.loads(harness.PROFILES.read_text())
        for case in profiles["cases"].values():
            self.profile["known_warnings"] = case["known_warnings"]
            result = self.assess(log=self.log +
                "Mechanical BC configuration: no nodes are fixed in Z-direction\n")
            self.assertEqual(result["status"], "fail")
            self.assertEqual(result["unexpected_warning_types"], ["mechanical_bc"])

    def test_new_warning_category_fails(self):
        result = self.assess(log=self.log + "Warning: new physical inconsistency\n")
        self.assertEqual(result["status"], "fail")
        self.assertEqual(result["unexpected_warning_types"], ["other"])

    def test_summary_rejects_nonfinite_missing_and_nonmonotone_rows(self):
        for data in ("No vectors found", "106 nan 20000 19276 2120000",
                     "106 220 20000 19276",
                     "106 220 20000 19276 2120000\n101 220 20000 19276 2020000"):
            with self.subTest(data=data), self.assertRaises(ValueError):
                harness.parse_summary(data)

    def test_pvd_times_and_known_triangle_area_volume(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            # Analytic triangle: area=1/2, width=2, volume=1.
            (root / "triangle.vtu").write_text('''<VTKFile><UnstructuredGrid><Piece>
              <Points><DataArray Name="Coordinates">0 0 0 1 0 0 0 1 0</DataArray></Points>
              <Cells><DataArray Name="connectivity">0 1 2</DataArray>
              <DataArray Name="offsets">3</DataArray></Cells>
              <CellData><DataArray Name="FractureWidth">2</DataArray></CellData>
              </Piece></UnstructuredGrid></VTKFile>''')
            pvd = root / "CASEFracture_on_W_nr.pvd"
            pvd.write_text('<VTKFile><Collection><DataSet timestep="172800" '
                           'file="triangle.vtu"/></Collection></VTKFile>')
            row = harness.geometry_history(root)[0]
            self.assertEqual(row["day"], 2)
            self.assertEqual(row["fracture_area_m2"], 0.5)
            self.assertEqual(row["fracture_volume_m3"], 1)
            pvd.write_text('<VTKFile><Collection><DataSet timestep="0" '
                           'file="missing.vtu"/></Collection></VTKFile>')
            with self.assertRaises(FileNotFoundError):
                harness.geometry_history(root)


if __name__ == "__main__":
    unittest.main()
