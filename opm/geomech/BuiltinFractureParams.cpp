#include <config.h>

#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/linalg/PropertyTree.hpp>

#include <stdexcept>
#include <string>

namespace {
const char* standardFractureParamJson()
{
    return R"opmjson(
{
    "hasfractures": "true",
    "verbosity": "1",
    "add_perfs_to_schedule": "true",
    "solver": {
        "method": "PostSolve",
        "implicit_flow": "false",
        "max_mech_it": "2"
    },
    "fractureparam": {
        "verbosity": "1",
        "method": {
            "iterate": "true",
            "max_it": "3",
            "tolerance": "0.0001"
        },
        "reduce_boundary": "false",
        "vem_stability_choice": "3",
        "mech_diagonal_scaling": "false",
        "patch_recovery_stress": "false",
        "smooth_force": "false",
        "smooth_force_length": "0",
        "addconnections": "true",
        "include_fracture_contributions": "false",
        "config": {
            "verbosity": "1",
            "type": "well_seed",
            "initial_fracture_width": "0.0001",
            "min_width": "0",
            "trires": "5",
            "gravity_off": "false",
            "scale_filtrate": "true"
        },
        "solver": {
            "verbosity": "2",
            "method": "if_propagate_trimesh",
            "target_cellcount": "300",
            "cellcount_threshold": "300",
            "numcell_threshold": "100000",
            "max_num_coarsening": "0",
            "max_iter_on_same_level": "100000000",
            "efac": "0.5",
            "rfac": "0.10000000000000001",
            "max_expand_iter": "20",
            "max_iter": "100",
            "tolerance": "9.9999999999999995e-07",
            "damping": "1",
            "min_width": "0.001",
            "max_width": "0.5",
            "max_change": "100000",
            "max_dwidth": "0.0050000000000000001",
            "max_dp": "1000000000",
            "remap_solution": "false",
            "linsolver": {
                "tol": "1e-10",
                "atol": "9.9999999999999995e-21",
                "max_iter": "1000",
                "verbosity": "0",
                "solver": "bicgstab",
                "preconditioner": {
                    "diag_mech": "false",
                    "diag_flow": "false",
                    "mech_first": "false",
                    "mech_press_coupling": "true",
                    "verbosity": "0",
                    "flow_solver": {
                        "solver": "umfpack",
                        "verbosity": "0",
                        "tol": "1",
                        "maxiter": "1",
                        "preconditioner": {
                            "type": "DILU",
                            "verbosity": "0"
                        }
                    }
                }
            },
            "area_change_fac": "3",
            "dt_limit": "0.10000000000000001",
            "damping_factor_perf": "2",
            "damping_factor_wi": "2",
            "failure_on_nonconvergence": "false",
            "force_limit": "0",
            "smooth_boundary": "false",
            "full_intersections": "false",
            "divide_wellidx": "false",
            "no_leakof_outercells": "false",
            "drop_fluid_mech_linearization": "true",
            "closed_cell_policy": "legacy",
            "coupling_legacy_first_damp": "true"
        },
        "reservoir": {
            "calculate_dist": "true",
            "dist": "1",
            "interpolate_stress": "false",
            "mobility": "0.0012999999999999999",
            "perm": "1e-13"
        },
        "control": {
            "type": "perf_pressure",
            "rate": "0.029000000000000001",
            "WI": "9.9999999999999995e-07"
        },
        "fractureWI": "9.9999999999999995e-07",
        "extended_fractures": "true",
        "KMax": "1000000",
        "write_pressure_system": "false",
        "write_fracture_system": "false",
        "pressuresolver": "umfpack",
        "fracturesolver": "notused"
    }
}
)opmjson";
}

const char* sequentialImplicitFractureParamJson()
{
    return R"opmjson(
{
    "__doc": "opm-flowgeomechanics fracture parameter file. method=SeqMechFrac (sequential-implicit). Grouped + annotated; see fracture_settings_REFERENCE.md",
    "hasfractures": "true",
    "add_perfs_to_schedule": "true",
    "verbosity": 0,
    "solver": {
        "__doc": "Top-level SeqMechFrac driver options (read from getGeoMechParam().solver). See fracture_settings_REFERENCE.md",
        "method": "SeqMechFrac",
        "implicit_flow": "true",
        "max_mech_it": "90",
        "max_frac_it": "90",
        "verbosity": 0
    },
    "fractureparam": {
        "__doc": "Fracture model parameters. Outer SeqMechFrac coupling knobs live here and under .solver",
        "verbosity": 0,
        "reduce_boundary": "false",
        "vem_stability_choice": "1",
        "vem_source": "true",
        "vem_stress": "true",
        "vem_force": "true",
        "stab_on_stress": "false",
        "smooth_force": "false",
        "addconnections": "true",
        "include_fracture_contributions": "false",
        "KMax": "1000000",
        "extended_fractures": "true",
        "fractureWI": "1e-6",
        "write_pressure_system": "false",
        "write_fracture_system": "false",
        "pressuresolver": "umfpack",
        "fracturesolver": "notused",
        "require_converged_fracture_for_wi_update": true,
        "topology_hysteresis": {
            "enable": false,
            "on_threshold": 1.0,
            "off_threshold": 0.3,
            "confirm_iterations": 1,
            "cooldown_iterations": 0,
            "emergency_threshold": 10.0
        },
        "method": {
            "iterate": "true",
            "max_it": "3",
            "tolerance": "0.0001"
        },
        "config": {
            "type": "well_seed",
            "initial_fracture_width": "0.0001",
            "min_width": "0",
            "trires": "10",
            "gravity_off": "false",
            "scale_filtrate": "true"
        },
        "solver": {
            "__doc": "Fracture inner-solver + outer-coupling options (read from fractureparam.solver). See fracture_settings_REFERENCE.md",
            "method": "if_propagate_trimesh",
            "__1_propagation_growth": "--- fracture grid growth / trimesh propagation (KEEP fidelity here) ---",
            "area_change_fac": 1.2,
            "dt_limit": 10,
            "efac": "0.5",
            "rfac": "0.1",
            "max_expand_iter": "20",
            "max_iter_on_same_level": "100000000",
            "target_cellcount": "400",
            "cellcount_threshold": "800",
            "numcell_threshold": "1500",
            "max_num_coarsening": "0",
            "smooth_boundary": "false",
            "remap_solution": "true",
            "full_intersections": "true",
            "divide_wellidx": "true",
            "no_leakof_outercells": "false",
            "__2_inner_nonlinear_solve": "--- inner width/pressure Newton (per fixed geometry) ---",
            "max_iter": "100",
            "tolerance": "1e-6",
            "damping": "1",
            "min_width": "0.001",
            "max_width": "0.5",
            "max_dwidth": "0.005",
            "max_dp": "1000000000",
            "max_change": "100000",
            "force_limit": "0.0",
            "failure_on_nonconvergence": "false",
            "__3_closed_cell_handling": "--- closed-cell (contact) logic; drives convergence AND growth ---",
            "closing_type": "org",
            "closed_cell_policy": "sticky",
            "close_force_tolerance": 0.0,
            "reopen_force_tolerance": 0.0,
            "reopen_width_tolerance": 0.0,
            "max_closed_cell_toggle_count": 1000000000,
            "max_closed_cell_toggle_fraction": 1.0,
            "__4_fracture_flow_mech_coupling": "--- coupling block C inside the fracture system ---",
            "drop_fluid_mech_linearization": true,
            "drop_tol_h": 100.0,
            "drop_tol_p": 1.0,
            "use_ad_pressure_assembly": true,
            "__5_state_guards": "--- robustness guards on the nonlinear state ---",
            "guard_state": true,
            "guard_max_abs_width": 1.0,
            "guard_max_abs_pressure": 100000000000.0,
            "__6_wellindex_coupling_stability": "--- OUTER coupling: when to update well CTF/perf-pressure from fracture, and how to damp/accelerate it ---",
            "ctf_change_threshold": "1e-15",
            "ctf_change_threshold_rel": 0.0,
            "ctf_relative_scale": 1.0,
            "perf_pressure_change_threshold": 1e+30,
            "perf_pressure_change_threshold_rel": 0.0,
            "perf_pressure_relative_scale": 1.0,
            "coupling_update_mode": "legacy",
            "damping_factor_perf": 0.0,
            "damping_factor_wi": 0.0,
            "enable_wi_coupling_update": false,
            "coupling_relax_init": 1.0,
            "coupling_relax_min": 0.05,
            "coupling_relax_max": 1.25,
            "coupling_acceptance_factor": 1.0,
            "coupling_fallback_legacy": true,
            "coupling_anderson_alpha_min": -1.0,
            "coupling_anderson_alpha_max": 2.0,
            "coupling_legacy_first_damp": false,
            "legacy_coupling_change_logic": false,
            "legacy_parent_setup_iteration": false,
            "__7_debug_output": "--- diagnostics (off for production) ---",
            "verbosity": "0",
            "dump_linear_system_on_failure": false,
            "write_coupled_linear_system": false,
            "linsolver": {
                "tol": "1e-10",
                "atol": "1e-20",
                "max_iter": "100",
                "verbosity": "0",
                "solver": "bicgstab",
                "preconditioner": {
                    "verbosity": 0,
                    "diag_mech": "false",
                    "diag_flow": "false",
                    "mech_first": false,
                    "fixed_stress": false,
                    "mode_policy": "manual",
                    "mode_switch_coupling_threshold": 0.25,
                    "mech_press_coupling": "true",
                    "update_mech_on_reuse": false,
                    "flow_solver": {
                        "solver": "umfpack",
                        "verbosity": "0",
                        "tol": "1",
                        "maxiter": 1,
                        "preconditioner": {
                            "type": "DILU",
                            "verbosity": 0
                        }
                    }
                }
            }
        },
        "reservoir": {
            "dist": "1",
            "calculate_dist": "true",
            "interpolate_stress": "false",
            "mobility": "0.0013",
            "perm": "1e-13"
        },
        "control": {
            "type": "perf_pressure",
            "rate": "0.029",
            "WI": "1e-6"
        }
    }
}
)opmjson";
}
} // namespace

namespace Opm::Detail {
bool isBuiltinFractureParamAlias(const std::string& filename)
{
    return filename == "standard" || filename == "sequential_implicit";
}

Opm::PropertyTree builtinFractureParam(const std::string& filename)
{
    if (filename == "standard") {
        return Opm::PropertyTree::fromJsonString(standardFractureParamJson());
    }

    if (filename == "sequential_implicit") {
        return Opm::PropertyTree::fromJsonString(sequentialImplicitFractureParamJson());
    }

    OPM_THROW(std::invalid_argument, "Unknown built-in fracture parameter alias: " + filename);
}
} // namespace Opm::Detail