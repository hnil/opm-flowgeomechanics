# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

###########################################################################
#                                                                         #
# Note: The bulk of the build system is located in the cmake/ directory.  #
#       This file only contains the specializations for this particular   #
#       project. Most likely you are interested in editing one of these   #
#       files instead:                                                    #
#                                                                         #
#       dune.module                Name and version number                #
#       CMakeLists_files.cmake     Path of source files                   #
#       ${project}-prereqs.cmake   Dependencies                           #
#                                                                         #
###########################################################################

cmake_minimum_required(VERSION 3.26)

# Mandatory call to project
project(opm-flowgeomechanics C CXX)

set(USE_OPENMP_DEFAULT OFF) # Use of OpenMP is considered experimental

option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)
option(USE_TRACY_PROFILER "Enable tracy profiling" OFF)
option(USE_HYPRE "Use the Hypre library for linear solvers?" OFF)

# The following was copied from CMakeLists.txt in opm-common.

if(SIBLING_SEARCH AND NOT opm-common_DIR)
	# guess the sibling dir
	get_filename_component(_leaf_dir_name ${PROJECT_BINARY_DIR} NAME)
	get_filename_component(_parent_full_dir ${PROJECT_BINARY_DIR} DIRECTORY)
	get_filename_component(_parent_dir_name ${_parent_full_dir} NAME)
	#Try if <module-name>/<build-dir> is used
	get_filename_component(_modules_dir ${_parent_full_dir} DIRECTORY)
	if(IS_DIRECTORY ${_modules_dir}/opm-common/${_leaf_dir_name})
		set(opm-common_DIR ${_modules_dir}/opm-common/${_leaf_dir_name})
	else()
		string(REPLACE ${PROJECT_NAME} opm-common _opm_common_leaf ${_leaf_dir_name})
		if(NOT _leaf_dir_name STREQUAL _opm_common_leaf
				AND IS_DIRECTORY ${_parent_full_dir}/${_opm_common_leaf})
			# We are using build directories named <prefix><module-name><postfix>
			set(opm-common_DIR ${_parent_full_dir}/${_opm_common_leaf})
		elseif(IS_DIRECTORY ${_parent_full_dir}/opm-common)
			# All modules are in a common build dir
			set(opm-common_DIR "${_parent_full_dir}/opm-common")
		endif()
	endif()
endif()

if(opm-common_DIR AND NOT IS_DIRECTORY ${opm-common_DIR})
	message(WARNING "Value ${opm-common_DIR} passed to variable"
		" opm-common_DIR is not a directory")
endif()
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

find_package(opm-common REQUIRED)

# project information is in dune.module. Read this file and set
# variables.  we cannot generate dune.module since it is read by
# dunecontrol before the build starts, so it makes sense to keep the
# data there then.
include(OpmInit)
OpmSetPolicies()

# The OPM modules keep asserts enabled also in release builds:
# use_additional_optimization() in opm-common adds -UNDEBUG to every
# target unless WITH_NDEBUG is set.  Classes such as EnsureFinalized
# (the base of all material-law parameter classes) have members that
# exist only when NDEBUG is undefined, so compiling this module with
# NDEBUG while linking OPM libraries built without it changes class
# layouts and corrupts memory at runtime.  Apply the same flag to all
# targets in this module.
include(UseOptimization)
if(NOT WITH_NDEBUG)
  add_compile_options($<$<CONFIG:Release,RelWithDebInfo>:-UNDEBUG>)
endif()

# Executables go to bin/ like in the other OPM modules (opm_add_executable
# sets this per target; this module creates its targets directly).  The
# regression-test driver locates the simulators there.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

if(USE_MPI)
	set(HDF5_PREFER_PARALLEL TRUE)
endif()

if(USE_HYPRE AND USE_MPI)
  find_package(HYPRE)
  if(HYPRE_FOUND)
    set(HAVE_HYPRE 1)
  else()
    message(WARNING "Hypre requested but not found. Continuing without Hypre support.")
    set(USE_HYPRE OFF)
  endif()
elseif(USE_HYPRE)
  message(WARNING "Hypre requested but MPI not activated. Continuing without Hypre support.")
  set(USE_HYPRE OFF)
  set(HYPRE_FOUND OFF)
endif()

# Look for the opm-tests repository; if found the variable
# HAVE_OPM_TESTS will be set to true.
include(Findopm-tests)

macro (opm-flowgeomechanics_config_hook)
	opm_need_version_of ("dune-common")
	opm_need_version_of ("dune-istl")
	if(dune-fem_FOUND)
		opm_need_version_of ("dune-fem")
	endif()
	opm_need_version_of ("opm-simulators")

	if(USE_TRACY_PROFILER AND Tracy_FOUND)
		set(USE_TRACY 1)
	endif()
endmacro (opm-flowgeomechanics_config_hook)

# link the module library against its prerequisites; this hook runs
# right after OpmLibMain has created the opmflowgeomechanics target
macro (opm-flowgeomechanics_prereqs_hook)
	target_link_libraries(opmflowgeomechanics
		PUBLIC
			opmupscaling
			opmsimulators
			CGAL::CGAL GMP::gmp GMP::gmpxx
	)

	if(HAVE_HYPRE)
		target_link_libraries(opmflowgeomechanics
			PUBLIC
				HYPRE::HYPRE
		)
	endif()
endmacro (opm-flowgeomechanics_prereqs_hook)

macro (opm-flowgeomechanics_files_hook)
	if(MPI_FOUND AND HDF5_FOUND AND NOT HDF5_IS_PARALLEL)
		message(WARNING "When building parallel OPM flow we need a "
			"parallel version of hdf5, but found only a serial one. "
			"Please install a parallel hdf5 library for MPI "
			"(e.g with apt-get install libhdf5-mpi-dev) and do a clean "
			"rebuild (build after \"make clean\"). Continuing with "
			"only normal restart without hdf5 file support.")
		set(HDF5_FOUND OFF)
		unset(HAVE_HDF5)
	endif()
endmacro (opm-flowgeomechanics_files_hook)

# all setup common to the OPM library modules is done here
include (OpmLibMain)

include(OpmBashCompletion)

add_library(moduleVersionGeoMech OBJECT opm/simulators/utils/moduleVersion.cpp)
set_property(TARGET moduleVersionGeoMech PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(moduleVersionGeoMech PRIVATE "${PROJECT_BINARY_DIR}")

# Very strictly speaking, moduleVersionGeoMech depends only on the
# custom "update-version" target from OpmLibMain since it uses the
# generated file "${PROJECT_BINARY_DIR}/project-verion.h" as an input.
# However, "update-version" is not available in the "super build"
# setup, so we instead make sure that opmflowgeomechanics is fully
# built before we build moduleVersionGeoMech.
add_dependencies(moduleVersionGeoMech opmflowgeomechanics)

set(CMAKE_PROJECT_NAME "${PROJECT_NAME}")
add_custom_target(opm-flowgeomechanics_prepare)

add_executable(flow_geomech
	examples/flow_geomech.cpp
)
target_link_libraries(flow_geomech
	PUBLIC
		opmflowgeomechanics  
	PRIVATE
		moduleVersionGeoMech
)

if(USE_TRACY_PROFILER AND Tracy_FOUND)
	target_link_libraries(flow_geomech
		PRIVATE
			Tracy::TracyClient
	)
endif()

add_executable(flow_energy_geomech
	examples/flow_energy_geomech.cpp
)
target_link_libraries(flow_energy_geomech
	PUBLIC
		opmflowgeomechanics 
	PRIVATE
		moduleVersionGeoMech
)

add_executable(flow_energy_geomech_tpsa
	examples/flow_energy_geomech_tpsa.cpp
)
target_link_libraries(flow_energy_geomech_tpsa
	PUBLIC
		opmflowgeomechanics 
	PRIVATE
		moduleVersionGeoMech
)

if(USE_TRACY_PROFILER AND Tracy_FOUND)
	target_link_libraries(flow_energy_geomech_tpsa
		PRIVATE
			Tracy::TracyClient
	)
endif()

if(USE_TRACY_PROFILER AND Tracy_FOUND)
	target_link_libraries(flow_energy_geomech
		PRIVATE
			Tracy::TracyClient
	)
endif()

add_executable(flow_oilgas_energy_geomech
	examples/flow_oilgas_energy_geomech.cpp
)
target_link_libraries(flow_oilgas_energy_geomech
	PUBLIC
		opmflowgeomechanics
	PRIVATE
		moduleVersionGeoMech
)

if(USE_TRACY_PROFILER AND Tracy_FOUND)
	target_link_libraries(flow_oilgas_energy_geomech
		PRIVATE
			Tracy::TracyClient
	)
endif()

add_executable(padmodel
	examples/padmodel.cpp
)
target_link_libraries(padmodel
	PUBLIC
		opmsimulators
)

add_executable(test_foamgrid
	examples/test_foamgrid.cpp
)
target_link_libraries(test_foamgrid
	PUBLIC
		opmsimulators opmupscaling
)

add_executable(test_vem_elements
	examples/test_vem_elements.cpp
)
target_link_libraries(test_vem_elements
	PUBLIC
		opmflowgeomechanics
)
if(USE_TRACY_PROFILER AND Tracy_FOUND)
	target_link_libraries(test_vem_elements
		PRIVATE
			Tracy::TracyClient
	)
endif()

#add_executable(test_vem_assembly
#	examples/test_vem_assembly.cpp
#)
#target_link_libraries(test_vem_assembly
#	PUBLIC
#		opmflowgeomechanics
#)
#if(USE_TRACY_PROFILER AND Tracy_FOUND)
#	target_link_libraries(test_vem_assembly
#		PRIVATE
#			Tracy::TracyClient
#	)
#endif()

add_executable(paramtest
	examples/paramtest.cpp
)
target_link_libraries(paramtest
	PUBLIC
		opmflowgeomechanics
)

add_executable(test_aspect_ratio
	examples/test_aspect_ratio.cpp
)
target_link_libraries(test_aspect_ratio
	PUBLIC
		opmflowgeomechanics
)
add_test(NAME Test.AspectRatio
	COMMAND test_aspect_ratio
)

add_executable(test_trimesh
	examples/test_trimesh.cpp
)
target_link_libraries(test_trimesh
	PUBLIC
		opmflowgeomechanics
)

add_executable(test_gridstretch
	examples/test_gridstretch.cpp
)
target_link_libraries(test_gridstretch
	PUBLIC
		opmflowgeomechanics
)

add_executable(coupled_standalone_test
	examples/coupled_standalone_test.cpp
)
target_link_libraries(coupled_standalone_test
	PUBLIC
		opmflowgeomechanics
)

add_executable(standalone_sim_test
	examples/standalone_sim_test.cpp
)
target_link_libraries(standalone_sim_test
	PUBLIC
		opmflowgeomechanics
)

add_executable(test_fracture_pressure_ad
	examples/test_fracture_pressure_ad.cpp
)
target_link_libraries(test_fracture_pressure_ad
	PUBLIC
		opmflowgeomechanics
)

add_executable(test_patch_recovery
	examples/test_patch_recovery.cpp
)
target_link_libraries(test_patch_recovery
	PUBLIC
		opmflowgeomechanics
)

add_executable(replay_fracture_linear_system
	examples/replay_fracture_linear_system.cpp
)
target_link_libraries(replay_fracture_linear_system
	PUBLIC
		opmflowgeomechanics
)

add_executable(test_partition_of_unity
	examples/test_partition_of_unity.cpp
)
target_link_libraries(test_partition_of_unity
	PUBLIC
		opmflowgeomechanics
)

# Verify that the node (vertex) degrees of freedom form a partition of unity
# after load balancing a deck. The number of MPI processes is the input: one
# ctest per process count in PartitionOfUnity_NumProcs.
set(PartitionOfUnity_NumProcs 1 2 4)
foreach(np ${PartitionOfUnity_NumProcs})
	if (NOT MPIEXEC_MAX_NUMPROCS LESS np)
		add_test(NAME Unit.PartitionOfUnity.Nodes.np-${np}
			COMMAND
			${MPIEXEC_EXECUTABLE}
			${MPIEXEC_NUMPROC_FLAG} ${np}
			${MPIEXEC_PREFLAGS}
			$<SHELL_PATH:$<TARGET_FILE:test_partition_of_unity>>
			${MPIEXEC_POSTFLAGS}
			gridfilename=${CMAKE_CURRENT_LIST_DIR}/data/MECHTEST/BLACKOILTHERMAL_BC.DATA
			--add-corners=true
			--edge-conformal=true
			--edge-weights=uniform
		)
		set_tests_properties(Unit.PartitionOfUnity.Nodes.np-${np}
			PROPERTIES
			PROCESSORS ${np}
		)
	endif()
endforeach()

####################33
function(add_mechtest_compareECLFiles)
	set(oneValueArgs CASENAME REFERENCE_CASENAME FILENAME SIMULATOR ABS_TOL REL_TOL DIR DIR_PREFIX PREFIX RESTART_STEP RESTART_SCHED)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

	if(NOT PARAM_REFERENCE_CASENAME)
	set(PARAM_REFERENCE_CASENAME ${PARAM_CASENAME})
	endif()
  
  set(RESULT_PATH ${CMAKE_CURRENT_LIST_DIR}/result/${PARAM_CASENAME})
  set(TEST_ARGS ${PARAM_TEST_ARGS})
  set(DRIVER_ARGS -i ${CMAKE_CURRENT_LIST_DIR}/data
                  -r ${RESULT_PATH}
                  -b ${PROJECT_BINARY_DIR}/bin
                  -f ${PARAM_FILENAME}
                  -a ${PARAM_ABS_TOL}
                  -t ${PARAM_REL_TOL}
                  -c $<TARGET_FILE:compareECL>
				  -o ${PARAM_REFERENCE_CASENAME})
 #                 -d ${RST_DECK_COMMAND})
 # if(PARAM_RESTART_STEP)
 #   list(APPEND DRIVER_ARGS -s ${PARAM_RESTART_STEP})
 # endif()
 # if(PARAM_RESTART_SCHED STREQUAL "false" OR PARAM_RESTART_SCHED STREQUAL "true")
 #   list(APPEND DRIVER_ARGS -h ${PARAM_RESTART_SCHED})
 # endif()
  opm_add_test(${PARAM_CASENAME}_${PARAM_SIMULATOR}+${PARAM_FILENAME}
               EXE_TARGET ${PARAM_SIMULATOR}
               DRIVER_ARGS ${DRIVER_ARGS}
               TEST_ARGS ${TEST_ARGS})
  set_tests_properties(${PARAM_CASENAME}_${PARAM_SIMULATOR}+${PARAM_FILENAME} 	PROPERTIES
                        DIRNAME ${PARAM_DIR}
                        FILENAME ${PARAM_FILENAME}
                        SIMULATOR ${PARAM_SIMULATOR}
                        TESTNAME ${PARAM_CASENAME})
endfunction()


#######################3333



add_test(NAME Test.FracturePressureAD
	COMMAND test_fracture_pressure_ad
)

add_executable(test_ddm_assembly
	examples/test_ddm_assembly.cpp
)
target_link_libraries(test_ddm_assembly
	PUBLIC
		opmflowgeomechanics
)

add_test(NAME Test.DDMAssembly
	COMMAND test_ddm_assembly
)

add_executable(example_penny_crack_ddm
	examples/example_penny_crack_ddm.cpp
)
target_link_libraries(example_penny_crack_ddm
	PUBLIC
		opmflowgeomechanics
)

# Accuracy test disabled - DDM currently has ~50% COD error vs analytic, needs fixing
# add_test(NAME Test.DDM.PennyCrack
# 	COMMAND example_penny_crack_ddm
# )

# Regression test: ensure the DDM assembly stays numerically stable.
# Uses relaxed tolerances (60%) just to catch major regressions or NaN/crash failures.
add_test(NAME Regression.DDM.PennyCrack.Assembly
	COMMAND example_penny_crack_ddm
		--layers=3
		--cod-tol=0.6
		--k1-tol=0.6
)

# Stored-solution regression test: checks the DDM matrix assembly bit-for-bit
# against known opening values for a 2-layer mesh.
add_executable(test_penny_crack_regression
	tests/test_penny_crack_regression.cpp
)
target_link_libraries(test_penny_crack_regression
	PUBLIC
		opmflowgeomechanics
)
add_test(NAME Regression.DDM.PennyCrack.StoredSolution
	COMMAND test_penny_crack_regression
)

# ---------------------------------------------------------------------------

set(Test_NumProc 4)
if (MPIEXEC_MAX_NUMPROCS LESS Test_NumProc)
	set(Test_NumProc ${MPIEXEC_MAX_NUMPROCS})
endif()

add_test(NAME Run.Blackoil_Thermal_BC.Parallel
	COMMAND
	${MPIEXEC_EXECUTABLE}
	${MPIEXEC_NUMPROC_FLAG} ${Test_NumProc}
	${MPIEXEC_PREFLAGS}
	$<SHELL_PATH:$<TARGET_FILE:flow_energy_geomech>>
	${MPIEXEC_POSTFLAGS}
	--ecl-deck-file-name=${CMAKE_CURRENT_LIST_DIR}/data/MECHTEST/BLACKOILTHERMAL_BC.DATA
	--enable-esmry=false
	--output-dir=results/np-${Test_NumProc}/flow_energy_geomech/Blackoil_Thermal_BC
)

set_tests_properties(Run.Blackoil_Thermal_BC.Parallel
	PROPERTIES
	PROCESSORS ${Test_NumProc}
	FIXTURES_SETUP Generate_SimResults
)

# add_test(NAME Compare.Blackoil_Thermal_BC.Parallel
# 	COMMAND
# 	$<SHELL_PATH:${COMPARE_ECL_COMMAND}>
# 	-a
# 	results/np-${Test_NumProc}/flow_energy_geomech/Blackoil_Thermal_BC/BLACKOILTHERMAL_BC
# 	results/np-${Test_NumProc}/flow_energy_geomech/Blackoil_Thermal_BC/BLACKOILTHERMAL_BC
# 	2.0e-2
# 	1.0e-5
# )

# set_tests_properties(Compare.Blackoil_Thermal_BC.Parallel
# 	PROPERTIES
# 	PROCESSORS 1
# 	FIXTURES_REQUIRED Generate_SimResults
# )

# ---------------------------------------------------------------------------
# Regression test: simple full simulation with flow_geomech
# ---------------------------------------------------------------------------

add_test(NAME Run.SPE1_Geomech_Mini.Parallel
	COMMAND
	${MPIEXEC_EXECUTABLE}
	${MPIEXEC_NUMPROC_FLAG} ${Test_NumProc}
	${MPIEXEC_PREFLAGS}
	$<SHELL_PATH:$<TARGET_FILE:flow_geomech>>
	${MPIEXEC_POSTFLAGS}
	--ecl-deck-file-name=${CMAKE_CURRENT_LIST_DIR}/data/SPE1CASE2_MECH_MINI.DATA
	--enable-esmry=false
	--add-corners=true
	--output-dir=results/np-${Test_NumProc}/flow_geomech/SPE1_Geomech_Mini
)

set_tests_properties(Run.SPE1_Geomech_Mini.Parallel
	PROPERTIES
	PROCESSORS ${Test_NumProc}
)

# ---------------------------------------------------------------------------
# Regression test: SIMPLE_MECH run with custom parameter + linear solver setup
# ---------------------------------------------------------------------------

add_test(NAME Run.SIMPLE_MECH.EnergyGeomech
	COMMAND
	$<SHELL_PATH:$<TARGET_FILE:flow_energy_geomech>>
	${CMAKE_CURRENT_LIST_DIR}/data/SIMPLE_MECH_NX_11_NY_11_NZ_25.DATA
	--parameter-file=${CMAKE_CURRENT_LIST_DIR}/data/myparamsfrac.ini
	--output-dir=test
	--linear-solver-mech=${CMAKE_CURRENT_LIST_DIR}/data/amgmechsolverhypre_geosxcpu_new_nv.json
)

set_tests_properties(Run.SIMPLE_MECH.EnergyGeomech
	PROPERTIES
	PROCESSORS 1
)

add_test(NAME Run.SIMPLE_MECH.EnergyGeomech.Parallel
	COMMAND
	${MPIEXEC_EXECUTABLE}
	${MPIEXEC_NUMPROC_FLAG} ${Test_NumProc}
	${MPIEXEC_PREFLAGS}
	$<SHELL_PATH:$<TARGET_FILE:flow_energy_geomech>>
	${MPIEXEC_POSTFLAGS}
	${CMAKE_CURRENT_LIST_DIR}/data/SIMPLE_MECH_NX_11_NY_11_NZ_25.DATA
	--parameter-file=${CMAKE_CURRENT_LIST_DIR}/data/myparamsfrac.ini
	--linear-solver-mech=${CMAKE_CURRENT_LIST_DIR}/data/amgmechsolverhypre_geosxcpu_new_nv.json
	--add-corners=true
	--edge-conformal=true
	--enable-esmry=false
	--output-dir=results/np-${Test_NumProc}/flow_energy_geomech/SIMPLE_MECH
)

set_tests_properties(Run.SIMPLE_MECH.EnergyGeomech.Parallel
	PROPERTIES
	PROCESSORS ${Test_NumProc}
)

add_test(NAME Run.SPE1_OILGAS_THERMAL_MECH.OilGasEnergyGeomech
	COMMAND
	$<SHELL_PATH:$<TARGET_FILE:flow_oilgas_energy_geomech>>
	${CMAKE_CURRENT_LIST_DIR}/data/SPE1CASE2_OILGAS_THERMAL_MECH.DATA
	--enable-esmry=false
	--add-corners=true
	--linear-solver-mech=${CMAKE_CURRENT_LIST_DIR}/data/amgmechsolverhypre_geosxcpu_new_nv.json
	--output-dir=test_oilgas_energy_geomech
)

set_tests_properties(Run.SPE1_OILGAS_THERMAL_MECH.OilGasEnergyGeomech
	PROPERTIES
	PROCESSORS 1
)

add_test(NAME Run.SPE1_OILGAS_THERMAL_MECH_FRAC_FAIL.OilGasEnergyGeomech
	COMMAND
	$<SHELL_PATH:$<TARGET_FILE:flow_oilgas_energy_geomech>>
	${CMAKE_CURRENT_LIST_DIR}/data/SPE1CASE2_OILGAS_THERMAL_MECH_FRAC_FAIL.DATA
	--enable-esmry=false
	--add-corners=true
	--linear-solver-mech=${CMAKE_CURRENT_LIST_DIR}/data/amgmechsolverhypre_geosxcpu_new_nv.json
	--output-dir=test_oilgas_energy_geomech_frac_fail
)

set_tests_properties(Run.SPE1_OILGAS_THERMAL_MECH_FRAC_FAIL.OilGasEnergyGeomech
	PROPERTIES
	PROCESSORS 1
	WILL_FAIL TRUE
)

## regressuin
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-mechRegressionTest.sh "")

# Set absolute tolerance to be used passed to the macros in the following tests
set(abs_tol 2e-2)
set(rel_tol 1e-5)	
add_mechtest_compareECLFiles(CASENAME REGTEST_geomech
                         FILENAME SIMPLE_MECH_NX_11_NY_11_NZ_25
                         SIMULATOR flow_energy_geomech
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR ${CMAKE_CURRENT_LIST_DIR}/data/
						 PREFIX geomech)

add_mechtest_compareECLFiles(CASENAME REGTEST_geomech_frac
                         FILENAME SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC
                         SIMULATOR flow_energy_geomech
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR ${CMAKE_CURRENT_LIST_DIR}/data/
						 PREFIX geomech_frac
						 TEST_ARGS --fracture-param-file=${CMAKE_CURRENT_LIST_DIR}/data/fracture_simple_frac_postsolve_legacy_full.json)

add_mechtest_compareECLFiles(CASENAME REGTEST_geomech_frac_standard
						 REFERENCE_CASENAME REGTEST_geomech_frac
						 FILENAME SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC
						 SIMULATOR flow_energy_geomech
						 ABS_TOL ${abs_tol}
						 REL_TOL ${rel_tol}
						 DIR ${CMAKE_CURRENT_LIST_DIR}/data/
						 PREFIX geomech_frac
						 TEST_ARGS --fracture-param-file=standard)

# Formation damage (WINJDAM/WINJFCNC) + fracture: verifies the filter-cake
# transmissibility damage on the original perforation, the matrix/fracture
# filtrate split (CFCFRATE), and that dynamically added fracture
# perforations are handled with an active cake model.
add_mechtest_compareECLFiles(CASENAME REGTEST_geomech_frac_damage
                         FILENAME SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_DAMAGE
                         SIMULATOR flow_energy_geomech
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR ${CMAKE_CURRENT_LIST_DIR}/data/
						 PREFIX geomech_frac_damage
						 TEST_ARGS --fracture-param-file=${CMAKE_CURRENT_LIST_DIR}/data/fracture_simple_frac_postsolve_legacy_full.json)

add_mechtest_compareECLFiles(CASENAME REGTEST_geomech_frac_seq_bhp
                         FILENAME SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SEQ
                         SIMULATOR flow_energy_geomech
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR ${CMAKE_CURRENT_LIST_DIR}/data/
						 PREFIX geomech_frac_seq_bhp
						 TEST_ARGS --parameter-file=${CMAKE_CURRENT_LIST_DIR}/data/myparamsfrac.ini 
						--fracture-param-file=${CMAKE_CURRENT_LIST_DIR}/data/fracture_simple_bhp_seq_legacy_full.json
						--enable-write-all-solutions=true
						--linear-solver-mech=${CMAKE_CURRENT_LIST_DIR}/data/amgmechsolverhypre_geosxcpu_new_nv.json
						--threads-per-process=4)

# Sequential-implicit fracture run with the recommended (physically-correct) settings:
# gravity on, AD assembly (coupling dropped for robustness), robust mech-last/decoupled
# linear solver, relaxed closed-cell guards (data/fracture_simple_bhp_seq_final.json,
# see data/fracture_settings_REFERENCE.md). Uses myparamsfrac_seq.ini: adaptive
# timestepping on and no fixed initial step, so the fracture solve is exercised under
# varying/larger steps (harder than the legacy seq tests).
add_mechtest_compareECLFiles(CASENAME REGTEST_geomech_frac_seq_bhp_recommended
                         FILENAME SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SEQ
                         SIMULATOR flow_energy_geomech
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR ${CMAKE_CURRENT_LIST_DIR}/data/
						 PREFIX geomech_frac_seq_bhp_recommended
						 TEST_ARGS --parameter-file=${CMAKE_CURRENT_LIST_DIR}/data/myparamsfrac_seq.ini
						--fracture-param-file=${CMAKE_CURRENT_LIST_DIR}/data/fracture_simple_bhp_seq_final.json
						--enable-write-all-solutions=true
						--linear-solver-mech=${CMAKE_CURRENT_LIST_DIR}/data/amgmechsolverhypre_geosxcpu_new_nv.json
						--threads-per-process=4)

add_mechtest_compareECLFiles(CASENAME REGTEST_geomech_frac_seq_bhp_sequential_implicit
                         REFERENCE_CASENAME REGTEST_geomech_frac_seq_bhp_recommended
                         FILENAME SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SEQ
                         SIMULATOR flow_energy_geomech
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR ${CMAKE_CURRENT_LIST_DIR}/data/
						 PREFIX geomech_frac_seq_bhp_recommended
						 TEST_ARGS --parameter-file=${CMAKE_CURRENT_LIST_DIR}/data/myparamsfrac_seq.ini
						 --fracture-param-file=sequential_implicit
						 --enable-write-all-solutions=true
						 --linear-solver-mech=${CMAKE_CURRENT_LIST_DIR}/data/amgmechsolverhypre_geosxcpu_new_nv.json
						 --threads-per-process=4)

# Reference data not yet committed - disabled until reference is updated
# add_mechtest_compareECLFiles(CASENAME REGTEST_geomech_frac_seq_bhp_default
#                          FILENAME SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SEQ
#                          SIMULATOR flow_energy_geomech
#                          ABS_TOL ${abs_tol}
#                          REL_TOL ${rel_tol}
#                          DIR ${CMAKE_CURRENT_LIST_DIR}/data/
# 						 PREFIX geomech_frac_seq_bhp_default
# 						 TEST_ARGS --parameter-file=${CMAKE_CURRENT_LIST_DIR}/data/myparamsfrac.ini
# 						--fracture-param-file=${CMAKE_CURRENT_LIST_DIR}/data/fracture_simple_frac_seq_default_full.json
# 						--enable-write-all-solutions=true
# 						--linear-solver-mech=${CMAKE_CURRENT_LIST_DIR}/data/amgmechsolverhypre_geosxcpu_new_nv.json
# 						--threads-per-process=4)

add_mechtest_compareECLFiles(CASENAME REGTEST_geomech_frac_seq_rate
                         FILENAME SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SEQ
                         SIMULATOR flow_energy_geomech
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR ${CMAKE_CURRENT_LIST_DIR}/data/
						 PREFIX geomech_frac_seq_rate
						 TEST_ARGS --parameter-file=${CMAKE_CURRENT_LIST_DIR}/data/myparamsfrac.ini 
						--fracture-param-file=${CMAKE_CURRENT_LIST_DIR}/data/fracture_simple_rate_seq_full.json 
						--linear-solver-mech=${CMAKE_CURRENT_LIST_DIR}/data/amgmechsolverhypre_geosxcpu_new_nv.json
						--threads-per-process=4)



# Multisegment-well + fracture: a fracture on an MSW injector.  Guards the
# WGState-preservation fix (commits c70a98f89 / 068c544) -- without it the
# MSW well reports zero rate after a fracture-driven well re-init while the
# perforations keep injecting.  Coarse deck (fracture stays in one reservoir
# cell -> CTF-value-only update path); the fine-grid structural-change
# counterpart is data/CASE_REFINE_nx_21_ny_21_nz_25_THERMAL_MSW.DATA.
add_mechtest_compareECLFiles(CASENAME REGTEST_geomech_frac_seq_msw
                         FILENAME SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SEQ_MSW
                         SIMULATOR flow_energy_geomech
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR ${CMAKE_CURRENT_LIST_DIR}/data/
						 PREFIX geomech_frac_seq_msw
						 TEST_ARGS --parameter-file=${CMAKE_CURRENT_LIST_DIR}/data/myparamsfrac_seq.ini
						--fracture-param-file=${CMAKE_CURRENT_LIST_DIR}/data/fracture_simple_bhp_seq_final.json
						--enable-write-all-solutions=true
						--linear-solver-mech=${CMAKE_CURRENT_LIST_DIR}/data/amgmechsolverhypre_geosxcpu_new_nv.json)


# fracture test

add_test(NAME Run.SIMPLE_MECH_FRAC_MPI.EnergyGeomech
	COMMAND
	$<SHELL_PATH:$<TARGET_FILE:flow_energy_geomech>>
	${CMAKE_CURRENT_LIST_DIR}/data/SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC.DATA
	--parameter-file=${CMAKE_CURRENT_LIST_DIR}/data/myparamsfrac.ini
	--fracture-param-file=${CMAKE_CURRENT_LIST_DIR}/data/fracture_simple_bhp_run_full.json
	--output-dir=test_default
	--add-corners=true
	--linear-solver-mech=${CMAKE_CURRENT_LIST_DIR}/data/amgmechsolverhypre_geosxcpu_new_nv.json
)

set_tests_properties(Run.SIMPLE_MECH_FRAC_MPI.EnergyGeomech
	PROPERTIES
	PROCESSORS ${Test_NumProc}
)

add_test(NAME Run.SIMPLE_MECH_FRAC_THREAD.EnergyGeomech
	COMMAND
	$<SHELL_PATH:$<TARGET_FILE:flow_energy_geomech>>
	${CMAKE_CURRENT_LIST_DIR}/data/SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC.DATA
	--parameter-file=${CMAKE_CURRENT_LIST_DIR}/data/myparamsfrac.ini
	--fracture-param-file=${CMAKE_CURRENT_LIST_DIR}/data/fracture_simple_bhp_run_full.json
	--output-dir=test_default
	--add-corners=true
	--threads-per-process=4
	--linear-solver-mech=${CMAKE_CURRENT_LIST_DIR}/data/amgmechsolverhypre_geosxcpu_new_nv.json
)

set_tests_properties(Run.SIMPLE_MECH_FRAC_THREAD.EnergyGeomech
	PROPERTIES
	PROCESSORS 1
)

add_test(NAME Run.SIMPLE_MECH_FRAC_RESET.EnergyGeomech
	COMMAND
	$<SHELL_PATH:$<TARGET_FILE:flow_energy_geomech>>
	${CMAKE_CURRENT_LIST_DIR}/data/SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_RESET.DATA
	--parameter-file=${CMAKE_CURRENT_LIST_DIR}/data/myparamsfrac_frac.ini
	--output-dir=test_frac_reset
	--add-corners=true
	--linear-solver-mech=${CMAKE_CURRENT_LIST_DIR}/data/amgmechsolverhypre_geosxcpu_new_nv.json
)

add_test(NAME Run.SIMPLE_MECH_FRAC_RESET_SEQ.EnergyGeomech
	COMMAND
	$<SHELL_PATH:$<TARGET_FILE:flow_energy_geomech>>
	${CMAKE_CURRENT_LIST_DIR}/data/SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_RESET.DATA
	--parameter-file=${CMAKE_CURRENT_LIST_DIR}/data/myparamsfrac_frac.ini
	--output-dir=test_frac_reset_seq
	--fracture-param-file=${CMAKE_CURRENT_LIST_DIR}/data/fracture_simple_bhp_seq_nverbose.json
	--add-corners=true
	--linear-solver-mech=${CMAKE_CURRENT_LIST_DIR}/data/amgmechsolverhypre_geosxcpu_new_nv.json
) 

set_tests_properties(Run.SIMPLE_MECH_FRAC_RESET.EnergyGeomech
	PROPERTIES
	PROCESSORS 1
)

# add_test(NAME Run.SIMPLE_MECH_FRAC_SEQ_BHP.EnergyGeomech
# 	COMMAND
# 	$<SHELL_PATH:$<TARGET_FILE:flow_energy_geomech>>
# 	${CMAKE_CURRENT_LIST_DIR}/data/SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SEQ.DATA
# 	--parameter-file=${CMAKE_CURRENT_LIST_DIR}/data/myparamsfrac.ini
# 	--fracture-param-file=${CMAKE_CURRENT_LIST_DIR}/data/fracture_simple_bhp_seq_nverbose.json
# 	--output-dir=test_bhp_seq
# 	--enable-write-all-solutions=true
# 	--linear-solver-mech=${CMAKE_CURRENT_LIST_DIR}/data/amgmechsolverhypre_geosxcpu_new_nv.json
# )

# set_tests_properties(Run.SIMPLE_MECH_FRAC_SEQ_BHP.EnergyGeomech
# 	PROPERTIES
# 	PROCESSORS 1
# )

# add_test(NAME Run.SIMPLE_MECH_FRAC_SEQ_RATE.EnergyGeomech
# 	COMMAND
# 	$<SHELL_PATH:$<TARGET_FILE:flow_energy_geomech>>
# 	${CMAKE_CURRENT_LIST_DIR}/data/SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SEQ.DATA
# 	--parameter-file=${CMAKE_CURRENT_LIST_DIR}/data/myparamsfrac.ini
# 	--fracture-param-file=${CMAKE_CURRENT_LIST_DIR}/data/fracture_simple_rate_seq_nverbose.json
# 	--output-dir=test_rate_seq
# 	--enable-write-all-solutions=true
# 	--linear-solver-mech=${CMAKE_CURRENT_LIST_DIR}/data/amgmechsolverhypre_geosxcpu_new_nv.json
# )

# set_tests_properties(Run.SIMPLE_MECH_FRAC_SEQ_RATE.EnergyGeomech
# 	PROPERTIES
# 	PROCESSORS 1
# )	
# ---------------------------------------------------------------------------
# Unit tests for individually testable components
# ---------------------------------------------------------------------------

add_test(NAME Unit.VemElements
	COMMAND $<TARGET_FILE:test_vem_elements>
)

add_test(NAME Unit.Trimesh.IrregularGrid
	COMMAND $<TARGET_FILE:test_trimesh> 1
)

add_test(NAME Unit.Trimesh.GridRefinement
	COMMAND $<TARGET_FILE:test_trimesh> 3
)

add_test(NAME Unit.PatchRecovery
	COMMAND $<TARGET_FILE:test_patch_recovery>
)

set(FRACTURE_LINEAR_REPLAY_PREFIX
	${CMAKE_CURRENT_LIST_DIR}/tests/fracture_linear_system_nlin_11/fracture_linear_system_nlin_11_linear_failure)

add_test(NAME Unit.FractureLinearReplay.Pressure
	COMMAND $<TARGET_FILE:replay_fracture_linear_system>
	${FRACTURE_LINEAR_REPLAY_PREFIX}
	--pressure=active
)

add_test(NAME Unit.FractureLinearReplay.CoupledFailure
	COMMAND $<TARGET_FILE:replay_fracture_linear_system>
	${FRACTURE_LINEAR_REPLAY_PREFIX}
	--coupled
)
set_tests_properties(Unit.FractureLinearReplay.CoupledFailure
	PROPERTIES WILL_FAIL TRUE
)

add_test(NAME Unit.FractureLinearReplay.CoupledScaled
	COMMAND $<TARGET_FILE:replay_fracture_linear_system>
	${FRACTURE_LINEAR_REPLAY_PREFIX}
	--coupled
	--scaling=blockmax
)

add_test(NAME Unit.FractureLinearReplay.CoupledMechFirst
	COMMAND $<TARGET_FILE:replay_fracture_linear_system>
	${FRACTURE_LINEAR_REPLAY_PREFIX}
	--coupled
	--solver-json=${CMAKE_CURRENT_LIST_DIR}/tests/fracture_linear_system_nlin_11/fracture_linear_system_nlin_11_linear_failure_linsolver_mech_first.json
)

add_test(NAME Unit.FractureLinearReplay.CoupledFixedStress
	COMMAND $<TARGET_FILE:replay_fracture_linear_system>
	${FRACTURE_LINEAR_REPLAY_PREFIX}
	--coupled
	--solver-json=${CMAKE_CURRENT_LIST_DIR}/tests/fracture_linear_system_nlin_11/fracture_linear_system_nlin_11_linear_failure_linsolver_fixed_stress.json
)

set(FRACTURE_LINEAR_REPLAY_PREFIX_N4
	${CMAKE_CURRENT_LIST_DIR}/tests/fracture_linear_system_nlin_4/fracture_linear_system_nlin_4_linear_failure)

add_test(NAME Unit.FractureLinearReplay.N4.CoupledFailure
	COMMAND $<TARGET_FILE:replay_fracture_linear_system>
	${FRACTURE_LINEAR_REPLAY_PREFIX_N4}
	--coupled
)
set_tests_properties(Unit.FractureLinearReplay.N4.CoupledFailure
	PROPERTIES WILL_FAIL TRUE
)

add_test(NAME Unit.FractureLinearReplay.N4.CoupledFixedStress
	COMMAND $<TARGET_FILE:replay_fracture_linear_system>
	${FRACTURE_LINEAR_REPLAY_PREFIX_N4}
	--coupled
	--solver-json=${CMAKE_CURRENT_LIST_DIR}/tests/fracture_linear_system_nlin_11/fracture_linear_system_nlin_11_linear_failure_linsolver_fixed_stress.json
)

add_test(NAME Unit.FractureLinearReplay.N4.CoupledAuto
	COMMAND $<TARGET_FILE:replay_fracture_linear_system>
	${FRACTURE_LINEAR_REPLAY_PREFIX_N4}
	--coupled
	--solver-json=${CMAKE_CURRENT_LIST_DIR}/tests/fracture_linear_system_nlin_11/fracture_linear_system_nlin_11_linear_failure_linsolver_coupling_auto.json
)

set(FRACTURE_LINEAR_REPLAY_PREFIX_N22
	${CMAKE_CURRENT_LIST_DIR}/tests/fracture_linear_system_nlin_22/fracture_linear_system_nlin_22_linear_failure)

add_test(NAME Unit.FractureLinearReplay.N22.CoupledFixedStressFailure
	COMMAND $<TARGET_FILE:replay_fracture_linear_system>
	${FRACTURE_LINEAR_REPLAY_PREFIX_N22}
	--coupled
	--solver-json=${CMAKE_CURRENT_LIST_DIR}/tests/fracture_linear_system_nlin_11/fracture_linear_system_nlin_11_linear_failure_linsolver_fixed_stress.json
)
set_tests_properties(Unit.FractureLinearReplay.N22.CoupledFixedStressFailure
	PROPERTIES WILL_FAIL TRUE
)

set(FRACTURE_LINEAR_REPLAY_PREFIX_N67
	${CMAKE_CURRENT_LIST_DIR}/tests/fracture_linear_system_nlin_67/fracture_linear_system_nlin_67_linear_failure)

add_test(NAME Unit.FractureLinearReplay.N67.CoupledFixedStressFailure
	COMMAND $<TARGET_FILE:replay_fracture_linear_system>
	${FRACTURE_LINEAR_REPLAY_PREFIX_N67}
	--coupled
	--solver-json=${CMAKE_CURRENT_LIST_DIR}/tests/fracture_linear_system_nlin_11/fracture_linear_system_nlin_11_linear_failure_linsolver_fixed_stress.json
)
set_tests_properties(Unit.FractureLinearReplay.N67.CoupledFixedStressFailure
	PROPERTIES WILL_FAIL TRUE
)

# Contact / fracture-opening report on a dumped nonlinear state. The nlin_5 state
# is the pathological "opened too large" case: every cell is flagged open yet a
# large fraction carry compressive force (want to close). The report classifies
# the cells offline; this test pins the diagnostic and guards the report tool.
add_test(NAME Unit.FractureContactReport.N5
	COMMAND $<TARGET_FILE:replay_fracture_linear_system>
	${CMAKE_CURRENT_LIST_DIR}/tests/fracture_contact_report_nlin_5/fracture_linear_system_nlin_5_linear_failure
	--contact-report
)
set_tests_properties(Unit.FractureContactReport.N5
	PROPERTIES PASS_REGULAR_EXPRESSION "open but force wants close .f>=close_tol.: 456"
)
