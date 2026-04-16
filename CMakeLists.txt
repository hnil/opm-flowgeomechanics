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
#find_package(CGAL REQUIRED)
#find_package(METIS REQUIRED)	
find_package(opm-common REQUIRED)

include(OpmInit)
OpmSetPolicies()

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

# not the same location as most of the other projects? this hook
# overrides
macro (dir_hook)
endmacro (dir_hook)

# project information is in dune.module. Read this file and set
# variables.  we cannot generate dune.module since it is read by
# dunecontrol before the build starts, so it makes sense to keep the
# data there then.
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)
include (OpmInit)

# Look for the opm-tests repository; if found the variable
# HAVE_OPM_TESTS will be set to true.
include(Findopm-tests)

# list of prerequisites for this particular project; this is in a
# separate file (in cmake/Modules sub-directory) because it is shared
# with the find module
include ("${project}-prereqs")

# read the list of components from this file (in the project
# directory); it should set various lists with the names of the files
# to include
include (CMakeLists_files.cmake)

macro (config_hook)
	opm_need_version_of ("dune-common")
	opm_need_version_of ("dune-istl")
	if(dune-fem_FOUND)
		opm_need_version_of ("dune-fem")
	endif()
	opm_need_version_of ("opm-simulators")

	if(USE_TRACY_PROFILER AND Tracy_FOUND)
		set(USE_TRACY 1)
	endif()

	include(UseDamaris)
endmacro (config_hook)

macro (prereqs_hook)
endmacro (prereqs_hook)

macro (fortran_hook)
endmacro (fortran_hook)

macro (files_hook)
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
	include (CMakeLists_files.cmake)
endmacro (files_hook)

macro (tests_hook)
endmacro (tests_hook)

# all setup common to the OPM library modules is done here
include (OpmLibMain)

include(OpmBashCompletion)

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

if(USE_TRACY_PROFILER AND Tracy_FOUND)
	target_link_libraries(flow_energy_geomech
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

####################33
function(add_mechtest_compareECLFiles)
  set(oneValueArgs CASENAME FILENAME SIMULATOR ABS_TOL REL_TOL DIR DIR_PREFIX PREFIX RESTART_STEP RESTART_SCHED)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  
  set(RESULT_PATH ${CMAKE_CURRENT_LIST_DIR}/result/${PARAM_CASENAME})
  set(TEST_ARGS ${PARAM_TEST_ARGS})
  set(DRIVER_ARGS -i ${CMAKE_CURRENT_LIST_DIR}/data
                  -r ${RESULT_PATH}
                  -b ${PROJECT_BINARY_DIR}/bin
                  -f ${PARAM_FILENAME}
                  -a ${PARAM_ABS_TOL}
                  -t ${PARAM_REL_TOL}
                  -c ${COMPARE_ECL_COMMAND}
				  -o ${PARAM_CASENAME})
 #                 -d ${RST_DECK_COMMAND})
 # if(PARAM_RESTART_STEP)
 #   list(APPEND DRIVER_ARGS -s ${PARAM_RESTART_STEP})
 # endif()
 # if(PARAM_RESTART_SCHED STREQUAL "false" OR PARAM_RESTART_SCHED STREQUAL "true")
 #   list(APPEND DRIVER_ARGS -h ${PARAM_RESTART_SCHED})
 # endif()
  opm_add_test(${PARAM_CASENAME}_${PARAM_SIMULATOR}+${PARAM_FILENAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
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
						 PREFIX geomech_frac)

add_mechtest_compareECLFiles(CASENAME REGTEST_geomech_frac_seq_bhp
                         FILENAME SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SEQ
                         SIMULATOR flow_energy_geomech
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR ${CMAKE_CURRENT_LIST_DIR}/data/
						 PREFIX geomech_frac_seq_bhp
						 TEST_ARGS --parameter-file=${CMAKE_CURRENT_LIST_DIR}/data/myparamsfrac.ini 
						--fracture-param-file=${CMAKE_CURRENT_LIST_DIR}/data/fracture_simple_bhp_seq_nverbose.json
						--enable-write-all-solutions=true
						--linear-solver-mech=${CMAKE_CURRENT_LIST_DIR}/data/amgmechsolverhypre_geosxcpu_new_nv.json
						--threads-per-process=4)

add_mechtest_compareECLFiles(CASENAME REGTEST_geomech_frac_seq_rate
                         FILENAME SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC_SEQ
                         SIMULATOR flow_energy_geomech
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR ${CMAKE_CURRENT_LIST_DIR}/data/
						 PREFIX geomech_frac_seq_rate
						 TEST_ARGS --parameter-file=${CMAKE_CURRENT_LIST_DIR}/data/myparamsfrac.ini 
						--fracture-param-file=${CMAKE_CURRENT_LIST_DIR}/data/fracture_simple_rate_seq_nverbose.json 
						--linear-solver-mech=${CMAKE_CURRENT_LIST_DIR}/data/amgmechsolverhypre_geosxcpu_new_nv.json
						--threads-per-process=4)




# fracture test

add_test(NAME Run.SIMPLE_MECH_FRAC_MPI.EnergyGeomech
	COMMAND
	$<SHELL_PATH:$<TARGET_FILE:flow_energy_geomech>>
	${CMAKE_CURRENT_LIST_DIR}/data/SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC.DATA
	--parameter-file=${CMAKE_CURRENT_LIST_DIR}/data/myparamsfrac.ini
	--output-dir=test_default
	--add-corners=true
	--linear-solver-mech=${CMAKE_CURRENT_LIST_DIR}/data/amgmechsolverhypre_geosxcpu_new_nv.json
)

set_tests_properties(Run.SIMPLE_MECH_FRAC_MPI.EnergyGeomech
	PROPERTIES
	PROCESSORS 4
)

add_test(NAME Run.SIMPLE_MECH_FRAC_THREAD.EnergyGeomech
	COMMAND
	$<SHELL_PATH:$<TARGET_FILE:flow_energy_geomech>>
	${CMAKE_CURRENT_LIST_DIR}/data/SIMPLE_MECH_NX_11_NY_11_NZ_25_FRAC.DATA
	--parameter-file=${CMAKE_CURRENT_LIST_DIR}/data/myparamsfrac.ini
	--output-dir=test_default
	--add-corners=true
	--threads-per-process=4
	--linear-solver-mech=${CMAKE_CURRENT_LIST_DIR}/data/amgmechsolverhypre_geosxcpu_new_nv.json
)

set_tests_properties(Run.SIMPLE_MECH_FRAC_THREAD.EnergyGeomech
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
