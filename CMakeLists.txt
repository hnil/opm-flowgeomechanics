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

# Mandatory call to project
project(opm-flowgeomechanics C CXX)

cmake_minimum_required(VERSION 3.15)

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
