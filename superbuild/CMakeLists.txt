# Builds opm-common, opm-grid, opm-simulators and opm-flowgeomechanics
# as one CMake project.  Copy this file into the directory holding the
# four checkouts; see docs/GETTING_STARTED.md.
cmake_minimum_required(VERSION 3.15)
project(opm-geomech-superbuild CXX C)

foreach(TARGET opm-common opm-grid opm-models opm-simulators)
  set(${TARGET}_DIR ${CMAKE_BINARY_DIR}/${TARGET})
endforeach()

set(SIBLING_SEARCH 0)
set_property(GLOBAL PROPERTY ALLOW_DUPLICATE_CUSTOM_TARGETS 1)

enable_testing()
add_subdirectory(opm-common)
add_subdirectory(opm-grid)
add_dependencies(opmgrid opmcommon)
set(ADD_DISABLED_CTESTS 0)
set(BUILD_TESTING 1)
add_subdirectory(opm-simulators)
add_dependencies(opmsimulators opmgrid)
add_subdirectory(opm-flowgeomechanics)
add_dependencies(opm-flowgeomechanics_prepare opmsimulators)
