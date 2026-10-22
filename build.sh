#!/bin/bash
# Build script for opm-flowgeomechanics and its OPM dependencies.
# Tested on Ubuntu 24.04 LTS.
#
# This script should be run from a workspace directory (NOT from inside
# the opm-flowgeomechanics repository).  It will clone every required
# repository as a sibling, create a "superbuild" CMakeLists.txt in the
# workspace root, and compile everything under a "build" subdirectory.
#
# Usage:
#   mkdir ~/opm-workspace && cd ~/opm-workspace
#   bash /path/to/build.sh        # or copy it here first

set -e

# ── configuration ──────────────────────────────────────────────────────
FORK="${FORK:-hnil}"
BRANCH="${BRANCH:-geomech_hypre}"
GEOMECH_BRANCH="${GEOMECH_BRANCH:-geomech_hypre_cgal}"
BUILD_TYPE="${BUILD_TYPE:-Release}"
JOBS="${JOBS:-$(nproc)}"

# ── install system dependencies ────────────────────────────────────────
echo "=== Installing system dependencies ==="
sudo apt-get update
sudo apt-get install -y \
    build-essential cmake git pkg-config \
    libboost-all-dev \
    libdune-common-dev libdune-grid-dev libdune-istl-dev \
    libdune-localfunctions-dev libdune-geometry-dev \
    libsuitesparse-dev libtrilinos-zoltan-dev \
    libopenmpi-dev mpi-default-dev \
    libhdf5-mpi-dev libsuperlu-dev \
    libcgal-dev libgmp-dev \
    libblas-dev liblapack-dev \
    libfmt-dev

# ── clone OPM modules from fork ───────────────────────────────────────
echo "=== Cloning OPM modules ==="
for module in common grid simulators upscaling; do
    if [ ! -d "opm-${module}" ]; then
        git clone -b "${BRANCH}" \
            "https://github.com/${FORK}/opm-${module}.git"
    else
        echo "opm-${module} already present – skipping clone"
    fi
done

# ── clone opm-tests (regression-test data used by ctest) ──────────────
echo "=== Cloning opm-tests ==="
if [ ! -d "opm-tests" ]; then
    git clone https://github.com/OPM/opm-tests.git
else
    echo "opm-tests already present – skipping clone"
fi

# ── clone dune-foamgrid (required, not in Ubuntu packages) ─────────────
echo "=== Cloning dune-foamgrid ==="
if [ ! -d "dune-foamgrid" ]; then
    git clone https://gitlab.dune-project.org/extensions/dune-foamgrid.git
else
    echo "dune-foamgrid already present – skipping clone"
fi

# ── clone opm-flowgeomechanics ─────────────────────────────────────────
echo "=== Cloning opm-flowgeomechanics ==="
if [ ! -d "opm-flowgeomechanics" ]; then
    git clone -b "${GEOMECH_BRANCH}" \
        "https://github.com/${FORK}/opm-flowgeomechanics.git"
else
    echo "opm-flowgeomechanics already present – skipping clone"
fi

# ── generate superbuild CMakeLists.txt ─────────────────────────────────
echo "=== Creating superbuild CMakeLists.txt ==="
cat > CMakeLists.txt << 'SUPERBUILD'
# Superbuild: compile all OPM modules (and opm-flowgeomechanics) together.
project(opm)
cmake_minimum_required(VERSION 3.26)

foreach(TARGET opm-common opm-grid opm-simulators opm-upscaling opm-flowgeomechanics)
  set(${TARGET}_DIR ${CMAKE_BINARY_DIR}/${TARGET})
endforeach()

set(SIBLING_SEARCH 0)
# Necessary – sadly this means no IDE project generation
set_property(GLOBAL PROPERTY ALLOW_DUPLICATE_CUSTOM_TARGETS 1)

enable_testing()
add_subdirectory(opm-common)
add_subdirectory(opm-grid)
add_dependencies(opmgrid opmcommon)
add_subdirectory(opm-simulators)
add_dependencies(opmsimulators opmgrid)
add_subdirectory(opm-upscaling)
add_dependencies(opmupscaling opmgrid)
add_subdirectory(opm-flowgeomechanics)
add_dependencies(opmflowgeomechanics opmsimulators opmupscaling)
SUPERBUILD

# ── build ──────────────────────────────────────────────────────────────
echo "=== Configuring and building ==="
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE="${BUILD_TYPE}"
make -j"${JOBS}"
echo "=== Build finished successfully ==="

# ── run tests ─────────────────────────────────────────────────────────
echo "=== Running ctests ==="
ctest -j"${JOBS}" --output-on-failure
echo "=== Tests finished ==="
