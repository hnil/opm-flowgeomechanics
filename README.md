# opm-flowgeomechanics

Simulator Extension to OPM Flow for Handling Geomechanics and Formation Damage

## Prerequisites

The build has been tested on **Ubuntu 24.04 LTS**. Install the required
system packages first:

```bash
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
```

`dune-foamgrid` is also required but is not available as an Ubuntu
package — the build script below clones it automatically.

For more background on OPM dependencies see the
[OPM project website](https://opm-project.org).

## Quick build (automated)

A helper script is provided that clones every dependency, generates a
*superbuild* `CMakeLists.txt`, and compiles everything:

```bash
mkdir ~/opm-workspace && cd ~/opm-workspace
bash /path/to/opm-flowgeomechanics/build.sh
```

The script accepts several environment variables:

| Variable          | Default               | Description                               |
|-------------------|-----------------------|-------------------------------------------|
| `FORK`            | `hnil`                | GitHub user/org that hosts the OPM forks  |
| `BRANCH`          | `geomech_hypre`       | Branch for OPM modules                    |
| `GEOMECH_BRANCH`  | `geomech_hypre_cgal`  | Branch for opm-flowgeomechanics           |
| `BUILD_TYPE`      | `Release`             | CMake build type                          |
| `JOBS`            | `$(nproc)`            | Parallel make jobs                        |

## Manual build (step by step)

If you prefer to run each step yourself:

```bash
# 1. Create a workspace
mkdir ~/opm-workspace && cd ~/opm-workspace

# 2. Clone the required OPM modules
for module in common grid simulators upscaling; do
    git clone -b geomech_hypre https://github.com/hnil/opm-${module}.git
done

# 3. Clone dune-foamgrid
git clone https://gitlab.dune-project.org/extensions/dune-foamgrid.git

# 4. Clone opm-flowgeomechanics
git clone -b geomech_hypre_cgal https://github.com/hnil/opm-flowgeomechanics.git

# 5. Create a superbuild CMakeLists.txt in the workspace root
#    (see build.sh for the full content, or copy it from
#     opm-flowgeomechanics/build.sh)
cat > CMakeLists.txt << 'EOF'
project(opm)
cmake_minimum_required(VERSION 3.26)
foreach(TARGET opm-common opm-grid opm-simulators opm-upscaling opm-flowgeomechanics)
  set(${TARGET}_DIR ${CMAKE_BINARY_DIR}/${TARGET})
endforeach()
set(SIBLING_SEARCH 0)
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
EOF

# 6. Configure and build
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

## License

Distributed under the GNU General Public License v3 or later (GPLv3+).
