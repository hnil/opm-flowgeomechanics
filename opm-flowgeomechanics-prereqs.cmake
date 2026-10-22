# defines that must be present in config.h for our headers
set (opm-flowgeomechanics_CONFIG_VAR
  HAVE_OPM_GRID
  HAVE_PTHREAD
  HAVE_MPI
  HAVE_PETSC
  HAVE_SUITESPARSE_UMFPACK_H
  HAVE_DUNE_ISTL
  DUNE_ISTL_VERSION_MAJOR
  DUNE_ISTL_VERSION_MINOR
  DUNE_ISTL_VERSION_REVISION
  HAVE_SUITESPARSE_UMFPACK
  HAVE_HDF5
  USE_TRACY
  HAVE_DUNE_PDELAB
  )

# dependencies; transitive requirements of the upstream OPM modules are
# resolved by their package configuration files.
find_package(Boost COMPONENTS date_time unit_test_framework REQUIRED)
find_package(dune-common REQUIRED)
find_package(dune-istl REQUIRED)
find_package(dune-foamgrid REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(CGAL REQUIRED)
find_package(GMP REQUIRED)
find_package(SuiteSparse COMPONENTS UMFPACK)
find_package(opm-grid REQUIRED)
find_package(opm-upscaling REQUIRED)
find_package(opm-simulators REQUIRED)

find_package(MPI)
find_package(HDF5)
find_package(fmt)
find_package(PETSc)
find_package(SuperLU)
find_package(dune-alugrid)
find_package(dune-polygongrid)
