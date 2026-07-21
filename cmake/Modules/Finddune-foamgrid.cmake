# - Find DUNE foamgrid library
#
# Sets HAVE_DUNE_FOAMGRID and makes the version number available in
# config.h.  Follows the pattern of the Finddune-* modules shipped with
# opm-common.
if(dune-foamgrid_FOUND)
  return()
endif()

if(dune-foamgrid_FIND_REQUIRED)
  find_package(dune-foamgrid CONFIG REQUIRED)
else()
  find_package(dune-foamgrid CONFIG)
endif()

if(dune-foamgrid_FOUND)
  target_compile_definitions(dunefoamgrid INTERFACE HAVE_DUNE_FOAMGRID=1)
  # make version number available in config.h
  include (UseDuneVer)
  find_dune_version ("dune" "foamgrid")
endif()
