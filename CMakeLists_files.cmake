# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# This file sets up five lists:
# MAIN_SOURCE_FILES     List of compilation units which will be included in
#                       the library. If it isn't on this list, it won't be
#                       part of the library. Please try to keep it sorted to
#                       maintain sanity.
#
# TEST_SOURCE_FILES     List of programs that will be run as unit tests.
#
# TEST_DATA_FILES       Files from the source three that should be made
#                       available in the corresponding location in the build
#                       tree in order to run tests there.
#
# EXAMPLE_SOURCE_FILES  Other programs that will be compiled as part of the
#                       build, but which is not part of the library nor is
#                       run as tests.
#
# PUBLIC_HEADER_FILES   List of public header files that should be
#                       distributed together with the library. The source
#                       files can of course include other files than these;
#                       you should only add to this list if the *user* of
#                       the library needs it.

list (APPEND MAIN_SOURCE_FILES
	opm/geomech/coupledsolver.cpp
	opm/geomech/CutDe.cpp
	opm/geomech/DiscreteDisplacement.cpp
	opm/geomech/FlexibleSolverMech.cpp
	opm/geomech/Fracture_fullSystemIteration.cpp
	opm/geomech/FractureMechanicsPreconditioner.cpp
	opm/geomech/Fracture.cpp
	opm/geomech/FractureModel.cpp
	opm/geomech/FractureWell.cpp
	opm/geomech/GeometryHelpers.cpp
	opm/geomech/GridStretcher.cpp
	opm/geomech/param_interior.cpp
	opm/geomech/RegularTrimesh.cpp
	opm/geomech/vem/vem.cpp
	opm/geomech/vem/vemutils.cpp
)

list (APPEND EXAMPLE_SOURCE_FILES
)

list (APPEND PROGRAM_SOURCE_FILES
)

list (APPEND PUBLIC_HEADER_FILES
	opm/geomech/BlackoilGeomechWellModel.hpp
	opm/geomech/BlackoilModelGeomech.hpp
	opm/geomech/boundaryutils.hh
	opm/geomech/convex_boundary.hpp
	opm/geomech/coupledsolver.hpp
	opm/geomech/CutDe.hpp
	opm/geomech/DiscreteDisplacement.hpp
	opm/geomech/DuneCommunicationHelpers.hpp
	opm/geomech/dune_utilities.hpp
	opm/geomech/eclgeomechmodel.hh
	opm/geomech/eclproblemgeomech.hh
	opm/geomech/elasticity_solver.hpp
	opm/geomech/elasticity_solver_impl.hpp
	opm/geomech/FlowGeomechLinearSolverParameters.hpp
	opm/geomech/Fracture.hpp
	opm/geomech/Fracture_impl.hpp
	opm/geomech/FractureModel.hpp
	opm/geomech/FractureModel_impl.hpp
	opm/geomech/FractureWell.hpp
	opm/geomech/GeometryHelpers.hpp
	opm/geomech/GridStretcher.hpp
	opm/geomech/Math.hpp
	opm/geomech/param_interior.hpp
	opm/geomech/RegularTrimesh.hpp
	opm/geomech/vem_elasticity_solver.hpp
	opm/geomech/vem_elasticity_solver_impl.hpp
	opm/geomech/vem/topology.hpp
	opm/geomech/vem/vem.hpp
	opm/geomech/vem/vemutils.hpp
	opm/geomech/vtkgeomechmodule.hh
)
