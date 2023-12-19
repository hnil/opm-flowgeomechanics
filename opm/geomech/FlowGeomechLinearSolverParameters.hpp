/*
  Copyright 2015, 2020 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2015 IRIS AS
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 NTNU
  Copyright 2015 Statoil AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_FLOW_GEOMECH_LINEARSOLVERPARAMETERS_HEADER_INCLUDED
#define OPM_FLOW_GEOMECH_LINEARSOLVERPARAMETERS_HEADER_INCLUDED

#include "opm/simulators/linalg/FlowLinearSolverParameters.hpp"
#include <opm/models/utils/basicproperties.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/models/utils/propertysystem.hh>


namespace Opm::Properties {

namespace TTag {
struct FlowGeomechIstlSolverParams {};
}

template<class TypeTag, class MyTypeTag>
struct LinearSolverMechReduction {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct LinearSolverMechMaxIter {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct LinearSolverMech {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct LinearSolverMechPrintJsonDefinition {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct LinearSolverMechVerbosity {
    using type = UndefinedProperty;
};
template<class TypeTag>
struct LinearSolverMechReduction<TypeTag, TTag::FlowGeomechIstlSolverParams> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-6;
};
template<class TypeTag>
struct LinearSolverMechMaxIter<TypeTag, TTag::FlowGeomechIstlSolverParams> {
    static constexpr int value = 200;
};
template<class TypeTag>
struct LinearSolverMechVerbosity<TypeTag, TTag::FlowGeomechIstlSolverParams> {
    static constexpr int value = 0;
};
template<class TypeTag>
struct LinearSolverMech<TypeTag, TTag::FlowGeomechIstlSolverParams> {
    static constexpr auto value = "ilu0";
};
template<class TypeTag>
struct LinearSolverMechPrintJsonDefinition<TypeTag, TTag::FlowGeomechIstlSolverParams> {
    static constexpr auto value = true;
};

} // namespace Opm::Properties

namespace Opm {

    /// This class carries all parameters for the NewtonIterationBlackoilInterleaved class.
    struct FlowLinearSolverParametersGeoMech : public FlowLinearSolverParameters
    {
        template <class TypeTag>
        void init()
        {
            // TODO: these parameters have undocumented non-trivial dependencies
            linear_solver_reduction_ = EWOMS_GET_PARAM(TypeTag, double, LinearSolverMechReduction);
            linear_solver_maxiter_ = EWOMS_GET_PARAM(TypeTag, int, LinearSolverMechMaxIter);
            linear_solver_verbosity_ = EWOMS_GET_PARAM(TypeTag, int, LinearSolverMechVerbosity);
            linsolver_ = EWOMS_GET_PARAM(TypeTag, std::string, LinearSolverMech);
            linear_solver_print_json_definition_ = EWOMS_GET_PARAM(TypeTag, bool, LinearSolverMechPrintJsonDefinition);
            linsolver_ = EWOMS_GET_PARAM(TypeTag, std::string, LinearSolverMech);
        }

        template <class TypeTag>
        static void registerParameters()
        {
            EWOMS_REGISTER_PARAM(TypeTag, double, LinearSolverMechReduction, "The minimum reduction of the residual which the linear solver must achieve");
            EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverMechMaxIter, "The maximum number of iterations of the linear solver");
            EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverMechVerbosity, "The verbosity level of the linear solver (0: off, 2: all)");
            EWOMS_REGISTER_PARAM(TypeTag, std::string, LinearSolverMech, "Configuration of solver. Valid options are: ilu0 (default), amg or umfpack. Alternatively, you can request a configuration to be read from a JSON file by giving the filename here, ending with '.json.'");
            EWOMS_REGISTER_PARAM(TypeTag, bool, LinearSolverMechPrintJsonDefinition, "Write the JSON definition of the linear solver setup to the DBG file.");
        }

        FlowLinearSolverParametersGeoMech() { reset(); }

        // set default values
        void reset()
        {
            linear_solver_reduction_  = 1e-6;
            linear_solver_maxiter_    = 200;
            linsolver_                = "ilu0";
        }
    };


} // namespace Opm

#endif // OPM_FLOW_GEOMECH_LINEARSOLVERPARAMETERS_HEADER_INCLUDED
