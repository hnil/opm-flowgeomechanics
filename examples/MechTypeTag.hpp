/*
  Copyright 2025 Equinor ASA.

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

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
//#include <opm/simulators/flow/Main.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
//#include <opm/flowexperimental/blackoilintensivequantitiessimple.hh>
#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/models/blackoil/blackoilconvectivemixingmodule.hh>
#include <opm/models/blackoil/blackoilenergymodules.hh>
#include <opm/geomech/GeoMechModel.hpp>
#include <opm/geomech/FlowProblemGeoMech.hpp>
#include <opm/geomech/BlackoilGeoMechWellModel.hpp>
#include <opm/geomech/NonlinearSystemBlackOilReservoirGeoMech.hpp>
namespace Opm
{
    namespace Properties
    {
        namespace TTag
        {
            struct FlowProblemMech {
	      using InheritsFrom = std::tuple<FlowProblem> ;//?? should it be blackoil
            };
        }

        // Set the problem class
        template <class TypeTag>
        struct Problem<TypeTag, TTag::FlowProblemMech> {
            using type = FlowProblemGeoMech<TypeTag>;
        };

        template <class TypeTag>
        struct NonlinearSystem<TypeTag, TTag::FlowProblemMech> {
              using type = NonlinearSystemBlackOilReservoirGeoMech<TypeTag>;
        };


        // template <class TypeTag>
        // struct WellModel<TypeTag, TTag::FlowProblemMech> {
        //     using type = BlackoilGeoMechWellModel<TypeTag>;
        // };

        // template<class TypeTag>
        // struct Model<TypeTag, TTag::FlowProblemMech> {
        //     using type = BlackOilModelFvLocal<TypeTag>;
        // };


        // template<class TypeTag>
        // struct EclWellModel<TypeTag, TTag::FlowProblemMech> {
        //     using type = BlackoilWellModelFvExtra<TypeTag>;
        // };

        // template<class TypeTag>
        // struct NewtonMethod<TypeTag, TTag::FlowProblemMech> {
        //     using type = EclNewtonMethodLinesearch<TypeTag>;
        // };
        template <class TypeTag>
        struct EnableMech<TypeTag, TTag::FlowProblemMech> {
            static constexpr bool value = true;
        };

        template <class TypeTag>
        struct EnableEnergy<TypeTag, TTag::FlowProblemMech> {
            static constexpr bool value = true;
        };

        template<class TypeTag>
        struct EnergyModuleType<TypeTag, TTag::FlowProblemMech>
        { static constexpr EnergyModules value = EnergyModules::FullyImplicitThermal; };

        // template <class TypeTag>
        // struct VtkWriteMoleFractions<TypeTag, TTag::FlowProblemMech> {
        //     static constexpr bool value = false;
        // };


        // template <class TypeTag>
        // struct EnableOpmRstFile<TypeTag, TTag::FlowProblemMech> {
        //     static constexpr bool value = true;
        // };

        // the default for the allowed volumetric error for oil per second

        // template<class TypeTag>
        // struct IntensiveQuantities<TypeTag, TTag::FlowProblemMech> {
        //     //using type = EclBlackOilIntensiveQuantities<TypeTag>;
        //     using type = BlackOilIntensiveQuantitiesSimple<TypeTag>;
        //     //using type = BlackOilIntensiveQuantities<TypeTag>;
        //     //using type = BlackOilIntensiveQuantitiesDryGas<TypeTag>;
        // };

        // Assemble per degree of freedom rather than per grid element, as every other
        // black-oil target does.  Beyond being the mainstream path, it is what lets a
        // degree of freedom that has no grid element -- a numerical aquifer or a fracture
        // represented as auxiliary cells -- carry the flow equations at all: the
        // element-driven linearizer walks straight past such a degree of freedom and
        // leaves its row empty.
        template<class TypeTag>
        struct Linearizer<TypeTag, TTag::FlowProblemMech> { using type = TpfaLinearizer<TypeTag>; };

        template<class TypeTag>
        struct LocalResidual<TypeTag, TTag::FlowProblemMech> { using type = BlackOilLocalResidualTPFA<TypeTag>; };

        template <class TypeTag>
        struct EnableDiffusion<TypeTag, TTag::FlowProblemMech> {
            static constexpr bool value = false;
        };

        template <class TypeTag>
        struct EnableDisgasInWater<TypeTag, TTag::FlowProblemMech> {
            static constexpr bool value = false;
        };


        template <class TypeTag>
        struct WellModel<TypeTag, TTag::FlowProblemMech> {
            using type = BlackoilGeoMechWellModel<TypeTag>;
          //using type = BlackoilWellModel<TypeTag>;
        };

        // static constexpr bool has_disgas_in_water = getPropValue<TypeTag, Properties::EnableDisgasInWater>();

        template <class TypeTag>
        struct Simulator<TypeTag, TTag::FlowProblemMech> {
            using type = Opm::Simulator<TypeTag>;
        };
        // simpler debugging

        // template <class TypeTag>
        // struct EnableAsyncEclOutput<TypeTag, TTag::FlowProblemMech> {
        //     static constexpr bool value = false;
        // };

    }
}
namespace Opm {
namespace Parameters {
    // template <class TypeTag>
    //     struct EnableVtkOutput<TypeTag, Properties::TTag::FlowProblemMech> {
    //         static constexpr bool value = true;
    //     };
    // template <class TypeTag>
    // struct ThreadsPerProcess<TypeTag, Properties::TTag::FlowProblemMech> {
    //     static constexpr int value = 1;
    // };

    // template <class TypeTag>
    // struct NewtonTolerance<TypeTag, Properties::TTag::FlowProblemMech> {
    //     using type = GetPropType<TypeTag, Properties::Scalar>;
    //     static constexpr type value = 1e-2;
    // };
    // template <class TypeTag>
    // struct EnableAsyncVtkOutput<TypeTag, Properties::TTag::FlowProblemMech> {
    //     static constexpr bool value = false;
    // };
}
}

// template<>
// class SupportsFaceTag<Dune::PolyhedralGrid<3, 3>>
//     : public std::bool_constant<true>
// {};
// template<>
// class SupportsFaceTag<Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>>
//     : public std::bool_constant<true>
// {};
