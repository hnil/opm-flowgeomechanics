#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
//#include <opm/simulators/flow/Main.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
//#include <opm/flowexperimental/blackoilintensivequantitiessimple.hh>
#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/geomech/GeoMechModel.hpp>
#include <opm/geomech/FlowProblemGeoMech.hpp>
#include <opm/geomech/BlackoilGeoMechWellModel.hpp>
#include <opm/geomech/BlackoilModelGeomech.hpp>
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
              using type = BlackoilModelGeomech<TypeTag>;
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

        // template<class TypeTag>
        // struct Linearizer<TypeTag, TTag::FlowProblemMech> { using type = TpfaLinearizer<TypeTag>; };

        // template<class TypeTag>
        // struct LocalResidual<TypeTag, TTag::FlowProblemMech> { using type = BlackOilLocalResidualTPFA<TypeTag>; };

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
