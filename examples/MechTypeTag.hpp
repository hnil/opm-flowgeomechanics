#include <opm/simulators/flow/FlowProblem.hpp>
//#include <opm/simulators/flow/Main.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
//#include <opm/flowexperimental/blackoilintensivequantitiessimple.hh>
#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/geomech/eclgeomechmodel.hh>
#include <opm/geomech/eclproblemgeomech.hh>

namespace Opm
{
    namespace Properties
    {
        namespace TTag
        {
            struct EclFlowProblemMech {
                using InheritsFrom = std::tuple<FlowProblem, VtkGeoMech, FlowGeomechIstlSolverParams>;
            };
        }

        // Set the problem class
        template <class TypeTag>
        struct Problem<TypeTag, TTag::EclFlowProblemMech> {
            using type = EclProblemGeoMech<TypeTag>;
        };


        // template <class TypeTag>
        // struct WellModel<TypeTag, TTag::EclFlowProblemMech> {
        //     using type = BlackoilGeomechWellModel<TypeTag>;
        // };

        // template<class TypeTag>
        // struct Model<TypeTag, TTag::EclFlowProblemMech> {
        //     using type = BlackOilModelFvLocal<TypeTag>;
        // };


        // template<class TypeTag>
        // struct EclWellModel<TypeTag, TTag::EclFlowProblemMech> {
        //     using type = BlackoilWellModelFvExtra<TypeTag>;
        // };

        // template<class TypeTag>
        // struct NewtonMethod<TypeTag, TTag::EclFlowProblemMech> {
        //     using type = EclNewtonMethodLinesearch<TypeTag>;
        // };
        template <class TypeTag>
        struct EnableMech<TypeTag, TTag::EclFlowProblemMech> {
            static constexpr bool value = true;
        };

        template <class TypeTag>
        struct EnableEnergy<TypeTag, TTag::EclFlowProblemMech> {
            static constexpr bool value = true;
        };


        template <class TypeTag>
        struct VtkWriteMoleFractions<TypeTag, TTag::EclFlowProblemMech> {
            static constexpr bool value = false;
        };


        template <class TypeTag>
        struct EnableOpmRstFile<TypeTag, TTag::EclFlowProblemMech> {
            static constexpr bool value = true;
        };

        // the default for the allowed volumetric error for oil per second

        // template<class TypeTag>
        // struct IntensiveQuantities<TypeTag, TTag::EclFlowProblemMech> {
        //     //using type = EclBlackOilIntensiveQuantities<TypeTag>;
        //     using type = BlackOilIntensiveQuantitiesSimple<TypeTag>;
        //     //using type = BlackOilIntensiveQuantities<TypeTag>;
        //     //using type = BlackOilIntensiveQuantitiesDryGas<TypeTag>;
        // };

        // template<class TypeTag>
        // struct Linearizer<TypeTag, TTag::EclFlowProblemMech> { using type = TpfaLinearizer<TypeTag>; };

        // template<class TypeTag>
        // struct LocalResidual<TypeTag, TTag::EclFlowProblemMech> { using type = BlackOilLocalResidualTPFA<TypeTag>; };

        template <class TypeTag>
        struct EnableDiffusion<TypeTag, TTag::EclFlowProblemMech> {
            static constexpr bool value = false;
        };

        template <class TypeTag>
        struct EnableDisgasInWater<TypeTag, TTag::EclFlowProblemMech> {
            static constexpr bool value = false;
        };

        // static constexpr bool has_disgas_in_water = getPropValue<TypeTag, Properties::EnableDisgasInWater>();

        template <class TypeTag>
        struct Simulator<TypeTag, TTag::EclFlowProblemMech> {
            using type = Opm::Simulator<TypeTag>;
        };
        // simpler debugging

        template <class TypeTag>
        struct EnableAsyncEclOutput<TypeTag, TTag::EclFlowProblemMech> {
            static constexpr bool value = false;
        };

    }
}
namespace Opm {
namespace Parameters {
    template <class TypeTag>
        struct EnableVtkOutput<TypeTag, Properties::TTag::EclFlowProblemMech> {
            static constexpr bool value = true;
        };
    template <class TypeTag>
    struct ThreadsPerProcess<TypeTag, Properties::TTag::EclFlowProblemMech> {
        static constexpr int value = 1;
    };

    template <class TypeTag>
    struct NewtonTolerance<TypeTag, Properties::TTag::EclFlowProblemMech> {
        using type = GetPropType<TypeTag, Properties::Scalar>;
        static constexpr type value = 1e-2;
    };
    template <class TypeTag>
    struct EnableAsyncVtkOutput<TypeTag, Properties::TTag::EclFlowProblemMech> {
        static constexpr bool value = false;
    };
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
