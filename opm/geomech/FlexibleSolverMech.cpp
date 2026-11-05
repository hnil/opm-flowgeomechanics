/*
  Copyright 2019, 2020 SINTEF Digital, Mathematics and Cybernetics.

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

#include "config.h"
#include <opm/common/TimingMacros.hpp>
#include <opm/simulators/linalg/FlexibleSolver_impl.hpp>
#include <opm/simulators/linalg/PreconditionerFactory_impl.hpp>
#include <opm/geomech/MechLinePlanePreconditioner.hpp>
template <int N>
using OBMMec = Dune::BCRSMatrix<Dune::FieldMatrix<double, N, N>>;
template <int N>
using BVMec = Dune::BlockVector<Dune::FieldVector<double, N>>;
;
using CommMec = Dune::OwnerOverlapCopyCommunication<int, int>;
// Sequential operators.
template <int N>
using SeqOpMec = Dune::MatrixAdapter<OBMMec<N>, BVMec<N>, BVMec<N>>;

template <int N>
using ParOpMec = Dune::OverlappingSchwarzOperator<OBMMec<N>, BVMec<N>, BVMec<N>, CommMec>;

INSTANTIATE_FLEXIBLESOLVER_OP(SeqOpMec<1>);
INSTANTIATE_FLEXIBLESOLVER_OP(ParOpMec<1>);

namespace Opm
{
void addMechLinePlanePreconditioner(const std::vector<double>& coords)
{
    using SeqInfo = Dune::Amg::SequentialInformation;
    using SeqPrec = MechLinePlanePreconditioner<OBMMec<1>, BVMec<1>, SeqInfo>;
    PreconditionerFactory<SeqOpMec<1>, SeqInfo>::addCreator(
        "mech_lines",
        [&coords](const SeqOpMec<1>& op, const PropertyTree& prm, const std::function<BVMec<1>()>&, std::size_t) {
            static const SeqInfo seq;
            auto inner = PreconditionerFactory<SeqOpMec<1>, SeqInfo>::create(op, prm.get_child("inner"));
            return std::make_shared<SeqPrec>(op.getmat(), coords, prm, inner, seq);
        });
    using ParPrec = MechLinePlanePreconditioner<OBMMec<1>, BVMec<1>, CommMec>;
    PreconditionerFactory<ParOpMec<1>, CommMec>::addCreator(
        "mech_lines",
        [&coords](const ParOpMec<1>& op,
                  const PropertyTree& prm,
                  const std::function<BVMec<1>()>&,
                  std::size_t,
                  const CommMec& comm) {
            auto inner = PreconditionerFactory<ParOpMec<1>, CommMec>::create(op, prm.get_child("inner"), comm);
            return std::make_shared<ParPrec>(op.getmat(), coords, prm, inner, comm);
        });
}
} // namespace Opm
