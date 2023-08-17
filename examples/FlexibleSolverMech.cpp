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
template <int N>
using OBMMec = Dune::BCRSMatrix<Dune::FieldMatrix<double, N, N>>;

// Sequential operators.
template <int N>
using SeqOpMec = Dune::MatrixAdapter<OBMMec<N>, BV<N>, BV<N>>;

template <int N>
using ParOpMec = Dune::OverlappingSchwarzOperator<OBMMec<N>, BV<N>, BV<N>, Comm>;

INSTANTIATE_FLEXIBLESOLVER_OP(SeqOpMec<1>);
INSTANTIATE_FLEXIBLESOLVER_OP(ParOpMec<1>);
