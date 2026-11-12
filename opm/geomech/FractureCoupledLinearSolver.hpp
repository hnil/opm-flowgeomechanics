/*
  Copyright 2026 Equinor ASA.

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
#pragma once

#include <functional>
#include <memory>
#include <string>

#include <dune/istl/operators.hh>
#include <dune/istl/solver.hh>

#include <opm/geomech/FractureMechanicsPreconditioner.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

namespace Opm
{
struct FractureSolveStats;

// Linear solve of one Newton step of the coupled fracture system
// [A I; C M] [dw; dp] = rhs: block scaling, Krylov solve and a rescue ladder.
// Keeps the preconditioner and Krylov solver between iterations of one grid solve.
class FractureCoupledLinearSolver
{
public:
    using SystemMatrix = FractureSystemMatrix;

    struct Options
    {
        std::string scaling{"none"};
        std::string solver;
        PropertyTree precond;
        bool updateMechOnReuse{false};
        double tol{0.0};
        double atol{0.0};
        int maxIter{0};
        int verbosity{0};
        int restart{100};
        bool gmresFallback{true};
        static Options read(const PropertyTree& prm);
    };

    // Writes the unscaled system for offline replay; returns the info file name or "".
    using DumpHook = std::function<std::string(const SystemMatrix& S,
                                               const char* tag,
                                               const Dune::InverseOperatorResult* result)>;

    struct Iteration
    {
        int index{0};
        int verbosity{0};
        int numClosed{0};
        bool newSystem{false}; // S was rebuilt: new operator, preconditioner and solver
        bool contactChanged{false};
        bool dump{false}; // call the hook with tag "manual" before solving
    };

    void setOptions(const Options& options) { opt_ = options; }

    // Scales S and rhs in place; solve() returns the unscaled update.
    void prepare(SystemMatrix& S, VectorHP& rhs);
    // Solves the prepared system; rhs is overwritten. Throws if no attempt converged.
    void solve(VectorHP& rhs, VectorHP& dx, const Iteration& it,
               FractureSolveStats& stats, const DumpHook& dump);
    bool scaled() const;

private:
    using Operator = Dune::MatrixAdapter<SystemMatrix, VectorHP, VectorHP>;
    using Solver = Dune::InverseOperator<VectorHP, VectorHP>;

    void scaleMatrix(int sign);
    std::unique_ptr<Solver> makeSolver(const std::string& type, FractureMechanicsPreconditioner& precond,
                                       double tol) const;

    Options opt_;
    SystemMatrix* S_{nullptr};
    double scaleMech_{1.0};
    double scalePressure_{1.0};
    std::unique_ptr<Operator> operator_;
    std::unique_ptr<FractureMechanicsPreconditioner> precond_;
    std::unique_ptr<Solver> solver_;
};

} // namespace Opm
