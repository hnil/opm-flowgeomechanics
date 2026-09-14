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
#ifndef OPM_WELL_LOCAL_SOLVER_HPP
#define OPM_WELL_LOCAL_SOLVER_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/models/discretization/common/tpfalinearizerstructs.hh>
#include <opm/simulators/flow/NewtonIterationContext.hpp>
#include <opm/simulators/linalg/extractMatrix.hpp>
#include <opm/simulators/wells/StandardWell.hpp>

#include <dune/istl/umfpack.hh>

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <set>
#include <string>
#include <vector>

namespace Opm {

/*!
 * \brief A local Newton solve of one well together with the cells it perforates.
 *
 * The domain is the well's fracture (auxiliary) cells and, optionally, its
 * perforated grid cells plus a number of rings of grid neighbours.  Everything
 * outside the domain is frozen at the current iterate: it enters only through
 * the off-diagonal columns the domain linearization already writes.  The
 * well's own equations enter through their Schur complement, which is why the
 * global run must carry the wells in the matrix
 * (--matrix-add-well-contributions=true).
 *
 * Built only from public entry points of the flow model: the linearizer's
 * domain-restricted linearization, the per-well assemble / contribute /
 * recover methods, the blackoil Newton update restricted to a DOF list, and
 * the sub-system extraction helpers.  Nothing in opm-simulators is modified.
 */
template <class TypeTag>
class WellLocalSolver
{
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using Matrix = typename SparseMatrixAdapter::IstlMatrix;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using StdWell = StandardWell<TypeTag>;

public:
    struct Settings {
        int ring = -1;          //!< -1: fracture cells only; >= 0: perforated cells + rings
        int maxIter = 20;
        double reduction = 1e-4; //!< stop when the local residual fell by this factor ...
        double absTol = 1e-9;    //!< ... or is below this (per-equation max, SI rates)
        int verbosity = 0;
    };

    struct Domain {
        std::string well;
        std::vector<int> cells; // sorted global DOF indices
        std::size_t numAux = 0;
        std::size_t numGrid = 0;
    };

    struct Report {
        bool converged = false;
        int iterations = 0;
        double residual0 = 0.0;
        double residual = 0.0;
    };

    explicit WellLocalSolver(Simulator& simulator)
        : simulator_(simulator)
    {}

    //! The domain of one well: its auxiliary perforations and, with ring >= 0,
    //! its grid perforations plus that many rings of grid neighbours.
    template <class AuxModule>
    Domain buildDomain(const std::string& wellName,
                       const AuxModule* aux,
                       const Settings& settings) const
    {
        Domain d;
        d.well = wellName;
        std::set<int> cells;

        if (aux != nullptr) {
            // every cell of the well's fractures, not only the perforated ones
            const auto fractureCells = aux->cellsOfWell(wellName);
            cells.insert(fractureCells.begin(), fractureCells.end());
            for (const auto& p : aux->wellPerforations(wellName)) {
                cells.insert(p.cell);
            }
        }
        d.numAux = cells.size();

        if (settings.ring >= 0) {
            const auto& model = simulator_.model();
            const auto numGridDof = static_cast<int>(model.numGridDof());
            const auto& well = simulator_.problem().wellModel().getWell(wellName);
            std::set<int> front;
            for (const int c : well.cells()) {
                if (c < numGridDof) {
                    front.insert(c);
                }
            }
            const auto& nbInfo = model.linearizer().getNeighborInfo();
            for (int r = 0; r < settings.ring; ++r) {
                std::set<int> next;
                for (const int c : front) {
                    for (const auto& nb : nbInfo[c]) {
                        const int n = static_cast<int>(nb.neighbor);
                        if (n < numGridDof && !front.count(n)) {
                            next.insert(n);
                        }
                    }
                }
                front.insert(next.begin(), next.end());
            }
            d.numGrid = front.size();
            cells.insert(front.begin(), front.end());
        }

        d.cells.assign(cells.begin(), cells.end());
        return d;
    }

    //! Solve the domain with the reservoir outside it frozen.
    Report solve(const Domain& domain, const double dt, const Settings& settings)
    {
        Report rep;
        if (domain.cells.empty()) {
            return rep;
        }
        auto& problem = simulator_.problem();
        auto& model = simulator_.model();
        auto& wellModel = problem.wellModel();
        auto& linearizer = model.linearizer();

        if (!wellModel.addMatrixContributions()) {
            OPM_THROW(std::runtime_error,
                      "WellLocalSolver needs the wells in the matrix "
                      "(--matrix-add-well-contributions=true)");
        }

        // the well container holds the live well objects; find ours
        WellInterface<TypeTag>* wellPtr = nullptr;
        for (const auto& w : wellModel.wellContainer()) {
            if (w->name() == domain.well) {
                wellPtr = w.get();
            }
        }
        if (wellPtr == nullptr) {
            return rep;
        }
        auto& well = *wellPtr;
        auto* stdWell = dynamic_cast<StdWell*>(wellPtr);

        // local iteration context: the domain linearization then resets only
        // the domain's rows, and the well skips its own inner iterations (the
        // local Newton iterates it)
        LocalContextGuard<std::remove_reference_t<decltype(problem)>> guard(problem);
        // the well code logs through the group-state helper's deferred logger,
        // which the well model only pushes around its own assembly
        auto loggerGuard = wellModel.groupStateHelper().pushLogger(/*do_mpi_gather=*/false);

        FullDomain<std::vector<int>> dom {domain.cells};
        auto& solution = model.solution(/*timeIdx=*/0);
        GlobalEqVector dxGlobal(model.numTotalDof());

        for (int it = 0; it < settings.maxIter; ++it) {
            // assemble the well at the current iterate, then the domain rows
            well.assembleWellEq(simulator_, dt, wellModel.groupStateHelper(), wellModel.wellState());
            wellModel.updateCellRates();
            linearizer.linearizeDomain(dom);

            // the well's Schur complement into the matrix and the residual, as
            // the global linearization does for every well
            auto& jacobian = linearizer.jacobian();
            auto& residual = linearizer.residual();
            well.addWellContributions(jacobian);
            {
                const auto& wcells = well.cells();
                GlobalEqVector rloc(wcells.size());
                for (std::size_t i = 0; i < wcells.size(); ++i) {
                    rloc[i] = residual[wcells[i]];
                }
                well.apply(rloc);
                for (std::size_t i = 0; i < wcells.size(); ++i) {
                    residual[wcells[i]] = rloc[i];
                }
            }

            const auto rDom = Details::extractVector(residual, domain.cells);
            double rmax = 0.0;
            for (const auto& blk : rDom) {
                for (const auto v : blk) {
                    rmax = std::max(rmax, std::abs(v));
                }
            }
            double wmax = 0.0;
            if (stdWell != nullptr) {
                for (const auto& blk : stdWell->linSys().residual()) {
                    for (const auto v : blk) {
                        wmax = std::max(wmax, std::abs(v));
                    }
                }
            }
            const double rnorm = std::max(rmax, wmax);
            if (it == 0) {
                rep.residual0 = rnorm;
            }
            rep.residual = rnorm;
            rep.iterations = it;
            if (settings.verbosity > 0) {
                OpmLog::info(fmt::format("WellLocalSolver {} it {} res {:.3e} (res {:.3e} well {:.3e})",
                                         domain.well, it, rnorm, rmax, wmax));
            }
            if (rnorm < settings.absTol || rnorm < settings.reduction * rep.residual0) {
                rep.converged = true;
                break;
            }

            // direct solve of the extracted block
            auto jDom = Details::extractMatrix(jacobian.istlMatrix(), domain.cells);
            GlobalEqVector xDom(domain.cells.size());
            xDom = 0.0;
            {
                Dune::UMFPack<Matrix> lu(jDom, 0, false);
                Dune::InverseOperatorResult res;
                auto rhs = rDom;
                lu.apply(xDom, rhs, res);
            }
            for (const auto& blk : xDom) {
                for (const auto v : blk) {
                    if (!std::isfinite(v)) {
                        OPM_THROW(NumericalProblem, "WellLocalSolver: non-finite update for well " + domain.well);
                    }
                }
            }

            // update the domain's primary variables (with the blackoil
            // limiters and variable switching), their cached intensive
            // quantities, and the well's own unknowns
            dxGlobal = 0.0;
            Details::setGlobal(xDom, domain.cells, dxGlobal);
            model.newtonMethod().update_(solution, solution, dxGlobal, dxGlobal, domain.cells);
            for (const int c : domain.cells) {
                IntensiveQuantities iq;
                iq.update(problem, solution[c], static_cast<unsigned>(c), /*timeIdx=*/0);
                model.updateCachedIntensiveQuantities(iq, static_cast<unsigned>(c), /*timeIdx=*/0);
            }
            well.recoverWellSolutionAndUpdateWellState(simulator_, dxGlobal,
                                                       wellModel.groupStateHelper(),
                                                       wellModel.wellState());
        }

        // leave the well's cell rates consistent with its final state
        well.assembleWellEq(simulator_, dt, wellModel.groupStateHelper(), wellModel.wellState());
        wellModel.updateCellRates();
        return rep;
    }

private:
    Simulator& simulator_;
};

} // namespace Opm

#endif // OPM_WELL_LOCAL_SOLVER_HPP
