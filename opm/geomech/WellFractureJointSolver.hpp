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
#ifndef OPM_WELL_FRACTURE_JOINT_SOLVER_HPP
#define OPM_WELL_FRACTURE_JOINT_SOLVER_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/geomech/FractureMechCoupling.hpp>
#include <opm/geomech/WellLocalSolver.hpp>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/umfpack.hh>

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace Opm {

/*!
 * \brief One Newton solve of a well, the pressures of all its fracture cells and
 *        the widths of all its fractures, with the reservoir frozen.
 *
 * Contact is Fischer-Burmeister, the same row as Fracture::fullSystemIteration.
 * Geometry is fixed. The well's connection factors are held fixed within the solve.
 */
template <class TypeTag>
class WellFractureJointSolver
{
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using Matrix = typename SparseMatrixAdapter::IstlMatrix;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using StdWell = StandardWell<TypeTag>;
    using MsWell = MultisegmentWell<TypeTag>;
    using LocalSolver = WellLocalSolver<TypeTag>;
    static constexpr int numEq = getPropValue<TypeTag, Properties::NumEq>();
    static constexpr bool enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>();
    using ScalarMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
    using ScalarVector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
    using WidthVector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
    using ADVector = Dune::FieldVector<Evaluation, numEq>;

public:
    using Domain = typename LocalSolver::Domain;
    using LocalSettings = typename LocalSolver::Settings;

    struct Settings {
        LocalSettings local {};
        double toleranceFb = 1e-6;   //!< scaled complementarity residual
        double toleranceAux = 1e-3;  //!< fracture cell mass residual relative to the well rate
        double fbMu = 0.0;           //!< Chen-Mangasarian smoothing
        double maxDwidth = 1e-2;     //!< largest width change per iteration [m]
        double relativeDwidth = 1.0; //!< largest width change relative to the width
        int fdColumns = 0;           //!< > 0: FD-check this many width columns once
    };

    struct Report {
        bool converged = false;
        int iterations = 0;
        double residual0 = 0.0;
        double residual = 0.0;
        double bhp = 0.0;
    };

    explicit WellFractureJointSolver(Simulator& simulator)
        : simulator_(simulator)
    {
    }

    template <class FractureModelT, class AuxModule>
    Report solve(const Domain& domain, FractureModelT& fractures, const AuxModule& aux,
                 const double dt, const Settings& settings)
    {
        Report rep;
        auto& problem = simulator_.problem();
        auto& model = simulator_.model();
        auto& wellModel = problem.wellModel();
        auto& linearizer = model.linearizer();
        if (wellModel.addMatrixContributions()) {
            OPM_THROW(std::runtime_error,
                      "solver.well_joint_solve needs --matrix-add-well-contributions=false");
        }

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
        auto* msWell = dynamic_cast<MsWell*>(wellPtr);
        if (stdWell == nullptr && msWell == nullptr) {
            OPM_THROW(std::runtime_error,
                      "solver.well_joint_solve: unsupported well type for " + domain.well);
        }
        auto withLinSys = [&](auto&& f) {
            if (stdWell != nullptr) {
                f(stdWell->linSys());
            } else {
                f(msWell->linSys());
            }
        };

        std::unordered_map<int, int> localOf;
        for (std::size_t i = 0; i < domain.cells.size(); ++i) {
            localOf[domain.cells[i]] = static_cast<int>(i);
        }

        const std::size_t nR = domain.cells.size() * numEq;
        std::vector<FractureBlock> blocks;
        std::size_t numWidths = 0;
        std::size_t fidx = 0;
        for (auto& wellFractures : fractures.wellFractures()) {
            for (auto& fracture : wellFractures) {
                if (fracture.wellInfo().name == domain.well) {
                    FractureBlock b;
                    b.fracture = &fracture;
                    b.index = fidx;
                    b.offset = numWidths;
                    b.nc = fracture.numCells();
                    b.dof = aux.cellDofs(fidx);
                    b.traction = fracture.jointNormalTraction();
                    numWidths += b.nc;
                    blocks.push_back(std::move(b));
                }
                ++fidx;
            }
        }
        if (blocks.empty()) {
            return rep;
        }

        std::optional<LocalContextGuard<std::remove_reference_t<decltype(problem)>>> guard;
        if (!settings.local.innerWellIterations) {
            guard.emplace(problem);
        }
        auto loggerGuard = wellModel.groupStateHelper().pushLogger(/*do_mpi_gather=*/false);
        FullDomain<std::vector<int>> dom {domain.cells};
        auto& solution = model.solution(/*timeIdx=*/0);
        GlobalEqVector dxGlobal(model.numTotalDof());

        auto refreshCache = [&]() {
            for (const int c : domain.cells) {
                IntensiveQuantities iq;
                iq.update(problem, solution[c], static_cast<unsigned>(c), /*timeIdx=*/0);
                model.updateCachedIntensiveQuantities(iq, static_cast<unsigned>(c), /*timeIdx=*/0);
            }
        };
        auto setWidths = [&](const std::vector<WidthVector>& widths) {
            for (std::size_t k = 0; k < blocks.size(); ++k) {
                blocks[k].fracture->setJointWidth(widths[k]);
            }
            problem.bindFractureAuxCells(/*allowTopologyChange=*/false, /*requireStableLayout=*/true);
            // the assembly reads its own copy of the transmissibilities
            linearizer.updateDiscretizationParameters();
        };
        auto currentWidths = [&]() {
            std::vector<WidthVector> w;
            for (const auto& b : blocks) {
                w.push_back(b.fracture->fractureWidth());
            }
            return w;
        };

        std::vector<typename std::remove_reference_t<decltype(solution)>::block_type> solution0;
        for (const int c : domain.cells) {
            solution0.push_back(solution[c]);
        }
        const auto wellState0 = wellModel.wellState();
        const auto widths0 = currentWidths();
        auto rollback = [&]() {
            for (std::size_t i = 0; i < domain.cells.size(); ++i) {
                solution[domain.cells[i]] = solution0[i];
            }
            setWidths(widths0);
            refreshCache();
            wellModel.wellState() = wellState0;
            well.updatePrimaryVariables(wellModel.groupStateHelper());
        };

        setWidths(widths0);

        double rmax = 0.0, wmax = 0.0, amax = 0.0, fbmax = 0.0, merit = 0.0;
        int worstAux = -1, worstAuxEq = -1;
        auto assembleAndMeasure = [&]() {
            well.assembleWellEq(simulator_, dt, wellModel.groupStateHelper(), wellModel.wellState());
            wellModel.updateCellRates();
            linearizer.linearizeDomain(dom);
            auto& residual = linearizer.residual();
            const auto& wcells = well.cells();
            GlobalEqVector rloc(wcells.size());
            for (std::size_t i = 0; i < wcells.size(); ++i) {
                rloc[i] = residual[wcells[i]];
            }
            well.apply(rloc);
            for (std::size_t i = 0; i < wcells.size(); ++i) {
                residual[wcells[i]] = rloc[i];
            }

            const auto& ws = wellModel.wellState().well(well.indexOfWell());
            double scale = 1.0e-6;
            for (const auto r : ws.surface_rates) {
                scale = std::max(scale, std::abs(r));
            }
            rmax = 0.0;
            amax = 0.0;
            worstAux = -1;
            worstAuxEq = -1;
            const auto numGridDof = model.numGridDof();
            for (const int c : domain.cells) {
                if (static_cast<unsigned>(c) >= numGridDof) {
                    // no CNV for a fracture cell: its volume is an aperture times an area
                    for (unsigned e = 0; e < residual[c].size(); ++e) {
                        const bool isEnergy = enableEnergy
                            && (static_cast<int>(e) == Indices::contiEnergyEqIdx);
                        if (!isEnergy) {
                            if (std::abs(residual[c][e]) / scale > amax) {
                                worstAux = c;
                                worstAuxEq = static_cast<int>(e);
                            }
                            amax = std::max(amax, std::abs(residual[c][e]) / scale);
                        }
                    }
                    continue;
                }
                const auto pv = problem.referencePorosity(c, /*timeIdx=*/0) * model.dofTotalVolume(c);
                if (!(pv > 0.0)) {
                    continue;
                }
                for (unsigned e = 0; e < residual[c].size(); ++e) {
                    const bool isEnergy = enableEnergy
                        && (static_cast<int>(e) == Indices::contiEnergyEqIdx);
                    const auto tol = isEnergy ? settings.local.toleranceCnvEnergy
                                              : settings.local.toleranceCnv;
                    rmax = std::max(rmax, double(std::abs(residual[c][e]) * dt / pv / tol));
                }
            }
            wmax = 0.0;
            double wmaxMass = 0.0;
            withLinSys([&](const auto& linSys) {
                for (const auto& blk : linSys.residual()) {
                    for (std::size_t e = 0; e < blk.size(); ++e) {
                        const double v = std::abs(double(blk[e])) / scale;
                        wmax = std::max(wmax, v);
                        // the last row is in pressure units, incomparable with the rates in a merit
                        if (e + 1 < blk.size()) {
                            wmaxMass = std::max(wmaxMass, v);
                        }
                    }
                }
            });
            fbmax = 0.0;
            // the scaled phi saturates at 2 on any shut cell under tension: useless for a step search
            double phiNorm2 = 0.0;
            for (auto& b : blocks) {
                this->evaluateContact_(b, model, settings.fbMu);
                for (std::size_t i = 0; i < b.nc; ++i) {
                    fbmax = std::max(fbmax, std::abs(b.phi[i]) / b.scale[i]);
                    phiNorm2 += std::pow(b.phi[i] / 1e6, 2);
                }
            }
            merit = std::sqrt(std::pow(rmax * settings.local.toleranceCnv, 2)
                              + wmaxMass * wmaxMass + amax * amax + phiNorm2);
            return std::max({rmax, wmax / settings.local.toleranceWell,
                             amax / settings.toleranceAux, fbmax / settings.toleranceFb});
        };

        double rnorm = assembleAndMeasure();
        rep.residual0 = rnorm;
        if (!std::isfinite(rnorm)) {
            rollback();
            return rep;
        }
        if (settings.fdColumns > 0) {
            this->checkWidthColumns_(blocks, localOf, dt, settings.fdColumns, setWidths, assembleAndMeasure);
            rnorm = assembleAndMeasure();
        }

        for (int it = 0; it < settings.local.maxIter; ++it) {
            rep.residual = rnorm;
            rep.iterations = it;
            if (settings.local.verbosity > 0) {
                OpmLog::info(fmt::format("WellFractureJoint {} it {} res {:.3e} "
                                         "(cnv {:.3e} well {:.3e} aux {:.3e} fb {:.3e})",
                                         domain.well, it, rnorm, rmax, wmax, amax, fbmax));
                if (settings.local.verbosity > 1 && worstAux >= 0) {
                    const auto& fs = model.intensiveQuantities(worstAux, 0).fluidState();
                    const auto oilIdx = FluidSystem::oilPhaseIdx;
                    const auto waterIdx = FluidSystem::waterPhaseIdx;
                    OpmLog::info(fmt::format("WellFractureJoint worst aux dof {} eq {}: res {:.3e} pw {:.2f} bar "
                                             "Sw {:.4f} So {:.4f} vol {:.3e} dt {:.3e}", worstAux, worstAuxEq,
                                             double(linearizer.residual()[worstAux][worstAuxEq]),
                                             getValue(fs.pressure(waterIdx)) / 1e5,
                                             getValue(fs.saturation(waterIdx)), getValue(fs.saturation(oilIdx)),
                                             double(model.dofTotalVolume(worstAux)), dt));
                    for (const auto& nb : linearizer.getNeighborInfo()[worstAux]) {
                        const auto& fsN = model.intensiveQuantities(nb.neighbor, 0).fluidState();
                        ADVector flux(0.0), darcy(0.0);
                        LocalResidual::computeFlux(flux, darcy, worstAux, nb.neighbor,
                                                   model.intensiveQuantities(worstAux, 0),
                                                   model.intensiveQuantities(nb.neighbor, 0),
                                                   nb.res_nbinfo, simulator_.problem().moduleParams());
                        OpmLog::info(fmt::format("WellFractureJoint   nbr {} trans {:.3e} dp {:.4e} bar T {:.2f}/{:.2f} "
                                                 "flux {:.3e}", nb.neighbor, double(nb.res_nbinfo.trans),
                                                 getValue(fs.pressure(waterIdx) - fsN.pressure(waterIdx)) / 1e5,
                                                 getValue(fs.temperature(waterIdx)), getValue(fsN.temperature(waterIdx)),
                                                 getValue(flux[worstAuxEq]) * nb.res_nbinfo.faceArea));
                    }
                }
            }
            if (rnorm < 1.0) {
                rep.converged = true;
                break;
            }

            std::vector<std::tuple<std::size_t, std::size_t, double>> entries;
            const auto rDom = Details::extractVector(linearizer.residual(), domain.cells);
            auto jDom = Details::extractMatrix(linearizer.jacobian().istlMatrix(), domain.cells);
            this->addFlowAndWell_(entries, jDom, domain.cells, well.cells(), localOf, withLinSys, nR);
            int nWellEq = 0;
            withLinSys([&](const auto& linSys) {
                nWellEq = static_cast<int>(linSys.getD().N() * linSys.getD()[0][0].N());
            });
            const std::size_t nFlow = nR + static_cast<std::size_t>(nWellEq);
            for (const auto& b : blocks) {
                this->addWidthBlocks_(entries, b, localOf, dt, nFlow);
            }

            const std::size_t n = nFlow + numWidths;
            ScalarVector rhs(n), x(n);
            rhs = 0.0;
            x = 0.0;
            for (std::size_t i = 0; i < rDom.size(); ++i) {
                for (int a = 0; a < numEq; ++a) {
                    rhs[i * numEq + a] = rDom[i][a];
                }
            }
            for (const auto& b : blocks) {
                for (std::size_t i = 0; i < b.nc; ++i) {
                    rhs[nFlow + b.offset + i] = b.phi[i];
                }
            }
            this->pinSaturationsOfWaterCells_(entries, rhs, domain.cells);
            this->solveSparse_(entries, n, rhs, x);
            bool finite = true;
            for (std::size_t i = 0; i < n; ++i) {
                finite = finite && std::isfinite(x[i][0]);
            }
            if (!finite) {
                break;
            }

            GlobalEqVector xDom(domain.cells.size());
            for (std::size_t i = 0; i < domain.cells.size(); ++i) {
                for (int a = 0; a < numEq; ++a) {
                    xDom[i][a] = x[i * numEq + a];
                }
            }
            if (settings.local.verbosity > 1) {
                double dpAux = 0.0, dpGrid = 0.0, dwMax = 0.0, wMax = 0.0, dwRel = 0.0;
                const auto numGridDof = model.numGridDof();
                for (std::size_t i = 0; i < domain.cells.size(); ++i) {
                    auto& d = (static_cast<unsigned>(domain.cells[i]) >= numGridDof) ? dpAux : dpGrid;
                    d = std::max(d, double(std::abs(xDom[i][Indices::pressureSwitchIdx])));
                }
                for (const auto& b : blocks) {
                    const auto& w = b.fracture->fractureWidth();
                    for (std::size_t i = 0; i < b.nc; ++i) {
                        const double dw = std::abs(x[nFlow + b.offset + i][0]);
                        dwMax = std::max(dwMax, dw);
                        wMax = std::max(wMax, double(w[i][0]));
                        dwRel = std::max(dwRel, dw / std::max(double(w[i][0]), 1e-6));
                    }
                }
                double dtAux = 0.0;
                if constexpr (enableEnergy) {
                    for (std::size_t i = 0; i < domain.cells.size(); ++i) {
                        if (static_cast<unsigned>(domain.cells[i]) >= numGridDof) {
                            dtAux = std::max(dtAux, double(std::abs(xDom[i][Indices::temperatureIdx])));
                        }
                    }
                }
                OpmLog::info(fmt::format("WellFractureJoint step: dp aux {:.3e} grid {:.3e} bar, dT aux {:.3e}, "
                                         "dw {:.3e} (max w {:.3e}, max dw/w {:.3e})",
                                         dpAux / 1e5, dpGrid / 1e5, dtAux, dwMax, wMax, dwRel));
            }
            std::vector<typename std::remove_reference_t<decltype(solution)>::block_type> solutionIt;
            for (const int c : domain.cells) {
                solutionIt.push_back(solution[c]);
            }
            const auto wellStateIt = wellModel.wellState();
            const auto widthsIt = currentWidths();
            const double meritPrev = merit;
            bool accepted = false;
            double alpha = 1.0;
            for (int ls = 0; ls < 6 && !accepted; ++ls, alpha *= 0.5) {
                if (ls > 0) {
                    for (std::size_t i = 0; i < domain.cells.size(); ++i) {
                        solution[domain.cells[i]] = solutionIt[i];
                    }
                    refreshCache();
                    wellModel.wellState() = wellStateIt;
                    well.updatePrimaryVariables(wellModel.groupStateHelper());
                    // the well update below must use equations of this point, not of the rejected one
                    setWidths(widthsIt);
                    well.assembleWellEq(simulator_, dt, wellModel.groupStateHelper(), wellModel.wellState());
                }
                auto widths = widthsIt;
                for (std::size_t k = 0; k < blocks.size(); ++k) {
                    for (std::size_t i = 0; i < blocks[k].nc; ++i) {
                        // relative limit: the conductivity goes with the cube of the width
                        const double limit = std::min(settings.maxDwidth,
                                                      std::max(settings.relativeDwidth * widthsIt[k][i][0],
                                                               blocks[k].fracture->cubicLawMinWidth()));
                        const double dw = std::clamp(alpha * x[nFlow + blocks[k].offset + i][0], -limit, limit);
                        widths[k][i][0] = std::max(0.0, widthsIt[k][i][0] - dw);
                    }
                }
                setWidths(widths);
                dxGlobal = 0.0;
                auto xStep = xDom;
                xStep *= alpha;
                Details::setGlobal(xStep, domain.cells, dxGlobal);
                model.newtonMethod().update_(solution, solution, dxGlobal, dxGlobal, domain.cells);
                refreshCache();
                // the well's B block is indexed by perforation, not by global DOF
                GlobalEqVector dxWell(well.cells().size());
                for (std::size_t i = 0; i < well.cells().size(); ++i) {
                    dxWell[i] = dxGlobal[well.cells()[i]];
                }
                well.recoverWellSolutionAndUpdateWellState(simulator_, dxWell,
                                                           wellModel.groupStateHelper(),
                                                           wellModel.wellState());
                rnorm = assembleAndMeasure();
                // a max over tolerance ratios would reject good steps: the rows start far apart
                accepted = std::isfinite(rnorm) && ((merit < meritPrev) || (ls == 5));
                if (settings.local.verbosity > 1) {
                    OpmLog::info(fmt::format("WellFractureJoint   alpha {:.4f}: merit {:.3e} -> {:.3e} "
                                             "(cnv {:.3e} well {:.3e} aux {:.3e} fb {:.3e})",
                                             alpha, meritPrev, merit, rmax, wmax, amax, fbmax));
                }
            }
            if (!accepted) {
                break;
            }
        }

        if (!rep.converged && !(rep.residual < rep.residual0)) {
            rollback();
            rep.residual = rep.residual0;
        }
        for (auto& b : blocks) {
            this->evaluateContact_(b, model, settings.fbMu);
            b.fracture->setJointContactState(b.pressure, b.isClosed);
        }
        well.assembleWellEq(simulator_, dt, wellModel.groupStateHelper(), wellModel.wellState());
        wellModel.updateCellRates();
        rep.bhp = wellModel.wellState().well(well.indexOfWell()).bhp;
        return rep;
    }

private:
    struct FractureBlock {
        Fracture* fracture = nullptr;
        std::size_t index = 0;
        std::size_t offset = 0;
        std::size_t nc = 0;
        std::vector<int> dof;
        WidthVector traction;
        std::vector<double> phi, da, db, c, scale, pressure;
        std::vector<int> isClosed;
    };

    //! phi_FB(c w, lambda) and its weights at the current width and flow pressure.
    template <class Model>
    void evaluateContact_(FractureBlock& b, const Model& model, const double mu) const
    {
        const auto waterIdx = FluidSystem::waterPhaseIdx;
        const auto& A = b.fracture->jointElasticMatrix();
        const auto& w = b.fracture->fractureWidth();
        const auto& pOwn = b.fracture->fracturePressure();
        b.phi.assign(b.nc, 0.0);
        b.da.assign(b.nc, 0.0);
        b.db.assign(b.nc, 0.0);
        b.c.assign(b.nc, 0.0);
        b.scale.assign(b.nc, 1.0);
        b.pressure.assign(b.nc, 0.0);
        b.isClosed.assign(b.nc, 0);
        constexpr double eps = 1e-30;
        for (std::size_t i = 0; i < b.nc; ++i) {
            const double p = (b.dof[i] >= 0)
                ? getValue(model.intensiveQuantities(b.dof[i], 0).fluidState().pressure(waterIdx))
                : pOwn[i][0];
            b.pressure[i] = p;
            double lambda = b.traction[i][0] - p;
            for (std::size_t j = 0; j < b.nc; ++j) {
                lambda -= A[i][j] * w[j][0];
            }
            const double ci = std::max(std::abs(double(A[i][i])), eps);
            const double a = ci * w[i][0];
            const double r = std::sqrt(a * a + lambda * lambda + 2.0 * mu) + eps;
            b.c[i] = ci;
            b.da[i] = 1.0 - a / r;
            b.db[i] = 1.0 - lambda / r;
            b.phi[i] = a + lambda - std::sqrt(a * a + lambda * lambda + 2.0 * mu);
            b.scale[i] = 1.0 + std::abs(a) + std::abs(lambda);
            b.isClosed[i] = (a < lambda) ? 1 : 0;
        }
    }

    //! dR_flow/dw of one fracture connection: the flux is linear in the transmissibility.
    template <class F>
    void forEachFluxDerivative_(const FractureBlock& b, const double, F&& f) const
    {
        auto& problem = simulator_.problem();
        const auto& model = simulator_.model();
        const auto& neighborInfo = model.linearizer().getNeighborInfo();
        const auto input = problem.fractureAuxCells()->couplingInput(*b.fracture, b.index, -1.0);
        using Eval = LocalAD<2>;
        for (const auto& [ci, cj, t1, t2] : input.htrans) {
            if (ci >= b.nc || cj >= b.nc || b.dof[ci] < 0 || b.dof[cj] < 0) {
                continue;
            }
            const unsigned gi = static_cast<unsigned>(b.dof[ci]);
            const unsigned gj = static_cast<unsigned>(b.dof[cj]);
            const Eval wi = Eval::variable(input.aperture[ci], 0);
            const Eval wj = Eval::variable(input.aperture[cj], 1);
            const Eval trans = detail::connectionTrans(wi, wj, t1, t2, input.cubic_law_min_width);
            if (!(trans.value > 0.0)) {
                continue;
            }
            for (const auto& nb : neighborInfo[gi]) {
                if (nb.neighbor != gj) {
                    continue;
                }
                ADVector flux(0.0), darcy(0.0);
                LocalResidual::computeFlux(flux, darcy, gi, gj,
                                           model.intensiveQuantities(gi, 0),
                                           model.intensiveQuantities(gj, 0),
                                           nb.res_nbinfo, problem.moduleParams());
                for (int e = 0; e < numEq; ++e) {
                    const double perTrans = getValue(flux[e]) * nb.res_nbinfo.faceArea / nb.res_nbinfo.trans;
                    const double dI = perTrans * trans.derivatives[0];
                    const double dJ = perTrans * trans.derivatives[1];
                    f(gi, e, ci, dI);
                    f(gi, e, cj, dJ);
                    f(gj, e, ci, -dI);
                    f(gj, e, cj, -dJ);
                }
                break;
            }
        }
    }

    //! dR_flow/dw of the storage: the cell's volume is its aperture times its area.
    template <class F>
    void forEachStorageDerivative_(const FractureBlock& b, const double dt, F&& f) const
    {
        auto& problem = simulator_.problem();
        const auto& model = simulator_.model();
        const auto input = problem.fractureAuxCells()->couplingInput(*b.fracture, b.index, dt);
        for (std::size_t c = 0; c < b.nc; ++c) {
            if (b.dof[c] < 0 || input.aperture[c] <= input.volume_min_width) {
                continue;
            }
            ADVector storage(0.0);
            LocalResidual::template computeStorage<Evaluation>(storage,
                                                               model.intensiveQuantities(b.dof[c], 0));
            for (int e = 0; e < numEq; ++e) {
                f(static_cast<unsigned>(b.dof[c]), e, c, getValue(storage[e]) * input.area[c] / dt);
            }
        }
    }

    template <class Entries, class WithLinSys>
    void addFlowAndWell_(Entries& entries, const Matrix& jDom, const std::vector<int>&,
                         const std::vector<int>& wellCells,
                         const std::unordered_map<int, int>& localOf,
                         WithLinSys&& withLinSys, const std::size_t nR) const
    {
        for (auto row = jDom.begin(); row != jDom.end(); ++row) {
            for (auto col = row->begin(); col != row->end(); ++col) {
                for (int a = 0; a < numEq; ++a) {
                    for (int b = 0; b < numEq; ++b) {
                        entries.emplace_back(row.index() * numEq + a, col.index() * numEq + b, (*col)[a][b]);
                    }
                }
            }
        }
        auto localOfPerf = [&](const std::size_t perf) {
            const auto it = localOf.find(wellCells[perf]);
            return (it == localOf.end()) ? -1 : it->second;
        };
        withLinSys([&](const auto& eq) {
            const auto& B = eq.getB();
            const auto& C = eq.getC();
            const auto& D = eq.getD();
            const int nW = static_cast<int>(D[0][0].N());
            for (auto row = B.begin(); row != B.end(); ++row) {
                for (auto col = row->begin(); col != row->end(); ++col) {
                    const int l = localOfPerf(col.index());
                    if (l < 0) {
                        continue;
                    }
                    for (int a = 0; a < nW; ++a) {
                        for (int b = 0; b < numEq; ++b) {
                            entries.emplace_back(nR + row.index() * nW + a, l * numEq + b, (*col)[a][b]);
                        }
                    }
                }
            }
            for (auto row = C.begin(); row != C.end(); ++row) {
                for (auto col = row->begin(); col != row->end(); ++col) {
                    const int l = localOfPerf(col.index());
                    if (l < 0) {
                        continue;
                    }
                    for (int a = 0; a < nW; ++a) {
                        for (int b = 0; b < numEq; ++b) {
                            entries.emplace_back(l * numEq + b, nR + row.index() * nW + a, (*col)[a][b]);
                        }
                    }
                }
            }
            for (auto row = D.begin(); row != D.end(); ++row) {
                for (auto col = row->begin(); col != row->end(); ++col) {
                    for (int a = 0; a < nW; ++a) {
                        for (int b = 0; b < nW; ++b) {
                            entries.emplace_back(nR + row.index() * nW + a, nR + col.index() * nW + b, (*col)[a][b]);
                        }
                    }
                }
            }
        });
    }

    template <class Entries>
    void addWidthBlocks_(Entries& entries, const FractureBlock& b,
                         const std::unordered_map<int, int>& localOf,
                         const double dt, const std::size_t nFlow) const
    {
        const auto& model = simulator_.model();
        const auto& A = b.fracture->jointElasticMatrix();
        const auto waterIdx = FluidSystem::waterPhaseIdx;
        const std::size_t first = nFlow + b.offset;
        for (std::size_t i = 0; i < b.nc; ++i) {
            for (std::size_t j = 0; j < b.nc; ++j) {
                const double v = -b.db[i] * A[i][j] + ((i == j) ? b.da[i] * b.c[i] : 0.0);
                entries.emplace_back(first + i, first + j, v);
            }
            if (b.dof[i] >= 0) {
                const auto loc = localOf.find(b.dof[i]);
                if (loc != localOf.end()) {
                    const auto& pw = model.intensiveQuantities(b.dof[i], 0).fluidState().pressure(waterIdx);
                    for (int pv = 0; pv < numEq; ++pv) {
                        entries.emplace_back(first + i, loc->second * numEq + pv,
                                             -b.db[i] * pw.derivative(pv));
                    }
                }
            }
        }
        auto add = [&](const unsigned dof, const int e, const std::size_t cell, const double v) {
            const auto loc = localOf.find(static_cast<int>(dof));
            if (loc != localOf.end()) {
                entries.emplace_back(loc->second * numEq + e, first + cell, v);
            }
        };
        this->forEachFluxDerivative_(b, dt, add);
        this->forEachStorageDerivative_(b, dt, add);
    }

    //! A fracture cell holding only water solves the water row for pressure. Its other mass
    //! rows have no storage to act on, and a saturation step there would create gas.
    template <class Entries>
    void pinSaturationsOfWaterCells_(Entries& entries, ScalarVector& rhs,
                                     const std::vector<int>& cells) const
    {
        const auto& model = simulator_.model();
        const auto numGridDof = model.numGridDof();
        const auto waterEq = Indices::conti0EqIdx
            + FluidSystem::canonicalToActiveCompIdx(FluidSystem::waterCompIdx);
        std::vector<std::pair<int, int>> pinned;
        for (int e = 0; e < numEq; ++e) {
            if (e >= Indices::conti0EqIdx && e < Indices::conti0EqIdx + Indices::numPhases
                && e != waterEq) {
                pinned.emplace_back(e, pinned.empty() ? Indices::waterSwitchIdx
                                                      : Indices::compositionSwitchIdx);
            }
        }
        std::vector<char> dropRow(rhs.size(), 0), dropCol(rhs.size(), 0);
        for (std::size_t l = 0; l < cells.size(); ++l) {
            if (static_cast<unsigned>(cells[l]) < numGridDof) {
                continue;
            }
            // Temperature is slow and left to the global Newton: in a nearly closed cell the
            // energy row hardly determines it, and limited steps there make the mass row flicker.
            if constexpr (enableEnergy) {
                dropRow[l * numEq + Indices::contiEnergyEqIdx] = 1;
                dropCol[l * numEq + Indices::temperatureIdx] = 1;
                rhs[l * numEq + Indices::contiEnergyEqIdx] = 0.0;
            }
            const auto sw = getValue(model.intensiveQuantities(cells[l], 0)
                                     .fluidState().saturation(FluidSystem::waterPhaseIdx));
            if (sw < 1.0 - 1e-6) {
                continue;
            }
            for (const auto& [e, pv] : pinned) {
                dropRow[l * numEq + e] = 1;
                dropCol[l * numEq + pv] = 1;
                rhs[l * numEq + e] = 0.0;
            }
        }
        entries.erase(std::remove_if(entries.begin(), entries.end(),
                                     [&](const auto& t) {
                                         return dropRow[std::get<0>(t)] || dropCol[std::get<1>(t)];
                                     }),
                      entries.end());
        for (std::size_t l = 0; l < cells.size(); ++l) {
            for (const auto& [e, pv] : pinned) {
                if (dropRow[l * numEq + e]) {
                    entries.emplace_back(l * numEq + e, l * numEq + pv, 1.0);
                }
            }
            if constexpr (enableEnergy) {
                if (dropRow[l * numEq + Indices::contiEnergyEqIdx]) {
                    entries.emplace_back(l * numEq + Indices::contiEnergyEqIdx,
                                         l * numEq + Indices::temperatureIdx, 1.0);
                }
            }
        }
    }

    template <class Entries>
    void solveSparse_(const Entries& entries, const std::size_t n,
                      ScalarVector& rhs, ScalarVector& x) const
    {
        Dune::MatrixIndexSet pattern(n, n);
        for (const auto& [r, c, v] : entries) {
            pattern.add(r, c);
        }
        ScalarMatrix S;
        pattern.exportIdx(S);
        S = 0.0;
        for (const auto& [r, c, v] : entries) {
            S[r][c] += v;
        }
        Dune::UMFPack<ScalarMatrix> lu(S, 0, false);
        Dune::InverseOperatorResult res;
        lu.apply(x, rhs, res);
    }

    //! Compares the assembled dR_flow/dw against central differences of the flow residual.
    template <class SetWidths, class Measure>
    void checkWidthColumns_(const std::vector<FractureBlock>& blocks,
                            const std::unordered_map<int, int>& localOf, const double dt,
                            const int columns, SetWidths&& setWidths, Measure&& measure) const
    {
        auto& residual = simulator_.model().linearizer().residual();
        std::vector<WidthVector> base;
        for (const auto& b : blocks) {
            base.push_back(b.fracture->fractureWidth());
        }
        double worst = 0.0;
        int compared = 0;
        for (std::size_t k = 0; k < blocks.size(); ++k) {
            const auto& b = blocks[k];
            std::vector<std::unordered_map<std::size_t, double>> analytic(b.nc);
            auto add = [&](const unsigned dof, const int e, const std::size_t cell, const double v) {
                const auto loc = localOf.find(static_cast<int>(dof));
                if (loc != localOf.end()) {
                    analytic[cell][loc->second * numEq + e] += v;
                }
            };
            measure();
            this->forEachFluxDerivative_(b, dt, add);
            this->forEachStorageDerivative_(b, dt, add);
            const std::size_t stride = std::max<std::size_t>(1, b.nc / std::max(columns, 1));
            for (std::size_t col = 0; col < b.nc; col += stride) {
                const double dw = 1e-6 * std::max(double(base[k][col][0]), 1e-4);
                auto evaluate = [&](const double sign) {
                    auto widths = base;
                    widths[k][col][0] += sign * dw;
                    setWidths(widths);
                    measure();
                    std::vector<double> r;
                    for (const auto& [idx, v] : analytic[col]) {
                        const auto cell = idx / numEq;
                        r.push_back(residual[domainCell_(localOf, cell)][idx % numEq]);
                    }
                    return r;
                };
                const auto plus = evaluate(1.0);
                const auto minus = evaluate(-1.0);
                std::size_t m = 0;
                for (const auto& [idx, v] : analytic[col]) {
                    const double fd = (plus[m] - minus[m]) / (2.0 * dw);
                    ++m;
                    const double err = std::abs(fd - v) / std::max({std::abs(fd), std::abs(v), 1e-12});
                    worst = std::max(worst, err);
                    ++compared;
                    OpmLog::info(fmt::format("WellFractureJoint FD fracture {} col {} row {}: ad {:.6e} fd {:.6e}",
                                             k, col, idx, v, fd));
                }
            }
        }
        setWidths(base);
        OpmLog::info(fmt::format("WellFractureJoint FD check: {} entries, worst relative error {:.3e}",
                                 compared, worst));
    }

    static int domainCell_(const std::unordered_map<int, int>& localOf, const std::size_t local)
    {
        for (const auto& [global, l] : localOf) {
            if (static_cast<std::size_t>(l) == local) {
                return global;
            }
        }
        return -1;
    }

    Simulator& simulator_;
};

} // namespace Opm

#endif // OPM_WELL_FRACTURE_JOINT_SOLVER_HPP
