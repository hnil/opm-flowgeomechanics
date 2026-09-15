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
#ifndef OPM_FRACTURE_MECH_COUPLING_HPP
#define OPM_FRACTURE_MECH_COUPLING_HPP

#include <opm/geomech/FracturePressureAssemblerAD.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace Opm
{

/*!
 * \file
 * \brief The two coupling blocks between an embedded fracture's flow system and
 *        its mechanics system, built with AD, plus a finite-difference check.
 *
 * When the fracture's cells are degrees of freedom of the flow problem, the
 * flow and the mechanics are coupled both ways:
 *
 *  - flow depends on the aperture: the cubic law makes the intra-fracture
 *    transmissibility go with \f$w^3\f$, and the cell's pore volume is its
 *    aperture times its area.  That is \f$\partial R_{flow}/\partial w\f$,
 *    the block this file calls the *flow-mech* coupling.  It is the non-trivial
 *    one and the reason the local well+fracture solve cannot hold the aperture
 *    fixed while it moves the pressure.
 *
 *  - mechanics depends on the fluid pressure, through the normal traction on
 *    an open cell.  That is \f$\partial R_{mech}/\partial p\f$, the *mech-flow*
 *    coupling: structurally an identity masked by the open/closed state.
 *
 * Both are assembled here from the same kernels the residual uses, so a
 * finite-difference check of the matrix against that residual is a check of
 * the production code and not of a separate copy of it.
 */

//! Everything the two coupling blocks are formed from.
struct MechCouplingInput
{
    //! Intra-fracture half transmissibilities, (i, j, t1, t2); the same list the
    //! fracture's own pressure assembly uses.
    std::vector<Htrans> htrans {};

    std::vector<double> aperture {};      //!< \f$w_i\f$ [m]
    std::vector<double> pressure {};      //!< \f$p_i\f$ [Pa], already a head if gravity is carried
    std::vector<double> area {};          //!< \f$A_i\f$ [m2], storage term only
    std::vector<double> face_mobility {}; //!< per-cell mobility, averaged onto the faces
    std::vector<double> density {};       //!< \f$\rho_i\f$, storage term only
    std::vector<double> volume_prev {};   //!< pore volume at the start of the step [m3]

    //! 1 where the cell is mechanically open, 0 where contact holds it shut.
    std::vector<int> open {};

    //! Floor under the aperture in the cubic law (the fracture's config.min_width).
    double cubic_law_min_width {1.0e-4};
    //! Floor under the aperture in the pore volume (solver.embedded_min_aperture).
    double volume_min_width {1.0e-3};
    //! Timestep [s]; <= 0 leaves the storage term out.
    double dt {-1.0};

    std::size_t numCells() const
    { return this->aperture.size(); }
};

namespace detail
{
    //! Lets one kernel serve both the AD assembly and the plain-double residual
    //! the finite-difference check evaluates.
    template <class Eval>
    struct EvalTraits
    {
        static Eval constant(const double v) { return Eval::constant(v); }
        static double value(const Eval& x) { return x.value; }
    };

    template <>
    struct EvalTraits<double>
    {
        static double constant(const double v) { return v; }
        static double value(const double x) { return x; }
    };

    //! Aperture entering the cubic law: floored, with a zero derivative below the
    //! floor -- a clamped cell's conductivity does not respond to its width.
    template <class Eval>
    Eval cubicAperture(const Eval& w, const double floor)
    {
        using T = EvalTraits<Eval>;
        return (T::value(w) > floor) ? w : T::constant(floor);
    }

    //! Transmissibility of one intra-fracture connection, cubic law over the two
    //! half transmissibilities -- the same combination FracturePressureAssemblerAD
    //! and FractureAuxCells::bind() form.
    template <class Eval>
    Eval connectionTrans(const Eval& wi, const Eval& wj,
                         const double t1, const double t2,
                         const double floor)
    {
        using T = EvalTraits<Eval>;
        const Eval hi = cubicAperture(wi, floor);
        const Eval hj = cubicAperture(wj, floor);
        const Eval inv = T::constant(12.0) / (hi * hi * hi * T::constant(t1))
            + T::constant(12.0) / (hj * hj * hj * T::constant(t2));
        return T::constant(1.0) / inv;
    }

    //! Pore volume of a cell, floored like the flow's auxiliary cell.
    template <class Eval>
    Eval poreVolume(const Eval& w, const double area, const double floor)
    {
        using T = EvalTraits<Eval>;
        const Eval h = (T::value(w) > floor) ? w : T::constant(floor);
        return h * T::constant(area);
    }
} // namespace detail

/*!
 * \brief The flow residual of the fracture cells as a function of the apertures.
 *
 * \f[ R_i = \sum_j T_{ij}(w)\,\lambda_{ij}\,(p_i - p_j)
 *          + \rho_i \frac{V_i(w_i) - V_i^{n}}{\Delta t} \f]
 *
 * Outflow positive.  Only the aperture varies; pressures, mobilities and the
 * leak-off connections are held at the values in \p input, which is what makes
 * this the residual whose derivative the flow-mech block is.
 */
inline std::vector<double>
flowResidualFromAperture(const MechCouplingInput& input,
                         const std::vector<double>& w)
{
    const std::size_t nc = input.numCells();
    std::vector<double> res(nc, 0.0);

    for (const auto& [i, j, t1, t2] : input.htrans) {
        if (i >= nc || j >= nc) {
            continue;
        }
        const double trans = detail::connectionTrans(w[i], w[j], t1, t2,
                                                     input.cubic_law_min_width);
        const double mob = 0.5 * (input.face_mobility[i] + input.face_mobility[j]);
        const double flux = trans * mob * (input.pressure[i] - input.pressure[j]);
        res[i] += flux;
        res[j] -= flux;
    }

    if (input.dt > 0.0) {
        for (std::size_t i = 0; i < nc; ++i) {
            const double vol = detail::poreVolume(w[i], input.area[i],
                                                  input.volume_min_width);
            const double vol0 = (i < input.volume_prev.size()) ? input.volume_prev[i] : vol;
            res[i] += input.density[i] * (vol - vol0) / input.dt;
        }
    }

    return res;
}

//! Sparsity of the flow-mech block: a cell couples to itself and to every cell
//! it shares a connection with.
inline std::unique_ptr<BCRSMatrix1x1>
buildFlowMechSparsity(const MechCouplingInput& input)
{
    const std::size_t nc = input.numCells();
    auto mat = std::make_unique<BCRSMatrix1x1>(nc, nc, 6, 0.4, BCRSMatrix1x1::implicit);

    for (std::size_t i = 0; i < nc; ++i) {
        mat->entry(i, i) = 0.0;
    }
    for (const auto& [i, j, t1, t2] : input.htrans) {
        static_cast<void>(t1);
        static_cast<void>(t2);
        if (i >= nc || j >= nc) {
            continue;
        }
        mat->entry(i, j) = 0.0;
        mat->entry(j, i) = 0.0;
    }
    mat->compress();
    return mat;
}

/*!
 * \brief \f$\partial R_{flow} / \partial w\f$, by automatic differentiation.
 *
 * The production path: every entry comes from the same kernels
 * flowResidualFromAperture() evaluates, so the matrix and the residual cannot
 * drift apart.
 */
inline std::unique_ptr<BCRSMatrix1x1>
buildFlowMechCoupling(const MechCouplingInput& input)
{
    using Eval = LocalAD<2>;
    constexpr int W_I = 0;
    constexpr int W_J = 1;

    const std::size_t nc = input.numCells();
    auto mat = buildFlowMechSparsity(input);
    *mat = 0.0;

    for (const auto& [i, j, t1, t2] : input.htrans) {
        if (i >= nc || j >= nc) {
            continue;
        }
        const Eval wi = Eval::variable(input.aperture[i], W_I);
        const Eval wj = Eval::variable(input.aperture[j], W_J);
        const Eval trans = detail::connectionTrans(wi, wj, t1, t2,
                                                   input.cubic_law_min_width);
        const double mob = 0.5 * (input.face_mobility[i] + input.face_mobility[j]);
        const double dp = input.pressure[i] - input.pressure[j];
        const Eval flux = trans * Eval::constant(mob * dp);

        (*mat)[i][i] += flux.derivatives[W_I];
        (*mat)[i][j] += flux.derivatives[W_J];
        (*mat)[j][i] -= flux.derivatives[W_I];
        (*mat)[j][j] -= flux.derivatives[W_J];
    }

    if (input.dt > 0.0) {
        for (std::size_t i = 0; i < nc; ++i) {
            const Eval wi = Eval::variable(input.aperture[i], W_I);
            const Eval vol = detail::poreVolume(wi, input.area[i],
                                                input.volume_min_width);
            (*mat)[i][i] += input.density[i] * vol.derivatives[W_I] / input.dt;
        }
    }

    return mat;
}

/*!
 * \brief \f$\partial R_{mech} / \partial p\f$: the identity masked by the
 *        open/closed state.
 *
 * A closed cell's mechanics row is replaced by \f$w_i = 0\f$ and carries no
 * pressure, so its column is empty; an open cell's normal traction takes the
 * fluid pressure with unit coefficient.  Structurally trivial -- the check
 * below guards the mask and the sparsity, not an algebraic expression.
 */
inline std::unique_ptr<BCRSMatrix1x1>
buildMechFlowCoupling(const MechCouplingInput& input)
{
    const std::size_t nc = input.numCells();
    auto mat = std::make_unique<BCRSMatrix1x1>(nc, nc, 2, 0.4, BCRSMatrix1x1::implicit);
    for (std::size_t i = 0; i < nc; ++i) {
        mat->entry(i, i) = 0.0;
    }
    mat->compress();
    *mat = 0.0;

    for (std::size_t i = 0; i < nc; ++i) {
        const bool isOpen = (i < input.open.size()) ? (input.open[i] != 0) : true;
        (*mat)[i][i] = isOpen ? 1.0 : 0.0;
    }
    return mat;
}

//! The mechanics residual's pressure part, the functional buildMechFlowCoupling()
//! differentiates: an open cell takes the fluid pressure, a closed one does not.
inline std::vector<double>
mechPressureTermFromPressure(const MechCouplingInput& input,
                             const std::vector<double>& p)
{
    const std::size_t nc = input.numCells();
    std::vector<double> res(nc, 0.0);
    for (std::size_t i = 0; i < nc; ++i) {
        const bool isOpen = (i < input.open.size()) ? (input.open[i] != 0) : true;
        res[i] = isOpen ? p[i] : 0.0;
    }
    return res;
}

//! What to check, and how hard.
struct CouplingCheckOptions
{
    bool enabled {false};        //!< master switch, so a production call is one line
    bool check_flow_mech {true}; //!< dR_flow/dw, the non-trivial block
    bool check_mech_flow {true}; //!< dR_mech/dp, the masked identity
    //! Columns to differentiate; <= 0 checks every one.  A smaller number is
    //! spread evenly over the matrix rather than taken from the front, so a
    //! cheap check still visits open and closed cells alike.
    int max_columns {-1};
    //! Relative step on the aperture, scaled by max(|w|, floor).
    double perturbation {1.0e-6};
    //! Relative tolerance; entries below abs_floor are compared absolutely.
    double tolerance {1.0e-5};
    double abs_floor {1.0e-12};
    int verbosity {0};
};

//! Outcome of a check.
struct CouplingCheckReport
{
    int checked {0};
    int failed {0};
    double max_rel_error {0.0};
    std::size_t worst_row {0};
    std::size_t worst_col {0};
    double worst_ad {0.0};
    double worst_fd {0.0};
    std::string block {};
    //! Largest aperture in the case, and the floor it is measured against: a
    //! check that compared nothing is not a passing check, it means every
    //! sampled cell sat at the floor where the block is identically zero.
    double max_aperture {0.0};
    double cubic_floor {0.0};

    bool ok() const
    { return this->failed == 0; }

    std::string summary() const
    {
        std::ostringstream os;
        os << "coupling check [" << this->block << "]: " << this->checked
           << " entries, " << this->failed << " failed, max relative error "
           << this->max_rel_error;
        if (this->checked == 0) {
            os << " (nothing to compare: max aperture " << this->max_aperture
               << " m against a cubic-law floor of " << this->cubic_floor
               << " m, so the block is identically zero)";
        }
        if (this->failed > 0) {
            os << " (worst at (" << this->worst_row << ',' << this->worst_col
               << "): AD " << this->worst_ad << " vs FD " << this->worst_fd << ')';
        }
        return os.str();
    }
};

namespace detail
{
    //! Columns to visit: all of the candidates, or maxColumns of them spread
    //! evenly over the list.
    inline std::vector<std::size_t>
    spreadOver(const std::vector<std::size_t>& candidates, const int maxColumns)
    {
        const std::size_t nc = candidates.size();
        if ((maxColumns <= 0) || (static_cast<std::size_t>(maxColumns) >= nc)) {
            return candidates;
        }
        const auto n = static_cast<std::size_t>(maxColumns);
        std::vector<std::size_t> cols;
        cols.reserve(n);
        for (std::size_t k = 0; k < n; ++k) {
            cols.push_back(candidates[(k * nc) / n]);
        }
        return cols;
    }

    inline std::vector<std::size_t>
    allColumns(const std::size_t nc)
    {
        std::vector<std::size_t> cols(nc);
        for (std::size_t k = 0; k < nc; ++k) {
            cols[k] = k;
        }
        return cols;
    }

    /*!
     * \brief Columns worth differentiating.
     *
     * A cell whose aperture is at the cubic-law floor has an identically zero
     * column -- the flow does not respond to its width at all -- so sampling
     * those proves nothing.  On a fracture that is mostly shut an even spread
     * hits almost only such cells, so prefer the ones above the floor when
     * there are any.
     */
    inline std::vector<std::size_t>
    columnsToCheck(const MechCouplingInput& input, const int maxColumns)
    {
        const std::size_t nc = input.numCells();
        std::vector<std::size_t> live;
        for (std::size_t k = 0; k < nc; ++k) {
            if (input.aperture[k] > input.cubic_law_min_width) {
                live.push_back(k);
            }
        }
        if (live.empty()) {
            return spreadOver(allColumns(nc), maxColumns);
        }
        return spreadOver(live, maxColumns);
    }

    inline void
    accumulate(CouplingCheckReport& rep, const std::size_t i, const std::size_t k,
               const double ad, const double fd, const CouplingCheckOptions& opt)
    {
        ++rep.checked;
        const double scale = std::max({std::abs(ad), std::abs(fd), opt.abs_floor});
        const double rel = std::abs(ad - fd) / scale;
        if (rel > rep.max_rel_error) {
            rep.max_rel_error = rel;
            rep.worst_row = i;
            rep.worst_col = k;
            rep.worst_ad = ad;
            rep.worst_fd = fd;
        }
        if (rel > opt.tolerance) {
            ++rep.failed;
        }
    }
} // namespace detail

/*!
 * \brief Check the flow-mech block against finite differences of the residual.
 *
 * Central differences where the aperture is safely above the cubic-law floor;
 * a one-sided step away from it otherwise, because the floor is a kink and the
 * matrix follows the clamped branch there by construction.
 */
inline CouplingCheckReport
checkFlowMechCouplingFD(const MechCouplingInput& input,
                        const CouplingCheckOptions& opt)
{
    CouplingCheckReport rep;
    rep.block = "dRflow/dw";
    rep.cubic_floor = input.cubic_law_min_width;
    rep.max_aperture = input.aperture.empty()
        ? 0.0 : *std::max_element(input.aperture.begin(), input.aperture.end());

    const std::size_t nc = input.numCells();
    const auto mat = buildFlowMechCoupling(input);
    const auto cols = detail::columnsToCheck(input, opt.max_columns);

    for (const auto k : cols) {
        const double floor = std::min(input.cubic_law_min_width, input.volume_min_width);
        const double scale = std::max(std::abs(input.aperture[k]), floor);
        const double eps = opt.perturbation * scale;

        // Stay on one branch of the floor: a step that crosses it would compare
        // a one-sided derivative against a two-sided difference.
        const double w0 = input.aperture[k];
        const bool clamped = (w0 <= std::max(input.cubic_law_min_width, input.volume_min_width));
        std::vector<double> wPlus = input.aperture;
        std::vector<double> wMinus = input.aperture;
        double denom = 2.0 * eps;
        if (clamped) {
            // below the kink the matrix says zero; differentiate downwards, where
            // the clamp stays active, so the difference must say zero as well
            wMinus[k] = w0 - eps;
            denom = eps;
        } else {
            wPlus[k] = w0 + eps;
            wMinus[k] = w0 - eps;
        }

        const auto rPlus = flowResidualFromAperture(input, wPlus);
        const auto rMinus = flowResidualFromAperture(input, wMinus);

        for (std::size_t i = 0; i < nc; ++i) {
            double ad = 0.0;
            const auto& row = (*mat)[i];
            const auto it = row.find(k);
            if (it != row.end()) {
                ad = (*it)[0][0];
            }
            const double fd = (rPlus[i] - rMinus[i]) / denom;
            if ((ad == 0.0) && (fd == 0.0)) {
                continue; // structural zero on both sides: nothing to compare
            }
            detail::accumulate(rep, i, k, ad, fd, opt);
        }
    }

    return rep;
}

//! Check the mech-flow block against finite differences of the pressure term.
inline CouplingCheckReport
checkMechFlowCouplingFD(const MechCouplingInput& input,
                        const CouplingCheckOptions& opt)
{
    CouplingCheckReport rep;
    rep.block = "dRmech/dp";
    rep.cubic_floor = input.cubic_law_min_width;
    rep.max_aperture = input.aperture.empty()
        ? 0.0 : *std::max_element(input.aperture.begin(), input.aperture.end());

    const std::size_t nc = input.numCells();
    const auto mat = buildMechFlowCoupling(input);
    const auto cols = detail::spreadOver(detail::allColumns(nc), opt.max_columns);

    for (const auto k : cols) {
        const double scale = std::max(std::abs(input.pressure[k]), 1.0);
        const double eps = opt.perturbation * scale;
        auto pPlus = input.pressure;
        auto pMinus = input.pressure;
        pPlus[k] += eps;
        pMinus[k] -= eps;

        const auto rPlus = mechPressureTermFromPressure(input, pPlus);
        const auto rMinus = mechPressureTermFromPressure(input, pMinus);

        for (std::size_t i = 0; i < nc; ++i) {
            double ad = 0.0;
            const auto& row = (*mat)[i];
            const auto it = row.find(k);
            if (it != row.end()) {
                ad = (*it)[0][0];
            }
            const double fd = (rPlus[i] - rMinus[i]) / (2.0 * eps);
            if ((ad == 0.0) && (fd == 0.0)) {
                continue;
            }
            detail::accumulate(rep, i, k, ad, fd, opt);
        }
    }

    return rep;
}

/*!
 * \brief Both checks, for a production call site.
 *
 * Returns false if either block failed.  Costs one residual evaluation per
 * checked column, so leave max_columns small when this runs inside a solve.
 */
inline bool
checkCouplingMatricesFD(const MechCouplingInput& input,
                        const CouplingCheckOptions& opt,
                        std::vector<CouplingCheckReport>* reports = nullptr)
{
    if (!opt.enabled) {
        return true;
    }
    bool ok = true;
    if (opt.check_flow_mech) {
        const auto rep = checkFlowMechCouplingFD(input, opt);
        ok = ok && rep.ok();
        if (reports != nullptr) {
            reports->push_back(rep);
        }
    }
    if (opt.check_mech_flow) {
        const auto rep = checkMechFlowCouplingFD(input, opt);
        ok = ok && rep.ok();
        if (reports != nullptr) {
            reports->push_back(rep);
        }
    }
    return ok;
}

} // namespace Opm

#endif // OPM_FRACTURE_MECH_COUPLING_HPP
