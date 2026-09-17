/*
  Finite-difference verification of the fracture <-> mechanics coupling blocks.

  The production code builds both blocks with AD (FractureMechCoupling.hpp);
  this test differentiates the same residual kernels numerically on a small
  fracture and compares.  It also exercises the "check only a few columns"
  path that a production run would use, and the cases the kink in the aperture
  floor creates.
*/
#include <config.h>

#include <opm/geomech/FractureMechCoupling.hpp>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace
{

//! A small fracture: `nc` cells in a chain, apertures spanning the floors, one
//! closed cell, a pressure gradient along the chain.
Opm::MechCouplingInput makeSmallCase(const std::size_t nc = 6)
{
    Opm::MechCouplingInput in;
    in.cubic_law_min_width = 1.0e-4;
    in.volume_min_width = 1.0e-3;
    in.dt = 86400.0;

    in.aperture.resize(nc);
    in.pressure.resize(nc);
    in.area.resize(nc);
    in.face_mobility.resize(nc);
    in.density.resize(nc);
    in.volume_prev.resize(nc);
    in.open.resize(nc);

    for (std::size_t i = 0; i < nc; ++i) {
        // a wide cell at the well, tapering to below both floors at the tip
        in.aperture[i] = 5.0e-3 * std::pow(0.6, static_cast<double>(i));
        in.pressure[i] = 300.0e5 - 2.0e5 * static_cast<double>(i);
        in.area[i] = 4.0 + 0.5 * static_cast<double>(i);
        in.face_mobility[i] = 1.0e-3 * (1.0 + 0.1 * static_cast<double>(i));
        in.density[i] = 1000.0;
        in.volume_prev[i] = in.area[i] * in.volume_min_width;
        in.open[i] = (i == nc - 2) ? 0 : 1; // one closed cell
    }

    for (std::size_t i = 0; i + 1 < nc; ++i) {
        in.htrans.emplace_back(i, i + 1, 3.0 + 0.2 * static_cast<double>(i), 2.5);
    }

    return in;
}

bool report(const Opm::CouplingCheckReport& rep, const std::string& what)
{
    std::cout << "  " << what << ": " << rep.summary() << '\n';
    return rep.ok();
}

bool testFullCheck()
{
    std::cout << "Test: every column, both blocks\n";
    const auto in = makeSmallCase();

    Opm::CouplingCheckOptions opt;
    opt.enabled = true;
    opt.tolerance = 1.0e-5;

    std::vector<Opm::CouplingCheckReport> reps;
    const bool ok = Opm::checkCouplingMatricesFD(in, opt, &reps);
    bool all = ok;
    for (const auto& r : reps) {
        all = report(r, "block") && all;
    }
    if (reps.size() != 2) {
        std::cerr << "  expected two blocks, got " << reps.size() << '\n';
        return false;
    }
    if (reps.front().checked == 0) {
        std::cerr << "  the flow-mech check compared nothing\n";
        return false;
    }
    return all;
}

//! Only a few columns: the mode a production run would switch on.
bool testFewColumns()
{
    std::cout << "Test: a few columns only\n";
    const auto in = makeSmallCase(12);

    Opm::CouplingCheckOptions opt;
    opt.enabled = true;
    opt.max_columns = 3;
    opt.check_mech_flow = false;

    std::vector<Opm::CouplingCheckReport> reps;
    const bool ok = Opm::checkCouplingMatricesFD(in, opt, &reps);
    if (!report(reps.front(), "flow-mech")) {
        return false;
    }

    Opm::CouplingCheckOptions all = opt;
    all.max_columns = -1;
    std::vector<Opm::CouplingCheckReport> repsAll;
    Opm::checkCouplingMatricesFD(in, all, &repsAll);

    if (!(reps.front().checked < repsAll.front().checked)) {
        std::cerr << "  limiting the columns did not reduce the work: "
                  << reps.front().checked << " vs " << repsAll.front().checked << '\n';
        return false;
    }
    std::cout << "  checked " << reps.front().checked << " entries instead of "
              << repsAll.front().checked << '\n';
    return ok;
}

//! Every cell below the floors: the clamped branch, where the matrix must be
//! zero in the cubic-law part and the volume part alike.
bool testAllClamped()
{
    std::cout << "Test: all apertures below the floors\n";
    auto in = makeSmallCase();
    for (auto& w : in.aperture) {
        w = 1.0e-6;
    }

    Opm::CouplingCheckOptions opt;
    opt.enabled = true;

    std::vector<Opm::CouplingCheckReport> reps;
    const bool ok = Opm::checkCouplingMatricesFD(in, opt, &reps);
    for (const auto& r : reps) {
        report(r, "block");
    }

    const auto mat = Opm::buildFlowMechCoupling(in);
    for (auto row = mat->begin(); row != mat->end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            if (std::abs((*col)[0][0]) > 0.0) {
                std::cerr << "  clamped cells still couple: (" << row.index() << ','
                          << col.index() << ") = " << (*col)[0][0] << '\n';
                return false;
            }
        }
    }
    return ok;
}

//! Without a timestep there is no storage term, so the diagonal is the flux
//! part alone -- and must still match finite differences.
bool testNoStorage()
{
    std::cout << "Test: no storage term\n";
    auto in = makeSmallCase();
    in.dt = -1.0;

    Opm::CouplingCheckOptions opt;
    opt.enabled = true;
    opt.check_mech_flow = false;

    std::vector<Opm::CouplingCheckReport> reps;
    const bool ok = Opm::checkCouplingMatricesFD(in, opt, &reps);
    return report(reps.front(), "flow-mech") && ok;
}

//! A negative control: the comparison must reject a matrix that does not belong
//! to the residual it is checked against.  Without this the test would pass just
//! as happily on a checker that compares nothing.
bool testCatchesAnError()
{
    std::cout << "Test: a mismatched matrix is rejected\n";
    const auto in = makeSmallCase();

    Opm::CouplingCheckOptions opt;
    opt.enabled = true;
    opt.tolerance = 1.0e-5;

    if (!Opm::checkFlowMechCouplingFD(in, opt).ok()) {
        std::cerr << "  the unmodified matrix already fails\n";
        return false;
    }

    // The AD matrix of `in` against finite differences of a case whose mobility
    // is twice as large: every flux entry must now be out by a factor of two.
    auto other = in;
    for (auto& m : other.face_mobility) {
        m *= 2.0;
    }
    const auto mat = Opm::buildFlowMechCoupling(in);

    int checked = 0;
    int flagged = 0;
    const std::size_t nc = in.numCells();
    for (std::size_t k = 0; k < nc; ++k) {
        const double eps = opt.perturbation * std::max(in.aperture[k], in.cubic_law_min_width);
        auto wPlus = other.aperture;
        auto wMinus = other.aperture;
        wPlus[k] += eps;
        wMinus[k] -= eps;
        const auto rPlus = Opm::flowResidualFromAperture(other, wPlus);
        const auto rMinus = Opm::flowResidualFromAperture(other, wMinus);
        for (std::size_t i = 0; i < nc; ++i) {
            double ad = 0.0;
            const auto it = (*mat)[i].find(k);
            if (it != (*mat)[i].end()) {
                ad = (*it)[0][0];
            }
            const double fd = (rPlus[i] - rMinus[i]) / (2.0 * eps);
            const double scale = std::max({std::abs(ad), std::abs(fd), 1.0e-12});
            if (scale <= 1.0e-12) {
                continue;
            }
            ++checked;
            if (std::abs(ad - fd) / scale > opt.tolerance) {
                ++flagged;
            }
        }
    }

    if (checked == 0) {
        std::cerr << "  nothing was compared\n";
        return false;
    }
    if (flagged == 0) {
        std::cerr << "  a matrix from a different case passed: the check does not "
                     "discriminate\n";
        return false;
    }
    std::cout << "  " << flagged << " of " << checked
              << " entries flagged, as they should be\n";
    return true;
}

} // anonymous namespace

int main()
{
    std::cout << "=== fracture <-> mechanics coupling, AD against finite differences ===\n";

    const bool ok = testFullCheck()
        && testFewColumns()
        && testAllClamped()
        && testNoStorage()
        && testCatchesAnError();

    std::cout << (ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
    return ok ? 0 : 1;
}
