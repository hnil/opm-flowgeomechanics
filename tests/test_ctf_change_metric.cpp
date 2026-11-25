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

#include <config.h>

#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/geomech/FractureCouplingOuterLoop.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Opm::Properties {
namespace TTag {
struct CtfMetricTest {};
}
template <class TypeTag>
struct EqVector<TypeTag, TTag::CtfMetricTest> { using type = std::vector<double>; };
template <class TypeTag>
struct SolutionVector<TypeTag, TTag::CtfMetricTest> { using type = std::vector<double>; };
} // namespace Opm::Properties

using namespace Opm;
using Perfs = std::vector<std::vector<RuntimePerforation>>;

namespace {
void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

struct MockProblem {
    PropertyTree prm;
    Perfs current;
    const PropertyTree& getFractureParam() const { return prm; }
    Perfs getAllExtraWellIndices() const { return current; }
};
struct MockSimulator {
    MockProblem p;
    const MockProblem& problem() const { return p; }
};

struct Loop : FractureCouplingOuterLoop<Properties::TTag::CtfMetricTest, Loop> {
    MockSimulator simulator_;
};

Perfs scaled(const std::vector<double>& ctf, double s)
{
    std::vector<RuntimePerforation> w;
    for (size_t i = 0; i < ctf.size(); ++i) {
        RuntimePerforation p;
        p.cell = static_cast<int>(i);
        p.ctf = ctf[i] * s;
        w.push_back(p);
    }
    return {w};
}

Loop loop(double threshold_well)
{
    Loop l;
    l.simulator_.p.prm.put("solver.verbosity", 0);
    if (threshold_well > 0.0) {
        l.simulator_.p.prm.put("solver.ctf_change_threshold_well", threshold_well);
    }
    return l;
}

// A 5 % change of the well's CTF sum counts as a change at any CTF magnitude.
void scaleFree()
{
    const std::vector<double> before {1.0, 2.0, 3.0};
    const std::vector<double> after {1.0, 2.0, 3.3}; // sum 6 -> 6.3
    for (const double s : {1e-22, 1e-12, 1.0, 1e6}) {
        Loop l = loop(1e-3);
        l.simulator_.p.current = scaled(after, s);
        const auto m = l.evaluateFractureCouplingChange(scaled(before, s));
        require(std::abs(m.max_ctf_change_well - 0.3 / 6.3) < 1e-12, "well metric is change sum over CTF sum");
        require(m.ctf_changed, "a 5 % change is a change at every scale");
        require(std::abs(m.composite_norm - (0.3 / 6.3) / 1e-3) < 1e-6, "composite norm uses the well metric");
    }
    std::cout << "scale free: ok for CTF scale 1e-22 .. 1e6\n";
}

// Negative control: the default absolute test misses the same growth on a tiny seed.
void absoluteDefaultMissesTinySeed()
{
    Loop l = loop(0.0);
    l.simulator_.p.current = scaled({1.0, 2.0, 3.3}, 1e-22);
    const auto m = l.evaluateFractureCouplingChange(scaled({1.0, 2.0, 3.0}, 1e-22));
    require(!m.ctf_changed, "negative control: absolute 1e-15 test sees no change at CTF 1e-22");
    std::cout << "absolute default: tiny-seed growth not seen (as documented)\n";
}

// Below the threshold the step counts as unchanged; the metric takes the worst well.
void thresholdAndWorstWell()
{
    Loop l = loop(1e-2);
    Perfs before = scaled({1.0, 1.0}, 1.0);
    Perfs after = scaled({1.0, 1.001}, 1.0);
    before.push_back(scaled({5.0}, 1e-20)[0]);
    after.push_back(scaled({5.5}, 1e-20)[0]); // second well: 0.5 / 5.5
    l.simulator_.p.current = after;
    const auto m = l.evaluateFractureCouplingChange(before);
    require(std::abs(m.max_ctf_change_well - 0.5 / 5.5) < 1e-12, "the worst well sets the metric");
    require(m.ctf_changed, "the small well's change is seen beside a large well");

    Loop q = loop(1e-2);
    q.simulator_.p.current = scaled({1.0, 1.001}, 1.0);
    require(!q.evaluateFractureCouplingChange(scaled({1.0, 1.0}, 1.0)).ctf_changed,
            "a 0.05 % change is below a 1 % threshold");
    std::cout << "threshold and worst well: ok\n";
}

// A changed connection cell is a structure change whatever the CTFs do.
void structureChange()
{
    Loop l = loop(1e-3);
    Perfs after = scaled({1.0, 2.0}, 1.0);
    after[0][1].cell = 7;
    l.simulator_.p.current = after;
    const auto m = l.evaluateFractureCouplingChange(scaled({1.0, 2.0}, 1.0));
    require(m.structure_changed && m.changed(), "a moved connection is a change");
    require(m.composite_norm >= 1.0, "a structure change sets the composite norm to at least 1");
    std::cout << "structure change: ok\n";
}
} // namespace

int main()
{
    try {
        scaleFree();
        absoluteDefaultMissesTinySeed();
        thresholdAndWorstWell();
        structureChange();
    } catch (const std::exception& e) {
        std::cerr << "FAILED: " << e.what() << "\n";
        return 1;
    }
    std::cout << "passed\n";
    return 0;
}
