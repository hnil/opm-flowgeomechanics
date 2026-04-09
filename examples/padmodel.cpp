/*
  Copyright 2018 Statoil ASA.

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

#include <opm/common/utility/FileSystem.hpp>

#include <opm/input/eclipse/EclipseState/Grid/SatfuncPropertyInitializers.hpp>
#include <opm/input/eclipse/EclipseState/IOConfig/IOConfig.hpp>
#include <opm/input/eclipse/EclipseState/InitConfig/InitConfig.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Deck/DeckSection.hpp>

#include <opm/input/eclipse/Parser/ErrorGuard.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>

#include <opm/input/eclipse/Parser/InputErrorAction.hpp>

#include <opm/input/eclipse/Parser/Parser.hpp>

#include <opm/input/eclipse/Parser/ParserKeywords/A.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/B.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/C.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/D.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/E.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/F.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/G.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/I.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/M.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/N.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/O.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/P.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/S.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/Z.hpp>

#include <opm/simulators/linalg/PropertyTree.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include <getopt.h>

namespace fs = std::filesystem;

namespace
{

struct GridSize
{
    int nx{};
    int ny{};
    int nz{};
};

struct ExtendParam
{
    int nz_upper{};
    double top_upper{};
    int nz_lower{};
    double bottom_lower{};
    double upper_poro{};
    double min_dist{};
    bool no_gap{};
    bool monotonic_zcorn{};
    bool vert_coord{};
    int upper_equilnum{};
    int nz_new{};
};

class PropsSectionDefaultValue
{
public:
    explicit PropsSectionDefaultValue(const Opm::Deck& deck);

    double operator()(const std::string& keyword) const;

private:
    std::unordered_map<std::string, double> epsDfltValue_{};

    static std::regex epsKwRegEx_;
    static std::unordered_map<std::string, double> fallbackDefaultValue_;
};

std::regex PropsSectionDefaultValue::epsKwRegEx_ { R"(I?([SKP].+)([XYZ]-?)?)" };

std::unordered_map<std::string, double>
PropsSectionDefaultValue::fallbackDefaultValue_ {
    std::pair { std::string { "SWATINIT" }, 1.0 },
};

PropsSectionDefaultValue::PropsSectionDefaultValue(const Opm::Deck& deck)
{
    using namespace std::string_literals;

    const auto tables = Opm::TableManager { deck };
    const auto rspec  = Opm::Runspec { deck };

    const auto rtep = Opm::satfunc::getRawTableEndpoints
        (tables, rspec.phases(),
         rspec.saturationFunctionControls()
         .minimumRelpermMobilityThreshold());

    const auto rfunc = Opm::satfunc::getRawFunctionValues
        (tables, rspec.phases(), rtep);

    auto cvrtPress = [&usys = deck.getActiveUnitSystem()](const double p)
    {
        return usys.from_si(Opm::UnitSystem::measure::pressure, p);
    };

    // Default values from region 1.  Both for drainage and imbibition.
    this->epsDfltValue_.insert({
        // Horizontally scaled end-points for gas.
        std::pair { "SGL"s  , rtep.connate.gas.front()  },
        std::pair { "SGLPC"s, rtep.connate.gas.front()  },
        std::pair { "SGCR"s , rtep.critical.gas.front() },
        std::pair { "SGU"s  , rtep.maximum.gas.front()  },

        // ------------------------------------------------------------
        // Horizontally scaled end-points for oil.
        std::pair { "SOGCR"s, rtep.critical.oil_in_gas  .front() },
        std::pair { "SOWCR"s, rtep.critical.oil_in_water.front() },

        // ------------------------------------------------------------
        // Horizontally scaled end-points for water.
        std::pair { "SWL"s  , rtep.connate.water.front()  },
        std::pair { "SWLPC"s, rtep.connate.water.front()  },
        std::pair { "SWCR"s , rtep.critical.water.front() },
        std::pair { "SWU"s  , rtep.maximum.water.front()  },

        // ============================================================

        // Vertically scaled end-points for gas.
        std::pair { "KRGR"s, rfunc.krg.r.front()           },
        std::pair { "KRG"s , rfunc.krg.max.front()         },
        std::pair { "PCG"s , cvrtPress(rfunc.pc.g.front()) },

        // ------------------------------------------------------------
        // Vertically scaled end-points for oil.
        std::pair { "KRORG"s, rfunc.kro.rg.front()  },
        std::pair { "KRORW"s, rfunc.kro.rw.front()  },
        std::pair { "KRO"s  , rfunc.kro.max.front() },

        // ------------------------------------------------------------
        // Vertically scaled end-points for water.
        std::pair { "KRWR"s, rfunc.krw.r.front()           },
        std::pair { "KRW"s , rfunc.krw.max.front()         },
        std::pair { "PCW"s , cvrtPress(rfunc.pc.w.front()) },
    });
}

double PropsSectionDefaultValue::operator()(const std::string& keyword) const
{
    if (auto m = std::smatch{}; std::regex_match(keyword, m, epsKwRegEx_)) {
        const auto valuePos = this->epsDfltValue_.find(m[1]);

        if (valuePos != this->epsDfltValue_.end()) {
            return valuePos->second;
        }
    }

    const auto valuePos = fallbackDefaultValue_.find(keyword);
    return (valuePos == fallbackDefaultValue_.end())
        ? 0.0 : valuePos->second;
}

std::string defaultJson()
{
    return {
        R"({
    "nz_upper": 2,
    "top_upper": 1000,
    "nz_lower": 3,
    "bottom_lower": 3600,
    "upper_poro": 0.1,
    "min_dist": 10,
    "no_gap": true,
    "vert_coord": true,
    "monotonic_zcorn": true
})" };
}

void populateExtendParameters(const Opm::PropertyTree& fracture_param,
                              ExtendParam&             extend_param)
{
    fracture_param.write_json(std::cout, true);

    extend_param.nz_upper = fracture_param.get<int>("nz_upper");
    extend_param.top_upper = fracture_param.get<double>("top_upper");
    extend_param.nz_lower = fracture_param.get<int>("nz_lower");
    extend_param.bottom_lower = fracture_param.get<double>("bottom_lower");
    extend_param.upper_poro = fracture_param.get<double>("upper_poro");
    extend_param.min_dist = fracture_param.get<double>("min_dist");
    extend_param.no_gap = fracture_param.get<bool>("no_gap");
    extend_param.monotonic_zcorn = fracture_param.get<bool>("monotonic_zcorn");
    extend_param.vert_coord = fracture_param.get<bool>("vert_coord");
}

void populateFallback(ExtendParam& extend_param)
{
    extend_param.nz_upper = 2;
    extend_param.top_upper = 1000;
    extend_param.nz_lower = 3;
    extend_param.bottom_lower = 3600;
    extend_param.upper_poro = 0.1;
    extend_param.min_dist = 10;
    extend_param.no_gap = true;
    extend_param.vert_coord = true;
    extend_param.monotonic_zcorn = true; // change results on norne?
}

void populateNoFracture(ExtendParam& extend_param)
{
    const auto dflt_json_fname = fs::temp_directory_path() /
        Opm::unique_path("extend-dflt-%%%%.json");

    std::ofstream dflt_json { dflt_json_fname };
    if (! dflt_json) {
        populateFallback(extend_param);
    }
    else {
        dflt_json << defaultJson() << '\n';
        dflt_json.close();

        const auto fracture_param = Opm::PropertyTree {
            dflt_json_fname.generic_string()
        };

        populateExtendParameters(fracture_param, extend_param);

        fs::remove(dflt_json_fname);
    }
}

void populateFromFile(std::string_view filename,
                      ExtendParam& extend_param)
{
    const auto fracture_param = Opm::PropertyTree {
        std::string { filename }
    };

    std::cout << "Fracture parameters from file '" << filename << "':\n";

    // set seed values
    populateExtendParameters(fracture_param, extend_param);
}

ExtendParam
getExtendParam(std::string_view filename)
{
    auto extend_param = ExtendParam {};

    try {
        populateFromFile(filename, extend_param);
    }
    catch (const std::exception&) {
        std::cerr << "Error reading fracture parameter file: '"
                  << filename << "'\n"
                  << "Defaulting to no fractures.\n";

        populateNoFracture(extend_param);
    }

    return extend_param;
}

GridSize
getDimens(const Opm::Deck& deck)
{
    GridSize grid_size{};

    if (deck.hasKeyword<Opm::ParserKeywords::DIMENS>()) {
        grid_size.nx = deck[Opm::ParserKeywords::DIMENS::keywordName]
            .back().getRecord(0).getItem("NX").get<int>(0);

        grid_size.ny = deck[Opm::ParserKeywords::DIMENS::keywordName]
            .back().getRecord(0).getItem("NY").get<int>(0);

        grid_size.nz = deck[Opm::ParserKeywords::DIMENS::keywordName]
            .back().getRecord(0).getItem("NZ").get<int>(0);
    }
    else {
        std::cerr << "No DIMENS keyword found in the deck\n";
    }

    return grid_size;
}

void
extendDimens(const GridSize& grid_size,
             ExtendParam& extend_param,
             Opm::Deck& deck)
{
    const int nz = grid_size.nz;

    if (deck.hasKeyword<Opm::ParserKeywords::DIMENS>()) {

        extend_param.nz_new = nz + extend_param.nz_upper + extend_param.nz_lower;
        {
            auto& NZit = const_cast<Opm::DeckItem&>
                (deck[Opm::ParserKeywords::DIMENS::keywordName].back().getRecord(0).getItem("NZ"));

            std::vector<int>& data = NZit.getData<int>();
            data[0] = extend_param.nz_new;
        }

        {
            auto& eqldims = const_cast<Opm::DeckItem&>
                (deck[Opm::ParserKeywords::EQLDIMS::keywordName].back().getRecord(0).getItem("NTEQUL"));

            std::vector<int>& data = eqldims.getData<int>();

            std::vector<Opm::value::status>& value_status
                = const_cast<std::vector<Opm::value::status>&>(eqldims.getValueStatus());

            if (value_status[0] == Opm::value::status::deck_value) {
                data[0] += 1;
            }
            else {
                value_status[0] = Opm::value::status::deck_value;
                data[0] = 2;
            }

            const int upper_equilnum = data[0];
            extend_param.upper_equilnum = upper_equilnum;
        }
    }
    else {
        std::cerr << "No DIMENS keyword found in the deck" << std::endl;
    }

    if (deck.hasKeyword<Opm::ParserKeywords::SPECGRID>()) {
        auto& NZit = const_cast<Opm::DeckItem&>
            (deck[Opm::ParserKeywords::SPECGRID::keywordName].back().getRecord(0).getItem("NZ"));

        std::vector<int>& data = NZit.getData<int>();

        data[0] = extend_param.nz_new;
    }
}

template <typename T>
void extendPropertyArray(const GridSize&      grid_size,
                         const ExtendParam&   extend_param,
                         const Opm::type_tag  type,
                         const T              default_value,
                         const Opm::DeckItem& property)
{
    if (property.getType() != type) {
        return;
    }

    const auto input_nc = static_cast<std::size_t>(grid_size.nx)
        * static_cast<std::size_t>(grid_size.ny)
        * static_cast<std::size_t>(grid_size.nz);

    auto& value = const_cast<Opm::DeckItem&>(property).template getData<T>();
    if (std::size(value) != input_nc) {
        return;
    }

    auto& status = const_cast<std::vector<Opm::value::status>&>
        (property.getValueStatus());

    const auto output_nc = static_cast<std::size_t>(grid_size.nx)
        * static_cast<std::size_t>(grid_size.ny)
        * static_cast<std::size_t>(extend_param.nz_upper +
                                   grid_size.nz          +
                                   extend_param.nz_lower);

    const auto value_offset = static_cast<std::size_t>(grid_size.nx)
        * static_cast<std::size_t>(grid_size.ny)
        * static_cast<std::size_t>(extend_param.nz_upper);

    auto remapped_value  = std::vector<T>(output_nc, default_value);
    auto remapped_status = std::vector<Opm::value::status>
        (output_nc, Opm::value::status::deck_value);

    std::copy(value .begin(), value .end(), remapped_value .begin() + value_offset);
    std::copy(status.begin(), status.end(), remapped_status.begin() + value_offset);

    value .swap(remapped_value);
    status.swap(remapped_status);
}

template <typename T>
void
extendGridSection(Opm::GRIDSection& gridsec,
                  const GridSize& grid_size,
                  const ExtendParam& extend_param,
                  const Opm::type_tag type)
{
    for (const auto& keyword : gridsec) {
        if ((keyword.size() != 1) || !keyword.isDataKeyword()) {
            continue;
        }

        T default_value {};
        if (keyword.name() == "PORO") {
            default_value = extend_param.upper_poro;
        }
        else if (keyword.name() == "NTG") {
            default_value = T{1};
        }

        extendPropertyArray(grid_size, extend_param,
                            type, default_value,
                            keyword.getRecord(0).getItem(0));
    }
}

void extendBCCon(Opm::GRIDSection& gridSect,
                 const ExtendParam& extend_param)
{
    // extend BCCON only downwards

    for (const auto& records : gridSect) {
        if (records.name() != "BCCON") {
            continue;
        }

        //std::cout << "Extending BCCON" << std::endl;

        for (const auto& record : records) {
            int& K1 = const_cast<int&>(record.getItem("K1").getData<int>()[0]);
            int& K2 = const_cast<int&>(record.getItem("K2").getData<int>()[0]);

            const auto direction = record.getItem("DIRECTION").getTrimmedString(0);
            //std::cout << direction << std::endl;

            if (direction == "Z+") {
                const int shift = extend_param.nz_upper + extend_param.nz_lower;
                K1 += shift;
                K2 += shift;
            }

            if (direction == "Z-") {
                assert(K1 == 1);
                assert(K2 == 1);
            }
        }
    }
}

void extendFaults(const ExtendParam& extend_param,
                  Opm::GRIDSection&  gridSect)
{
    for (const auto* flt : gridSect.getKeywordList<Opm::ParserKeywords::FAULTS>()) {
        for (const auto& face : *flt) {
            const_cast<int&>(face.getItem<Opm::ParserKeywords::FAULTS::IZ1>()
                             .getData<int>().front()) += extend_param.nz_upper;

            const_cast<int&>(face.getItem<Opm::ParserKeywords::FAULTS::IZ2>()
                             .getData<int>().front()) += extend_param.nz_upper;
        }
    }
}

void extendOperators(const ExtendParam& extend_param,
                     Opm::Deck&         deck)
{
    for (const auto& opKw : {
            Opm::ParserKeywords::ADD     ::keywordName,
            Opm::ParserKeywords::BOX     ::keywordName,
            Opm::ParserKeywords::EQUALS  ::keywordName,
            Opm::ParserKeywords::MAXVALUE::keywordName,
            Opm::ParserKeywords::MINVALUE::keywordName,
            Opm::ParserKeywords::MULTIPLY::keywordName,
            Opm::ParserKeywords::OPERATE ::keywordName,
        })
    {
        for (const auto& operation : deck.getKeywordList(opKw)) {
            for (const auto& record : *operation) {
                for (const auto* kname : { "K1", "K2", }) {
                    const auto& k = record.getItem(kname);
                    if (! k.defaultApplied(0)) {
                        const_cast<Opm::DeckItem&>(k)
                            .getData<int>().front()
                            += extend_param.nz_upper;
                    }
                }
            }
        }
    }
}

void extendSpecialGridProperties(Opm::GRIDSection& gridSect,
                                 const ExtendParam& extend_param)
{
    if (! gridSect.hasKeyword<Opm::ParserKeywords::ACTNUM>()) {
        // All cells active.
        return;
    }

    auto& actnum = const_cast<std::vector<int>&>
        (gridSect.get<Opm::ParserKeywords::ACTNUM>().back().getIntData());

    auto* const poro = gridSect.hasKeyword<Opm::ParserKeywords::PORO>()
        ? const_cast<std::vector<double>*>
        (&gridSect.get<Opm::ParserKeywords::PORO>().back().getRawDoubleData())
        : nullptr;

    auto* const permx = gridSect.hasKeyword<Opm::ParserKeywords::PERMX>()
        ? const_cast<std::vector<double>*>
        (&gridSect.get<Opm::ParserKeywords::PERMX>().back().getRawDoubleData())
        : nullptr;

    auto* const permy = gridSect.hasKeyword<Opm::ParserKeywords::PERMY>()
        ? const_cast<std::vector<double>*>
        (&gridSect.get<Opm::ParserKeywords::PERMY>().back().getRawDoubleData())
        : nullptr;

    auto* const permz = gridSect.hasKeyword<Opm::ParserKeywords::PERMZ>()
        ? const_cast<std::vector<double>*>
        (&gridSect.get<Opm::ParserKeywords::PERMZ>().back().getRawDoubleData())
        : nullptr;

    auto* const ntg = gridSect.hasKeyword<Opm::ParserKeywords::NTG>()
        ? const_cast<std::vector<double>*>
        (&gridSect.get<Opm::ParserKeywords::NTG>().back().getRawDoubleData())
        : nullptr;

    {
        auto& actnum_status = const_cast<std::vector<Opm::value::status>&>
            (gridSect.get<Opm::ParserKeywords::ACTNUM>().back().getValueStatus());

        actnum_status.assign(actnum.size(), Opm::value::status::deck_value);
    }

    for (auto i = 0*actnum.size(); i != actnum.size(); ++i) {
        if (actnum[i] == 1) {
            // Input cell is active.  Nothing to do.
            continue;
        }

        // If we get here, the input cell is inactive.  Reset to active and
        // set permeability to zero, NTG to one, while porosity gets the
        // 'upper_poro' runtime parameter.

        actnum[i] = 1;

        // maybe one should have used values from deck if resonable and
        // check if deck_value
        if (poro != nullptr) {
            (*poro)[i] = extend_param.upper_poro;
        }

        if (permx != nullptr) {
            (*permx)[i] = 0.0;
        }

        if (permy != nullptr) {
            (*permy)[i] = 0.0;
        }

        if (permz != nullptr) {
            (*permz)[i] = 0.0;
        }

        if (ntg != nullptr) {
            (*ntg)[i] = 1.0;
        }
    }
}

void
extendGRDECL(Opm::GRIDSection& gridsec,
             const GridSize& grid_size,
             const ExtendParam& extend_param)
{
    const int nx = grid_size.nx;
    const int ny = grid_size.ny;
    const int nz = grid_size.nz;
    const int nz_upper = extend_param.nz_upper;
    const int nz_lower = extend_param.nz_lower;
    const double top_upper = extend_param.top_upper;
    const double bottom_lower = extend_param.bottom_lower;
    const double min_dist = extend_param.min_dist;
    const double no_gap = extend_param.no_gap;
    const int nz_new = extend_param.nz_new;

    std::vector<double>& coord = const_cast<std::vector<double>&>
        (gridsec.get<Opm::ParserKeywords::COORD>().back().getRawDoubleData());

    std::vector<double>& zcorn = const_cast<std::vector<double>&>
        (gridsec.get<Opm::ParserKeywords::ZCORN>().back().getRawDoubleData());

    std::vector<Opm::value::status>& zcorn_status = const_cast<std::vector<Opm::value::status>&>
        (gridsec.get<Opm::ParserKeywords::ZCORN>().back().getValueStatus());

    // include all cells

    // make pilars straight
    if (extend_param.vert_coord) {
        const int ncoord = coord.size();
        const int nxny = ncoord / 6;

        assert((nx + 1) * (ny + 1) == nxny);
        assert(ncoord % 6 == 0);

        for (int i = 0; i < nxny; ++i) {
            coord[6 * i + 3] = coord[6 * i + 0];
            coord[6 * i + 4] = coord[6 * i + 1];
        }
    }

    if (extend_param.monotonic_zcorn) {
        auto zcorn_status_new = std::vector<Opm::value::status>((2 * nx) * (2 * ny) * (2 * nz_new),
                                                                Opm::value::status::deck_value);

        auto zcorn_new = std::vector<double>((2 * nx) * (2 * ny) * (2 * nz_new), 0.0);

        // extend, make monotonic, fill gaps zcorn
        for (int i = 0; i < 2 * nx; ++i) {
            for (int j = 0; j < 2 * ny; ++j) {
                double minz = 1e20;
                double maxz = -1e20;

                for (int k = 0; k < 2 * nz; ++k) {
                    const int index = i + j * (2 * nx) + k * (2 * nx) * (2 * ny);
                    const double z = zcorn[index];
                    if (z > 0.0) {
                        minz = std::min(minz, z);
                        maxz = std::max(maxz, z);
                    }
                }

                const double dz_upper = std::max(0.0, (minz - top_upper) / nz_upper);
                const double dz_lower = std::max(0.0, (bottom_lower - maxz) / nz_lower);

                std::vector<double> zcornvert(2 * nz_new);
                if (nz_upper > 0) {
                    zcornvert[0] = std::min(top_upper, minz);
                } else {
                    zcornvert[0] = minz;
                }

                double z_prev = zcornvert[0];
                for (int k = 1; k < 2 * nz_new; k++) {
                    const int k_old = k - 2 * nz_upper;
                    const int ind_old = i + j * (2 * nx) + k_old * (2 * nx) * (2 * ny);
                    if (k_old >= 0 && k_old < 2 * nz) {
                        double z = zcorn[ind_old];
                        double dz = z - z_prev;
                        if (dz < min_dist) {
                            dz = 0;
                        }

                        if (z == 0.0) {
                            dz = 0;
                        }

                        if (no_gap) {
                            // if we jup to new logical cell wee should not have gaps
                            if (k % 2 == 0) {
                                dz = 0;
                            }
                        }

                        zcornvert[k] = zcornvert[k - 1] + dz;
                        z_prev = zcornvert[k];
                    } else {
                        double dz_cell = dz_lower;
                        if (k < 2 * nz_upper) {
                            dz_cell = dz_upper;
                        }

                        if (k % 2 == 1) {
                            zcornvert[k] = zcornvert[k - 1] + dz_cell;
                        } else {
                            zcornvert[k] = zcornvert[k - 1];
                        }

                        z_prev = zcornvert[k];
                    }
                }

                for (int k = 0; k < 2 * nz_new; k++) {
                    int ind = i + j * (2 * nx) + k * (2 * nx) * (2 * ny);
                    zcorn_new[ind] = zcornvert[k];
                }
            }
        }

        zcorn_status = zcorn_status_new;
        zcorn = zcorn_new;
    }
}

void
extendProps(const GridSize&    grid_size,
            const ExtendParam& extend_param,
            Opm::Deck&         deck)
{
    const auto defaultValue = PropsSectionDefaultValue { deck };

    for (const auto& keyword : Opm::PROPSSection { deck }) {
        if ((keyword.size() != 1) || !keyword.isDataKeyword()) {
            continue;
        }

        extendPropertyArray(grid_size, extend_param,
                            Opm::type_tag::fdouble,
                            defaultValue(keyword.name()),
                            keyword.getRecord(0).getItem(0));
    }
}

void
extendRegions(Opm::REGIONSSection& regions,
              const GridSize& grid_size,
              const ExtendParam& extend_param,
              const std::vector<int>& actnum_old)
{
    const auto input_nc = static_cast<std::size_t>(grid_size.nx)
        * static_cast<std::size_t>(grid_size.ny)
        * static_cast<std::size_t>(grid_size.nz);

    const auto output_nc = static_cast<std::size_t>(grid_size.nx)
        * static_cast<std::size_t>(grid_size.ny)
        * static_cast<std::size_t>(extend_param.nz_upper +
                                   grid_size.nz          +
                                   extend_param.nz_lower);

    const auto value_offset = static_cast<std::size_t>(grid_size.nx)
        * static_cast<std::size_t>(grid_size.ny)
        * static_cast<std::size_t>(extend_param.nz_upper);

    if (regions.hasKeyword<Opm::ParserKeywords::EQLNUM>()) {
        std::cout << "EQLNUM keyword found in the deck" << std::endl;
    }
    else {
        std::cout << "No EQLNUM keyword found in the deck" << std::endl;
    }

    for (const auto& keyword : regions) {
        if ((keyword.size() != 1) || !keyword.isDataKeyword()) {
            continue;
        }

        const auto& record = keyword.getRecord(0);
        if (record.size() != 1) {
            continue;
        }

        const auto& item = record.getItem(0);
        if ((item.getType() != Opm::type_tag::integer) ||
            (item.data_size() != input_nc))
        {
            continue;
        }

        auto& value = const_cast<Opm::DeckItem&>(item).getData<int>();

        const auto default_value =
            [&value, &extend_param, is_eqlnum = keyword.name() == "EQLNUM"]()
        {
            if (is_eqlnum) { return extend_param.upper_equilnum; }

            return std::max(1, *std::min_element(value.begin(), value.end()));
        }();

        auto& status = const_cast<std::vector<Opm::value::status>&>
            (item.getValueStatus());

        auto remapped_value  = std::vector<int>(output_nc, default_value);
        auto remapped_status = std::vector<Opm::value::status>
            (output_nc, Opm::value::status::deck_value);

        for (auto i = 0*input_nc; i != input_nc; ++i) {
            if (actnum_old[i] == 1) {
                remapped_value [value_offset + i] = value [i];
                remapped_status[value_offset + i] = status[i];
            }
        }

        value .swap(remapped_value);
        status.swap(remapped_status);
    }
}

void extendEQUIL(Opm::DeckKeyword& equil)
{
    using Kw = Opm::ParserKeywords::EQUIL;

    // copy first record
    assert(! equil.empty());

    // Equilibration data for extended region (above/below formation).
    // Based on first equilibration region.
    auto record = equil.getRecord(0);

    double& goc = const_cast<double&>(record.getItem<Kw::GOC>().getData<double>()[0]);
    double& woc = const_cast<double&>(record.getItem<Kw::OWC>().getData<double>()[0]);

    goc = woc = 0.0;

    equil.addRecord(std::move(record));
}

void extendDepthTable(const ExtendParam& extend_param,
                      Opm::DeckKeyword&  propertyVD)
{
    assert(! propertyVD.empty());

    auto newPropertyVD = propertyVD.emptyStructuralCopy();

    for (const auto& table : propertyVD) {
        auto i = std::vector { table.getItem(0).emptyStructuralCopy() };

        const auto& values = table.getItem(0).getData<double>();

        // Constant extrapolation above formation.
        if (extend_param.top_upper < values.front()) {
            i.front().push_back(extend_param.top_upper);
            i.front().push_back(values[0 + 1]);
        }

        // Input table unchanged within formation's depth range.
        for (const auto& value : values) {
            i.front().push_back(value);
        }

        // Constant extrapolation below formation.
        if (extend_param.bottom_lower > values[values.size() - 2]) {
            i.front().push_back(extend_param.bottom_lower);
            i.front().push_back(values.back());
        }

        newPropertyVD.addRecord(Opm::DeckRecord { std::move(i) });
    }

    // Add table corresponding to "extended" region (propertyVD.size() + 1).
    {
        auto first = newPropertyVD.getRecord(0);
        newPropertyVD.addRecord(std::move(first));
    }

    propertyVD = newPropertyVD;
}

void
extendSolution(const ExtendParam&    extend_param,
               Opm::SOLUTIONSection& solution)
{
    for (const auto& keyword : solution) {
        if (keyword.name() == "EQUIL") {
            extendEQUIL(const_cast<Opm::DeckKeyword&>(keyword));
        }
        else if ((keyword.name() == "RSVD") ||
                 (keyword.name() == "RVVD") ||
                 (keyword.name() == "RTEMPVD") ||
                 (keyword.name() == "PBVD") ||
                 (keyword.name() == "PDVD"))
        {
            extendDepthTable(extend_param,
                             const_cast<Opm::DeckKeyword&>(keyword));
        }
    }
}

void
extendSchedule(Opm::SCHEDULESection& schedule,
               const ExtendParam& extend_param)
{
    // fix kz for schedule i.e. COMPDAT
    const int nz_upper = extend_param.nz_upper;

    for (const auto& records : schedule) {
        if ((records.name() != "COMPDAT") && (records.name() != "WSEED")
            && (records.name() != "COMPSEGS")) {
            continue;
        }

        //std::cout << "Extending schedule record: " << records.name() << std::endl;

        if (records.name() == "WSEED") {
            for (const auto& record : records) {
                int& K = const_cast<int&>(record.getItem("K").getData<int>()[0]);
                K += nz_upper;
            }
        } else if (records.name() == "COMPSEGS") {
            auto is_first = true;

            for (const auto& record : records) {
                if (is_first) {
                    is_first = false;
                    continue;
                }

                int& K = const_cast<int&>(record.getItem("K").getData<int>()[0]);
                K += nz_upper;
            }
        } else {
            for (const auto& record : records) {
                int& K1 = const_cast<int&>(record.getItem("K1").getData<int>()[0]);
                int& K2 = const_cast<int&>(record.getItem("K2").getData<int>()[0]);
                K1 += nz_upper;
                K2 += nz_upper;
            }
        }
    }
}

Opm::Deck
manipulate_deck(std::string_view deck_file,
                std::string_view param_file,
                std::ostream& os)
{
    auto parseContext = Opm::ParseContext(Opm::InputErrorAction::WARN);
    Opm::ErrorGuard errors;

    // Pad vertically
    auto deck = Opm::Parser{}.parseFile(std::string {deck_file}, parseContext, errors);
    auto extend_param = getExtendParam(param_file);

    std::cout << "Reading deck file: " << deck_file << std::endl;



    // write out json file

    std::cout << "Extend parameters:\n"
              << "nz_upper: "          << extend_param.nz_upper        << '\n'
              << "top_upper: "         << extend_param.top_upper       << '\n'
              << "nz_lower: "          << extend_param.nz_lower        << '\n'
              << "bottom_lower: "      << extend_param.bottom_lower    << '\n'
              << "upper_poro: "        << extend_param.upper_poro      << '\n'
              << "min_dist: "          << extend_param.min_dist        << '\n'
              << "no_gap: "            << extend_param.no_gap          << '\n'
              << "monotonic_zcorn: "   << extend_param.monotonic_zcorn << '\n'
              << "vert_coord: "        << extend_param.vert_coord
              << std::endl;

    const auto grid_size = getDimens(deck);

    // int nz_new,upper_equilnum;
    auto gridsec = Opm::GRIDSection {deck}; // fixing CARTDIMS? + ZCORN + double valued arrays.

    {
        extendDimens(grid_size, extend_param, deck);

        extendBCCon(gridsec, extend_param);

        extendGridSection<double>(gridsec, grid_size, extend_param, Opm::type_tag::fdouble);
        extendGridSection<int>(gridsec, grid_size, extend_param, Opm::type_tag::integer);

        extendSpecialGridProperties(gridsec, extend_param);

        extendFaults(extend_param, gridsec);

        extendGRDECL(gridsec, grid_size, extend_param);
    }

    if (Opm::DeckSection::hasPROPS(deck)) {
        // Yes, we intentionally operate on the entire Deck here rather than
        // just the PROPSSection.  We need the saturation function tables
        // along with the run's active phases to properly remap scaled
        // end-points for the saturation functions.
        extendProps(grid_size, extend_param, deck);
    }

    if (Opm::DeckSection::hasREGIONS(deck)) {
        const auto actnum_old = gridsec.hasKeyword<Opm::ParserKeywords::ACTNUM>()
            ? gridsec.get<Opm::ParserKeywords::ACTNUM>().back().getIntData()
            : std::vector<int>(grid_size.nx * grid_size.ny * grid_size.nz, 1);

        // NB NB need to always add EQUILNUM
        auto regions = Opm::REGIONSSection {deck}; // extend and add for equilnum

        extendRegions(regions, grid_size, extend_param, actnum_old);
    }

    if (Opm::DeckSection::hasSOLUTION(deck)) {
        auto solution = Opm::SOLUTIONSection {deck}; // extend and add for equilnum

        extendSolution(extend_param, solution);
    }

    if (Opm::DeckSection::hasSCHEDULE(deck)) {
        auto schedule = Opm::SCHEDULESection {deck}; // fixing COMPDAT

        extendSchedule(schedule, extend_param);
    }

    extendOperators(extend_param, deck);

    if (deck.hasKeyword<Opm::ParserKeywords::EQLDIMS>()) {
        std::cout << "EQLDIMS keyword found in the deck" << std::endl;
    }
    else {
        std::cout << "No EQLDIMS keyword found in the deck" << std::endl;
    }

    os << deck;

    return deck;
}

void print_help()
{
    std::cerr << R"(
The padmodel program will load a deck, resolve all include files and then
print it out again on stdout. All comments will be stripped and the value
types will be validated.

By passing the option -o you can redirect the output to a file or a
directory.

It will also pad the grid vertically according to parameters specified in a
JSON file (default "fracture_param.json").

Print on stdout:

    padmodel -p fracture_param.json /path/to/case/CASE.DATA

Print MY_CASE.DATA in /tmp:

    padmodel -p fracture_param.json -o /tmp /path/to/MY_CASE.DATA

Print NEW_CASE in cwd:

    padmodel -o NEW_CASE.DATA -p fracture_param.json path/to/MY_CASE.DATA

As an alternative to the -o option you can use -c. This flag is nominally
equivalent to "-o -", but restart and import files referred to in the deck
are also copied to the output stream. The -o and -c options are mutually
exclusive.
)";
}

void print_example_json()
{
    const auto dflt_json_fname = fs::temp_directory_path() /
        Opm::unique_path("extend-dflt-%%%%.json");

    std::ofstream dflt_json { dflt_json_fname };

    std::cerr << "\n============================================================"
        "\n\nExample JSON parameters (tool's built-in defaults):\n\n";

    if (! dflt_json) {
        std::cerr << defaultJson() << '\n';
    }
    else {
        dflt_json << defaultJson() << '\n';
        dflt_json.close();

        const auto fracture_param = Opm::PropertyTree {
            dflt_json_fname.generic_string()
        };

        fs::remove(dflt_json_fname);

        fracture_param.write_json(std::cerr, true);
    }

    std::cerr << "\n============================================================\n\n";
}

void
copy_file(const fs::path& source_dir, fs::path fname, const fs::path& target_dir)
{
    if (fname.is_absolute()) {
        const auto prefix_len = fs::canonical(source_dir).string().size();

        fname = fs::canonical(fname);
        fname = fs::path(fname.string().substr(prefix_len + 1));
    }

    const auto source_file = source_dir / fname;
    const auto target_file = target_dir / fname;
    {
        const auto& parent_path = target_file.parent_path();
        if (!parent_path.empty() && !fs::is_directory(parent_path)) {
            fs::create_directories(parent_path);
        }
    }

    fs::copy_file(source_file, target_file, fs::copy_options::overwrite_existing);

    std::cerr << "Copying file " << source_file.string() << " -> " << target_file.string() << std::endl;
}

} // Anonymous namespace

int main(int argc, char** argv)
{
    bool stdout_output = true;
    bool copy_binary = false;
    bool help = argc < 2;
    auto coutput_arg = std::string_view{};
    auto param_file_arg = std::string_view{"fracture_param.json"};

    while (true) {
        const auto c = getopt(argc, argv, "c:ho:p:");
        if (c == -1) {
            break;
        }

        switch (c) {
        case 'c':
            stdout_output = false;
            copy_binary = true;
            coutput_arg = optarg;
            break;

        case 'h':
            help = true;
            break;

        case 'o':
            stdout_output = false;
            coutput_arg = optarg;
            break;

        case 'p':
            param_file_arg = optarg;
            break;
        }
    }

    const auto arg_offset = optind;
    if (help || (arg_offset >= argc)) {
        print_help();

        if (help) {
            print_example_json();
            return EXIT_SUCCESS;
        }

        return EXIT_FAILURE;
    }

    if (stdout_output) {
        manipulate_deck(argv[arg_offset], param_file_arg, std::cout);
    }
    else {
        const auto input_arg = fs::path { argv[arg_offset] };

        auto output_arg = fs::path { coutput_arg };
        auto output_dir = fs::path { coutput_arg };

        std::ofstream os{};
        if (fs::is_directory(output_arg)) {
            os.open(output_arg / input_arg.filename());
        }
        else {
            os.open(output_arg);
            output_dir = output_arg.parent_path();
        }

        const auto deck = manipulate_deck(argv[arg_offset], param_file_arg, os);

        if (copy_binary) {
            if (const auto init_config = Opm::InitConfig { deck, Opm::Runspec {deck}.phases() };
                init_config.restartRequested())
            {
                const auto io_config = Opm::IOConfig { deck };

                const auto restart_file = fs::path {
                    io_config.getRestartFileName(init_config.getRestartRootName(),
                                                 init_config.getRestartStep(), false)
                };

                copy_file(input_arg.parent_path(), restart_file, output_dir);
            }

            using IMPORT = Opm::ParserKeywords::IMPORT;
            for (const auto& import_keyword : deck.get<IMPORT>()) {
                const auto& fname = import_keyword.getRecord(0)
                    .getItem<IMPORT::FILE>().get<std::string>(0);

                copy_file(input_arg.parent_path(), fname, output_dir);
            }

            using PYACTION = Opm::ParserKeywords::PYACTION;
            for (const auto& pyaction_keyword : deck.get<PYACTION>()) {
                const auto& fname = pyaction_keyword.getRecord(1)
                    .getItem<PYACTION::FILENAME>().get<std::string>(0);

                copy_file(input_arg.parent_path(), fname, output_dir);
            }

            using GDFILE = Opm::ParserKeywords::GDFILE;
            if (deck.hasKeyword<GDFILE>()) {
                const auto& gdfile_keyword = deck.get<GDFILE>().back();
                const auto& fname = gdfile_keyword.getRecord(0)
                    .getItem<GDFILE::filename>().get<std::string>(0);

                copy_file(input_arg.parent_path(), fname, output_dir);
            }
        }
    }
}
