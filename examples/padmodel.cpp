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

#include <opm/input/eclipse/EclipseState/IOConfig/IOConfig.hpp>
#include <opm/input/eclipse/EclipseState/InitConfig/InitConfig.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Deck/DeckSection.hpp>

#include <opm/input/eclipse/Parser/ErrorGuard.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>

#include <opm/input/eclipse/Parser/InputErrorAction.hpp>

#include <opm/input/eclipse/Parser/Parser.hpp>

#include <opm/input/eclipse/Parser/ParserKeywords/A.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/C.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/D.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/E.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/G.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/I.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/N.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/P.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/S.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/Z.hpp>

#include <opm/simulators/linalg/PropertyTree.hpp>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>
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

void populateNoFracture(ExtendParam& extend_param)
{
    extend_param.nz_upper = 2;
    extend_param.top_upper = 1000;
    extend_param.no_gap = true;
    extend_param.vert_coord = true;
    extend_param.monotonic_zcorn = true; // change results on norne?
    extend_param.min_dist = 10;
    extend_param.upper_poro = 0.1;
    extend_param.nz_lower = 3;
    extend_param.bottom_lower = 3600;

    std::string default_json(R"(
{
    "nz_upper": 2,
    "top_upper": 1000,
    "nz_lower": 3,
    "bottom_lower": 3600,
    "upper_poro": 0.1,
    "min_dist": 10,
    "no_gap": true,
    "monotonic_zcorn": true,
    "vert_coord": true
}
)");
    std::cout << "Using default fracture parameters:\n"
              << default_json << std::endl;
}

void populateFromFile(std::string_view filename,
                      ExtendParam& extend_param)
{
    const auto fracture_param = Opm::PropertyTree {
        std::string { filename }
    };
    std::cout << "Fracture parameters from file '" << filename << "':\n";
    fracture_param.write_json(std::cout,true);
    // set seed values
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

ExtendParam
getExtendParam(std::string_view filename)
{
    auto extend_param = ExtendParam {};

    try {
        populateFromFile(filename, extend_param);
    }
    catch (const std::exception& e) {
        std::cout << "Error reading fracture parameter file: '" << filename << "'\n"
                  << "Defaulting to no fractures." << std::endl;
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
void
extendGridSection(Opm::GRIDSection& gridsec,
                  const GridSize& grid_size,
                  const ExtendParam& extend_param,
                  const Opm::type_tag type)
{
    const int nx = grid_size.nx;
    const int ny = grid_size.ny;
    const int nz = grid_size.nz;
    const int nz_upper = extend_param.nz_upper;
    const int nz_new = extend_param.nz_new;
    const int nc_new = nx * ny * nz_new;
    const int nc = nx * ny * nz;

    for (const auto& records : gridsec) {
        if (records.size() != 1) {
            continue;
        }

        const auto& record = records.getRecord(0);
        if (record.size() != 1) {
            continue;
        }

        const auto& item = record.getItem(0);
        if (item.getType() != type) {
            continue;
        }

        auto& val = const_cast<std::vector<T>&>(item.getData<T>());
        if (static_cast<int>(val.size()) != nc) {
            continue;
        }

        T default_value = 0;
        if (records.name() == "PORO") {
            default_value = extend_param.upper_poro;
        }
        else if (records.name() == "NTG") {
            default_value = T{1};
        }

        //std::cout << "Extend data_size: " << records.name() << std::endl;

        auto& val_status = const_cast<std::vector<Opm::value::status>&>(item.getValueStatus());

        std::vector<T> val_new(nc_new);
        std::vector<Opm::value::status> val_status_new(nc_new, Opm::value::status::deck_value);

        const int nz_bottom = nz_upper + nz;
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int k = 0; k < nz_new; ++k) {
                    const int ind_new = i + nx * j + nx * ny * k;
                    const int ind_old = i + nx * j + nx * ny * (k - nz_upper);

                    if ((k >= nz_upper) && (k < nz_bottom)) {
                        val_new[ind_new] = val[ind_old];
                        val_status_new[ind_new] = val_status[ind_old];
                    } else {
                        val_status_new[ind_new] = Opm::value::status::deck_value;
                        val_new[ind_new] = default_value;
                    }
                }
            }
        }

        val_status = val_status_new;
        val = val_new;
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

void extendSpecialGridProperties(Opm::GRIDSection& gridSect,
                                 const ExtendParam& extend_param)
{
    auto& actnum = const_cast<std::vector<int>&>
        (gridSect.get<Opm::ParserKeywords::ACTNUM>().back().getIntData());

    auto* poro = gridSect.hasKeyword<Opm::ParserKeywords::PORO>()
        ? const_cast<std::vector<double>*>
        (&gridSect.get<Opm::ParserKeywords::PORO>().back().getRawDoubleData())
        : nullptr;

    auto* permx = gridSect.hasKeyword<Opm::ParserKeywords::PERMX>()
        ? const_cast<std::vector<double>*>
        (&gridSect.get<Opm::ParserKeywords::PERMX>().back().getRawDoubleData())
        : nullptr;

    auto* permy = gridSect.hasKeyword<Opm::ParserKeywords::PERMY>()
        ? const_cast<std::vector<double>*>
        (&gridSect.get<Opm::ParserKeywords::PERMY>().back().getRawDoubleData())
        : nullptr;

    auto* permz = gridSect.hasKeyword<Opm::ParserKeywords::PERMZ>()
        ? const_cast<std::vector<double>*>
        (&gridSect.get<Opm::ParserKeywords::PERMZ>().back().getRawDoubleData())
        : nullptr;

    auto* ntg = gridSect.hasKeyword<Opm::ParserKeywords::NTG>()
        ? const_cast<std::vector<double>*>
        (&gridSect.get<Opm::ParserKeywords::NTG>().back().getRawDoubleData())
        : nullptr;

    {
        auto& actnum_status = const_cast<std::vector<Opm::value::status>&>
            (gridSect.get<Opm::ParserKeywords::ACTNUM>().back().getValueStatus());

        actnum_status.assign(actnum.size(), Opm::value::status::deck_value);
    }

    for (auto i = 0*actnum.size(); i < actnum.size(); ++i) {
        if (actnum[i] == 1) {
            continue;
        }

        actnum[i] = 1;

        // maybe one sould have used values from deck if resonable and
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
extendRegions(Opm::REGIONSSection& regions,
              const GridSize& grid_size,
              const ExtendParam& extend_param,
              const std::vector<int>& actnum_old)
{
    const int nx = grid_size.nx;
    const int ny = grid_size.ny;
    const int nz = grid_size.nz;
    const int nz_upper = extend_param.nz_upper;
    const int nz_new = extend_param.nz_new;
    const int nc_new = nx * ny * nz_new;
    const int nc = nx * ny * nz;
    const int upper_equilnum = extend_param.upper_equilnum;

    if (regions.hasKeyword<Opm::ParserKeywords::EQLNUM>()) {
        std::cout << "EQLNUM keyword found in the deck" << std::endl;
    }
    else {
        std::cout << "No EQLNUM keyword found in the deck" << std::endl;
    }

    for (const auto& records : regions) {
        if (records.size() != 1) {
            continue;
        }

        const auto& record = records.getRecord(0);
        if (record.size() != 1) {
            continue;
        }

        const auto& item = record.getItem(0);
        if (item.getType() != Opm::type_tag::integer) {
            continue;
        }

        std::vector<int>& val = const_cast<std::vector<int>&>(item.getData<int>());
        if (static_cast<int>(val.size()) != nc) {
            continue;
        }

        const int min_val = *std::min_element(val.begin(), val.end());
        int exteded_val = min_val;
        if (records.name() == "EQLNUM") {
            exteded_val = upper_equilnum;
        }

        //std::cout << "Extend regions: " << records.name() << std::endl;
        auto& val_status = const_cast<std::vector<Opm::value::status>&>(item.getValueStatus());

        std::vector<int> val_new(nc_new);
        std::vector<Opm::value::status> val_status_new(nc_new, Opm::value::status::deck_value);

        const int nz_bottom = nz_upper + nz;
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int k = 0; k < nz_new; ++k) {
                    const int k_old = k - nz_upper;
                    const int ind_new = i + nx * j + nx * ny * k;
                    const int ind_old = i + nx * j + nx * ny * (k_old);

                    if (k >= nz_upper && k < nz_bottom && actnum_old[ind_old] == 1) {
                        val_new[ind_new] = val[ind_old];
                        val_status_new[ind_new] = val_status[ind_old];
                    } else {
                        val_status_new[ind_new] = Opm::value::status::deck_value;
                        val_new[ind_new] = exteded_val;
                    }
                }
            }
        }

        val_status = val_status_new;
        val = val_new;
    }
}

void
extendSolution(Opm::SOLUTIONSection& solution,
               [[maybe_unused]] const ExtendParam& extend_param)
{
    for (const auto& records : solution) {
        if (records.name() == "EQUIL") {

            // copy first record
            auto& deckkeyword = const_cast<Opm::DeckKeyword&>(records);
            assert(records.size() > 0);

            auto record = records.getRecord(0);

            double& datum_pressure
                = const_cast<double&>(record.getItem("DATUM_PRESSURE").getData<double>()[0]);

            double& datum_depth
                = const_cast<double&>(record.getItem("DATUM_DEPTH").getData<double>()[0]);

            datum_pressure = 1.0;
            datum_depth = 0.0;

            double& woc = const_cast<double&>(record.getItem("OWC").getData<double>()[0]);
            double& goc = const_cast<double&>(record.getItem("GOC").getData<double>()[0]);
            woc = 0.0;
            goc = 0.0;

            deckkeyword.addRecord(std::move(record));
        }

        if ((records.name() == "RSVD") || records.name() == "RTEMPVD") {
            assert(records.size() > 0);

            auto& deckkeyword = const_cast<Opm::DeckKeyword&>(records);
            auto record = records.getRecord(0);

            deckkeyword.addRecord(std::move(record));
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
    const auto actnum_old = gridsec.get<Opm::ParserKeywords::ACTNUM>().back().getIntData();

    extendDimens(grid_size, extend_param, deck);

    extendBCCon(gridsec, extend_param);

    extendGridSection<double>(gridsec, grid_size, extend_param, Opm::type_tag::fdouble);
    extendGridSection<int>(gridsec, grid_size, extend_param, Opm::type_tag::integer);

    extendSpecialGridProperties(gridsec, extend_param);

    extendGRDECL(gridsec, grid_size, extend_param);

    {
        // NB NB need to always add EQUILNUM
        auto regions = Opm::REGIONSSection {deck}; // extend and add for equilnum

        extendRegions(regions, grid_size, extend_param, actnum_old);
    }

    {
        auto solution = Opm::SOLUTIONSection {deck}; // extend and add for equilnum

        extendSolution(solution, extend_param);
    }

    {
        auto schedule = Opm::SCHEDULESection {deck}; // fixing COMPDAT

        extendSchedule(schedule, extend_param);
    }

    for (const auto& records : deck) {
        // DeckView
        if ((records.name() != "EQUALS") &&
            (records.name() != "MULTIPLY") &&
            (records.name() != "ADD") &&
            (records.name() != "COPY"))
        {
            continue;
        }

        //std::cout << "Exending '" << records.name() << '\'' << std::endl;

        const int nz_upper = extend_param.nz_upper;

        for (const auto& record : records) {
            //std::cout << "Record: " << record << std::endl;

            int& K1 = const_cast<int&>(record.getItem("K1").getData<int>()[0]);
            int& K2 = const_cast<int&>(record.getItem("K2").getData<int>()[0]);
            K1 += nz_upper;
            K2 += nz_upper;
        }
    }

    if (deck.hasKeyword<Opm::ParserKeywords::EQLDIMS>()) {
        std::cout << "EQLDIMS keyword found in the deck" << std::endl;
    }
    else {
        std::cout << "No EQLDIMS keyword found in the deck" << std::endl;
    }

    os << deck;

    return deck;
}

void
print_help_and_exit()
{
    std::cerr << R"(
The manipulatedeck program will load a deck, resolve all include
files and then print it out again on stdout. All comments
will be stripped and the value types will be validated.

By passing the option -o you can redirect the output to a file
or a directory.

It will also pad the grid vertically according to parameters
specified in a json file (default fracture_param.json).

Print on stdout:

   manipulatedeck  /path/to/case/CASE.DATA -p fracture_param.json


Print MY_CASE.DATA in /tmp:

    manipulatedeck -o /tmp /path/to/MY_CASE.DATA -p fracture_param.json


Print NEW_CASE in cwd:

    opmpack -o NEW_CASE.DATA path/to/MY_CASE.DATA -p fracture_param.json

As an alternative to the -o option you can use -c; that is equivalent to -o -
but restart and import files referred to in the deck are also copied. The -o and
-c options are mutually exclusive.
)";

    std::exit(EXIT_FAILURE);
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
    auto coutput_arg = std::string_view{};
    auto param_file_arg = std::string_view{"fracture_param.json"};

    while (true) {
        const auto c = getopt(argc, argv, "c:o:p:");
        if (c == -1) {
            break;
        }

        switch (c) {
        case 'c':
            stdout_output = false;
            copy_binary = true;
            coutput_arg = optarg;
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
    if (arg_offset >= argc) {
        print_help_and_exit();
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
