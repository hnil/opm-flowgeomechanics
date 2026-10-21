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
#ifndef OPM_MECH_GRID_CONTEXT_HPP
#define OPM_MECH_GRID_CONTEXT_HPP

#include <opm/geomech/MechFlowMap.hpp>
#include <opm/geomech/MechGridPartition.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/common/CartesianIndexMapper.hpp>
#include <opm/grid/cpgrid/coarsening/CornerPointCoarsening.hpp>
#include <opm/grid/cpgrid/GeometricCheck.hpp>
#include <opm/grid/cpgrid/RetainedCornerPointInput.hpp>
#include <opm/grid/cpgrid/refinement/GridStateWriter.hpp>
#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>
#include <opm/grid/cpgrid/refinement/RefinementRequest.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/input/eclipse/EclipseState/Grid/LgrCollection.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <dune/grid/common/mcmgmapper.hh>

#include <algorithm>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <type_traits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace Opm::Parameters {
    struct MechGridCheck {
        inline static std::string value{"error"};
    };
}

namespace Opm
{

/// The mechanics grid when it is not the flow grid: a coarsening of the deck's
/// grid, built with the same corner points, plus the map to the flow cells.
///
/// In parallel the flow grid must have been partitioned by
/// MechPartition::flowPartition, which puts every flow cell on the rank of its
/// coarse cell; the mechanics grid is then distributed the same way, so a
/// mechanics cell and all of its flow cells share a rank.
class MechGridContext
{
public:
    static std::vector<Coarsening::CoarsenRequest> readRecords(const std::string& filename)
    { return MechPartition::readRecords(filename); }

    /// Build the mechanics grid from the deck's grid and the records, and the
    /// map from the flow grid's cells to its own.
    /// `processed` is what the flow grid was built from (MINPV and PINCH
    /// included), and exists on rank 0 only.
    enum class Method { Auto, Grdecl, Merge, Collapse };

    static Method method(const std::string& name)
    {
        if (name == "auto")   { return Method::Auto; }
        if (name == "grdecl") { return Method::Grdecl; }
        if (name == "merge")  { return Method::Merge; }
        if (name == "collapse") { return Method::Collapse; }
        OPM_THROW(std::runtime_error,
                  "Unknown mechanics coarsening method '" + name
                  + "'; use grdecl, merge, collapse or auto");
    }

    template <class FlowGrid, class CartesianMapper>
    MechGridContext(const RetainedCornerPointInput* processed,
                    const std::vector<Coarsening::CoarsenRequest>& requests,
                    const FlowGrid& flowGrid,
                    const CartesianMapper& flowCartesianMapper,
                    const Method wanted = Method::Auto)
    {
        if (flowGrid.comm().size() > 1) {
            // The balancer runs on rank 0 only, so ask it there and let every
            // rank throw together: a throw on one rank alone would leave the
            // others in a collective call.
            int haveIt = (flowGrid.comm().rank() != 0) || !MechPartition::Registry::empty();
            haveIt = flowGrid.comm().min(haveIt);
            if (!haveIt) {
                OPM_THROW_NOLOG(std::runtime_error,
                                "In parallel the flow grid must be partitioned by the "
                                "mechanics coarsening: install "
                                "MechPartition::installFlowPartition in main() before the "
                                "simulator is created.");
            }
        }
        // Only rank 0 holds the processed description, so it decides for all.
        const bool isRoot = flowGrid.comm().rank() == 0;
        int problem = 0;
        if (isRoot) {
            problem = (processed == nullptr) ? 1 : 0;
        }
        problem = flowGrid.comm().max(problem);
        if (problem != 0) {
            OPM_THROW_NOLOG(std::runtime_error,
                            "The mechanics grid needs the description the flow grid was built "
                            "from: call RetainCornerPointInput::enable() before the grid is "
                            "created (MechPartition::installFlowPartition does).");
        }

        const auto size = flowGrid.logicalCartesianSize();
        const std::array<int,3> fineDims{size[0], size[1], size[2]};

        // Corner-point processing is a rank-0 job with collective steps, so
        // every rank calls it and only rank 0 passes the description.
        // auto merges cells with one face between coarse cells (collapse): the
        // corner-point coarsening where that exists, and also part of a column.
        merged_ = wanted != Method::Grdecl;
        const bool collapse = wanted == Method::Collapse || wanted == Method::Auto;

        const auto cornerPoint = [&]() {
            Coarsening::Grdecl fine;
            fine.dims = fineDims;
            fine.coord = processed->coord;
            fine.zcorn = processed->zcorn;
            fine.actnum = processed->actnum;

            Coarsening::Options options;
            // Mechanics wants rock everywhere the block has volume, so cells
            // that are inactive or that MINPV removed are absorbed instead of
            // left as holes.
            options.activity = Coarsening::Activity::FillHoles;
            // Without edge-conformal processing, MINPV leaves the removed
            // cells' volume as a gap between their neighbours; the mechanics
            // body takes it back as rock.
            options.allowVerticalGaps = !processed->edgeConformal;
            const auto result = Coarsening::coarsenCornerPoint(fine, requests, options);
            for (const auto& note : result.report.notes) {
                OpmLog::info("Mechanics grid: " + note);
            }
            if (result.report.absorbedGapVolume > 0.0) {
                OpmLog::warning("Mechanics grid: took "
                                + std::to_string(result.report.absorbedGapVolume)
                                + " m3 of gaps left by cell removal back as rock. Run with "
                                "--edge-conformal=true so the removal merges the cells "
                                "geometrically; gaps under cells the mechanics does not "
                                "coarsen stay holes in the body.");
            }
            return result;
        };

        // Where the records are a corner-point coarsening, auto keeps its checks
        // (faults, face counts, nested layer groupings).
        if (wanted == Method::Auto) {
            std::string refused;
            if (isRoot) {
                try {
                    Coarsening::cartesianMap(fineDims, requests);
                    try {
                        cornerPoint();
                    } catch (const std::exception& e) {
                        refused = e.what();
                    }
                } catch (const std::invalid_argument& e) {
                    OpmLog::info(std::string("Mechanics grid: not a corner-point coarsening (")
                                 + e.what() + "); merging cells of the flow grid.");
                }
            }
            if (flowGrid.comm().max(static_cast<int>(!refused.empty())) != 0) {
                OPM_THROW_NOLOG(std::runtime_error,
                                "Mechanics coarsening refused" + (refused.empty() ? std::string{}
                                                                  : ": " + refused));
            }
        }

        // The index side needs no geometry, so every rank can have it.
        Coarsening::CartesianMap cartesian;
        if (!merged_) {
            cartesian = Coarsening::cartesianMap(fineDims, requests);
        }

        std::unique_ptr<EclipseGrid> coarseGrid;
        if (processed != nullptr && !merged_) {
            const auto result = cornerPoint();
            coarseGrid = std::make_unique<EclipseGrid>(result.grid.dims, result.grid.coord,
                                                       result.grid.zcorn,
                                                       result.grid.actnum.data());
        }

        grid_ = std::make_unique<Dune::CpGrid>(flowGrid.comm());
        if (merged_) {
            grdecl input{};
            std::vector<int> allActive;
            Coarsening::BlockLayout layout;
            if (processed != nullptr) {
                input.dims[0] = fineDims[0];
                input.dims[1] = fineDims[1];
                input.dims[2] = fineDims[2];
                input.coord = processed->coord.data();
                input.zcorn = processed->zcorn.data();
                // Rock wherever there is volume, as the grdecl route's FillHoles:
                // inactive cells are mechanics cells too. Zero-thickness cells stay
                // inactive, so the cells above and below them still touch.
                allActive = withVolume(fineDims, processed->zcorn);
                input.actnum = allActive.data();
                layout = Coarsening::blockLayout(fineDims, requests);
                std::ostringstream os;
                os << "Mechanics grid: merging " << fineDims[0]*fineDims[1]*fineDims[2]
                   << " cells into " << layout.boxes.size() << " blocks";
                OpmLog::info(os.str());
            }
            grid_->processEclipseFormatCoarsened(input, layout.blockOfCartesian, layout.boxes,
                                                 /*edge_conformal*/ true,
                                                 /*collapse_coarse_faces*/ collapse);
        } else {
            grid_->processEclipseFormat(coarseGrid.get(), /*ecl_state*/ nullptr,
                                        /*periodic_extension*/ false, /*turn_normals*/ false,
                                        /*clip_z*/ false, /*pinchActive*/ false,
                                        /*edge_conformal*/ true);
        }

        if (flowGrid.comm().size() > 1) {
            // The balancer gave every Cartesian cell a rank; the merge keeps
            // that index space, the grdecl route has its own.
            const auto& partOfFine = MechPartition::Registry::get();
            std::vector<int> partOfMech;
            if (merged_) {
                partOfMech = partOfFine;
            } else if (flowGrid.comm().rank() == 0) {
                partOfMech.assign(static_cast<std::size_t>(cartesian.coarseDims[0])
                                  * cartesian.coarseDims[1]*cartesian.coarseDims[2], 0);
                for (std::size_t f = 0; f < cartesian.fineToCoarse.size(); ++f) {
                    if (cartesian.fineToCoarse[f] >= 0) {
                        partOfMech[cartesian.fineToCoarse[f]] = partOfFine[f];
                    }
                }
            }
            distribute(partOfMech);
        }
        if (merged_) {
            // Anchor index of each cell's block, in the fine index space.
            const auto layout = Coarsening::blockLayout(fineDims, requests);
            std::vector<int> anchorOfCartesian(layout.blockOfCartesian.size(), -1);
            for (std::size_t c = 0; c < layout.blockOfCartesian.size(); ++c) {
                const auto& box = layout.boxes[layout.blockOfCartesian[c]];
                anchorOfCartesian[c] = box[0] + fineDims[0]*(box[1] + fineDims[1]*box[2]);
            }
            buildMap(anchorOfCartesian, fineDims, flowGrid, flowCartesianMapper);
        } else {
            buildMap(cartesian.fineToCoarse, cartesian.coarseDims, flowGrid, flowCartesianMapper);
        }
    }

    /// The flow grid's corner points processed edge-conformal, with cells
    /// thinner than `thickness` (0: none) collapsed into the cell below (above
    /// at the bottom), inactive cells with volume filled the same way when
    /// `fillHoles`, coarsened by `requests` (each block one polyhedral cell
    /// with `merge`, else a corner-point coarsening), and the deck's LGRs
    /// refined again. Flow is untouched. Serial only.
    template <class FlowGrid, class CartesianMapper>
    MechGridContext(const RetainedCornerPointInput* processed,
                    const double thickness,
                    const bool fillHoles,
                    const std::vector<Coarsening::CoarsenRequest>& requests,
                    const bool merge,
                    const LgrCollection& lgrs,
                    const FlowGrid& flowGrid,
                    const CartesianMapper& flowCartesianMapper)
    {
        if (flowGrid.comm().size() > 1) {
            OPM_THROW(std::runtime_error, "--mech-grid=separate runs serially only");
        }
        if (processed == nullptr) {
            OPM_THROW(std::runtime_error,
                      "The mechanics grid needs the description the flow grid was built from.");
        }
        // A copy: processing a grid may replace what `processed` points to.
        auto input = std::make_shared<RetainedCornerPointInput>(*processed);
        if (!input->edgeConformal) {
            closeVerticalGaps_(*input);
            input->edgeConformal = true;
        }
        Coarsening::BlockLayout layout;
        std::vector<char> merged;
        if (merge && !requests.empty()) {
            layout = mergeLayout_(*input, requests, lgrs, thickness);
            merged.resize(layout.blockOfCartesian.size());
            for (std::size_t c = 0; c < merged.size(); ++c) {
                const auto& box = layout.boxes[layout.blockOfCartesian[c]];
                merged[c] = box[0] != box[3] || box[1] != box[4] || box[2] != box[5];
            }
        }
        // A sliver inside a merged block is part of one cell anyway.
        const auto absorber = collapseThinCells_(*input, thickness, fillHoles, lgrs, merged);
        fineDims_ = input->dims;
        fromCornerPoints_ = true;

        if (merge && !requests.empty()) {
            buildMerged_(input, layout, absorber, lgrs, flowGrid.comm());
        } else {
            auto mechInput = input;
            std::vector<int> fineToCoarse(absorber.size());
            std::iota(fineToCoarse.begin(), fineToCoarse.end(), 0);
            if (!requests.empty()) {
                Coarsening::Grdecl fine{input->dims, input->coord, input->zcorn, input->actnum};
                Coarsening::Options options;
                options.activity = Coarsening::Activity::FillHoles;
                options.allowVerticalGaps = true;
                auto result = Coarsening::coarsenCornerPoint(fine, requests, options);
                for (const auto& note : result.report.notes) {
                    OpmLog::info("Mechanics grid: " + note);
                }
                mechInput = std::make_shared<RetainedCornerPointInput>();
                mechInput->dims = result.grid.dims;
                mechInput->coord = std::move(result.grid.coord);
                mechInput->zcorn = std::move(result.grid.zcorn);
                mechInput->actnum = std::move(result.grid.actnum);
                mechInput->edgeConformal = true;
                mechInput->pinchActive = true;
                fineToCoarse = std::move(result.fineToCoarse);
                blocks_ = std::move(result.blocks);
            }

            grid_ = std::make_unique<Dune::CpGrid>(flowGrid.comm());
            grid_->setLgrBackend(Refinement::Backend::Conforming);
            const grdecl raw{ { mechInput->dims[0], mechInput->dims[1], mechInput->dims[2] },
                              mechInput->coord.data(), mechInput->zcorn.data(),
                              mechInput->actnum.empty() ? nullptr : mechInput->actnum.data() };
            // Edge-conformal processing joins cells across zero-thickness ones as
            // a pinch would; the retained flag is the deck's.
            grid_->processEclipseFormat(raw, mechInput->pinchNnc, /*pinchActive*/ true,
                                        /*edge_conformal*/ true);
            Refinement::GridStateWriter::setRetainedCornerPointInput(*grid_->currentData().front(),
                                                                     mechInput);
            if (lgrs.size() > 0) {
                grid_->addLgrsUpdateLeafView(lgrRequests_(lgrs, input->dims, mechInput->dims,
                                                          fineToCoarse));
            }
            fineToMech_.resize(absorber.size());
            for (std::size_t f = 0; f < absorber.size(); ++f) {
                fineToMech_[f] = fineToCoarse[absorber[f]];
            }
        }
        buildLeafMap_(fineToMech_, flowGrid, flowCartesianMapper);
    }

    const Dune::CpGrid& grid() const { return *grid_; }
    /// Built from the flow grid's corner points: BCCON boxes, in the flow
    /// grid's Cartesian space, are found with cellOfFineFace().
    bool fromCornerPoints() const { return fromCornerPoints_; }
    const std::array<int,3>& fineDims() const { return fineDims_; }

    /// The leaf grid as VTK hexahedra (8 corners per cell, unshared), with
    /// each cell's level and number of flow cells; for looking at the grid.
    void writeVtu(const std::string& filename) const
    {
        const auto& gv = grid_->leafGridView();
        std::vector<int> flowCells(gv.size(0), 0);
        for (std::size_t f = 0; f < map_.numFlowCells(); ++f) {
            if (const int m = map_.mechCell(f); m >= 0) {
                ++flowCells[m];
            }
        }
        std::ofstream os(filename);
        os.precision(12);
        const auto n = gv.size(0);
        os << "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\""
           << " byte_order=\"LittleEndian\">\n"
           << "<UnstructuredGrid><Piece NumberOfPoints=\"" << 8*n << "\" NumberOfCells=\"" << n
           << "\">\n";
        // ResInsight's reader needs a displacement field.
        os << "<PointData Vectors=\"disp\">\n<DataArray type=\"Float32\" Name=\"disp\" "
           << "NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (std::size_t p = 0; p < 8*n; ++p) {
            os << "0 0 0\n";
        }
        os << "</DataArray>\n</PointData>\n<Points><DataArray type=\"Float64\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        // Dune's corner order to VTK's hexahedron order.
        constexpr int toVtk[8] = {0, 1, 3, 2, 4, 5, 7, 6};
        std::vector<int> level(n), index(n);
        for (const auto& cell : elements(gv)) {
            const int c = gv.indexSet().index(cell);
            level[c] = cell.level();
            for (const int k : toVtk) {
                os << cell.geometry().corner(k) << "\n";
            }
        }
        os << "</DataArray></Points>\n<Cells><DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
        for (std::size_t p = 0; p < 8*n; ++p) {
            os << p << ((p % 8 == 7) ? "\n" : " ");
        }
        os << "</DataArray><DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
        for (std::size_t c = 1; c <= n; ++c) {
            os << 8*c << "\n";
        }
        os << "</DataArray><DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        for (std::size_t c = 0; c < n; ++c) {
            os << "12\n";
        }
        os << "</DataArray></Cells>\n<CellData>\n<DataArray type=\"Float32\" Name=\"LEVEL\" format=\"ascii\">\n";
        for (const int v : level) {
            os << v << "\n";
        }
        os << "</DataArray>\n<DataArray type=\"Float32\" Name=\"FLOWCELLS\" format=\"ascii\">\n";
        for (const int v : flowCells) {
            os << v << "\n";
        }
        // Logical (parent) IJK, for slicing in ResInsight.
        os << "</DataArray>\n<DataArray type=\"Float32\" Name=\"ijk\" NumberOfComponents=\"3\" "
           << "format=\"ascii\">\n";
        const Dune::CartesianIndexMapper<Dune::CpGrid> mapper(*grid_);
        for (std::size_t c = 0; c < n; ++c) {
            std::array<int,3> ijk{};
            mapper.cartesianCoordinate(static_cast<int>(c), ijk);
            os << ijk[0] << " " << ijk[1] << " " << ijk[2] << "\n";
        }
        os << "</DataArray>\n</CellData>\n</Piece></UnstructuredGrid>\n</VTKFile>\n";
    }

    /// The unrefined mechanics cell whose face carries the given face of flow
    /// Cartesian cell `fineCart`, or -1 if that face is inside a coarse cell.
    int cellOfFineFace(const int fineCart, const FaceDir::DirEnum dir) const
    {
        const int m = mechCellOfMechCart_[fineToMech_[fineCart]];
        if (m < 0 || blocks_.empty()) {
            return m;
        }
        const int nx = fineDims_[0], ny = fineDims_[1];
        const std::array<int,3> ijk{fineCart % nx, (fineCart/nx) % ny, fineCart/(nx*ny)};
        const auto& block = blocks_[fineToMech_[fineCart]];
        switch (dir) {
        case FaceDir::XMinus: return ijk[0] == block.startIJK[0] ? m : -1;
        case FaceDir::XPlus:  return ijk[0] == block.endIJK[0] - 1 ? m : -1;
        case FaceDir::YMinus: return ijk[1] == block.startIJK[1] ? m : -1;
        case FaceDir::YPlus:  return ijk[1] == block.endIJK[1] - 1 ? m : -1;
        case FaceDir::ZMinus: return ijk[2] == block.startIJK[2] ? m : -1;
        case FaceDir::ZPlus:  return ijk[2] == block.endIJK[2] - 1 ? m : -1;
        default: return -1;
        }
    }
    Dune::CpGrid& grid() { return *grid_; }
    const MechFlowMap& map() const { return map_; }

    /// Fill the overlap cells of a mechanics-cell vector from their owners.
    /// Restriction only fills cells whose flow children are local.
    template <class Vector>
    void communicate(Vector& mechCellData) const
    {
        if (grid_->comm().size() == 1) {
            return;
        }
        const auto& indexSet = grid_->leafGridView().indexSet();
        VectorHandle<Vector, std::decay_t<decltype(indexSet)>> handle(mechCellData, indexSet);
        grid_->communicate(handle, Dune::InteriorBorder_All_Interface,
                           Dune::ForwardCommunication);
    }

private:
    /// Copies one value per cell, owner to copy.
    template <class Vector, class IndexSet>
    class VectorHandle
    {
    public:
        VectorHandle(Vector& data, const IndexSet& indexSet)
            : data_(data), index_set_(indexSet) {}

        using DataType = double;
        bool contains(int /*dim*/, int codim) const { return codim == 0; }
        bool fixedSize(int /*dim*/, int /*codim*/) const { return true; }
        template <class Entity> std::size_t size(const Entity&) const { return 1; }

        template <class Buffer, class Entity>
        void gather(Buffer& buffer, const Entity& e) const
        { buffer.write(value(data_[index_set_.index(e)])); }

        template <class Buffer, class Entity>
        void scatter(Buffer& buffer, const Entity& e, std::size_t)
        {
            double v;
            buffer.read(v);
            data_[index_set_.index(e)] = v;
        }

    private:
        template <class T>
        static double value(const T& v)
        {
            if constexpr (std::is_arithmetic_v<T>) {
                return v;
            } else {
                return v[0];
            }
        }

        Vector& data_;
        const IndexSet& index_set_;
    };

    /// Each block of `requests` one cell, the union of its fine cells with all
    /// their outer faces: a local coarsening, the inverse of an LGR. Keeps the
    /// fine Cartesian space; a block is known by its first cell.
    /// The blocks of `requests`, less those that cannot be merged: a corner
    /// column without a cell (no eight corners), touching an LGR box (which
    /// needs fine neighbours to refine against), or holding a cell thinner
    /// than `thickness` (merging a sliver leaves hanging edges; it is
    /// collapsed instead), or crossed by a fault.
    static Coarsening::BlockLayout mergeLayout_(const RetainedCornerPointInput& input,
                                                const std::vector<Coarsening::CoarsenRequest>& requests,
                                                const LgrCollection& lgrs,
                                                const double thickness)
    {
        const auto& d = input.dims;
        auto layout = Coarsening::blockLayout(d, requests);
        const auto cart = [&d](int i, int j, int k) { return i + d[0]*(j + d[1]*k); };
        const auto present = [&](int i, int j, int k) {
            if (!input.actnum.empty() && input.actnum[cart(i, j, k)] == 0) {
                return false;
            }
            for (int c = 0; c < 4; ++c) {
                const std::size_t top = 2*i + (c & 1) + 2*d[0]*(2*j + (c >> 1))
                    + 4*static_cast<std::size_t>(d[0])*d[1]*(2*k);
                if (input.zcorn[top + 4*static_cast<std::size_t>(d[0])*d[1]] != input.zcorn[top]) {
                    return true;
                }
            }
            return false;
        };
        const auto touchesLgr = [&lgrs](const std::array<int,6>& box) {
            for (std::size_t l = 0; l < lgrs.size(); ++l) {
                const auto& b = lgrs.getLgr(l);
                if (box[0] <= b.I2() + 1 && box[3] >= b.I1() - 1 && box[1] <= b.J2() + 1
                    && box[4] >= b.J1() - 1 && box[2] <= b.K2() + 1 && box[5] >= b.K1() - 1) {
                    return true;
                }
            }
            return false;
        };
        const auto thin = [&](int i, int j, int k) {
            double mean = 0.0;
            for (int c = 0; c < 4; ++c) {
                const std::size_t top = 2*i + (c & 1) + 2*d[0]*(2*j + (c >> 1))
                    + 4*static_cast<std::size_t>(d[0])*d[1]*(2*k);
                mean += 0.25*(input.zcorn[top + 4*static_cast<std::size_t>(d[0])*d[1]] - input.zcorn[top]);
            }
            return mean > 0.0 && mean < thickness;
        };
        const auto zc = [&](int i, int j, int k, int di, int dj, int dk) {
            return input.zcorn[2*i + di + 2*d[0]*(2*j + dj)
                               + 4*static_cast<std::size_t>(d[0])*d[1]*(2*k + dk)];
        };
        // A throw between two columns of the block.
        const auto faulted = [&](const std::array<int,6>& box) {
            for (int k = box[2]; k <= box[5]; ++k) {
                for (int j = box[1]; j <= box[4]; ++j) {
                    for (int i = box[0]; i <= box[3]; ++i) {
                        for (int a = 0; a < 2; ++a) {
                            for (int dk = 0; dk < 2; ++dk) {
                                if (i < box[3] && zc(i, j, k, 1, a, dk) != zc(i + 1, j, k, 0, a, dk)) {
                                    return true;
                                }
                                if (j < box[4] && zc(i, j, k, a, 1, dk) != zc(i, j + 1, k, a, 0, dk)) {
                                    return true;
                                }
                            }
                        }
                    }
                }
            }
            return false;
        };
        int noCorner = 0, nearLgr = 0, sliver = 0, fault = 0;
        const std::size_t numBoxes = layout.boxes.size();
        std::vector<char> broken(numBoxes, 0);
        for (std::size_t b = 0; b < numBoxes; ++b) {
            const auto box = layout.boxes[b];
            if (box[0] == box[3] && box[1] == box[4] && box[2] == box[5]) {
                continue;
            }
            if (touchesLgr(box)) {
                broken[b] = 1;
                ++nearLgr;
                continue;
            }
            for (int k = box[2]; k <= box[5] && !broken[b]; ++k) {
                for (int j = box[1]; j <= box[4] && !broken[b]; ++j) {
                    for (int i = box[0]; i <= box[3] && !broken[b]; ++i) {
                        broken[b] = thin(i, j, k);
                    }
                }
            }
            if (broken[b]) {
                ++sliver;
                continue;
            }
            if (faulted(box)) {
                broken[b] = 1;
                ++fault;
                continue;
            }
            for (int corner = 0; corner < 4 && !broken[b]; ++corner) {
                const int i = (corner & 1) ? box[3] : box[0];
                const int j = (corner & 2) ? box[4] : box[1];
                bool any = false;
                for (int k = box[2]; k <= box[5] && !any; ++k) {
                    any = present(i, j, k);
                }
                broken[b] = !any;
            }
            noCorner += broken[b];
        }
        for (std::size_t c = 0; c < layout.blockOfCartesian.size(); ++c) {
            if (broken[layout.blockOfCartesian[c]]) {
                const int i = c % d[0], j = (c/d[0]) % d[1], k = c/(d[0]*d[1]);
                layout.blockOfCartesian[c] = static_cast<int>(layout.boxes.size());
                layout.boxes.push_back({i, j, k, i, j, k});
            }
        }
        if (noCorner + nearLgr + sliver + fault > 0) {
            OpmLog::info("Mechanics grid: left blocks unmerged: " + std::to_string(nearLgr)
                         + " next to an LGR, " + std::to_string(sliver) + " with a thin cell, "
                         + std::to_string(fault) + " crossed by a fault, "
                         + std::to_string(noCorner) + " with an empty corner column");
        }
        return layout;
    }

    /// Each block of `layout` one cell, the union of its fine cells with all
    /// their outer faces: a local coarsening, the inverse of an LGR. Keeps the
    /// fine Cartesian space; a block is known by its first cell.
    template <class Comm>
    void buildMerged_(const std::shared_ptr<RetainedCornerPointInput>& input,
                      const Coarsening::BlockLayout& layout,
                      const std::vector<int>& absorber, const LgrCollection& lgrs,
                      const Comm& comm)
    {
        const auto& d = input->dims;
        const auto cart = [&d](int i, int j, int k) { return i + d[0]*(j + d[1]*k); };
        blocks_.assign(layout.blockOfCartesian.size(), Coarsening::Block{});
        std::vector<int> anchor(layout.blockOfCartesian.size());
        for (std::size_t c = 0; c < anchor.size(); ++c) {
            const auto& box = layout.boxes[layout.blockOfCartesian[c]];
            anchor[c] = cart(box[0], box[1], box[2]);
            blocks_[anchor[c]] = { {box[0], box[1], box[2]}, {box[3] + 1, box[4] + 1, box[5] + 1},
                                   anchor[c] };
        }
        for (std::size_t l = 0; l < lgrs.size(); ++l) {
            const auto& b = lgrs.getLgr(l);
            for (int k = b.K1(); k <= b.K2(); ++k) {
                for (int j = b.J1(); j <= b.J2(); ++j) {
                    for (int i = b.I1(); i <= b.I2(); ++i) {
                        if (anchor[cart(i, j, k)] != cart(i, j, k)
                            || blocks_[cart(i, j, k)].endIJK[0] - i != 1
                            || blocks_[cart(i, j, k)].endIJK[1] - j != 1
                            || blocks_[cart(i, j, k)].endIJK[2] - k != 1) {
                            OPM_THROW(std::runtime_error, "Mechanics coarsening merges cells of LGR "
                                      + b.NAME() + "; keep the records out of the box");
                        }
                    }
                }
            }
        }
        grid_ = std::make_unique<Dune::CpGrid>(comm);
        grid_->setLgrBackend(Refinement::Backend::Conforming);
        const grdecl raw{ { d[0], d[1], d[2] }, input->coord.data(), input->zcorn.data(),
                          input->actnum.empty() ? nullptr : input->actnum.data() };
        grid_->processEclipseFormatCoarsened(raw, layout.blockOfCartesian, layout.boxes,
                                             /*edge_conformal*/ true,
                                             /*collapse_coarse_faces*/ false);
        Refinement::GridStateWriter::setRetainedCornerPointInput(*grid_->currentData().front(),
                                                                 input);
        std::ostringstream os;
        os << "Mechanics grid: merged " << anchor.size() << " cells into "
           << layout.boxes.size() << " blocks";
        OpmLog::info(os.str());
        if (lgrs.size() > 0) {
            std::vector<int> identity(anchor.size());
            std::iota(identity.begin(), identity.end(), 0);
            grid_->addLgrsUpdateLeafView(lgrRequests_(lgrs, d, d, identity));
        }
        fineToMech_.resize(absorber.size());
        for (std::size_t f = 0; f < absorber.size(); ++f) {
            fineToMech_[f] = anchor[absorber[f]];
        }
    }

    /// Without edge-conformal processing MINPV leaves its removed cells'
    /// volume as a gap in the column; the next cell below takes it as rock.
    static void closeVerticalGaps_(RetainedCornerPointInput& input)
    {
        const int nx = input.dims[0], ny = input.dims[1], nz = input.dims[2];
        auto& z = input.zcorn;
        const auto zi = [nx, ny](int i, int j, int k, int c, int dk) {
            return static_cast<std::size_t>((2*k + dk)*4*nx*ny + (2*j + (c >> 1))*2*nx + 2*i + (c & 1));
        };
        const auto solid = [&](int i, int j, int k) {
            if (!input.actnum.empty() && input.actnum[i + nx*(j + ny*k)] == 0) {
                return false;
            }
            for (int c = 0; c < 4; ++c) {
                if (z[zi(i, j, k, c, 1)] != z[zi(i, j, k, c, 0)]) {
                    return true;
                }
            }
            return false;
        };
        double gapVolumeProxy = 0.0;
        int closed = 0;
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int above = -1;
                for (int k = 0; k < nz; ++k) {
                    if (!solid(i, j, k)) {
                        continue;
                    }
                    // Only gaps with nothing of volume in between: active
                    // cells with volume and ACTNUM=0 holes are left alone.
                    bool empty = above >= 0;
                    for (int kk = above + 1; empty && kk < k; ++kk) {
                        for (int c = 0; c < 4; ++c) {
                            empty = empty && z[zi(i, j, kk, c, 1)] == z[zi(i, j, kk, c, 0)];
                        }
                    }
                    bool moved = false;
                    for (int c = 0; empty && c < 4; ++c) {
                        const double bottom = z[zi(i, j, above, c, 1)];
                        if (z[zi(i, j, k, c, 0)] > bottom) {
                            gapVolumeProxy += z[zi(i, j, k, c, 0)] - bottom;
                            z[zi(i, j, k, c, 0)] = bottom;
                            moved = true;
                        }
                        for (int kk = above + 1; kk < k; ++kk) {
                            z[zi(i, j, kk, c, 0)] = z[zi(i, j, kk, c, 1)] = bottom;
                        }
                    }
                    closed += moved;
                    above = k;
                }
            }
        }
        if (closed > 0) {
            OpmLog::info("Mechanics grid: closed " + std::to_string(closed)
                         + " vertical gaps left by removed cells (sum of corner gaps "
                         + std::to_string(gapVolumeProxy) + " m)");
        }
    }

    /// Collapse thin active cells, and with `fillHoles` inactive cells that
    /// have volume, into the touching active cell below (above if none);
    /// returns, for every Cartesian cell, the cell that took it over.
    static std::vector<int> collapseThinCells_(RetainedCornerPointInput& input,
                                               const double thickness,
                                               const bool fillHoles,
                                               const LgrCollection& lgrs,
                                               const std::vector<char>& keep = {})
    {
        const int nx = input.dims[0], ny = input.dims[1], nz = input.dims[2];
        auto& z = input.zcorn;
        const auto zi = [nx, ny](int i, int j, int k, int c, int dk) {
            return static_cast<std::size_t>((2*k + dk)*4*nx*ny + (2*j + (c >> 1))*2*nx + 2*i + (c & 1));
        };
        const auto active = [&input, nx, ny](int i, int j, int k) {
            return input.actnum.empty() || input.actnum[i + nx*(j + ny*k)] != 0;
        };
        const auto inLgr = [&lgrs](int i, int j, int k) {
            for (std::size_t l = 0; l < lgrs.size(); ++l) {
                const auto& b = lgrs.getLgr(l);
                if (i >= b.I1() && i <= b.I2() && j >= b.J1() && j <= b.J2()
                    && k >= b.K1() && k <= b.K2()) {
                    return true;
                }
            }
            return false;
        };
        const auto meanThickness = [&](int i, int j, int k) {
            double mean = 0.0;
            for (int c = 0; c < 4; ++c) {
                mean += 0.25*(z[zi(i, j, k, c, 1)] - z[zi(i, j, k, c, 0)]);
            }
            return mean;
        };
        // Bottom of `upper` on top of `lower` at all four corners.
        const auto touches = [&](int i, int j, int upper, int lower) {
            for (int c = 0; c < 4; ++c) {
                if (z[zi(i, j, upper, c, 1)] != z[zi(i, j, lower, c, 0)]) {
                    return false;
                }
            }
            return true;
        };

        std::vector<int> absorber(static_cast<std::size_t>(nx)*ny*nz);
        std::iota(absorber.begin(), absorber.end(), 0);
        int thin = 0, holes = 0, skipped = 0;
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                const auto absorbable = [&](int k) {
                    const double t = meanThickness(i, j, k);
                    return t > 0.0 && (active(i, j, k) ? t < thickness : fillHoles);
                };
                for (int k = 0; k < nz; ++k) {
                    if (!absorbable(k)) {
                        continue;
                    }
                    const int self = i + nx*(j + ny*k);
                    if (!keep.empty() && keep[self]) {
                        continue;
                    }
                    if (!inLgr(i, j, k)) {
                        // Down, also into a hole that is absorbed in its turn.
                        if (k + 1 < nz && meanThickness(i, j, k + 1) > 0.0
                            && (active(i, j, k + 1) || fillHoles)
                            && !inLgr(i, j, k + 1) && touches(i, j, k, k + 1)) {
                            for (int c = 0; c < 4; ++c) {
                                const double top = z[zi(i, j, k, c, 0)];
                                z[zi(i, j, k, c, 1)] = top;
                                z[zi(i, j, k + 1, c, 0)] = top;
                            }
                            absorber[self] = self + nx*ny;
                            (active(i, j, k) ? thin : holes) += 1;
                            continue;
                        }
                        // Up, past cells already collapsed.
                        int a = k - 1;
                        while (a >= 0 && meanThickness(i, j, a) <= 0.0) {
                            --a;
                        }
                        if (a >= 0 && active(i, j, a) && !inLgr(i, j, a) && touches(i, j, a, k)) {
                            for (int c = 0; c < 4; ++c) {
                                const double bottom = z[zi(i, j, k, c, 1)];
                                for (int kk = a + 1; kk <= k; ++kk) {
                                    z[zi(i, j, kk, c, 0)] = z[zi(i, j, kk, c, 1)] = bottom;
                                }
                                z[zi(i, j, a, c, 1)] = bottom;
                            }
                            absorber[self] = i + nx*(j + ny*a);
                            (active(i, j, k) ? thin : holes) += 1;
                            continue;
                        }
                    }
                    ++skipped;
                }
            }
        }
        // A run of collapsed cells ends in the last one's absorber.
        for (auto& a : absorber) {
            while (absorber[a] != a) {
                a = absorber[a];
            }
        }
        if (thin + holes + skipped == 0) {
            return absorber;
        }
        std::ostringstream os;
        os << "Mechanics grid: merged " << thin << " cells thinner than " << thickness
           << " m and " << holes << " inactive cells into a vertical neighbour";
        if (skipped > 0) {
            os << "; kept " << skipped << " (no touching active neighbour, or in an LGR box)";
        }
        OpmLog::info(os.str());
        return absorber;
    }

    /// The deck's CARFIN boxes, as the flow grid's vanguard requests them.
    static std::vector<Refinement::BlockRefinement>
    lgrRequests_(const LgrCollection& lgrs, const std::array<int,3>& fineDims,
                 const std::array<int,3>& mechDims, const std::vector<int>& fineToMech)
    {
        std::vector<Refinement::BlockRefinement> requests;
        for (std::size_t l = 0; l < lgrs.size(); ++l) {
            const auto& b = lgrs.getLgr(l);
            if (b.PARENT_NAME() != "GLOBAL") {
                OPM_THROW(std::runtime_error,
                          "--mech-thin-cell-thickness does not support nested LGRs");
            }
            Refinement::BlockRefinement r;
            r.name = b.NAME();
            r.parentGridName = b.PARENT_NAME();
            r.cellsPerDim = { b.NX()/(b.I2() + 1 - b.I1()), b.NY()/(b.J2() + 1 - b.J1()),
                              b.NZ()/(b.K2() + 1 - b.K1()) };
            // The box in the mechanics grid's index space; coarsening must
            // leave its cells one to one.
            const auto fine = [&fineDims](int i, int j, int k) {
                return i + fineDims[0]*(j + fineDims[1]*k);
            };
            const int c0 = fineToMech[fine(b.I1(), b.J1(), b.K1())];
            const std::array<int,3> m0{ c0 % mechDims[0], (c0/mechDims[0]) % mechDims[1],
                                        c0/(mechDims[0]*mechDims[1]) };
            for (int k = b.K1(); k <= b.K2(); ++k) {
                for (int j = b.J1(); j <= b.J2(); ++j) {
                    for (int i = b.I1(); i <= b.I2(); ++i) {
                        const int want = (m0[0] + i - b.I1())
                            + mechDims[0]*((m0[1] + j - b.J1()) + mechDims[1]*(m0[2] + k - b.K1()));
                        if (fineToMech[fine(i, j, k)] != want) {
                            OPM_THROW(std::runtime_error, "Mechanics coarsening merges cells of LGR "
                                      + b.NAME() + " or shifts it unevenly; keep the box's pillar "
                                      "lines and layers out of the coarsening records");
                        }
                    }
                }
            }
            r.startIJK = m0;
            r.endIJK = { m0[0] + b.I2() + 1 - b.I1(), m0[1] + b.J2() + 1 - b.J1(),
                         m0[2] + b.K2() + 1 - b.K1() };
            if (b.isGraded()) {
                for (std::size_t d = 0; d < 3; ++d) {
                    auto columns = b.refinedColumns(d);
                    r.subdivision[d] = Refinement::AxisSubdivision{ std::move(columns.parentOffset),
                                                                    std::move(columns.fracLo),
                                                                    std::move(columns.fracHi) };
                }
            }
            r.minpvRemoved = b.minpvRemoved();
            r.pillarsFromBoxLayer = lgrs.pillarsFromBoxLayer();
            requests.push_back(std::move(r));
        }
        return requests;
    }

    /// Flow leaf cell -> mechanics leaf cell. Unrefined cells go by Cartesian
    /// index (through the absorber), refined ones to the nearest refined
    /// mechanics cell of the same parent.
    template <class FlowGrid, class CartesianMapper>
    void buildLeafMap_(const std::vector<int>& absorber, const FlowGrid& flowGrid,
                       const CartesianMapper& flowCartesianMapper)
    {
        const Dune::CartesianIndexMapper<Dune::CpGrid> mechMapper(*grid_);
        const auto& mechView = grid_->leafGridView();
        std::vector<int> coarseMech(absorber.size(), -1);
        std::multimap<int, int> refinedMech;
        std::vector<Dune::FieldVector<double,3>> mechCentre(mechView.size(0));
        for (const auto& cell : elements(mechView)) {
            const int m = mechView.indexSet().index(cell);
            const int cart = mechMapper.cartesianIndex(m);
            mechCentre[m] = cell.geometry().center();
            if (cell.level() > 0) {
                refinedMech.emplace(cart, m);
            } else {
                coarseMech[cart] = m;
            }
        }

        const auto& flowView = flowGrid.leafGridView();
        const std::size_t numFlowCells = flowView.size(0);
        std::vector<int> flowToMech(numFlowCells, -1);
        std::vector<double> weight(numFlowCells, 0.0);
        double worstMiss = 0.0;
        for (const auto& cell : elements(flowView)) {
            const int f = flowView.indexSet().index(cell);
            weight[f] = cell.geometry().volume();
            const int cart = flowCartesianMapper.cartesianIndex(f);
            if (cell.level() == 0) {
                flowToMech[f] = coarseMech[absorber[cart]];
                continue;
            }
            const auto centre = cell.geometry().center();
            double best = std::numeric_limits<double>::max();
            const auto [first, last] = refinedMech.equal_range(absorber[cart]);
            for (auto it = first; it != last; ++it) {
                const double d = (mechCentre[it->second] - centre).two_norm();
                if (d < best) {
                    best = d;
                    flowToMech[f] = it->second;
                }
            }
            if (flowToMech[f] >= 0) {
                worstMiss = std::max(worstMiss, best);
            }
        }
        const auto unmapped = std::count(flowToMech.begin(), flowToMech.end(), -1);
        if (unmapped > 0) {
            OPM_THROW(std::runtime_error, std::to_string(unmapped)
                      + " flow cells have no cell in the thin-collapsed mechanics grid");
        }
        map_ = MechFlowMap(std::move(flowToMech), std::move(weight), mechView.size(0));
        mechCellOfMechCart_ = std::move(coarseMech);

        std::ostringstream os;
        os << "Mechanics grid: " << mechView.size(0) << " cells for " << numFlowCells
           << " flow cells; refined cells matched within " << worstMiss << " m";
        OpmLog::info(os.str());
    }

    /// 1 for every cell whose top and bottom differ at some corner.
    static std::vector<int> withVolume(const std::array<int,3>& dims,
                                       const std::vector<double>& zcorn)
    {
        const std::size_t nx = dims[0], ny = dims[1], nz = dims[2];
        std::vector<int> active(nx*ny*nz, 0);
        for (std::size_t k = 0; k < nz; ++k) {
            for (std::size_t j = 0; j < ny; ++j) {
                for (std::size_t i = 0; i < nx; ++i) {
                    for (std::size_t c = 0; c < 4; ++c) {
                        const std::size_t top = 2*i + (c & 1) + 2*nx*(2*j + (c >> 1)) + 4*nx*ny*(2*k);
                        if (zcorn[top + 4*nx*ny] != zcorn[top]) {
                            active[i + nx*(j + ny*k)] = 1;
                            break;
                        }
                    }
                }
            }
        }
        return active;
    }

    /// `partOfMechCartesian` gives the rank of each cell of the mechanics
    /// grid's own Cartesian space.
    void distribute(const std::vector<int>& partOfMechCartesian)
    {
        std::vector<int> parts;
        if (grid_->comm().rank() == 0) {
            const Dune::CartesianIndexMapper<Dune::CpGrid> mapper(*grid_);
            parts.resize(grid_->leafGridView().size(0));
            for (std::size_t c = 0; c < parts.size(); ++c) {
                parts[c] = partOfMechCartesian[mapper.cartesianIndex(static_cast<int>(c))];
            }
        }
        grid_->loadBalance(parts, /*ownersFirst*/ false, /*addCornerCells*/ true,
                           /*overlapLayers*/ 1);
    }

    /// `toMech` takes a fine Cartesian index to the mechanics grid's own
    /// Cartesian index: the coarse index for the grdecl route, the block's
    /// anchor for the merge, which keeps the fine index space.
    template <class FlowGrid, class CartesianMapper, class ToMech>
    void buildMap(const ToMech& toMech, const std::array<int,3>& mechDims,
                  const FlowGrid& flowGrid, const CartesianMapper& flowCartesianMapper)
    {
        const Dune::CartesianIndexMapper<Dune::CpGrid> mechMapper(*grid_);
        const auto& mechView = grid_->leafGridView();
        std::vector<int> cartesianToMech(static_cast<std::size_t>(mechDims[0])
                                         * mechDims[1]*mechDims[2], -1);
        for (const auto& cell : elements(mechView)) {
            const int idx = mechView.indexSet().index(cell);
            cartesianToMech[mechMapper.cartesianIndex(idx)] = idx;
        }

        const auto& flowView = flowGrid.leafGridView();

        const std::size_t numFlowCells = flowView.size(0);
        std::vector<int> flowToMech(numFlowCells, -1);
        std::vector<double> weight(numFlowCells, 0.0);

        std::size_t unmappedInterior = 0, unmappedOverlap = 0;
        for (const auto& cell : elements(flowView)) {
            const int f = flowView.indexSet().index(cell);
            weight[f] = cell.geometry().volume();
            const int coarseCartesian = toMech[flowCartesianMapper.cartesianIndex(f)];
            const int m = (coarseCartesian < 0) ? -1 : cartesianToMech[coarseCartesian];
            flowToMech[f] = m;
            if (m < 0) {
                if (cell.partitionType() == Dune::InteriorEntity) {
                    ++unmappedInterior;
                } else {
                    ++unmappedOverlap;
                }
            }
        }
        if (unmappedInterior > 0) {
            // In parallel this means the partitions disagree: an interior flow
            // cell must always be on the rank that owns its mechanics cell.
            OPM_THROW(std::runtime_error,
                      std::to_string(unmappedInterior) + " interior flow cells have no "
                      "mechanics cell on this rank; the mechanics grid must cover the flow "
                      "grid, and the two partitions must correspond.");
        }
        if (unmappedOverlap > 0) {
            OpmLog::info("Mechanics grid: " + std::to_string(unmappedOverlap)
                         + " overlap flow cells have no local mechanics cell; they take no "
                         "part in the mechanics.");
        }

        map_ = MechFlowMap(std::move(flowToMech), std::move(weight), mechView.size(0));

        std::ostringstream os;
        os << "Mechanics grid: " << mechView.size(0) << " cells for " << numFlowCells
           << " flow cells on this rank";
        OpmLog::info(os.str());
    }

    std::unique_ptr<Dune::CpGrid> grid_;
    MechFlowMap map_;
    bool merged_{false};
    bool fromCornerPoints_{false};
    std::array<int,3> fineDims_{};
    std::vector<int> fineToMech_;
    std::vector<int> mechCellOfMechCart_;
    std::vector<Coarsening::Block> blocks_;
};

} // namespace Opm

#endif // OPM_MECH_GRID_CONTEXT_HPP
