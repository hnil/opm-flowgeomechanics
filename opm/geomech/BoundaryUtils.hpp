// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
#ifndef OPM_BOUNDARY_UTILS_HPP
#define OPM_BOUNDARY_UTILS_HPP
#include <cassert>
#include <cstring>
#include <iostream>
#include <set>
#include <stdexcept>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <opm/common/ErrorMacros.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/input/eclipse/Schedule/BCState.hpp>
#include <opm/geomech/NodeConstraint.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <type_traits>
#include <vector>
namespace Opm{
    namespace Elasticity{
        template<int dimension>
        unsigned cartesianIndex(const std::array<int,dimension>& coords,
                                const std::array<int, dimension>& cartesianDimensions
            )
        {
            unsigned cartIndex = coords[0];
            int factor = cartesianDimensions[0];
            for (unsigned i = 1; i < dimension; ++i) {
                cartIndex += coords[i]*factor;
                factor *= cartesianDimensions[i];
            }

            return cartIndex;
        }

        inline Opm::FaceDir::DirEnum faceToFaceDir(int insideFaceIdx){
            Opm::FaceDir::DirEnum faceDirection;
            switch (insideFaceIdx) {
            case 0:
                faceDirection = FaceDir::XMinus;
                break;
            case 1:
                faceDirection = FaceDir::XPlus;
                break;
            case 2:
                faceDirection = FaceDir::YMinus;
                break;
            case 3:
                faceDirection = FaceDir::YPlus;
                break;
            case 4:
                faceDirection = FaceDir::ZMinus;
                break;
            case 5:
                faceDirection = FaceDir::ZPlus;
                break;
            default:
                OPM_THROW(std::logic_error,
                          "Internal error in initialization of aquifer.");
            }
            return faceDirection;
        }

        inline int faceDirToFace(Opm::FaceDir::DirEnum faceDirection){
            int face;
            switch (faceDirection) {
                case FaceDir::XMinus:
                    face = 0;
                    break;
                case FaceDir::XPlus:
                    face = 1;
                    break;
                case FaceDir::YMinus:
                    face = 2;
                    break;
                case FaceDir::YPlus:
                    face = 3;
                    break;
                case FaceDir::ZMinus:
                    face = 4;
                    break;
                case FaceDir::ZPlus:
                    face = 5;
                    break;
                default:
                    OPM_THROW(std::logic_error,
                        "Internal error in initialization of aquifer.");
            }
            return face;
        }

        inline std::array<int,4> faceDirToNodes(Opm::FaceDir::DirEnum faceDirection){
            std::array<int,4> face;
            switch (faceDirection) {
                case FaceDir::XMinus:
                    face = {0,4,6,2};
                    break;
                case FaceDir::XPlus:
                    face = {1,3,7,5};
                    break;
                case FaceDir::YMinus:
                    face = {0,1,5,4};
                    break;
                case FaceDir::YPlus:
                    face = {3,2,6,7};
                    break;
                case FaceDir::ZMinus:
                    face = {0,2,3,1};
                    break;
                case FaceDir::ZPlus:
                    face = {4,5,7,6};
                    break;
                default:
                    OPM_THROW(std::logic_error,
                        "Internal error in initialization of aquifer.");
            }
            return face;
        }


        /// Constrain the nodes on the grid's outer boundary. Used for a
        /// mechanics grid of its own, whose boundary is not the one BCCON
        /// describes. Rollers still remove every rigid mode of a box, and
        /// leaving the top free lets the overburden move vertically.
        enum class OuterBoundary {
            Fixed,          ///< all three components on every outer face
            Roller,         ///< the normal component on every outer face
            RollerFreeTop,  ///< as Roller, but the Z- face is left free
        };

        /// Boundary faces at a node whose normals differ by less than this
        /// angle form one patch and give one roller direction.
        inline constexpr double rollerPatchAngle = 30.0 * M_PI / 180.0;

        /// Roller directions at one node from the area normals of its boundary
        /// faces: one per patch, largest patch first, orthonormalised.
        inline NodeConstraint rollerConstraint(std::size_t node,
                                               const std::vector<std::array<double,3>>& normals)
        {
            using Vec = std::array<double,3>;
            const auto dot = [](const Vec& a, const Vec& b) {
                return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
            };
            std::vector<Vec> patches;
            for (const auto& n : normals) {
                const double len = std::sqrt(dot(n, n));
                if (len == 0.0) {
                    continue;
                }
                auto p = std::find_if(patches.begin(), patches.end(), [&](const Vec& q) {
                    return dot(n, q) > std::cos(rollerPatchAngle) * len * std::sqrt(dot(q, q));
                });
                if (p == patches.end()) {
                    patches.push_back(n);
                } else {
                    for (int d = 0; d < 3; ++d) {
                        (*p)[d] += n[d];
                    }
                }
            }
            std::stable_sort(patches.begin(), patches.end(),
                             [&dot](const Vec& a, const Vec& b) { return dot(a, a) > dot(b, b); });
            std::vector<Vec> dirs;
            for (auto d : patches) {
                const double len = std::sqrt(dot(d, d));
                for (auto& x : d) {
                    x /= len;
                }
                for (const auto& q : dirs) {
                    const double proj = dot(d, q);
                    for (int k = 0; k < 3; ++k) {
                        d[k] -= proj * q[k];
                    }
                }
                const double rest = std::sqrt(dot(d, d));
                // a patch nearly in the span of the larger ones adds no direction
                if (rest < std::sin(rollerPatchAngle)) {
                    continue;
                }
                for (auto& x : d) {
                    x /= rest;
                    if (std::abs(x) < 1e-12) {
                        x = 0.0; // round-off tilt of an axis-aligned face
                    }
                }
                const double len2 = std::sqrt(dot(d, d));
                for (auto& x : d) {
                    x /= len2;
                }
                dirs.push_back(d);
                if (dirs.size() == 3) {
                    break;
                }
            }

            NodeConstraint c;
            c.node = node;
            const auto axisOf = [](const Vec& d) {
                int nonzero = 0, axis = -1;
                for (int k = 0; k < 3; ++k) {
                    if (d[k] != 0.0) {
                        ++nonzero;
                        axis = k;
                    }
                }
                return (nonzero == 1 && std::abs(d[axis]) == 1.0) ? axis : -1;
            };
            // a full set or unit axes need no rotation (roller values are zero)
            if (dirs.size() == 3) {
                c.fixed = {true, true, true};
                return c;
            }
            if (std::all_of(dirs.begin(), dirs.end(), [&](const Vec& d) { return axisOf(d) >= 0; })) {
                for (const auto& d : dirs) {
                    c.fixed[axisOf(d)] = true;
                }
                return c;
            }
            c.frame[0] = dirs[0];
            if (dirs.size() == 1) {
                int least = 0;
                for (int k = 1; k < 3; ++k) {
                    if (std::abs(dirs[0][k]) < std::abs(dirs[0][least])) {
                        least = k;
                    }
                }
                Vec t{};
                t[least] = 1.0;
                const double proj = dot(t, dirs[0]);
                for (int k = 0; k < 3; ++k) {
                    t[k] -= proj * dirs[0][k];
                }
                const double len = std::sqrt(dot(t, t));
                for (auto& x : t) {
                    x /= len;
                }
                c.frame[1] = t;
            } else {
                c.frame[1] = dirs[1];
            }
            const auto& a = c.frame[0];
            const auto& b = c.frame[1];
            c.frame[2] = {a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]};
            for (std::size_t i = 0; i < dirs.size(); ++i) {
                c.fixed[i] = true;
            }
            return c;
        }

        template<class GridView>
        struct NormalsHandle {
            const GridView& gv;
            const std::vector<std::vector<double>>& own;
            std::vector<std::vector<double>>& recv;
            using DataType = double;
            bool contains(int, int codim) const { return codim == 3; }
            bool fixedSize(int, int) const { return false; }
            template <class E> std::size_t size(const E& e) const
            { return own[gv.indexSet().index(e)].size(); }
            template <class B, class E> void gather(B& b, const E& e) const
            { for (const double x : own[gv.indexSet().index(e)]) b.write(x); }
            template <class B, class E> void scatter(B& b, const E& e, std::size_t n)
            {
                auto& r = recv[gv.indexSet().index(e)];
                for (std::size_t k = 0; k < n; ++k) { double x; b.read(x); r.push_back(x); }
            }
        };

        /// Boundary faces are found from interior cells only, so a vertex shared
        /// between ranks gets its normals from every rank; all ranks must then
        /// constrain it alike, as the elimination is symmetric.
        template<class Grid>
        void exchangeBoundaryNormals(const Grid& grid,
            std::vector<std::pair<std::size_t, std::array<double,3>>>& nodeNormals)
        {
            const auto gv = grid.leafGridView();
            using GridView = std::decay_t<decltype(gv)>;
            std::vector<std::vector<double>> own(gv.size(3)), recv(gv.size(3));
            for (const auto& [node, n] : nodeNormals) {
                own[node].insert(own[node].end(), n.begin(), n.end());
            }
            NormalsHandle<GridView> handle{gv, own, recv};
            grid.communicate(handle, Dune::All_All_Interface, Dune::ForwardCommunication);
            for (std::size_t v = 0; v < recv.size(); ++v) {
                for (std::size_t k = 0; k + 2 < recv[v].size(); k += 3) {
                    nodeNormals.emplace_back(v, std::array<double,3>{recv[v][k], recv[v][k+1], recv[v][k+2]});
                }
            }
        }

        /// Every node of a boundary face, not just the cell's eight corners: a
        /// cell merged from several others has faces with more nodes than
        /// that, and leaving them free would make the boundary soft.
        template<class Grid>
        void nodesAtOuterBoundary(std::vector<NodeConstraint>& constraints,
                                  const Grid& grid,
                                  const OuterBoundary mode)
        {
            const auto& gv = grid.leafGridView();
            std::vector<std::pair<std::size_t, std::array<double,3>>> nodeNormals;
            for (const auto& cell : elements(gv)) {
                // In parallel the outermost overlap cells have faces to cells
                // this rank does not hold, which look like domain boundary and
                // would constrain the interior.
                if (cell.partitionType() != Dune::InteriorEntity) {
                    continue;
                }
                const int cellIdx = gv.indexSet().index(cell);
                for (int f = 0; f < grid.numCellFaces(cellIdx); ++f) {
                    const int face = grid.cellFace(cellIdx, f);
                    const int other = (grid.faceCell(face, 0) == cellIdx)
                        ? grid.faceCell(face, 1) : grid.faceCell(face, 0);
                    if (other >= 0) {
                        continue;                      // not the domain boundary
                    }
                    const auto normal = grid.faceAreaNormalEcl(face);
                    int axis = 0;
                    for (int d = 1; d < 3; ++d) {
                        if (std::abs(normal[d]) > std::abs(normal[axis])) {
                            axis = d;
                        }
                    }
                    if (mode == OuterBoundary::RollerFreeTop && axis == 2
                        && grid.faceCentroid(face)[2] < grid.cellCentroid(cellIdx)[2]) {
                        continue;                      // the top of the model
                    }
                    for (int v = 0; v < grid.numFaceVertices(face); ++v) {
                        nodeNormals.emplace_back(grid.faceVertex(face, v),
                                                 std::array<double,3>{normal[0], normal[1], normal[2]});
                    }
                }
            }
            const bool parallel = gv.comm().size() > 1;
            if (parallel) {
                exchangeBoundaryNormals(grid, nodeNormals);
            }
            std::stable_sort(nodeNormals.begin(), nodeNormals.end(),
                             [](const auto& a, const auto& b) { return a.first < b.first; });
            constraints.clear();
            std::vector<std::array<double,3>> normals;
            for (std::size_t i = 0; i < nodeNormals.size();) {
                const auto node = nodeNormals[i].first;
                normals.clear();
                for (; i < nodeNormals.size() && nodeNormals[i].first == node; ++i) {
                    normals.push_back(nodeNormals[i].second);
                }
                if (parallel) {
                    // received normals arrive in rank order; the roller frame must not
                    std::sort(normals.begin(), normals.end());
                }
                if (mode == OuterBoundary::Fixed) {
                    NodeConstraint c;
                    c.node = node;
                    c.fixed = {true, true, true};
                    constraints.push_back(c);
                } else {
                    constraints.push_back(rollerConstraint(node, normals));
                }
            }
        }

        template<class BCConfig, class GvType, class CartMapperType>
        void nodesAtBoundary(std::vector<std::tuple<size_t,MechBCValue>>& bc_nodes,
                             const BCConfig& bcconfigs,
                             const BCState& bcprops,
                             const GvType& gv,
                             const CartMapperType& cartesianIndexMapper
            ){
            if (bcprops.size() == 0) {
                return;
            }
            size_t numCartDof = cartesianIndexMapper.cartesianSize();
            unsigned numElems = gv.size(/*codim=*/0);
            std::vector<int> cartesianToCompressedElemIdx(numCartDof, -1);
            for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx){
                cartesianToCompressedElemIdx[cartesianIndexMapper.cartesianIndex(elemIdx)] = elemIdx;
            }
            nodesAtBoundary(bc_nodes, bcconfigs, bcprops, gv,
                            cartesianIndexMapper.cartesianDimensions(),
                            [&cartesianToCompressedElemIdx](int cartindex, Opm::FaceDir::DirEnum) {
                                return cartesianToCompressedElemIdx[cartindex];
                            });
        }

        //! BCCON boxes in a Cartesian space of `cartdim`; `cellOfFace` gives the
        //! grid cell carrying that cell's face in the BC direction, or -1.
        template<class BCConfig, class GvType, class CellOfFace>
        void nodesAtBoundary(std::vector<std::tuple<size_t,MechBCValue>>& bc_nodes,
                             const BCConfig& bcconfigs,
                             const BCState& bcprops,
                             const GvType& gv,
                             const std::array<int, 3>& cartdim,
                             const CellOfFace& cellOfFace
            ){
            static const int dim = 3;
            if (bcprops.size() > 0) {
                for (const auto& bcconfig : bcconfigs) {
                    for (const auto& bcprop : bcprops) {


                        if(bcprop.index == bcconfig.index){
                            // double search since structure is strange
                            const auto& bcface = bcconfig;

                            if((bcface.i1 < 0) || (bcface.j1<0) || (bcface.k1<0)){
                                throw std::logic_error("Lower range of BC wrong");
                            }
                            if( (bcface.i2 > cartdim[0]) || (bcface.j2> cartdim[1]) || (bcface.k2 > cartdim[2])){
                                throw std::logic_error("Upper range of BC wrong");
                            }

                            const auto& type = bcprop.bcmechtype;
                            if (type == Opm::BCMECHType::FREE) {
                                // do nothing
                            }else if (type == Opm::BCMECHType::FIXED) {
                                std::set<size_t> effected_cells;
                                for (int i = bcface.i1; i <= bcface.i2; ++i) {
                                    for (int j = bcface.j1; j <= bcface.j2; ++j) {
                                        for (int k = bcface.k1; k <= bcface.k2; ++k) {


                                            std::array<int, 3> tmp = {i,j,k};
                                            int cartindex =
                                                cartesianIndex<3>(tmp,cartdim);
                                            auto elemIdx = cellOfFace(cartindex, bcface.dir);
                                            if (elemIdx>-1){
                                                effected_cells.insert(elemIdx);
                                            }
                                        }
                                    }
                                }

                                MechBCValue bcval = bcprop.mechbcvalue;
                                for(const auto& cell:elements(gv)){
                                    auto index = gv.indexSet().index(cell);
                                    auto it = effected_cells.find(index);
                                    if(!(it == effected_cells.end())){
                                        // fix all noted for now
                                        std::array<int,4> nodes = faceDirToNodes(bcface.dir);
                                        for(auto nind: nodes){
                                            auto global_ind = gv.indexSet().subIndex(cell,nind,dim);
                                            bc_nodes.push_back(
                                                std::tuple<size_t,MechBCValue>(global_ind,bcval)
                                                );
                                        }
                                    }
                                }
                            } else {
                                throw std::logic_error("invalid mechanical type for BC. Use FREE or FIXED in BCMECH");
                            }
                        }
                    }
                }
            }else{
                return;
            }
            auto compare = [](std::tuple<size_t,MechBCValue> const &t1,
                              std::tuple<size_t,MechBCValue> const &t2)
            {
                return std::get<0>(t1) < std::get<0>(t2);
            };
            auto isequal = [](std::tuple<size_t,MechBCValue> const &t1,
                              std::tuple<size_t,MechBCValue> const &t2)
            {
                return std::get<0>(t1) == std::get<0>(t2);
            };
            if (bc_nodes.empty()) {
                // In parallel a rank may hold no cell on a BCCON face.
                return;
            }
            //TODO: make unique with masks removing boundary set on same dofs mayr romve masks
            const std::vector<std::tuple<size_t,MechBCValue>> org_bc_nodes = bc_nodes;
            std::sort(bc_nodes.begin(), bc_nodes.end(), compare); // {1 1 2 3 4 4 5}
            auto last = std::unique(bc_nodes.begin(), bc_nodes.end(), isequal);
            // v now holds {1 2 3 4 5 x x}, where 'x' is indeterminate
            bc_nodes.erase(last, bc_nodes.end());

            // make index
            auto max_elem = std::max_element(bc_nodes.begin(), bc_nodes.end(),
                                          [](const auto& a, const auto& b){
                                              return std::get<0>(a) < std::get<0>(b);
                                          });
            int maxnode = std::get<0>(*max_elem);
            std::vector<int> mymap(maxnode+1,-1);
            for(int i=0; i< int(bc_nodes.size()); ++i){
                mymap[std::get<0>(bc_nodes[i])] = i;
            }

            for(int i=0; i< int(org_bc_nodes.size()); ++i){
                auto node = std::get<0>(org_bc_nodes[i]);
                assert(node<mymap.size());
                int index = mymap[node];
                assert(index>-1);
                auto mask = std::get<1>(org_bc_nodes[i]).fixeddir;
                for(int mdim=0; mdim<3; mdim++){
                    if(mask[mdim]){
                        assert(index < int(bc_nodes.size()));
                        auto& bc_node = bc_nodes[index];
                        assert(node == std::get<0>(bc_node));
                        std::get<1>(bc_node).fixeddir[mdim] = true;
                        assert(std::get<1>(bc_nodes[index]).fixeddir[mdim] == true);
                    }
                }
            }
        }
    }
}

#endif
