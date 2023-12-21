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
#ifndef BOUNDARYUTILS_HH
#define BOUNDARYUTILS_HH
#include <cstring>
#include <iostream>
#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>
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

        Opm::FaceDir::DirEnum faceToFaceDir(int insideFaceIdx){
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

        int faceDirToFace(Opm::FaceDir::DirEnum faceDirection){
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

        std::array<int,4> faceDirToNodes(Opm::FaceDir::DirEnum faceDirection){
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


        template<class BCConfig, class GvType, class CartMapperType>
        void nodesAtBoundary(std::vector<std::tuple<size_t,MechBCValue>>& bc_nodes,
                             const BCConfig& bcconfigs,
                             const BCProp& bcprops,
                             const GvType& gv,
                             const CartMapperType& cartesianIndexMapper
            ){
            static const int dim = 3;
            if (bcprops.size() > 0) {
                //nonTrivialBoundaryConditions_ = true;
                //auto gv = grid.leafGridView();
                size_t numCartDof = cartesianIndexMapper.cartesianSize();
                unsigned numElems = gv.size(/*codim=*/0);
                std::vector<int> cartesianToCompressedElemIdx(numCartDof, -1);
                for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx){
                    cartesianToCompressedElemIdx[cartesianIndexMapper.cartesianIndex(elemIdx)] = elemIdx;
                }
                std::array<int, 3> cartdim = cartesianIndexMapper.cartesianDimensions();
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
                                            auto elemIdx = cartesianToCompressedElemIdx[cartindex];
                                            if (elemIdx>-1){
                                                effected_cells.insert(elemIdx);
                                            }
                                        }
                                    }
                                }

                                MechBCValue bcval = *bcprop.mechbcvalue;
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
                                throw std::logic_error("invalid type for BC. Use FREE or RATE");
                            }
                        }
                    }
                }
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
            std::sort(bc_nodes.begin(), bc_nodes.end(), compare); // {1 1 2 3 4 4 5}
            auto last = std::unique(bc_nodes.begin(), bc_nodes.end(), isequal);
            // v now holds {1 2 3 4 5 x x}, where 'x' is indeterminate
            bc_nodes.erase(last, bc_nodes.end());
        }
    }
}

#endif
