/*
  Copyright 2025 Equinor ASA.

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

#pragma once
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/version.hh>

#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>



#include <opm/common/OpmLog/OpmLog.hpp>

#include <cstring>
#include <iostream>
#include <sstream>
#include <unistd.h>
#include <unistd.h>

#include "DuneUtilities.hpp"
namespace Opm {
  void cellCellCommunication(const Dune::CpGrid& grid, Dune::MPIHelper& mpihelper){
    Dune::Communication<MPI_Comm> world_comm = Dune::MPIHelper::getCommunication();
    const auto& gv = grid.leafGridView();
    auto indexset = gv.indexSet();
    auto gidSet = gv.grid().globalIdSet();
    using AttributeSet = Dune::OwnerOverlapCopyAttributeSet;
    using ParallelIndexSet = Dune::OwnerOverlapCopyCommunication<int, int>::ParallelIndexSet;
    using RemoteIndices = Dune::OwnerOverlapCopyCommunication<int, int>::RemoteIndices;
    ParallelIndexSet cell_indexset;
    RemoteIndices remote_indices(cell_indexset, cell_indexset, world_comm);
    cell_indexset.beginResize();
    for (auto& elem : Dune::elements(gv)) {
      auto index = elem.index();
      auto gid = gidSet.id(elem);
      if (elem.partitionType() == Dune::PartitionType::InteriorEntity) {
	cell_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::owner, true));
      } else if (elem.partitionType() == Dune::PartitionType::OverlapEntity) {
	cell_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::copy, true));
      } else if (elem.partitionType() == Dune::PartitionType::BorderEntity) {
      } else if (elem.partitionType() == Dune::PartitionType::FrontEntity) {
      } else if (elem.partitionType() == Dune::PartitionType::GhostEntity) {
      } else {
	std::cout << "Unknown partition type" << std::endl;
      }
    }

    cell_indexset.endResize();
    remote_indices.rebuild<false>();
    std::vector<int> myrank(gv.size(0), world_comm.rank());
        

    Dune::Interface cominterface;
    using OwnerSet = Dune::OwnerOverlapCopyCommunication<int, int>::OwnerSet;
    using AllSet =  Dune::OwnerOverlapCopyCommunication<int, int>::AllSet;
    OwnerSet soureFlags;
    AllSet destFlags;


    cominterface.build(remote_indices, soureFlags, destFlags);
    Dune::BufferedCommunicator cell_cell_comm;
    using Vector = std::vector<int>;
    cell_cell_comm.template build<Vector>(cominterface);
    cell_cell_comm.template forward<Dune::CopyGatherScatter<Vector>>(myrank, myrank);
    std::vector<int> numpros(gv.size(0), 1);

    // all to all
    Dune::Interface all_cominterface;
    all_cominterface.build(remote_indices, destFlags, destFlags);
    Dune::BufferedCommunicator all_cell_cell_comm;
    all_cell_cell_comm.template build<Vector>(all_cominterface);
    all_cell_cell_comm.template forward<Dune::FreeAddGatherScatter<Vector>>(numpros, numpros);

    using VarVector = std::vector<std::vector<int>>;
    VarVector allranks(gv.size(0));
    for(int i=0; i<gv.size(0); ++i){
      allranks[i].resize(numpros[i]+1,9999);
      allranks[i][0]=world_comm.rank();
    }

    world_comm.barrier();
    Dune::BufferedCommunicator varvec_cell_cell_comm;
    varvec_cell_cell_comm.build(allranks, allranks, all_cominterface);
    varvec_cell_cell_comm.template forward<Dune::VariableVectorAdderGatherScatter>(allranks, allranks);

    cell_cell_comm.free();

    world_comm.barrier();
        
    for (int rank = 0; rank < mpihelper.size(); ++rank) {
      if (rank == mpihelper.rank()) {
	std::cout << "Grid size: " << gv.size(0) << " on rank " << mpihelper.rank() << std::endl;  
	for (auto& elem : Dune::elements(gv)) {
	  auto lindex = elem.index();
	  auto gid = gidSet.id(elem);
	  std::cout << "Element global: id " << gid << " local " << lindex << " type " << elem.partitionType()
		    << " owner rank " << myrank[lindex] 
		    << " num pros " << numpros[lindex] 
		    << " all ranks ";
	  for(auto& r:allranks[lindex]){
	    std::cout << r << " ";
	  } 
	  std::cout << std::endl;
	}
      }
      // for bliding communicator see ParallelIstlInformation. (or Dune::OwnerOverlapCopyCommunication)
      usleep(100);
      world_comm.barrier();
    }
  }

  void vertexVertexCommunication(const Dune::CpGrid& grid, Dune::MPIHelper& mpihelper){
    Dune::Communication<MPI_Comm> world_comm = Dune::MPIHelper::getCommunication();
    const auto& gv = grid.leafGridView();
    auto indexset = gv.indexSet();
    auto gidSet = gv.grid().globalIdSet();
    using AttributeSet = Dune::OwnerOverlapCopyAttributeSet;
    using ParallelIndexSet = Dune::OwnerOverlapCopyCommunication<int, int>::ParallelIndexSet;
    using RemoteIndices = Dune::OwnerOverlapCopyCommunication<int, int>::RemoteIndices;
    ParallelIndexSet vertex_indexset;
    RemoteIndices remote_indices(vertex_indexset, vertex_indexset, world_comm);
    vertex_indexset.beginResize();
    for (auto& vertex : Dune::vertices(gv)) {
      auto index = vertex.index();
      auto gid = gidSet.id(vertex);
      if (vertex.partitionType() == Dune::PartitionType::InteriorEntity) {
	vertex_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::owner, true));
      } else if (vertex.partitionType() == Dune::PartitionType::OverlapEntity) {
	vertex_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::copy, true));
      } else if (vertex.partitionType() == Dune::PartitionType::BorderEntity) {
      } else if (vertex.partitionType() == Dune::PartitionType::FrontEntity) {
      } else if (vertex.partitionType() == Dune::PartitionType::GhostEntity) {
      } else {
	std::cout << "Unknown partition type" << std::endl;
      }
    }

    vertex_indexset.endResize();
    remote_indices.rebuild<false>();
    std::vector<int> myrank(gv.size(3), world_comm.rank());
        

    Dune::Interface cominterface;
    using OwnerSet = Dune::OwnerOverlapCopyCommunication<int, int>::OwnerSet;
    using AllSet =  Dune::OwnerOverlapCopyCommunication<int, int>::AllSet;
    OwnerSet soureFlags;
    AllSet destFlags;


    cominterface.build(remote_indices, soureFlags, destFlags);
    Dune::BufferedCommunicator vertex_vertex_comm;
    using Vector = std::vector<int>;
    vertex_vertex_comm.template build<Vector>(cominterface);
    vertex_vertex_comm.template forward<Dune::CopyGatherScatter<Vector>>(myrank, myrank);
    std::vector<int> numpros(gv.size(3), 1);

    // all to all
    Dune::Interface all_cominterface;
    all_cominterface.build(remote_indices, destFlags, destFlags);
    Dune::BufferedCommunicator all_vertex_vertex_comm;
    all_vertex_vertex_comm.template build<Vector>(all_cominterface);
    all_vertex_vertex_comm.template forward<Dune::FreeAddGatherScatter<Vector>>(numpros, numpros);

    using VarVector = std::vector<std::vector<int>>;
    VarVector allranks(gv.size(3));
    for(int i=0; i<gv.size(3); ++i){
      allranks[i].resize(numpros[i]+1,9999);
      allranks[i][0]=world_comm.rank();
    }

    world_comm.barrier();
    Dune::BufferedCommunicator varvec_vertex_vertex_comm;
    varvec_vertex_vertex_comm.build(allranks, allranks, all_cominterface);
    varvec_vertex_vertex_comm.template forward<Dune::VariableVectorAdderGatherScatter>(allranks, allranks);

    vertex_vertex_comm.free();

    world_comm.barrier();
        
    for (int rank = 0; rank < mpihelper.size(); ++rank) {
      if (rank == mpihelper.rank()) {
	std::cout << "Grid size: " << gv.size(0) << " on rank " << mpihelper.rank() << std::endl;  
	for (auto& elem : Dune::vertices(gv)) {
	  auto lindex = elem.index();
	  auto gid = gidSet.id(elem);
	  std::cout << "Vertex global: id " << gid << " local " << lindex << " type " << elem.partitionType()
		    << " owner rank " << myrank[lindex] 
		    << " num pros " << numpros[lindex] 
		    << " all ranks ";
	  for(auto& r:allranks[lindex]){
	    std::cout << r << " ";
	  } 
	  std::cout << std::endl;
	}
      }
      // for bliding communicator see ParallelIstlInformation. (or Dune::OwnerOverlapCopyCommunication)
      usleep(100);
      world_comm.barrier();
    }
  }

  template <int codim>
  void
  entityEntityCommunication(const Dune::CpGrid& grid, Dune::MPIHelper& mpihelper)
  {
      Dune::Codim<codim> mycodim;
      Dune::Communication<MPI_Comm> world_comm = Dune::MPIHelper::getCommunication();
      const auto& gv = grid.leafGridView();
      auto indexset = gv.indexSet();
      auto gidSet = gv.grid().globalIdSet();
      using AttributeSet = Dune::OwnerOverlapCopyAttributeSet;
      using ParallelIndexSet = Dune::OwnerOverlapCopyCommunication<int, int>::ParallelIndexSet;
      using RemoteIndices = Dune::OwnerOverlapCopyCommunication<int, int>::RemoteIndices;
      ParallelIndexSet entity_indexset;
      RemoteIndices remote_indices(entity_indexset, entity_indexset, world_comm);
      entity_indexset.beginResize();
      for (auto& entity : Dune::entities(gv, mycodim)) {
          auto index = entity.index();
          auto gid = gidSet.id(entity);
          if (entity.partitionType() == Dune::PartitionType::InteriorEntity) {
              entity_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::owner, true));
          } else if (entity.partitionType() == Dune::PartitionType::OverlapEntity) {
              entity_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::copy, true));
          } else if (entity.partitionType() == Dune::PartitionType::BorderEntity) {
            entity_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::copy, true));
          } else if (entity.partitionType() == Dune::PartitionType::FrontEntity) {
            entity_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::copy, true));
          } else if (entity.partitionType() == Dune::PartitionType::GhostEntity) {
            entity_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::copy, true));
          } else {
              std::cout << "Unknown partition type" << std::endl;
          }
      }

      entity_indexset.endResize();
      remote_indices.rebuild<false>();
      std::vector<int> myrank(gv.size(codim), world_comm.rank());


      Dune::Interface cominterface;
      using OwnerSet = Dune::OwnerOverlapCopyCommunication<int, int>::
          OwnerSet;
      using AllSet = Dune::OwnerOverlapCopyCommunication<int, int>::AllSet;
      OwnerSet soureFlags;
      AllSet destFlags;


      cominterface.build(remote_indices, soureFlags, destFlags);
      Dune::BufferedCommunicator entity_entity_comm;
      using Vector = std::vector<int>;
      entity_entity_comm.template build<Vector>(cominterface);
      entity_entity_comm.template forward<Dune::CopyGatherScatter<Vector>>(myrank, myrank);
      std::vector<int> numpros(gv.size(codim), 1);

      // all to all
      Dune::Interface all_cominterface;
      all_cominterface.build(remote_indices, destFlags, destFlags);
      Dune::BufferedCommunicator all_entity_entity_comm;
      all_entity_entity_comm.template build<Vector>(all_cominterface);
      all_entity_entity_comm.template forward<Dune::FreeAddGatherScatter<Vector>>(numpros, numpros);

      using VarVector = std::vector<std::vector<int>>;
      VarVector allranks(gv.size(3));
      for (int i = 0; i < gv.size(codim); ++i) {
          allranks[i].resize(numpros[i] + 1, 9999);
          allranks[i][0] = world_comm.rank();
      }

      world_comm.barrier();
      Dune::BufferedCommunicator varvec_entity_entity_comm;
      varvec_entity_entity_comm.build(allranks, allranks, all_cominterface);
      varvec_entity_entity_comm.template forward<Dune::VariableVectorAdderGatherScatter>(allranks, allranks);

      entity_entity_comm.free();

      world_comm.barrier();

      for (int rank = 0; rank < mpihelper.size(); ++rank) {
          if (rank == mpihelper.rank()) {
              std::cout << "Grid size: " << gv.size(0) << " on rank " << mpihelper.rank() << std::endl;


              for (auto& elem : Dune::entities(gv, mycodim)) {
                  auto lindex = elem.index();
                  auto gid = gidSet.id(elem);
                  std::cout << "Entity<" << codim << "> global: id " << gid << " local " << lindex << " type "
                            << elem.partitionType() << " owner rank " << myrank[lindex] << " num pros "
                            << numpros[lindex] << " all ranks ";
                  for (auto& r : allranks[lindex]) {
                      std::cout << r << " ";
                  }
                  std::cout << std::endl;
              }
          }
          // for bliding communicator see ParallelIstlInformation. (or Dune::OwnerOverlapCopyCommunication)
          usleep(100);
          world_comm.barrier();
      }
  }

  template <int codim>
  Dune::OwnerOverlapCopyCommunication<int, int>::ParallelIndexSet
  makeEntityEntityCommunication(const Dune::CpGrid& grid, bool verbose = false);

  template <int codim, class Grid>
  Dune::OwnerOverlapCopyCommunication<int, int>::ParallelIndexSet makeEntityEntityCommunication(const Grid& /*grid*/, bool verbose = false){
    assert(false);// dummy to get polygrid which is not prallel to compile
    Dune::OwnerOverlapCopyCommunication<int, int>::ParallelIndexSet entity_indexset;
    return entity_indexset;
  }


  template <int codim>
  Dune::OwnerOverlapCopyCommunication<int, int>::ParallelIndexSet
  makeEntityEntityCommunication(const Dune::CpGrid& grid, bool verbose) // = false)
  {
     //first find maximum rank of entity to ensure unique owner
    Dune::Codim<codim> mycodim;
    using Vector = std::vector<int>;
    Dune::Communication<MPI_Comm> world_comm = grid.comm();
    const auto& gv = grid.leafGridView();
    auto indexset = gv.indexSet();
    auto gidSet = gv.grid().globalIdSet();
      
    Vector maxranks_tmp_ib(gv.size(codim), world_comm.rank());
    using GridType = Dune::CpGrid;
    using GridView =  typename GridType::LeafGridView;
    Dune::MaxEntityVectorVectorDataHandle<GridView,Vector> datahandle_ib(maxranks_tmp_ib, gv, codim);
    gv.communicate(datahandle_ib, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);
    const auto& allranks_ib = datahandle_ib.all_data();
    int myrank = world_comm.rank();
    std::vector<int> maxrank_ib(gv.size(codim), myrank);
    for(size_t j=0; j<allranks_ib.size(); ++j){
         const auto& all = allranks_ib[j];  
        for(size_t i=0; i<all.size(); ++i){
            maxrank_ib[j] = std::max(maxrank_ib[j], all[i]);
        }
    }
    Vector maxranks_tmp_all(gv.size(codim), world_comm.rank());
    Dune::MaxEntityVectorVectorDataHandle<GridView,Vector> datahandle_all(maxranks_tmp_all, gv, codim);
    gv.communicate(datahandle_all, Dune::All_All_Interface, Dune::ForwardCommunication);
    const auto& allranks_all = datahandle_all.all_data();


    world_comm.barrier();
      
      using AttributeSet = Dune::OwnerOverlapCopyAttributeSet;
      using ParallelIndexSet = Dune::OwnerOverlapCopyCommunication<int, int>::ParallelIndexSet;
      using RemoteIndices = Dune::OwnerOverlapCopyCommunication<int, int>::RemoteIndices;
      ParallelIndexSet entity_indexset;
      RemoteIndices remote_indices(entity_indexset, entity_indexset, world_comm);
      entity_indexset.beginResize();
      for (auto& entity : Dune::entities(gv, mycodim)) {
          auto index = entity.index();
          auto gid = gidSet.id(entity);
          if (entity.partitionType() == Dune::PartitionType::InteriorEntity) {
              entity_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::owner, true));
          } else if (entity.partitionType() == Dune::PartitionType::OverlapEntity) {
              entity_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::copy, true));
          } else if (entity.partitionType() == Dune::PartitionType::BorderEntity) {
            if(myrank == maxrank_ib[index]){
                entity_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::owner, true));
            }else{
                entity_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::copy, true));
            }
          } else if (entity.partitionType() == Dune::PartitionType::FrontEntity) {
            entity_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::copy, true));
          } else if (entity.partitionType() == Dune::PartitionType::GhostEntity) {
                entity_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::copy, true));
          } else {
              OpmLog::warning("Unknown partition type");
          }
      }
      entity_indexset.endResize();
      remote_indices.rebuild<false>();

      // check if unique owner and that all indices has an owner
      using AllSet = Dune::OwnerOverlapCopyCommunication<int, int>::AllSet;
      AllSet allSet;
      
      using Vector = std::vector<int>;
      Vector numpros(gv.size(codim), 0);
      for(const auto& ind:entity_indexset){
          if(ind.local().attribute() == Dune::OwnerOverlapCopyAttributeSet::owner){
              numpros[ind.local().local()] = 1;
          }
      }
      Dune::Interface all_cominterface;
      all_cominterface.build(remote_indices, allSet, allSet);
      Dune::BufferedCommunicator all_entity_entity_comm;
      all_entity_entity_comm.template build<Vector>(all_cominterface);
      all_entity_entity_comm.template forward<Dune::FreeAddGatherScatter<Vector>>(numpros, numpros);
      std::map<int,int> local_global;
      std::map<int,int> global_local;
      for (auto& entity : Dune::entities(gv, mycodim)) {
        auto index = entity.index();
        auto gid = gidSet.id(entity);
        local_global[index] = gid;
        global_local[gid] = index;
      }
      bool ok = true;
      for (auto& entity : Dune::entities(gv, mycodim)) {
          auto index = entity.index();
          auto gid = gidSet.id(entity);
          auto ind = entity_indexset.at(gid);
          assert(index == ind.local().local());
          if(numpros[ind.local().local()] != 1) {
            ok = false;
           if(!verbose || true){
              std::ostringstream msg;
              msg << "Num procs " << numpros[ind.local().local()];
              msg << " Global " << ind.global();
              msg << " Local " << ind.local().local() << " Attribute " << ind.local().attribute();
              msg << " rank " << world_comm.rank();
              msg << " enity is " << entity.partitionType();
              msg << " max_rank " << maxrank_ib[ind.local().local()];
              msg << " all ranks ib ";
              for(auto other:allranks_ib[ind.local().local()]){
                  msg << other << " ";
              }
              msg << " all ranks all ";
              for(auto other:allranks_all[ind.local().local()]){
                  msg << other << " ";
              }
              msg << '\n';
              {
                Dune::Codim<0> cellcodim;
                for (auto& cell : Dune::entities(gv, cellcodim)) {
                  auto cgid = gidSet.id(cell);
                  auto cind = cell.index();
                  auto faces = grid.cellFaceRow(cind);
                  auto correct_cell =false;
                  for(int lface=0; lface < faces.size(); ++lface){
                    int face = faces[lface].index();
                    int nv = grid.numFaceVertices(face);
                    for(int lnode=0; lnode < nv; ++lnode){
                      auto node = grid.faceVertex(face, lnode);
                      if( node == index){
                        correct_cell = true;
                      }
                    }
                  }
                  if(correct_cell == true){
                   // Cell/face/node dump using only the public grid interface;
                   // the internal partition-type indicators of faces and points
                   // are not accessible from here.
                   msg << "Rank " << world_comm.rank() << " Cell global " << cgid  << " local " << cind << " partition " << cell.partitionType() << '\n';
                   for(int lface=0; lface < faces.size(); ++lface){
                    int face = faces[lface].index();
                    int nv = grid.numFaceVertices(face);
                    msg << "Face " << face << ":" ;
                    for(int lnode=0; lnode < nv; ++lnode){
                      auto node = grid.faceVertex(face, lnode);
                      auto gnode = local_global[node];
                      msg << "{" << gnode << "," << node << "}";
                    }
                    msg << '\n';
                  }
                  }

                }
              }
              OpmLog::debug(msg.str());
              ok = false;
           } else{
            ok = false;
           }
          }
      }
      world_comm.barrier();
      int not_ok = !ok;
      not_ok = world_comm.sum(not_ok);
      if(not_ok > 0){
          OpmLog::error("Owner is not a partition of unity");
          DUNE_THROW(Dune::Exception, "Owner is not a partition of unity");
      }
      if(verbose){
      for (int rank = 0; rank < world_comm.size(); ++rank) {
          if (rank == world_comm.rank()) {
             std::ostringstream msg;
             msg << "Grid size: " << gv.size(0) << " on rank " << world_comm.rank() << '\n';
              for (const auto& ind : entity_indexset) {
                msg << "Global " << ind.global();
                msg << " Local " << ind.local().local() << " Attribute " << ind.local().attribute();
                msg << " rank " << world_comm.rank();
                msg << " max_rank ib " << maxrank_ib[ind.local().local()];
                  if (ind.local().attribute() == Dune::OwnerOverlapCopyAttributeSet::owner) {
                      if (numpros[ind.local().local()] != 1) {
                          msg << " Owner entity  " << numpros[ind.local().local()] <<" owners "<< '\n';
                      }else{
                          msg << '\n';
                      }
                  } else {
                      if(numpros[ind.local().local()] != 1) {
                          msg << " Owner is not unique " << numpros[ind.local().local()] << '\n';
                      }else{
                          msg << '\n';
                      }
                  }
              }
              OpmLog::debug(msg.str());
          }
          usleep(100);
          world_comm.barrier();
      }
      }
      return entity_indexset;
  }
  
  Dune::OwnerOverlapCopyCommunication<int, int>::ParallelIndexSet
  entityToDofIndexSet(const Dune::OwnerOverlapCopyCommunication<int, int>::ParallelIndexSet& entity_indexset, int ndof, bool verbose = false)
  {
    
    using ParallelIndexSet = Dune::OwnerOverlapCopyCommunication<int, int>::ParallelIndexSet;
    ParallelIndexSet dof_indexset;
    dof_indexset.beginResize();
    for(const auto& ind: entity_indexset){
      auto lind = ind.local().local();
      auto at = ind.local().attribute();
      auto gid = ind.global();
      for(int dof=0; dof<ndof; ++dof){
        dof_indexset.add(ndof*gid+dof, ParallelIndexSet::LocalIndex(ndof*lind+dof, at, true));
      }
    }
    dof_indexset.endResize();
    if(verbose){
     std::ostringstream msg;
     for (const auto& ind : dof_indexset) {
          msg << " Global " << ind.global();
          msg << " Local " << ind.local().local() << " Attribute " << ind.local().attribute();
      }
     OpmLog::debug(msg.str());
    }
    return dof_indexset;
  }

}
