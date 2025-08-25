#pragma once
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/version.hh>

#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>


//#include <dune/grid/utility/globalindexset.hh>

#include <cstring>
#include <iostream>
#include <unistd.h>
#include <unistd.h>

#include "dune_utilities.hpp"
namespace Opm {
  void cellCellCommunication(Dune::CpGrid grid, Dune::MPIHelper& mpihelper){
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
	// cell_indexset.add(gidset.id(elem),
	//            ParallelIndexSet::LocalIndex(index, AttributeSet::owner, true));
      } else if (elem.partitionType() == Dune::PartitionType::BorderEntity) {
	// cell_indexset.add(gidset.id(elem),
	//             ParallelIndexSet::LocalIndex(index, Dune::AttributeSet::overlap, true));
      } else if (elem.partitionType() == Dune::PartitionType::FrontEntity) {
	// cell_indexset.add(gidset.id(elem),
	//            ParallelIndexSet::LocalIndex(index, AttributeSet::border, true));
      } else if (elem.partitionType() == Dune::PartitionType::GhostEntity) {
	// cell_indexset.add(gidset.id(elem),
	//             ParallelIndexSet::LocalIndex(index, AttributeSet::ghost, true));
      } else {
	std::cout << "Unknown partition type" << std::endl;
      }
    }

    cell_indexset.endResize();
    remote_indices.rebuild<false>();
    //Dune::OwnerOverlapCopyCommunication<int, int> cell_cell_comm(world_comm, Dune::SolverCategory::overlapping);
    //cell_cell_comm.remoteIndices() = remote_indices;
    //cell_cell_comm.indexSet() = cell_indexset;
    std::vector<int> myrank(gv.size(0), world_comm.rank());
    // //cell_cell_comm.copyOwnerToAll(myrank, myrank);
        

    Dune::Interface cominterface;
    using OwnerSet = Dune::OwnerOverlapCopyCommunication<int, int>::OwnerSet;//Dune::EnumItem<AttributeSet,Dune::OwnerOverlapCopyAttributeSet::owner>;
    using AllSet =  Dune::OwnerOverlapCopyCommunication<int, int>::AllSet;//Dune::AllSet<AttributeSet>;
    OwnerSet soureFlags;
    //AllSet destFlags;
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
      //allranks[i].resize(4,1000);
      allranks[i][0]=world_comm.rank();
    }

    world_comm.barrier();
    //OwnerSet ownerflags;
    //Dune::NegateSet<OwnerSet> not_ownerflags;
    //typename Dune::RedistributeInformation< Dune::OwnerOverlapCopyCommunication<int, int>  > ri;
    //ri.getInterface().free();
    //ri.getInterface().build(remote_indices, ownerflags, not_ownerflags);
    //ri.template redistribute<Dune::VariableVectorAdderGatherScatter>(allranks, allranks);
    Dune::BufferedCommunicator varvec_cell_cell_comm;
    varvec_cell_cell_comm.build(allranks, allranks, all_cominterface);
    varvec_cell_cell_comm.template forward<Dune::VariableVectorAdderGatherScatter>(allranks, allranks);
    //varvec_cell_cell_comm.free();

    cell_cell_comm.free();
    //remote_indices.free();

    world_comm.barrier();
        
    for (int rank = 0; rank < mpihelper.size(); ++rank) {
      if (rank == mpihelper.rank()) {
	std::cout << "Grid size: " << gv.size(0) << " on rank " << mpihelper.rank() << std::endl;  
	for (auto& elem : Dune::elements(gv)) {
	  auto lindex = elem.index();
	  // auto gindex = gindexset.index(elem);
	  auto gid = gidSet.id(elem);
	  std::cout << "Element global: id " << gid << " local " << lindex << " type " << elem.partitionType()
		    << " owner rank " << myrank[lindex] 
		    << " num pros " << numpros[lindex] 
		    << " all ranks ";
	  for(auto& r:allranks[lindex]){
	    //if(r!=1000){
	    std::cout << r << " ";
	    //}
	  } 
	  std::cout << std::endl;
	}
      }
      // for bliding communicator see ParallelIstlInformation. (or Dune::OwnerOverlapCopyCommunication)
      usleep(100);
      world_comm.barrier();
    }
  }

  void vertexVertexCommunication(Dune::CpGrid grid, Dune::MPIHelper& mpihelper){
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
	// vertex_indexset.add(gidset.id(elem),
	//            ParallelIndexSet::LocalIndex(index, AttributeSet::owner, true));
      } else if (vertex.partitionType() == Dune::PartitionType::BorderEntity) {
	// vertex_indexset.add(gidset.id(elem),
	//             ParallelIndexSet::LocalIndex(index, Dune::AttributeSet::overlap, true));
      } else if (vertex.partitionType() == Dune::PartitionType::FrontEntity) {
	// vertex_indexset.add(gidset.id(elem),
	//            ParallelIndexSet::LocalIndex(index, AttributeSet::border, true));
      } else if (vertex.partitionType() == Dune::PartitionType::GhostEntity) {
	// vertex_indexset.add(gidset.id(elem),
	//             ParallelIndexSet::LocalIndex(index, AttributeSet::ghost, true));
      } else {
	std::cout << "Unknown partition type" << std::endl;
      }
    }

    vertex_indexset.endResize();
    remote_indices.rebuild<false>();
    //Dune::OwnerOverlapCopyCommunication<int, int> vertex_vertex_comm(world_comm, Dune::SolverCategory::overlapping);
    //vertex_vertex_comm.remoteIndices() = remote_indices;
    //vertex_vertex_comm.indexSet() = vertex_indexset;
    std::vector<int> myrank(gv.size(3), world_comm.rank());
    // //vertex_vertex_comm.copyOwnerToAll(myrank, myrank);
        

    Dune::Interface cominterface;
    using OwnerSet = Dune::OwnerOverlapCopyCommunication<int, int>::OwnerSet;//Dune::EnumItem<AttributeSet,Dune::OwnerOverlapCopyAttributeSet::owner>;
    using AllSet =  Dune::OwnerOverlapCopyCommunication<int, int>::AllSet;//Dune::AllSet<AttributeSet>;
    OwnerSet soureFlags;
    //AllSet destFlags;
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
      //allranks[i].resize(4,1000);
      allranks[i][0]=world_comm.rank();
    }

    world_comm.barrier();
    //OwnerSet ownerflags;
    //Dune::NegateSet<OwnerSet> not_ownerflags;
    //typename Dune::RedistributeInformation< Dune::OwnerOverlapCopyCommunication<int, int>  > ri;
    //ri.getInterface().free();
    //ri.getInterface().build(remote_indices, ownerflags, not_ownerflags);
    //ri.template redistribute<Dune::VariableVectorAdderGatherScatter>(allranks, allranks);
    Dune::BufferedCommunicator varvec_vertex_vertex_comm;
    varvec_vertex_vertex_comm.build(allranks, allranks, all_cominterface);
    varvec_vertex_vertex_comm.template forward<Dune::VariableVectorAdderGatherScatter>(allranks, allranks);
    //varvec_vertex_vertex_comm.free();

    vertex_vertex_comm.free();
    //remote_indices.free();

    world_comm.barrier();
        
    for (int rank = 0; rank < mpihelper.size(); ++rank) {
      if (rank == mpihelper.rank()) {
	std::cout << "Grid size: " << gv.size(0) << " on rank " << mpihelper.rank() << std::endl;  
	for (auto& elem : Dune::vertices(gv)) {
	  auto lindex = elem.index();
	  // auto gindex = gindexset.index(elem);
	  auto gid = gidSet.id(elem);
	  std::cout << "Vertex global: id " << gid << " local " << lindex << " type " << elem.partitionType()
		    << " owner rank " << myrank[lindex] 
		    << " num pros " << numpros[lindex] 
		    << " all ranks ";
	  for(auto& r:allranks[lindex]){
	    //if(r!=1000){
	    std::cout << r << " ";
	    //}
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
  entityEntityCommunication(Dune::CpGrid grid, Dune::MPIHelper& mpihelper)
  {
      Dune::Codim<codim> mycodim;
      // Dune::PartitionType::AllPartition
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
              // entity_indexset.add(gidset.id(elem),
              //            ParallelIndexSet::LocalIndex(index, AttributeSet::owner, true));
          } else if (entity.partitionType() == Dune::PartitionType::BorderEntity) {
            entity_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::copy, true));
              // entity_indexset.add(gidset.id(elem),
              //             ParallelIndexSet::LocalIndex(index, Dune::AttributeSet::overlap, true));
          } else if (entity.partitionType() == Dune::PartitionType::FrontEntity) {
            entity_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::copy, true));
              // entity_indexset.add(gidset.id(elem),
              //            ParallelIndexSet::LocalIndex(index, AttributeSet::border, true));
          } else if (entity.partitionType() == Dune::PartitionType::GhostEntity) {
            entity_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::copy, true));
              // entity_indexset.add(gidset.id(elem),
              //             ParallelIndexSet::LocalIndex(index, AttributeSet::ghost, true));
          } else {
              std::cout << "Unknown partition type" << std::endl;
          }
      }

      entity_indexset.endResize();
      remote_indices.rebuild<false>();
      // Dune::OwnerOverlapCopyCommunication<int, int> entity_entity_comm(world_comm,
      // Dune::SolverCategory::overlapping); entity_entity_comm.remoteIndices() = remote_indices;
      // entity_entity_comm.indexSet() = entity_indexset;
      std::vector<int> myrank(gv.size(codim), world_comm.rank());
      // //entity_entity_comm.copyOwnerToAll(myrank, myrank);


      Dune::Interface cominterface;
      using OwnerSet = Dune::OwnerOverlapCopyCommunication<int, int>::
          OwnerSet; // Dune::EnumItem<AttributeSet,Dune::OwnerOverlapCopyAttributeSet::owner>;
      using AllSet = Dune::OwnerOverlapCopyCommunication<int, int>::AllSet; // Dune::AllSet<AttributeSet>;
      OwnerSet soureFlags;
      // AllSet destFlags;
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
          // allranks[i].resize(4,1000);
          allranks[i][0] = world_comm.rank();
      }

      world_comm.barrier();
      // OwnerSet ownerflags;
      // Dune::NegateSet<OwnerSet> not_ownerflags;
      // typename Dune::RedistributeInformation< Dune::OwnerOverlapCopyCommunication<int, int>  > ri;
      // ri.getInterface().free();
      // ri.getInterface().build(remote_indices, ownerflags, not_ownerflags);
      // ri.template redistribute<Dune::VariableVectorAdderGatherScatter>(allranks, allranks);
      Dune::BufferedCommunicator varvec_entity_entity_comm;
      varvec_entity_entity_comm.build(allranks, allranks, all_cominterface);
      varvec_entity_entity_comm.template forward<Dune::VariableVectorAdderGatherScatter>(allranks, allranks);
      // varvec_entity_entity_comm.free();

      entity_entity_comm.free();
      // remote_indices.free();

      world_comm.barrier();

      for (int rank = 0; rank < mpihelper.size(); ++rank) {
          if (rank == mpihelper.rank()) {
              std::cout << "Grid size: " << gv.size(0) << " on rank " << mpihelper.rank() << std::endl;


              for (auto& elem : Dune::entities(gv, mycodim)) {
                  auto lindex = elem.index();
                  // auto gindex = gindexset.index(elem);
                  auto gid = gidSet.id(elem);
                  std::cout << "Entity<" << codim << "> global: id " << gid << " local " << lindex << " type "
                            << elem.partitionType() << " owner rank " << myrank[lindex] << " num pros "
                            << numpros[lindex] << " all ranks ";
                  for (auto& r : allranks[lindex]) {
                      // if(r!=1000){
                      std::cout << r << " ";
                      //}
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
  Dune::OwnerOverlapCopyCommunication<int, int>::ParallelIndexSet makeEntityEntityCommunication(const Grid& grid, bool verbose = false){
    //static_assert(false);
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
    Dune::Communication<MPI_Comm> world_comm = grid.comm();//Dune::MPIHelper::getCommunication();
    const auto& gv = grid.leafGridView();
    auto indexset = gv.indexSet();
    auto gidSet = gv.grid().globalIdSet();
      
    Vector maxranks_tmp_ib(gv.size(codim), world_comm.rank());
    using GridType = Dune::CpGrid;
    using GridView =  typename GridType::LeafGridView;
    Dune::MaxEntityVectorVectorDataHandle<GridView,Vector> datahandle_ib(maxranks_tmp_ib, gv, codim);
    gv.communicate(datahandle_ib, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);
    //gv.communicate(datahandle, Dune::All_All_Interface, Dune::ForwardCommunication);
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
    // //gv.communicate(datahandle, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);
    gv.communicate(datahandle_all, Dune::All_All_Interface, Dune::ForwardCommunication);
    const auto& allranks_all = datahandle_all.all_data();


    //std::cout << "Getting max rank of node " << std::endl;
    world_comm.barrier();
    //std::cout << "Finnish max rank of node " << std::endl;
      
      // Dune::PartitionType::AllPartition
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
            // if(allranks_all[index].size() > 1){
            //   std::cout << "Error My rank" << world_comm.rank() << " parition type " << entity.partitionType() << " global id " << gid << " local id " << index;
            //   std::cout << "all other ranks ";
            //   for (auto rr : allranks_all[index]){
            //     std::cout << rr << " ";
            //   }
            //   std::cout << std::endl;

            // }
            // if(myrank == maxrank_ib[index]){
              entity_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::owner, true));
            //  }else{
            //    //assert(false);
            //    std::cout<< "Warning!!!: Owner do not seem to be unique ?" << std::endl;
            //    entity_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::copy, true));
            //  }
          } else if (entity.partitionType() == Dune::PartitionType::OverlapEntity) {
              entity_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::copy, true));
              // entity_indexset.add(gidset.id(elem),
              //            ParallelIndexSet::LocalIndex(index, AttributeSet::owner, true));
          } else if (entity.partitionType() == Dune::PartitionType::BorderEntity) {
            if(myrank == maxrank_ib[index]){
                entity_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::owner, true));
            }else{
                entity_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::copy, true));
            }
              // entity_indexset.add(gidset.id(elem),
              //             ParallelIndexSet::LocalIndex(index, Dune::AttributeSet::overlap, true));
          } else if (entity.partitionType() == Dune::PartitionType::FrontEntity) {
            entity_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::copy, true));
              // entity_indexset.add(gidset.id(elem),
              //            ParallelIndexSet::LocalIndex(index, AttributeSet::border, true));
          } else if (entity.partitionType() == Dune::PartitionType::GhostEntity) {
                entity_indexset.add(gid, ParallelIndexSet::LocalIndex(index, AttributeSet::copy, true));
              // entity_indexset.add(gidset.id(elem),
              //             ParallelIndexSet::LocalIndex(index, AttributeSet::ghost, true));
          } else {
              std::cout << "Unknown partition type" << std::endl;
          }
      }
      entity_indexset.endResize();
      remote_indices.rebuild<false>();

      // check if unique owner and that all indices has an owner
      using AllSet = Dune::OwnerOverlapCopyCommunication<int, int>::AllSet; // Dune::AllSet<AttributeSet>;
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
      //for (const auto& ind : entity_indexset) {
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
           if(!verbose || true){
              std::cout << "Num procs " << numpros[ind.local().local()]; 
              std::cout << " Global " << ind.global(); 
              std::cout << " Local " << ind.local().local() << " Attribute " << ind.local().attribute();// << std::endl;
              std::cout << " rank " << world_comm.rank();
              std::cout << " enity is " << entity.partitionType();
              std::cout << " max_rank " << maxrank_ib[ind.local().local()];
              std::cout << " all ranks ib ";
              for(auto other:allranks_ib[ind.local().local()]){
                  std::cout << other << " ";
              }
              std::cout << " all ranks all ";
              for(auto other:allranks_all[ind.local().local()]){
                  std::cout << other << " ";
              }
              std::cout << std::endl;
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
                   auto& cv =  grid.current_view_data_;
                   auto& pi = cv->partition_type_indicator_;
                   std::cout << "Rank " << world_comm.rank() << " Cell global " << cgid  << " local " << cind << " parition " << Dune::PartitionType(pi->cell_indicator_[cind]) << std::endl; 
                   for(int lface=0; lface < faces.size(); ++lface){
                    int face = faces[lface].index();
                    int nv = grid.numFaceVertices(face);
                    std::cout << "Face " << face << " " << pi->getFacePartitionType(face) << ":" ;
                    for(int lnode=0; lnode < nv; ++lnode){
                      auto node = grid.faceVertex(face, lnode);
                      auto gnode = local_global[node];
                      //assert(index == node );
                      std::cout << "{" << gnode << "," << node <<"," << Dune::PartitionType(pi->point_indicator_[node]) << "}";
                    }
                    std::cout << std::endl;
                  } 
                  }

                }
              }
              ok = false;
              //assert(false);
              //DUNE_THROW(Dune::Exception, "Owner is not a partition of unity");
           } else{
            ok = false;
           }
          }
      }
      world_comm.barrier();
      int not_ok = !ok;
      not_ok = world_comm.sum(not_ok);
      if(not_ok > 0){
          std::cout << "Owner is not a partition of unity" << std::endl;
          //DUNE_THROW(Dune::Exception, "Owner is not a partition of unity");
      }
      if(verbose){
      for (int rank = 0; rank < world_comm.size(); ++rank) {
          if (rank == world_comm.rank()) {
             std::cout << "Grid size: " << gv.size(0) << " on rank " << world_comm.rank() << std::endl;
              for (const auto& ind : entity_indexset) {
                std::cout << "Global " << ind.global(); 
                std::cout << " Local " << ind.local().local() << " Attribute " << ind.local().attribute();// << std::endl;
                std::cout << " rank " << world_comm.rank();
                std::cout << " max_rank ib " << maxrank_ib[ind.local().local()];
                //std::cout << " max_rank all " << maxrank_all[ind.local().local()];
                  if (ind.local().attribute() == Dune::OwnerOverlapCopyAttributeSet::owner) {
                      if (numpros[ind.local().local()] != 1) {
                          //DUNE_THROW(Dune::Exception, "Owner entity has more than one owner");
                          std::cout << " Owner entity  " << numpros[ind.local().local()] <<" owners "<< std::endl;
                      }else{
                          std::cout << std::endl;
                      }
                  } else {
                      if(numpros[ind.local().local()] != 1) {
                          //DUNE_THROW(Dune::Exception, " Owner is not unique");
                          std::cout << " Owner is not unique " << numpros[ind.local().local()] << std::endl;
                      }else{
                          std::cout << std::endl;
                      }
                  }
              }
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
     for (const auto& ind : dof_indexset) { 
          std::cout << " Global " << ind.global(); 
          std::cout << " Local " << ind.local().local() << " Attribute " << ind.local().attribute();// << std::endl;
      }
    }
    return dof_indexset;
  }

}
