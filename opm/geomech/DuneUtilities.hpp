/*
  Copyright 2024 Equinor ASA.

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
#include <dune/grid/common/datahandleif.hh>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <sstream>
namespace Dune{
    template<class GridView, class Vector>
    class EntityVectorVectorDataHandle
  :  public Dune::CommDataHandleIF< EntityVectorVectorDataHandle<GridView,Vector>, typename Vector::value_type>
{
public:

  /// \brief the data type we send
  using DataType = typename Vector::value_type;

  /// \brief Constructor
  /// \param data The vector of data vectors
  /// \param gridView The gridview the data is attached to.
  EntityVectorVectorDataHandle(Vector& data, const GridView& gridView, int codim)
    : data_(data), gridView_(gridView), codim_(codim)
  {}

  bool contains(int /* dim */, int codim) const
  {
    return codim == codim_;
  }

  // Must be spelled fixedSize (Dune >= 2.8): otherwise the CRTP dispatch in
  // CommDataHandleIF::fixedSize calls itself (undefined behaviour).
  bool fixedSize(int /* dim */, int /* codim */) const
  {
    return true;
  }
  template<class EntityType>
  std::size_t size(const EntityType /* entity */) const
  {
    return 1;
  }


  template<class BufferType, class EntityType>
  void gather(BufferType& buffer, const EntityType& e) const
  {

      buffer.write(data_[gridView_.indexSet().index(e)]);
  }

  template<class BufferType, class EntityType>
  void scatter(BufferType& buffer, const EntityType& e,
               [[maybe_unused]] std::size_t n)
  {
    assert(n == 1);
    buffer.read(data_[gridView_.indexSet().index(e)]);
  }
private:
  Vector& data_;
  const GridView& gridView_;
  int codim_;
};


 template<class GridView, class Vector>
    class MaxEntityVectorVectorDataHandle
  :  public Dune::CommDataHandleIF< MaxEntityVectorVectorDataHandle<GridView,Vector>, typename Vector::value_type>
{
public:

  /// \brief the data type we send
  using DataType = typename Vector::value_type;

  /// \brief Constructor
  /// \param data The vector of data vectors
  /// \param gridView The gridview the data is attached to.
  MaxEntityVectorVectorDataHandle(Vector& data, const GridView& gridView, int codim)
    : data_(data), gridView_(gridView), codim_(codim)
  {
    all_data_.resize(gridView.size(codim));
  }

  bool contains(int /* dim */, int codim) const
  {
    return codim == codim_;
  }

  // Must be spelled fixedSize (Dune >= 2.8): otherwise the CRTP dispatch in
  // CommDataHandleIF::fixedSize calls itself (undefined behaviour).
  bool fixedSize(int /* dim */, int /* codim */) const
  {
    return true;
  }
  template<class EntityType>
  std::size_t size(const EntityType /* entity */) const
  {
    return 1;
  }


  template<class BufferType, class EntityType>
  void gather(BufferType& buffer, const EntityType& e) const
  {

      buffer.write(data_[gridView_.indexSet().index(e)]);
  }

  template<class BufferType, class EntityType>
  void scatter(BufferType& buffer, const EntityType& e,
               [[maybe_unused]] std::size_t n)
  {
    assert(n == 1);
    DataType value;
    value = 5;
    buffer.read(value);
    data_[gridView_.indexSet().index(e)] = std::max(data_[gridView_.indexSet().index(e)],value);
    all_data_[gridView_.indexSet().index(e)].push_back(value);
  }
  const std::vector<std::vector<int>>& all_data() const
  {
    return all_data_;
  }
private:
  Vector& data_;
  std::vector<std::vector<int>> all_data_;
  const GridView& gridView_;
  int codim_;
};




    template<typename T>
    struct MaxGatherScatter
    {
      typedef typename CommPolicy<T>::IndexedType V;

      static V gather(const T& a, std::size_t i)
      {
        return a[i];
      }

      static void scatter(T& a, V v, std::size_t i)
      {
        a[i] = std::max(a[i],v);
      }
    };
    template <>
    struct CommPolicy< std::vector< std::vector<int> > >
    {
        typedef std::vector<std::vector<int>> Type;
        typedef int IndexedType;
        typedef VariableSize IndexedTypeFlag;
        static const void* getAddress(const Type& v, int index){
            return &(v[index][0]);
        }
        static size_t getSize(const Type&v, int index){
            return v[index].size();
        }
    };

    struct VariableVectorAdderGatherScatter
    {
        typedef std::vector< std::vector< int > > Container;
        typedef typename CommPolicy<Container>::IndexedType IndexedType;
        static const IndexedType gather(const Container& cont, std::size_t i, std::size_t j)
        {
            return cont[i][j];
        }
        static void scatter(Container& cont, const IndexedType& data, std::size_t i, std::size_t j)
        {
            std::set<int> tmp(cont[i].begin(),cont[i].end());
            size_t tmp_size = tmp.size();
            tmp.insert(data);
            if(tmp_size != tmp.size()){
                if(not(tmp.size() <= cont[i].size())){
                    std::ostringstream msg;
                    msg << " size error\ntmp\n";
                    for(auto& t:tmp){
                        msg << t << " ";
                    }
                    msg << "\ncont\n";
                    for(auto& t:cont[i]){
                        msg << t << " ";
                    }
                    Opm::OpmLog::error(msg.str());

                    assert(false);
                }
                std::copy(tmp.begin(),tmp.end(),cont[i].begin());
            }
        }

    };
        /** \brief gather/scatter callback for communcation */
    template<typename T>
    struct FreeCopyGatherScatter
    {
        typedef typename CommPolicy<T>::IndexedType V;

        static V gather(const T& a, std::size_t i)
        {
            return a[i];
        }

        static void scatter(T& a, V v, std::size_t i)
        {
            a[i] = v;
        }
    };
    template<typename T>
    struct FreeAddGatherScatter
    {
        typedef typename CommPolicy<T>::IndexedType V;

        static V gather(const T& a, std::size_t i)
	{
            return a[i];
	}

        static void scatter(T& a, V v, std::size_t i)
        {
            a[i] += v;
        }
    };

}