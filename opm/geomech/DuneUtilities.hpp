#pragma once
#include <dune/grid/common/datahandleif.hh>
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

#if DUNE_VERSION_LT(DUNE_GRID, 2, 8)
  bool fixedsize(int /* dim */, int /* codim */) const
  {
    return true;
  }
#else

  bool fixedSize(int /* dim */, int /* codim */) const
  {
    return true;
  }
#endif
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

#if DUNE_VERSION_LT(DUNE_GRID, 2, 8)
  bool fixedsize(int /* dim */, int /* codim */) const
  {
    return true;
  }
#else

  bool fixedSize(int /* dim */, int /* codim */) const
  {
    return true;
  }
#endif
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
    //std::cout << "value " << value << std::endl;
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
                    std::cout << " size error" << std::endl;
                    std::cout << "tmp" << std::endl;
                    for(auto& t:tmp){
                        std::cout << t << " ";
                    }
                    std::cout << std::endl;
                    std::cout << "cont" << std::endl;
                    for(auto& t:cont[i]){
                        std::cout << t << " ";
                    }

                    /*std::stringstream ss;
                    ss << " tmp ";
                    //print(ss,tmp);
                    ss << std::endl;
                    LOG_INFO(ss.str());
                    LOG_INFO("i j data "<< i << " " << j << " " << data);
                    LOG_INFO("size of precalculated set error " << int(tmp.size()) << " " << int(cont[i].size()));*/
                    assert(false);
                }
                std::copy(tmp.begin(),tmp.end(),cont[i].begin());
            }
            //cont[i].push_back(data);                                                                                                                                                        
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