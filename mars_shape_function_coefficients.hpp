#ifndef MARS_SHAPE_FUNCTION_COEFFICIENT_HPP
#define MARS_REFERENCE_MAP_HPP

#include "mars_base_elementfunctionspace.hpp"
#include "mars_vector.hpp"
#include "mars_array.hpp"

#include "mars_tuple_utilities.hpp"
#include "mars_shape_function.hpp"
#include "mars_general_form.hpp"
namespace mars{

template<typename Elem,Integer FEFamily,Integer Order>
class SingleShapeFunctionCoefficientsCollection;

template<typename Elem>
class SignedNormal;

template<typename Elem,Integer FEFamily>
class ShapeFunctionCoefficientsCollectionSingleType;


template<Integer Dim,Integer ManifoldDim>
class ShapeFunctionCoefficientsCollectionSingleType<Simplex<Dim,ManifoldDim>,LagrangeFE>
{
 public:
 using type=std::tuple<>;
 ShapeFunctionCoefficientsCollectionSingleType(){}
};

template<Integer Dim,Integer ManifoldDim>
class ShapeFunctionCoefficientsCollectionSingleType<Simplex<Dim,ManifoldDim>,RaviartThomasFE>
{
 public:
 using type=SignedNormal<Simplex<Dim,ManifoldDim>>;
 ShapeFunctionCoefficientsCollectionSingleType(){}

 void init(const Mesh<Dim,ManifoldDim>& mesh)
 {signed_normal_.init(mesh);}

 const type& operator()()const{return signed_normal_;}

 private:
 type signed_normal_;
};











template<typename Tuple,Integer Nmax,Integer N>
class ShapeFunctionGlobalCoefficientTupleTypeHelper2;

template<typename Tuple,Integer Nmax>
class ShapeFunctionGlobalCoefficientTupleTypeHelper2<Tuple,Nmax,Nmax>
{
public:
    using Elem=GetType<Tuple,Nmax,0>;
    static constexpr Integer FEFamily=GetType<Tuple,Nmax,1>::value;
    using type=std::tuple<typename ShapeFunctionCoefficientsCollectionSingleType<Elem,FEFamily>::type>;
};

template<typename Tuple,Integer Nmax,Integer N>
class ShapeFunctionGlobalCoefficientTupleTypeHelper2
{
public:
    using Elem=GetType<Tuple,N,0>;
    static constexpr Integer FEFamily=GetType<Tuple,N,1>::value;
    using single_type=std::tuple<typename ShapeFunctionCoefficientsCollectionSingleType<Elem,FEFamily>::type>;
    using type=decltype(std::tuple_cat(std::declval<single_type>(),
      std::declval<typename ShapeFunctionGlobalCoefficientTupleTypeHelper2<Tuple,Nmax,N+1>::type>()));
};

template<typename UniqueBaseFunctionSpaces>
using ShapeFunctionGlobalCoefficientTupleType2=
typename ShapeFunctionGlobalCoefficientTupleTypeHelper2<UniqueBaseFunctionSpaces,TupleTypeSize<UniqueBaseFunctionSpaces>::value-1,0>::type;






















template<typename Tuple,Integer Nmax,Integer N>
class ShapeFunctionGlobalCoefficientTupleTypeHelper;

template<typename Tuple,Integer Nmax>
class ShapeFunctionGlobalCoefficientTupleTypeHelper<Tuple,Nmax,Nmax>
{
public:
    // using Elem=GetType<Tuple,Nmax,0>;
    // using BaseFunctionSpace=GetType<Tuple,Nmax,1>;
    using Space=GetType<Tuple,Nmax>;
    using Elem=typename Space::Elem;
    using BaseFunctionSpace=Elem2FunctionSpace<Space>;
    using type=std::tuple<typename ShapeFunctionCoefficientsCollectionSingleType<Elem,BaseFunctionSpace::FEFamily>::type>;
};

template<typename Tuple,Integer Nmax,Integer N>
class ShapeFunctionGlobalCoefficientTupleTypeHelper
{
public:
    // using Elem=GetType<Tuple,N,0>;
    // using BaseFunctionSpace=GetType<Tuple,N,1>;
    using Space=GetType<Tuple,N>;
    using Elem=typename Space::Elem;
    using BaseFunctionSpace=Elem2FunctionSpace<Space>;
    using single_type=std::tuple<typename ShapeFunctionCoefficientsCollectionSingleType<Elem,BaseFunctionSpace::FEFamily>::type>;
    using type=decltype(std::tuple_cat(std::declval<single_type>(),
      std::declval<typename ShapeFunctionGlobalCoefficientTupleTypeHelper<Tuple,Nmax,N+1>::type>()));
};

template<typename UniqueBaseFunctionSpaces>
using ShapeFunctionGlobalCoefficientTupleType=
typename ShapeFunctionGlobalCoefficientTupleTypeHelper<UniqueBaseFunctionSpaces,TupleTypeSize<UniqueBaseFunctionSpaces>::value-1,0>::type;







template<typename Tuple,Integer Nmax,Integer N>
class ShapeFunctionLocalCoefficientTupleTypeHelper;

template<typename Tuple,Integer Nmax>
class ShapeFunctionLocalCoefficientTupleTypeHelper<Tuple,Nmax,Nmax>
{
public:
    // using Elem=GetType<Tuple,Nmax,0>;
    // using BaseFunctionSpace=GetType<Tuple,Nmax,1>;
    using Space=GetType<Tuple,Nmax>;
    using Elem=typename Space::Elem;
    using BaseFunctionSpace=Elem2FunctionSpace<Space>;
    static constexpr Integer NComponents=BaseFunctionSpace::NComponents;
    static constexpr Integer Ntot=FunctionSpaceDofsPerElem<ElemFunctionSpace<Elem,BaseFunctionSpace>>::value;
    static constexpr Integer Ndofs=Ntot/NComponents;
    using type=typename std::conditional<BaseFunctionSpace::FEFamily==LagrangeFE,std::tuple<std::tuple<>>, std::tuple<Array<Real,Ndofs> > >::type;
};

template<typename Tuple,Integer Nmax,Integer N>
class ShapeFunctionLocalCoefficientTupleTypeHelper
{
public:
    // using Elem=GetType<Tuple,N,0>;
    // using BaseFunctionSpace=GetType<Tuple,N,1>;
    using Space=GetType<Tuple,N>;
    using Elem=typename Space::Elem;
    using BaseFunctionSpace=Elem2FunctionSpace<Space>;
    static constexpr Integer NComponents=BaseFunctionSpace::NComponents;
    static constexpr Integer Ntot=FunctionSpaceDofsPerElem<ElemFunctionSpace<Elem,BaseFunctionSpace>>::value;
    static constexpr Integer Ndofs=Ntot/NComponents;
    using single_type=typename std::conditional<BaseFunctionSpace::FEFamily==LagrangeFE,std::tuple<std::tuple<>>, std::tuple<Array<Real,Ndofs> > >::type;
    using type=decltype(std::tuple_cat(std::declval<single_type>(),
      std::declval<typename ShapeFunctionLocalCoefficientTupleTypeHelper<Tuple,Nmax,N+1>::type>()));
};

template<typename Tuple>
using ShapeFunctionLocalCoefficientTupleType=typename ShapeFunctionLocalCoefficientTupleTypeHelper<Tuple,TupleTypeSize<Tuple>::value-1,0>::type;


template<typename Elem,Integer FEFamily,Integer Order, typename ConstInput, typename ShapeFunctionCoefficientsCollection>
void shape_function_coefficients_init(const ConstInput& mesh_ptr,ShapeFunctionCoefficientsCollection& coeff);


template<typename GeneralForm,typename ...GeneralForms>
class ShapeFunctionCoefficientsCollection;

template<typename Form_,typename...Forms_>
class ShapeFunctionCoefficientsCollection<GeneralForm<Form_>,GeneralForm<Forms_>...>
{
public:
  using Form=MultipleAddition<typename GeneralForm<Form_>::Form,typename GeneralForm<Forms_>::Form...>;
  using FunctionSpace=typename GeneralForm<Form_>::FunctionSpace;
  using MeshT=typename FunctionSpace::MeshT;
  using Elem=typename FunctionSpace::Elem;
  using NormalType= typename SignedNormal<Elem>::NormalType;
  using UniqueElementFunctionSpacesTupleType=typename FunctionSpace::UniqueElementFunctionSpacesTupleType;  
  using SpacesToUniqueFEFamily=typename FunctionSpace::SpacesToUniqueFEFamily;
  using GlobalTuple=ShapeFunctionGlobalCoefficientTupleType2<UniqueElementFEFamily2<UniqueElementFunctionSpacesTupleType>>;
  using LocalTuple=ShapeFunctionLocalCoefficientTupleType<UniqueElementFunctionSpacesTupleType>; 
  static constexpr Integer Nglobal=TupleTypeSize<GlobalTuple>::value-1;
  static constexpr Integer Nlocal=TupleTypeSize<LocalTuple>::value-1;


  ShapeFunctionCoefficientsCollection(MeshT& mesh):
  mesh_(mesh)
  {}
  // Initialize the global shape function coefficients tuple with the mesh
  // RT elements have different coefficients for different order, 
  // but all of them are based on the face "orientation" of the element
  // so we distinguihs between global coefficient and the local ones
  // template<Integer Dim, Integer ManifoldDim, typename T>
  // typename std::enable_if_t< IsSame<T,std::tuple<>>::value, void >
  // init_aux_aux(Mesh<Dim,ManifoldDim>& mesh, T& t)
  // {}

  // template<Integer Dim, Integer ManifoldDim, typename T>
  // typename std::enable_if_t< IsDifferent<T,std::tuple<>>::value, void >
  // init_aux_aux(Mesh<Dim,ManifoldDim>& mesh, T& t)
  // {
  //   t.init(mesh);
  // }
  

  // template<Integer M=0, Integer Dim,Integer ManifoldDim>
  // typename std::enable_if_t< (M>Nglobal), void >
  // init_aux(Mesh<Dim,ManifoldDim>& mesh)
  // {}

  // template<Integer M=0, Integer Dim,Integer ManifoldDim>
  // typename std::enable_if_t< (M<=Nglobal), void >
  // init_aux(Mesh<Dim,ManifoldDim>& mesh)
  // {
  //   init_aux_aux(mesh,tuple_get<M>(global_tuple_));
  //   init_aux<M+1>(mesh);
  // }


  // template<Integer Dim,Integer ManifoldDim>
  // void init(Mesh<Dim,ManifoldDim>& mesh)
  // {
  //  init_aux(mesh);
  // }


 
  template<Integer M=0>
  typename std::enable_if_t< (M>Nglobal), void >
  init_aux()
  {}

  template<Integer M=0>
  typename std::enable_if_t< (M<=Nglobal), void >
  init_aux()
  {
    mesh_.init_signed_normal();
    init_aux<M+1>();
  }

  void init()
  {
   init_aux();
  }

  // Initialize the local shape function coefficients tuple with local informations (elem id)
  template<typename Elem,Integer FEFamily,Integer Order,typename S,typename T>
  typename std::enable_if_t< IsSame<T,std::tuple<>>::value, void >
  init_aux_aux(const Integer elem_id,const S& s, T& t)
  {}

  template<typename Elem,Integer FEFamily,Integer Order,typename S,typename T>
  typename std::enable_if_t< IsDifferent<T,std::tuple<>>::value, void >
  init_aux_aux(const Integer elem_id,const S& s, T& t)
  {
    // shape_function_coefficients_init<Elem,FEFamily,Order>(s.sign(elem_id),t);
    SingleShapeFunctionCoefficientsCollection<Elem,FEFamily,Order>::apply(s.sign(elem_id),t);
  }


  template<Integer M=0>
  typename std::enable_if_t< (M>Nlocal), void >
  init_aux(const Integer elem_id)
  {}

  template<Integer M=0>
  typename std::enable_if_t< (M<=Nlocal), void >
  init_aux(const Integer elem_id)
  {
    constexpr Integer N=GetType<SpacesToUniqueFEFamily,M>::value;
    using FunctionSpace=GetType<UniqueElementFunctionSpacesTupleType,M>;
    // using Elem=GetType<FunctionSpace,0>;
    // using BaseFunctionSpace=GetType<FunctionSpace,1>;
    using Elem=typename FunctionSpace::Elem;
    using BaseFunctionSpace=Elem2FunctionSpace<FunctionSpace>;
    constexpr Integer FEFamily=BaseFunctionSpace::FEFamily;
    constexpr Integer Order=BaseFunctionSpace::Order;
    // std::cout<<"local, global=="<<M<<", "<<N<<std::endl;
    // std::cout<<"FEFamily=="<<FEFamily<<std::endl;
    

    auto& signed_normal=mesh_.signed_normal();
    // init_aux_aux<Elem,FEFamily,Order>(elem_id,tuple_get<N>(global_tuple_),tuple_get<M>(local_tuple_));
    init_aux_aux<Elem,FEFamily,Order>(elem_id,signed_normal,tuple_get<M>(local_tuple_));
    init_aux<M+1>(elem_id);

  }

  void init(const Integer elem_id)
  {
   init_aux(elem_id);

   // const auto& stampa=tuple_get<0>(global_tuple_).sign();
   //  std::cout<<".  STAMPA "<<std::endl;

   // for(Integer i=0;i<stampa.size();i++)
   // {
   //  std::cout<<stampa[i]<<std::endl;
   // }
  }

  
   const auto& operator()()const
  {
    return local_tuple_;
  }

  // template<Integer M>
  //  const auto& value()const
  // {
  //   return tuple_get<GetType<SpacesToUniqueFEFamily,M>::value>(local_tuple_);
  // }

private:
    GlobalTuple global_tuple_;
    LocalTuple local_tuple_;
    NormalType normals_;
    MeshT& mesh_;
};


template<typename MeshT,typename ConstFormReference,typename...ConstFormReferences>
constexpr auto shape_function_coefficients(MeshT& mesh,const ConstFormReference& form,const ConstFormReferences&...forms)
{return ShapeFunctionCoefficientsCollection<ConstFormReference,ConstFormReferences...>(mesh) ; }



}
#endif //MARS_REFERENCE_MAP_HPP