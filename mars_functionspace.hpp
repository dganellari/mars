/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Written by Gabriele Rovi (April 2019)                                                                             ////////
////// We define the FunctionSpace class:                                                                                ////////
////// 1) It takes a mesh and 1 or more FunctionSpaces (Lagrange1<2>, RT0<1>...)                                         ////////                                        
////// 2) It builds the dofmap: a vector (long n_elements), whose component is the array of all the dofs of the element  ////////
////// 2) dofmap(space_id,elem_id) returns the dofs of element elem_id corresponding to the space space_id               ////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MARS_FunctionSpace_HPP
#define MARS_FunctionSpace_HPP

#include "mars_base.hpp"
#include "mars_elementfunctionspace.hpp"
#include "mars_functionspace_dofmap.hpp"
#include "mars_shape_function.hpp"
namespace mars{


// template<Integer Nelem_dofs>
// class IFunctionSpace
// {
//   public:
//   virtual~IFunctionSpace(){};
//   virtual std::shared_ptr<IFunctionSpace> get(const Integer i)=0;
//   template<Integer Nelem_dofs>
//   virtual const std::vector<std::array<Integer, Nelem_dofs>>& dofmap(){};//return dofmap_;};
// };


class IFunctionSpace
{
  public:
  virtual~IFunctionSpace(){};
  virtual std::shared_ptr<IFunctionSpace> get(const Integer i)=0;

  // virtual const std::vector<std::array<Integer, Nelem_dofs>>& dofmap(){};//return dofmap_;};
};


// class CollectorFunctionSpace: public IFunctionSpace
// {
// public: 

//     template<typename T,typename...Args>
//     typename std::enable_if<0==sizeof...(Args), void>::type
//     add(const T& t,Args...more)
//     {
//      children.push_back(std::make_shared<T>(t));
//     }

//     template<typename T,typename...Args>
//     typename std::enable_if<0<sizeof...(Args), void>::type
//     add(const T& t,Args...more)
//     {
//       children.push_back(std::make_shared<T>(t));
//       add(more...);
//     }

//     template<typename...Parameters>
//     CollectorFunctionSpace(const Parameters&...params){add(params...);}

//     virtual std::shared_ptr<IFunctionSpace> get(const Integer i){return children[i];}
// private:
//   std::vector<std::shared_ptr<IFunctionSpace>> children;

// };
template <int N, typename... Ts>
struct get;

template <int N, typename T, typename... Ts>
struct get<N, std::tuple<T, Ts...>>
{
    using type = typename get<N - 1, std::tuple<Ts...>>::type;
};

template <typename T, typename... Ts>
struct get<0, std::tuple<T, Ts...>>
{
    using type = T;
};

template<typename...Args>
std::tuple<Args...> add_costant (const std::tuple<Args...>& t1,const Integer& value);


template<std::size_t Nelem_dofs>
std::vector<std::array<Integer, Nelem_dofs>> add_costant
(const std::vector<std::array<Integer, Nelem_dofs>>& t1,const Integer& value)
{
  std::vector<std::array<Integer, Nelem_dofs>> t2;
  const auto& t1size=t1.size();
  t2.resize(t1size);

  std::cout<<"value=="<<value<<std::endl;
  for(Integer ii=0;ii<t1size;ii++)
    for(Integer jj=0;jj<t1[ii].size();jj++)
       t2[ii][jj]=t1[ii][jj]+value;
  return t2;
}

template<Integer N,typename...Args>
typename std::enable_if<sizeof...(Args)==N+1, void>::type
 add_costant_tuple_loop(const std::tuple<Args...>& t1, std::tuple<Args...>& t2,const Integer& value)
{
std::get<N>(t2)=add_costant(std::get<N>(t1),value);
}

template<Integer N,typename...Args>
typename std::enable_if< N+1<sizeof...(Args), void>::type
 add_costant_tuple_loop(const std::tuple<Args...>& t1, std::tuple<Args...>& t2,const Integer& value)
{
  std::get<N>(t2)=add_costant(std::get<N>(t1),value);
  add_costant_tuple_loop<N+1,Args...>(t1,t2,value);
}


template<typename...Args>
std::tuple<Args...> add_costant (const std::tuple<Args...>& t1,const Integer& value)
{
  std::tuple<Args...> t2=t1;

  add_costant_tuple_loop<0,Args...>(t1,t2,value);

  return t2;
}


template<typename MeshT, typename BaseFunctionSpace,typename...BaseFunctionSpaces>
class FunctionSpace //: public IFunctionSpace
{
public:


      // virtual std::shared_ptr<IFunctionSpace> get(const Integer i){return std::make_shared<FunctionSpace>(*this);}

      using Elem= typename MeshT::Elem;

      static constexpr Integer Nsubspaces=1+sizeof...(BaseFunctionSpaces);

      static constexpr Integer Nelem_dofs=DofsPerElemNums<Elem,BaseFunctionSpace,BaseFunctionSpaces...>::value;


      using type_dofmap=std::vector<std::array<Integer, Nelem_dofs>>;
      using type_offset=std::array<std::vector<Integer>, Nsubspaces>;
      using type_space_dofs=std::array<std::vector<std::vector<Integer>>, Nsubspaces>;
      using type_space_infos=std::array<std::array<Integer,4>,Nsubspaces>;
      using type_base_function_space=std::tuple<BaseFunctionSpace,BaseFunctionSpaces...>;
      inline Integer n_subspaces()const{return Nsubspaces;};

      inline const Integer& components (const Integer& space_id)const{return space_infos_[space_id][3];};

      inline const Integer& n_elem_dofs()const{return Nelem_dofs;};

      inline Integer n_elem_dofs(const Integer& space_id)const{
                                  const auto& os=offset_[space_id];
                                  const auto size=os[os.size()-1]-os[0];
                                  return size;}

      inline Integer n_elem_dofs(const Integer& space_id,const Integer& component_id)const{
                                  const auto& size=n_elem_dofs(space_id);
                                  return (size/space_infos_[space_id][3]);}


      inline const Integer& n_dofs()const{return n_dofs_;};

      inline const Integer n_dofs(const Integer& space_id,const Integer& component_id)const
                                 {return space_dofs_[space_id][component_id].size(); };


      inline void set_new_start(const Integer& new_start)
      {
         dof_count_start_=new_start;
         const auto& n_elements=dofmap_.size();
         dofmap_new_start_.resize(n_elements);
         for(Integer elem_id=0;elem_id<n_elements;elem_id++)
           for(Integer nn=0;nn<Nelem_dofs;nn++)
              dofmap_new_start_[elem_id][nn]=dofmap_[elem_id][nn]+dof_count_start_;

         for(Integer nn=0;nn<Nsubspaces;nn++)
          {
            const auto& size1=space_dofs_[nn].size();
            space_dofs_new_start_[nn].resize(size1);
            for(Integer ii=0;ii<size1;ii++)
            { 
              const auto& size2=space_dofs_[nn][ii].size();
              space_dofs_new_start_[nn][ii].resize(size2);
              for(Integer jj=0;jj<size2;jj++)
                {
                  space_dofs_new_start_[nn][ii][jj]=space_dofs_[nn][ii][jj]+new_start;
                }
            }
          }
      }

      inline const type_dofmap& dofmap()const{return dofmap_;};

      inline const type_dofmap& dofmap_new_start()const{return dofmap_new_start_;};

      inline void  dofmap(const type_dofmap& dm)const{dm=dofmap_;};

      inline void  dofmap_new_start(const type_dofmap& dm)const{dm=dofmap_new_start_;};


      inline const std::array<Integer, Nelem_dofs>& dofmap(const Integer& elem_id)const
                         {return dofmap_[elem_id];};

      inline const std::array<Integer, Nelem_dofs>& dofmap_new_start(const Integer& elem_id)const
                         {return dofmap_new_start_[elem_id];};


      inline void  dofmap(const Integer& elem_id, const std::array<Integer, Nelem_dofs> & elem_dm)const
                         {elem_dm=dofmap_[elem_id];};

      inline void  dofmap_new_start(const Integer& elem_id, const std::array<Integer, Nelem_dofs> & elem_dm)const
                         {elem_dm=dofmap_new_start_[elem_id];};

      inline std::vector<Integer> 
                   dofmap(const Integer& space_id,const Integer& elem_id)const{
                        const auto& os=offset_[space_id];
                        const auto& size=n_elem_dofs(space_id);
                        std::vector<Integer> output(size);
                        for(Integer nn=0;nn<size;nn++)
                             output[nn]=dofmap_[elem_id][nn+os[0]];
                        return output;};

      inline std::vector<Integer> dofmap_new_start(const Integer& space_id,const Integer& elem_id)const{
                        const auto& os=offset_[space_id];
                        const auto& size=n_elem_dofs(space_id);
                        std::vector<Integer> output(size);
                        for(Integer nn=0;nn<size;nn++)
                             output[nn]=dofmap_new_start_[elem_id][nn+os[0]];
                        return output;};


      inline std::vector<Integer> 
                   dofmap(const Integer& space_id,const Integer& component_id,const Integer& elem_id)const{
                        const auto& os=offset_[space_id];
                        const auto& size= n_elem_dofs(space_id);
                        std::vector<Integer> output(size);
                        const auto& comp=components(space_id);
                        space_infos_[space_id][3];
                        for(Integer nn=component_id;nn<size;nn=nn+comp)
                             output[nn]=dofmap_[elem_id][nn+os[0]];
                        return output;};

      inline std::vector<Integer> dofmap_new_start(const Integer& space_id,const Integer& component_id,const Integer& elem_id)const{
                        const auto& os=offset_[space_id];
                        const auto& size= n_elem_dofs(space_id);
                        std::vector<Integer> output(size);
                        const auto& comp=components(space_id);
                        space_infos_[space_id][3];
                        for(Integer nn=component_id;nn<size;nn=nn+comp)
                             output[nn]=dofmap_new_start_[elem_id][nn+os[0]];
                        return output;};

      inline void dofmap(const Integer& space_id,const Integer& elem_id, const std::vector<Integer>& elem_space_dm)const
                       {
                        const auto& os=offset_[space_id];
                        const auto& size=n_elem_dofs(space_id);
                        elem_space_dm.resize(size);
                        for(Integer nn=0;nn<size;nn++)
                             elem_space_dm[nn]=dofmap_[elem_id][nn+os[0]];
                       };

      inline void dofmap_new_start(const Integer& space_id,const Integer& elem_id, const std::vector<Integer>& elem_space_dm)const
                       {
                        const auto& os=offset_[space_id];
                        const auto& size=n_elem_dofs(space_id);
                        elem_space_dm.resize(size);
                        for(Integer nn=0;nn<size;nn++)
                             elem_space_dm[nn]=dofmap_new_start_[elem_id][nn+os[0]];
                       };


      inline const std::array<std::vector<Integer>, Nsubspaces>& offset() const {return offset_;};

      inline void  offset(const std::array<std::vector<Integer>, Nsubspaces> &os)const {os=offset_;};


      inline const std::vector<Integer>& offset(const Integer& space_id)const{return offset_[space_id];};

      inline void offset(Integer space_id, const std::vector<Integer>& space_os)const {space_os=offset_[space_id];};


      inline const std::vector<Integer>& space_dofs(const Integer& space_id,const Integer& component_id) const
                                         {return space_dofs_[space_id][component_id];};


      inline const std::vector<Integer>& space_dofs_new_start(const Integer& space_id,const Integer& component_id) const
                                         {
                                          return space_dofs_new_start_[space_id][component_id];};

      inline void space_dofs(const Integer& space_id, const Integer& component_id,std::vector<Integer>& spacedofs)const
                            {spacedofs.resize(n_dofs(space_id,component_id));
                             spacedofs=space_dofs_[space_id][component_id];};
      inline void space_dofs_new_start(const Integer& space_id, const Integer& component_id,std::vector<Integer>& spacedofs)const
                            {spacedofs.resize(n_dofs(space_id,component_id));
                             spacedofs=space_dofs_[space_id][component_id]+dof_count_start_;};

      inline const type_space_infos& space_info()const{return space_infos_;};

      inline std::shared_ptr< MeshT > mesh()const {return mesh_;};

      inline void set_dof_count_start(const Integer&value){dof_count_start_=value;};

      
      FunctionSpace(const MeshT& mesh,const Integer dof_count_start=0):
      dof_count_start_(dof_count_start),
      mesh_(std::make_shared< MeshT >(mesh)),
      is_dof_count_start_set_(dof_count_start)
      {
      function_space_info<0,Nsubspaces,BaseFunctionSpace,BaseFunctionSpaces...>(space_infos_);
      dofmap_fespace<BaseFunctionSpace,BaseFunctionSpaces...>(mesh,dofmap_,offset_,n_dofs_,space_infos_,space_dofs_);     
      };

      template<typename...QuadratureRules>
      FunctionSpace(const MeshT& mesh,const Integer dof_count_start=0):
      dof_count_start_(dof_count_start),
      mesh_(std::make_shared< MeshT >(mesh)),
      is_dof_count_start_set_(dof_count_start)
      {
      function_space_info<0,Nsubspaces,BaseFunctionSpace,BaseFunctionSpaces...>(space_infos_);
      dofmap_fespace<BaseFunctionSpace,BaseFunctionSpaces...>(mesh,dofmap_,offset_,n_dofs_,space_infos_,space_dofs_);     
      // set_quadrature_rule<0,QuadratureRules...>();
      // if(is_dof_count_start_set_==true)
      //    set_new_start(dof_count_start);
      };

private:
   
      std::shared_ptr< MeshT > mesh_;
      Integer dof_count_start_;
      bool is_dof_count_start_set_;
      Integer n_dofs_;
      type_dofmap dofmap_;
      type_dofmap dofmap_new_start_;
      type_offset offset_;
      type_space_dofs space_dofs_;
      type_space_dofs space_dofs_new_start_;
      type_space_infos space_infos_;
      type_base_function_space shape_functions_;
};




template<typename FunctionSpaceType, typename...QuadratureRules>
class FiniteElementSpace
{
 public:
      using type_base_function_space=typename FunctionSpaceType::type_base_function_space;


      template<Integer N,typename QR,typename...QRs>
      struct tuple_shape_function_type
      {    

            using BFS=typename get<N, type_base_function_space>::type;
            using rest = typename tuple_shape_function_type<N+1,QRs...>::type; 
            using tuple_ens=std::tuple<ShapeFunctionOperator<QR,typename QR::Elem, BFS>>;
            using type = decltype( std::tuple_cat( std::declval< tuple_ens >(), 
                                                   std::declval< rest >() ) );
      };


      template<Integer N,typename QR>
      struct tuple_shape_function_type<N,QR>
      {
           using BFS=typename get<N, type_base_function_space>::type;
           using type = typename std::tuple<ShapeFunctionOperator<QR,typename QR::Elem, BFS>>;
      };




      template<Integer N=0,typename QR,typename...QRs>
      typename std::enable_if< 0==sizeof...(QRs),
                               typename tuple_shape_function_type<N, QR>::type
                             >::type 
      set_quadrature_rule()
      {       
       using BFS=typename get<N, type_base_function_space>::type;
       ShapeFunctionOperator<QR,typename QR::Elem,typename get<N, type_base_function_space>::type> shape_function;
       return std::tuple<decltype(shape_function)> (shape_function);
      }

      template<Integer N=0,typename QR,typename...QRs>
      typename std::enable_if< 0<sizeof...(QRs),
                               typename tuple_shape_function_type<N,QR,QRs...>::type
                             >::type 
      set_quadrature_rule()
      {       
       using BFS=typename get<N, type_base_function_space>::type;
       ShapeFunctionOperator<QR,typename QR::Elem, BFS> shape_function;
       return std::tuple_cat(std::make_tuple(shape_function),set_quadrature_rule<N+1,QRs...>()  );
      }


      FiniteElementSpace(const FunctionSpaceType& fs): 
      function_space_(std::make_shared<FunctionSpaceType>(fs)),
      shape_function_(set_quadrature_rule<0,QuadratureRules...>())
      {

      }


      const FunctionSpaceType& function_space()const{return *function_space_;};
      const std::shared_ptr<FunctionSpaceType>& function_space_ptr()const{return function_space_;};



 private:
 std::shared_ptr<FunctionSpaceType> function_space_;
 typename tuple_shape_function_type<0,QuadratureRules...>::type shape_function_;

};






template<typename MeshT,typename ...Parameters >
class FEProva
{
  public:
     using Elem=typename MeshT::Elem;
     

      template<typename Arg1,typename Arg2,typename...Args>
      struct tuple_space_and_quadrature_rule_types
      {    
       public:
            using rest1 = typename tuple_space_and_quadrature_rule_types<Args...>::type1; 
            using rest2 = typename tuple_space_and_quadrature_rule_types<Args...>::type2; 

            using single_type1=std::tuple<Arg1>;
            using single_type2=std::tuple< typename Arg2:: template rule<Elem> >;

            using type1 = decltype( std::tuple_cat( std::declval< single_type1 >(), 
                                                    std::declval< rest1 >() ) );
            using type2 = decltype( std::tuple_cat( std::declval< single_type2 >(), 
                                                    std::declval< rest2 >() ) );           
      };


      template<typename Arg1,typename Arg2>
      struct tuple_space_and_quadrature_rule_types<Arg1,Arg2>
      {
      public:
           using type1 = typename std::tuple<Arg1>;
           using type2 = typename std::tuple< typename Arg2:: template rule<Elem> >;
      };

      using BaseFunctionSpaces=typename tuple_space_and_quadrature_rule_types<Parameters...>::type1;
      using QuadratureRules=typename tuple_space_and_quadrature_rule_types<Parameters...>::type2;
      static constexpr Integer Nsubspaces=std::tuple_size<BaseFunctionSpaces>::value;


      template<Integer N,Integer M=0>
      class ShapeFunctionTupleType
      {    

      public:
            using rest = typename ShapeFunctionTupleType<N+1>::type; 

            using BaseFunctionSpace=typename get<N, BaseFunctionSpaces>::type;
            using QuadratureRule=typename get<N, QuadratureRules>::type;

            using single_type=std::tuple<ShapeFunctionOperator<QuadratureRule,Elem, BaseFunctionSpace>>;
            using type = decltype( std::tuple_cat( std::declval< single_type >(), 
                                                   std::declval< rest >() ) );
      };


      template<Integer M>
      class ShapeFunctionTupleType<Nsubspaces-1,M>
      {
       public:
           static constexpr Integer N=Nsubspaces-1;
           using BaseFunctionSpace=typename get<N, BaseFunctionSpaces>::type;
           using QuadratureRule=typename get<N, QuadratureRules>::type;
           using single_type=std::tuple<ShapeFunctionOperator<QuadratureRule,Elem, BaseFunctionSpace>>;
           using type = typename std::tuple<ShapeFunctionOperator<QuadratureRule,Elem, BaseFunctionSpace>>;
      };

      
     using ShapeFunctionType=typename ShapeFunctionTupleType<0>::type;


      template<Integer N=0>
      typename std::enable_if< N==Nsubspaces-1,typename ShapeFunctionTupleType<N>::type>::type 
      set_quadrature_rule()
      {   
           using BaseFunctionSpace=typename get<N, BaseFunctionSpaces>::type;
           using QuadratureRule=typename get<N, QuadratureRules>::type;    
           ShapeFunctionOperator<QuadratureRule,Elem,BaseFunctionSpace> shape_function;
           return std::tuple<decltype(shape_function)> (shape_function);
          }

      template<Integer N=0>
      typename std::enable_if< N<Nsubspaces-1,typename ShapeFunctionTupleType<N>::type>::type 
      set_quadrature_rule()
      {       
           using BaseFunctionSpace=typename get<N, BaseFunctionSpaces>::type;
           using QuadratureRule=typename get<N, QuadratureRules>::type;    
           ShapeFunctionOperator<QuadratureRule,Elem,BaseFunctionSpace> shape_function;
           return std::tuple_cat(std::make_tuple(shape_function),set_quadrature_rule<N+1>()  );
      }

      FEProva(const MeshT& mesh):
      shape_functions_(set_quadrature_rule())
      {}



      const ShapeFunctionType& shape_functions()const{return shape_functions_;}
            ShapeFunctionType  shape_functions()     {return shape_functions_;}

  private:
  ShapeFunctionType shape_functions_;


};









template<typename...Args>
class MixedSpace
{
public:
using type_dofmap=std::tuple<typename Args::type_dofmap...>;

inline const Integer& n_dofs()const{return n_dofs_;}; 

inline const type_dofmap& dofmap()const{return dofmap_;};

template<typename OtherArg,typename...OtherArgs>
class tot_subspaces
{ public: static constexpr Integer value= OtherArg::Nsubspaces+tot_subspaces<OtherArgs...>::value;};

template<typename OtherArg>
class tot_subspaces<OtherArg>
{ public:  static constexpr Integer value= OtherArg::Nsubspaces;};


template<Integer N>
typename std::enable_if< 0==N, Integer >::type
tot_n_dofs()
{ const auto& tmp=std::get<N>(spaces_);
  return tmp->n_dofs();}


template<Integer N>
typename std::enable_if< 0<N, Integer >::type
tot_n_dofs()
{ 
static_assert(N>0, " tuple cannot have negative components");
const auto& tmp=std::get<N>(spaces_);
return tmp->n_dofs()+tot_n_dofs<N-1>();}



template<typename OtherArg,typename...OtherArgs>
struct tuple_type
{
      using rest = typename tuple_type<OtherArgs...>::type; 
      using tuple_ens=std::tuple<OtherArg>;
      using type = decltype( std::tuple_cat( std::declval< tuple_ens >(), std::declval< rest >() ) );
};


template<typename OtherArg>
struct tuple_type<OtherArg>
{
     using type = typename std::tuple<OtherArg>;
};


template<Integer N,typename OtherArg, typename...OtherArgs>
typename std::enable_if< 0==sizeof...(OtherArgs),typename tuple_type<OtherArg,OtherArgs...>::type >::type  
tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
{ std::cout<<"N=="<<N<<". "<<tot_n_dofs<N-1>()<<std::endl;
  return std::tuple<OtherArg>(add_costant(otherarg,tot_n_dofs<N-1>()));}


template<Integer N,typename OtherArg, typename...OtherArgs>
typename std::enable_if< 0<N && 0<sizeof...(OtherArgs),typename tuple_type<OtherArg,OtherArgs...>::type >::type  
tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
{std::cout<<"N=="<<N<<". "<<tot_n_dofs<N-1>()<<std::endl;
  return std::tuple_cat(std::tuple<OtherArg>( add_costant(otherarg,tot_n_dofs<N-1>()) ),
                       tuple_make<N+1,OtherArgs...>(otherargs...));}

template<Integer N,typename OtherArg, typename...OtherArgs>
typename std::enable_if< 0==N && 0<sizeof...(OtherArgs),typename tuple_type<OtherArg,OtherArgs...>::type >::type  
tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
{std::cout<<"N=="<<N<<". "<<0<<std::endl;
  return std::tuple_cat(std::tuple<OtherArg>( otherarg ),
                       tuple_make<N+1,OtherArgs...>(otherargs...));}

MixedSpace(const Args&...args):
spaces_(std::make_tuple(std::make_shared<Args>(args)...)),
n_dofs_(tot_n_dofs<sizeof...(Args)-1>()),
dofmap_(tuple_make<0,typename Args::type_dofmap...>(args.dofmap()...))
{}

static constexpr Integer Nsubspaces=tot_subspaces<Args...>::value;
private:
      std::tuple<std::shared_ptr<Args>...> spaces_;
      Integer n_dofs_;
      type_dofmap dofmap_;

      // std::array<std::vector<Integer>, Nsubspaces> offset_;
      // std::array<std::vector<std::vector<Integer>>, Nsubspaces> space_dofs_;
      // std::array<std::array<Integer,4>,Nsubspaces> space_infos_;


};

template<typename...Args>
MixedSpace<Args...> MixedFunctionSpace(const Args&...args){return MixedSpace<Args...>(args...);};



template<typename...Args>
class TestSpace
{
public:
 TestSpace(const MixedSpace<Args...>& spaces):
 spaces_(spaces)
 {}

 private:
   MixedSpace<Args...> spaces_;
};


template<typename...Args>
TestSpace<Args...> Test(const MixedSpace<Args...>& spaces){return TestSpace<Args...>(spaces);};














class Leaf;

 class Tree
{
public:
      Tree():
      leaf_(0)
      {};

      Tree(Integer i):
      leaf_(i)
      {};
// void print(){
//   for(Integer nn=0;nn< vec_.size();nn++)
//       vec_[nn]->print();}


//virtual void add(Tree& ifs) {vec_.push_back(std::make_shared<Tree>(ifs));};
// virtual void add(Tree& ifs) {vec_.push_back(ifs);};
// void add(Leaf& ifs) {children.push_back(ifs);};
     void add(const Tree& c)
    {
        c.print();
        children.push_back(std::make_shared<Tree>(c));
    } 

    //  void add(std::shared_ptr<Leaf> c)
    // {
    //     children.push_back(c);
    // }


// std::shared_ptr<Tree> operator[] (int x) {
//           return vec_[x];}
void print()const {
  for(Integer nn=0;nn<children.size();nn++)
      {//std::cout<<"----------------TREE == "<<leaf_<<std::endl;
        children[nn]->print();}
}
//void print(){std::cout<<"----------------LEAF == "<<leaf_<<std::endl;}
Tree& operator[] (int x) {
          return *children[x];}
private:
      std::vector<std::shared_ptr<Tree>> children; 
      Integer leaf_;
      //std::vector<Tree> vec_; 
};


class Leaf : public Tree
{
public:
      Leaf(Integer i):
      leaf_(i)
      {children.push_back(std::make_shared<Tree>(Tree(i)));};

      void print(){std::cout<<"----------------LEAF == "<<leaf_<<std::endl;}
      Integer val(){return leaf_;};
      Leaf operator[] (int x) {
                return *this;
            }

private:
      Integer leaf_;
      std::vector<std::shared_ptr<Tree>> children; 
};











class Component
{
  public:
    virtual void traverse() = 0;
    virtual~Component();
    //virtual Component& operator[] (int x);
      
};
class Primitive: public Component
{
    int value;
  public:
    Primitive(int val)
    {
        value = val;
    }
    void traverse() override
    {
        std::cout << " primitive=="<<value<<" ";
    }

   // Component& operator[] (int x) override {
   //        return *this;
   //    }


};

class Composite: public Component
{
    std::vector <std::shared_ptr< Component > > children;
    int value;
  public:
    Composite(int val)
    {
        value = val;
    }


     void add(const Component&);

     void add(const Primitive& c)
    {
        children.push_back(std::make_shared<Primitive>(c));
    } 

     void add(const Composite& c)
    {
        children.push_back(std::make_shared<Composite>(c));
    }   


    void add(std::shared_ptr<Component> c)
    {
        children.push_back(c);
    }

    //template<typename Type, typename...Types>

     template<typename...Args>
    typename std::enable_if<0==sizeof...(Args), void>::type
    add(const Primitive& t,Args...more)
    {
     children.push_back(std::make_shared<Primitive>(t));
    };

    template<typename...Args>
    typename std::enable_if<0==sizeof...(Args), void>::type
    add(const Composite& t,Args...more)
    {
     children.push_back(std::make_shared<Composite>(t));
    };


    template<typename...Args>
    typename std::enable_if<0<sizeof...(Args), void>::type
    add(const Primitive& t,Args...more)
    {
      children.push_back(std::make_shared<Primitive>(t));
      add(more...);
    };


    template<typename...Args>
    typename std::enable_if<0<sizeof...(Args), void>::type
    add(const Composite& t,Args...more)
    {
      children.push_back(std::make_shared<Composite>(t));
      add(more...);
    };






    void traverse() override
    {
        std::cout << " composite=="<< value;
        for (int i = 0; i < children.size(); i++)
          children[i]->traverse();
    }


//    // Component& operator[] (int x) {
//    //        return *children[x];
//    //    }

   std::shared_ptr< Component > operator[] (int x) {
          return children[x];
      }

};






// template< typename VectorT, typename FS,typename...FSs>
// typename std::enable_if<0<sizeof...(FSs),void>::type
// init_functionspacesystem(VectorT& vec,FS fs, FSs... fss)
// {
// vec.push_back(std::make_shared<FS>(fs));  
// init_functionspacesystem<FSs...>(vec,fss...);

// };

// template< typename VectorT, typename FS,typename...FSs,Integer>
// typename std::enable_if<0==sizeof...(FSs),void>::type
// init_functionspacesystem(VectorT& vec, FS fs, FSs... fss)
// {
// vec.push_back(std::make_shared<FS>(fs));  
// };

// template<typename FS,typename...FSs>
// FunctionSpaceSystem(FunctionSpace fs, FunctionSpaces... fss)
// {
// init_functionspacesystem<FunctionSpace,FunctionSpaces...>(nature_,fs,fss...);
// };




}


#endif