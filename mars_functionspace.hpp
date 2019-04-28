#ifndef MARS_FunctionSpace_HPP
#define MARS_FunctionSpace_HPP


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Written by Gabriele Rovi (April 2019)                                                                             ////////
////// We define the FunctionSpace class:                                                                                ////////
////// 1) It takes a mesh and 1 or more FunctionSpaces (Lagrange1<2>, RT0<1>...)                                         ////////                                        
////// 2) It builds the dofmap: a vector (long n_elements), whose component is the array of all the dofs of the element  ////////
////// 2) dofmap(space_id,elem_id) returns the dofs of element elem_id corresponding to the space space_id               ////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "mars_base.hpp"
#include "mars_elementfunctionspace.hpp"
#include "mars_functionspace_dofmap.hpp"

namespace mars{




// template<typename MeshT, typename BaseFunctionSpace,typename...BaseFunctionSpaces>
// class FunctionSpace
// {
// public:
// static constexpr Integer Nsubspaces=1+sizeof...(BaseFunctionSpaces);
// static constexpr Integer Nelem_dofs=DofsPerElemNums< typename MeshT::Elem,BaseFunctionSpace,BaseFunctionSpaces...>::value;

// inline const Integer n_subspaces(){return Nsubspaces;};
// inline const Integer n_elem_dofs(){return Nelem_dofs;};
// inline const Integer n_elem_dofs(Integer space_id){
//                             const auto& os=offset_[space_id];
//                            const auto size=os[os.size()-1]-os[0];
//                            return size;}

// inline const Integer n_dofs(){return n_dofs_;};
// inline const std::vector<std::array<Integer, Nelem_dofs>> dofmap(){return dofmap_;};
// inline const std::array<Integer, Nelem_dofs> dofmap(Integer elem_id){return dofmap_[elem_id];};

// inline const std::vector<Integer> dofmap(Integer space_id,Integer elem_id){
//             const auto& os=offset_[space_id];
//             const auto& size=n_elem_dofs(space_id);//os[os.size()-1]-os[0];
//             std::vector<Integer> output(size);
//             for(Integer nn=0;nn<size;nn++)
//                  output[nn]=dofmap_[elem_id][nn+os[0]];

//             return output;
//             };



// inline const std::array<std::vector<Integer>, Nsubspaces> offset(){return offset_;};
// inline const std::vector<Integer> offset(Integer space_id){return offset_[space_id];};
// inline const std::vector<Integer> space_dofs(Integer space_id){return space_dofs_[space_id];};
// inline void  space_dof(Integer space_id,std::vector<Integer>& spacedofs){spacedofs=space_dofs_[space_id];};


// FunctionSpace(const MeshT& mesh):
// mesh_(std::make_shared< MeshT >(mesh))
// {
// dofmap_fespace<BaseFunctionSpace,BaseFunctionSpaces...>(mesh,dofmap_,offset_,n_dofs_,space_dofs_);
// };

// private:
// std::shared_ptr< MeshT > mesh_;
// Integer n_dofs_;
// std::vector<std::array<Integer, Nelem_dofs>> dofmap_;
// std::array<std::vector<Integer>, Nsubspaces> offset_;
// std::array<std::vector<Integer>, Nsubspaces> space_dofs_;

// };



template<typename MeshT, typename BaseFunctionSpace,typename...BaseFunctionSpaces>
class FunctionSpace
{
public:
static constexpr Integer Nsubspaces=1+sizeof...(BaseFunctionSpaces);
static constexpr Integer Nelem_dofs=DofsPerElemNums< typename MeshT::Elem,BaseFunctionSpace,BaseFunctionSpaces...>::value;

inline const Integer n_subspaces(){return Nsubspaces;};
inline const Integer n_elem_dofs(){return Nelem_dofs;};
inline const Integer n_elem_dofs(Integer space_id){
                            const auto& os=offset_[space_id];
                           const auto size=os[os.size()-1]-os[0];
                           return size;}

inline const Integer n_dofs(){return n_dofs_;};
inline void n_dofs(const Integer& ndofs){ndofs=n_dofs_;};

inline const std::vector<std::array<Integer, Nelem_dofs>> dofmap(){return dofmap_;};
inline void  dofmap(const std::vector<std::array<Integer, Nelem_dofs>>& dm){dm=dofmap_;};


inline const std::array<Integer, Nelem_dofs> dofmap(Integer elem_id){return dofmap_[elem_id];};
inline void  dofmap(Integer elem_id, const std::array<Integer, Nelem_dofs> & elem_dm){elem_dm=dofmap_[elem_id];};


inline const std::vector<Integer> dofmap(Integer space_id,Integer elem_id){
            const auto& os=offset_[space_id];
            const auto& size=n_elem_dofs(space_id);
            std::vector<Integer> output(size);
            for(Integer nn=0;nn<size;nn++)
                 output[nn]=dofmap_[elem_id][nn+os[0]];
            return output;
            };

inline void dofmap(Integer space_id,Integer elem_id, const std::vector<Integer>& elem_space_dm)
                 {
                  const auto& os=offset_[space_id];
                  const auto& size=n_elem_dofs(space_id);
                  elem_space_dm.resize(size);
                  for(Integer nn=0;nn<size;nn++)
                       elem_space_dm[nn]=dofmap_[elem_id][nn+os[0]];
                 };

inline const std::array<std::vector<Integer>, Nsubspaces> offset(){return offset_;};
inline void  offset(const std::array<std::vector<Integer>, Nsubspaces> &os){os=offset_;};


inline const std::vector<Integer> offset(Integer space_id){return offset_[space_id];};
inline void offset(Integer space_id, const std::vector<Integer>& space_os){space_os=offset_[space_id];};


inline const std::vector<Integer> space_dofs(Integer space_id){return space_dofs_[space_id];};
inline void  space_dof(Integer space_id,std::vector<Integer>& spacedofs){spacedofs=space_dofs_[space_id];};


FunctionSpace(const MeshT& mesh):
mesh_(std::make_shared< MeshT >(mesh))
{
dofmap_fespace<BaseFunctionSpace,BaseFunctionSpaces...>(mesh,dofmap_,offset_,n_dofs_,space_dofs_);
};

private:
std::shared_ptr< MeshT > mesh_;
Integer n_dofs_;
std::vector<std::array<Integer, Nelem_dofs>> dofmap_;
std::array<std::vector<Integer>, Nsubspaces> offset_;
std::array<std::vector<Integer>, Nsubspaces> space_dofs_;

};




}


#endif