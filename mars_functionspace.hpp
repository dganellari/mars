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

namespace mars{



class InterfaceFunctionSpace
{
public:
 virtual void add(InterfaceFunctionSpace *ifs) {};
 virtual ~InterfaceFunctionSpace(){};
private:
};

template<typename MeshT, typename BaseFunctionSpace,typename...BaseFunctionSpaces>
class FunctionSpace : public InterfaceFunctionSpace
{
public:
      using Elem= typename MeshT::Elem;

      static constexpr Integer Nsubspaces=1+sizeof...(BaseFunctionSpaces);

      static constexpr Integer Nelem_dofs=DofsPerElemNums< typename MeshT::Elem,BaseFunctionSpace,BaseFunctionSpaces...>::value;

      inline const Integer n_subspaces()const{return Nsubspaces;};

      inline const Integer components (const Integer& space_id)const{return space_infos_[space_id][3];};

      inline const Integer n_elem_dofs()const{return Nelem_dofs;};

      inline const Integer n_elem_dofs(const Integer& space_id)const{
                                  const auto& os=offset_[space_id];
                                  const auto size=os[os.size()-1]-os[0];
                                  return size;}

      inline const Integer n_elem_dofs(const Integer& space_id,const Integer& component_id)const{
                                  const auto& size=n_elem_dofs(space_id);
                                  return (size/space_infos_[space_id][3]);}


      inline const Integer n_dofs()const{return n_dofs_;};

      inline const Integer n_dofs(const Integer& space_id,const Integer& component_id)const
                                 {return space_dofs_[space_id][component_id].size(); };


      inline const std::vector<std::array<Integer, Nelem_dofs>> dofmap()const{return dofmap_;};

      inline void  dofmap(const std::vector<std::array<Integer, Nelem_dofs>>& dm)const{dm=dofmap_;};


      inline const std::array<Integer, Nelem_dofs> dofmap(const Integer& elem_id)const
                         {return dofmap_[elem_id];};

      inline void  dofmap(const Integer& elem_id, const std::array<Integer, Nelem_dofs> & elem_dm)const
                         {elem_dm=dofmap_[elem_id];};


      inline const std::vector<Integer> 
                   dofmap(const Integer& space_id,const Integer& elem_id)const{
                        const auto& os=offset_[space_id];
                        const auto& size=n_elem_dofs(space_id);
                        std::vector<Integer> output(size);
                        for(Integer nn=0;nn<size;nn++)
                             output[nn]=dofmap_[elem_id][nn+os[0]];
                        return output;};

      inline const std::vector<Integer> 
                   dofmap(const Integer& space_id,const Integer& component_id,const Integer& elem_id)const{
                        const auto& os=offset_[space_id];
                        const auto& size= n_elem_dofs(space_id);
                        std::vector<Integer> output(size);
                        const auto& comp=components(space_id);
                        space_infos_[space_id][3];
                        for(Integer nn=component_id;nn<size;nn=nn+comp)
                             output[nn]=dofmap_[elem_id][nn+os[0]];
                        return output;};


      inline void dofmap(const Integer& space_id,const Integer& elem_id, const std::vector<Integer>& elem_space_dm)const
                       {
                        const auto& os=offset_[space_id];
                        const auto& size=n_elem_dofs(space_id);
                        elem_space_dm.resize(size);
                        for(Integer nn=0;nn<size;nn++)
                             elem_space_dm[nn]=dofmap_[elem_id][nn+os[0]];
                       };

      inline const std::array<std::vector<Integer>, Nsubspaces> offset() const {return offset_;};

      inline void  offset(const std::array<std::vector<Integer>, Nsubspaces> &os)const {os=offset_;};


      inline const std::vector<Integer> offset(const Integer& space_id)const{return offset_[space_id];};

      inline void offset(Integer space_id, const std::vector<Integer>& space_os)const {space_os=offset_[space_id];};


      inline const std::vector<Integer>& space_dofs(const Integer& space_id,const Integer& component_id) const
                                         {return space_dofs_[space_id][component_id];};

      inline void space_dofs(const Integer& space_id, const Integer& component_id,std::vector<Integer>& spacedofs)const
                            {spacedofs.resize(n_dofs(space_id,component_id));
                             spacedofs=space_dofs_[space_id][component_id];};

      inline const std::array<std::array<Integer,4>,Nsubspaces> space_info()const{return space_infos_;};

      inline const std::shared_ptr< MeshT > mesh()const {return mesh_;};

      FunctionSpace(const MeshT& mesh):
      mesh_(std::make_shared< MeshT >(mesh))
      {
      function_space_info<0,Nsubspaces,BaseFunctionSpace,BaseFunctionSpaces...>(space_infos_);
      dofmap_fespace<BaseFunctionSpace,BaseFunctionSpaces...>(mesh,dofmap_,offset_,n_dofs_,space_infos_,space_dofs_);     
      };

private:
      std::shared_ptr< MeshT > mesh_;
      Integer n_dofs_;
      std::vector<std::array<Integer, Nelem_dofs>> dofmap_;
      std::array<std::vector<Integer>, Nsubspaces> offset_;
      std::array<std::vector<std::vector<Integer>>, Nsubspaces> space_dofs_;
      std::array<std::array<Integer,4>,Nsubspaces> space_infos_;
};











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