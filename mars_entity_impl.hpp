#ifndef MARS_ENTITY_IMPL_HPP
#define MARS_ENTITY_IMPL_HPP

#include "mars_entity.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Written by Gabriele Rovi (April 2019)                                                                                            ////////
////// We define the following class:                                                                                                   ////////
////// Entity<Integer Dim, Integer ManifoldDim, Integer EntityDim>                                                                      ////////
////// Given a mesh, with elements whose dimension is ManifoldDim, defined in a Dim-dimensional space                                   ////////
////// We want to define the entities of dimension 0<=EntityDim<=ManifoldDim                                                            ////////
////// The class has two fields:                                                                                                        ////////                                                        
//////  1) entity_2_elem_: for each entity id, returns a 2-D array A such that:                                                         ////////
//////  A[0]= id of one of the element to wich entity belongs                                                                           ////////
//////  A[1]= local number of the entity inside the element A[0]                                                                        ////////
//////  2) elem_2_entity: for each element it returns a entity_combinations-D array containing all the ids of the type of the entity    ////////                                    ////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace mars {


// template<Integer Dim,Integer ManifoldDim, Integer EntityDim> 
// void ElemEntity<Simplex<Dim,ManifoldDim>,EntityDim>::add_bisection(const Bisection<Mesh<Dim,ManifoldDim>>& bisection)
// {bisection_ptr_=std::make_shared<Bisection<Mesh<Dim,ManifoldDim>>>(bisection);}

// template<Integer Dim,Integer ManifoldDim> 
// void ElemEntity<Simplex<Dim,ManifoldDim>,0>::add_bisection(const Bisection<Mesh<Dim,ManifoldDim>>& bisection)
// {bisection_ptr_=std::make_shared<Bisection<Mesh<Dim,ManifoldDim>>>(bisection);}

// template<Integer Dim,Integer ManifoldDim> 
// void ElemEntity<Simplex<Dim,ManifoldDim>,ManifoldDim>::add_bisection(const Bisection<Mesh<Dim,ManifoldDim>>& bisection)
// {bisection_ptr_=std::make_shared<Bisection<Mesh<Dim,ManifoldDim>>>(bisection);}

// template<Integer Dim,Integer ManifoldDim, Integer EntityDim> 
// void ElemEntity<Simplex<Dim,ManifoldDim>,EntityDim>::add_bisection(const std::shared_ptr<Bisection<Mesh<Dim,ManifoldDim>>> bisection_ptr)
// {bisection_ptr_=bisection_ptr;}


template<Integer Dim,Integer ManifoldDim, Integer EntityDim> 
void ElemEntity<Simplex<Dim,ManifoldDim>,EntityDim>::init(const Mesh<Dim,ManifoldDim>& mesh)
{mesh_ptr_=std::make_shared<Mesh<Dim,ManifoldDim>>(mesh);} 

template<Integer Dim,Integer ManifoldDim, Integer EntityDim> 
ElemEntity<Simplex<Dim,ManifoldDim>,EntityDim>::ElemEntity(const std::shared_ptr<Mesh<Dim,ManifoldDim>>& mesh_ptr):
mesh_ptr_(mesh_ptr)
{} 

 template<Integer Dim,Integer ManifoldDim,Integer EntityDim>
  ElemEntity<Simplex<Dim,ManifoldDim>,EntityDim>::ElemEntity
  (const Mesh<Dim,ManifoldDim>& mesh,const std::vector< std::vector<Integer> >& node_2_element,const Integer level):
       mesh_ptr_(std::make_shared<Mesh<Dim,ManifoldDim>>(mesh))
      {
        init_elem_entity(mesh,node_2_element,level); 
      };


 template<Integer Dim,Integer ManifoldDim,Integer EntityDim>
  ElemEntity<Simplex<Dim,ManifoldDim>,EntityDim>::ElemEntity
  (const std::shared_ptr<Mesh<Dim,ManifoldDim>>& mesh_ptr,const std::vector< std::vector<Integer> >& node_2_element,const Integer level):
       mesh_ptr_(mesh_ptr)
      {
        init_elem_entity(*mesh_ptr,node_2_element,level); 
      };










template<Integer Dim,Integer ManifoldDim,Integer EntityDim>
   void ElemEntity<Simplex<Dim,ManifoldDim>,EntityDim>::init_elem_entity
   (const Mesh<Dim,ManifoldDim>& mesh,
    const std::vector< std::vector<Integer> >& node_2_element,
    const Integer level)
  {



    static_assert(Dim>=0 , " the space dimension must be non negative ");
    static_assert(ManifoldDim>=0 , " the manifold dimension must be non negative ");
    static_assert(EntityDim>=0 , " the entity dimension must be non negative ");    
    static_assert(Dim>=ManifoldDim , " the space dimension must be greater than or equal to the manifold dimension ");
    static_assert(ManifoldDim>=EntityDim , " the manifold dimension must be greater than or equal to the entity dimension ");


    const auto& n_elements = mesh.n_elements();
    Integer entity_number=0;
    //std::array<Integer, EntityDim+1> tmp;
    std::array<Integer, 2> tmp;
    Integer entity_e1[EntityDim+1];
    Integer entity_e2[EntityDim+1];
    Integer entity_e3[EntityDim+1];
    std::array<Integer, EntityDim+1> entity_nodes1;
    std::array<Integer, EntityDim+1> entity_nodes2;
    std::array<Integer, EntityDim+1> entity_nodes3;
    elem_2_entity_.resize(n_elements,std::vector<Integer>());
    // loop on all the elements
    for(Integer elem_iter1 = 0; elem_iter1 < n_elements; ++elem_iter1) 
    {
      // if(!elem_belongs_to_level(mesh_ptr_,elem_iter1,level,bisection_ptr_)) continue;
      // if(mesh.is_active(elem_iter1))
      {

        elem_2_entity_[elem_iter1].resize(entity_combinations());

       const auto &el1 = mesh.elem(elem_iter1);
       const auto &el_nodes1=el1.nodes;
       std::vector<Integer> el_neighs;
        // find all the elements touching the element=elem_iter1
       for(Integer ii=0;ii<ManifoldDim+1;ii++)
         for(Integer jj=0;jj<node_2_element[el_nodes1[ii]].size();jj++)
         {
          el_neighs.push_back(node_2_element[el_nodes1[ii]][jj]);
        } 

        // remove all the elements larger than  elem_iter1
        el_neighs.erase(std::remove_if(el_neighs.begin(), el_neighs.end(), [elem_iter1](int index_tmp) { return index_tmp>=elem_iter1; }), el_neighs.end());       
        // sort and make unique the neighborhood elements
        std::sort( el_neighs.begin(), el_neighs.end() );
        el_neighs.erase( std::unique( el_neighs.begin(), el_neighs.end() ), el_neighs.end() );
        
        // loop on all the entities of the actual element e1
        for(Integer iter_entity_e1=0;iter_entity_e1< entity_combinations();iter_entity_e1++)
        {     
         Combinations<ManifoldDim + 1, EntityDim+1>::generate(iter_entity_e1,entity_e1);

         // build the face nodes of iter_entity_e1
         for(int nn=0;nn<EntityDim+1;nn++)
         {
           entity_nodes1[nn]=el_nodes1[entity_e1[nn]];
         }
         
         std::sort(std::begin(entity_nodes1), std::end(entity_nodes1)); 

         bool found_entity=false; 


           // loop on all the neighborhood elements
         for(Integer elem_iter2 = 0; elem_iter2 < el_neighs.size(); ++elem_iter2) 
         {   
          const auto &el2_neigh=el_neighs[elem_iter2];
          const auto &el2 = mesh.elem(el2_neigh);
          const auto &el_nodes2=el2.nodes;

            // loop on all the entity of the element e2
          for(Integer iter_entity_e2=0;iter_entity_e2<entity_combinations();iter_entity_e2++)
          {

            Combinations<ManifoldDim + 1, EntityDim+1>::generate(iter_entity_e2,entity_e2);
            
            for(int nn=0;nn<EntityDim+1;nn++)
              entity_nodes2[nn]=el_nodes2[entity_e2[nn]];   

            std::sort(std::begin(entity_nodes2), std::end(entity_nodes2));  
            
            if (std::equal(std::begin(entity_nodes1), std::end(entity_nodes1), std::begin(entity_nodes2)))
            {
                       // all the entities already computed on elem=elem_iter2
             const auto &elem2entity=elem_2_entity_[el2_neigh];
                       // loop on each of these, to find the one that equals iter_entity_e1
             for(auto i: elem2entity)
             {
                          // we generate the combination of dimension=EntityDim with local numbering=entity_2_elem_[i][1]
              Combinations<ManifoldDim + 1, EntityDim+1>::generate(entity_2_elem_[i][1],entity_e3);
                          // then we find the sorted nodes
              for(int nn=0;nn<EntityDim+1;nn++)
                entity_nodes3[nn]=mesh.elem(entity_2_elem_[i][0]).nodes[entity_e3[nn]]; 
              std::sort(std::begin(entity_nodes3), std::end(entity_nodes3));  

                          //auto nodes_i=entity_2_elem_[i][0];
              if(std::equal(std::begin(entity_nodes3), std::end(entity_nodes3), std::begin(entity_nodes1)))
              {
               elem_2_entity_[elem_iter1][iter_entity_e1]=i;
               break;
             }
           }
           found_entity=true;
           break;
         }         
       }  

       if(found_entity==true)
         break;

     }
          // if no neighbborhood element share the same entity, then we must add this entity   
     if(found_entity==false)
     {
                   //std::copy(std::begin(entity_nodes1), std::end(entity_nodes1), std::begin(tmp));
                   //entity_.push_back( tmp);
       tmp[0]=elem_iter1; 
       tmp[1]=iter_entity_e1;
       entity_2_elem_.push_back(tmp);
       elem_2_entity_[elem_iter1][iter_entity_e1]=entity_number;                   
       entity_number++;
     }   
   }          
 }    
}  
size_=entity_number;  
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////      ENTITY == NODE                                                                            ////////////
////////////      We consider mesh nodes and we do not renumber them                                        ////////////
////////////      This is why we need a specialization for nodes                                            ////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<Integer Dim,Integer ManifoldDim> 
    ElemEntity<Simplex<Dim,ManifoldDim>,0>::ElemEntity(const Mesh<Dim,ManifoldDim>& mesh):
    mesh_ptr_(std::make_shared<Mesh<Dim,ManifoldDim>>(mesh))
    {
     // init_elem_entity(mesh);//,entity_2_elem_,elem_2_entity_,size_);
    }

template<Integer Dim,Integer ManifoldDim> 
void ElemEntity<Simplex<Dim,ManifoldDim>,0>::init(const Mesh<Dim,ManifoldDim>& mesh)
{mesh_ptr_=std::make_shared<Mesh<Dim,ManifoldDim>>(mesh);} 


template<Integer Dim,Integer ManifoldDim>
ElemEntity<Simplex<Dim,ManifoldDim>,0>::ElemEntity
(const Mesh<Dim,ManifoldDim>& mesh,
 const std::vector< std::vector<Integer> >& node_2_element,
 const Integer level):
     mesh_ptr_(std::make_shared<Mesh<Dim,ManifoldDim>>(mesh))
    { 
      //init_elem_entity(mesh,node_2_element,entity_2_elem_,elem_2_entity_,size_);
      // init_elem_entity(mesh,entity_2_elem_,elem_2_entity_,size_); 
      init_elem_entity(mesh,node_2_element,level); 
    };

template<Integer Dim,Integer ManifoldDim>
void ElemEntity<Simplex<Dim,ManifoldDim>,0>::init_elem_entity
(const Integer level)
{
  static_assert(Dim>=0 , " the space dimension must be non negative ");
  static_assert(ManifoldDim>=0 , " the manifold dimension must be non negative ");
  static_assert(Dim>=ManifoldDim , " the space dimension dimension must be greater than or equal to the manifold dimension ");
  const auto& n_elements = mesh_ptr_->n_elements();
  const auto& n_nodes = mesh_ptr_->n_nodes();

  elem_2_entity_.resize(n_elements, std::vector<Integer>());
  entity_2_elem_.resize(n_nodes);
  for(Integer elem_iter = 0; elem_iter < n_elements; ++elem_iter) 
  {
    // if(!elem_belongs_to_level(mesh_ptr_,elem_iter,level,bisection_ptr_)) continue;
    auto& nodes=mesh_ptr_->elem(elem_iter).nodes;
    elem_2_entity_[elem_iter].resize(entity_combinations());
    // elem_2_entity_[elem_iter] = nodes;
    for(Integer i=0;i<elem_2_entity_[elem_iter].size();i++)
      elem_2_entity_[elem_iter][i] = nodes[i];

    for(Integer node_iter=0;node_iter<elem_2_entity_[elem_iter].size();node_iter++)
    {
      entity_2_elem_[nodes[node_iter]][0]=elem_iter;
      entity_2_elem_[nodes[node_iter]][1]=node_iter;
    }
  }
  size_=n_nodes;

};


template<Integer Dim,Integer ManifoldDim>
void ElemEntity<Simplex<Dim,ManifoldDim>,0>::init_elem_entity
 (const Mesh<Dim,ManifoldDim>& mesh,
  const std::vector< std::vector<Integer> >& node_2_element,
  const Integer level)
{
  static_assert(Dim>=0 , " the space dimension must be non negative ");
  static_assert(ManifoldDim>=0 , " the manifold dimension must be non negative ");
  static_assert(Dim>=ManifoldDim , " the space dimension dimension must be greater than or equal to the manifold dimension ");
  const auto& n_elements = mesh.n_elements();
  const auto& n_nodes = mesh.n_nodes();

  if(bisection_ptr_==nullptr && level>-1)
  {
   throw std::invalid_argument("ElemConnectivity ERROR: if we consider levels, then a bisection pointer must be provided");
  }

  elem_2_entity_.resize(n_elements, std::vector<Integer>());
  entity_2_elem_.resize(n_nodes);
  for(Integer elem_iter = 0; elem_iter < n_elements; ++elem_iter) 
  {
    // std::cout<<"belongs = "<<elem_belongs_to_level(mesh_ptr_,elem_iter,level,bisection_ptr_)<<std::endl;
    // if(!elem_belongs_to_level(mesh_ptr_,elem_iter,level,bisection_ptr_)) continue;
    // std::cout<<elem_iter<<std::endl;
    // std::cout<<"elem ="<<elem_iter<<" "<<elem_belongs_to_level(mesh_ptr_,elem_iter,level,bisection_ptr_)<<std::endl;
    auto& nodes=mesh.elem(elem_iter).nodes;
    elem_2_entity_[elem_iter].resize(entity_combinations());

    // elem_2_entity_[elem_iter] = nodes;
    for(Integer i=0;i<elem_2_entity_[elem_iter].size();i++)
      elem_2_entity_[elem_iter][i] = nodes[i];

    for(Integer node_iter=0;node_iter<elem_2_entity_[elem_iter].size();node_iter++)
    {
      entity_2_elem_[nodes[node_iter]][0]=elem_iter;
      entity_2_elem_[nodes[node_iter]][1]=node_iter;
    }
  }
  size_=n_nodes;

};





 













////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////      ENTITY == ELEM                                                                            ////////////
////////////      We consider mesh elems and we do not renumber them                                        ////////////
////////////      This is why we need a specialization for elems                                            ////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<Integer Dim,Integer ManifoldDim>
ElemEntity<Simplex<Dim,ManifoldDim>,ManifoldDim>::ElemEntity
(const Mesh<Dim,ManifoldDim>& mesh,const std::vector< std::vector<Integer> >& node_2_element,const Integer level):
       mesh_ptr_(std::make_shared<Mesh<Dim,ManifoldDim>>(mesh))
      {init_elem_entity(mesh,node_2_element,level);}//,entity_2_elem_,elem_2_entity_,size_); };


template<Integer Dim,Integer ManifoldDim>
void ElemEntity<Simplex<Dim,ManifoldDim>,ManifoldDim>::init_elem_entity
 (const Mesh<Dim,ManifoldDim>& mesh,
  const std::vector< std::vector<Integer> >& node_2_element,
  const Integer level)
{

  static_assert(Dim>=0 , " the space dimension must be non negative ");
  static_assert(ManifoldDim>=0 , " the manifold dimension must be non negative ");
  static_assert(Dim>=ManifoldDim , " the space dimension dimension must be greater than or equal to the manifold dimension ");
  const auto& n_elements = mesh.n_elements();

  
  elem_2_entity_.resize(n_elements,std::vector<Integer>());
  entity_2_elem_.resize(n_elements);
  for(Integer elem_iter = 0; elem_iter < n_elements; ++elem_iter) 
  {

    // if(!elem_belongs_to_level(mesh_ptr_,elem_iter,level,bisection_ptr_)) continue;
    auto& nodes=mesh.elem(elem_iter).nodes;
    elem_2_entity_[elem_iter].resize(entity_combinations());
    elem_2_entity_[elem_iter][0] = elem_iter;
    entity_2_elem_[elem_iter][0]=elem_iter;
    entity_2_elem_[elem_iter][1]=0;

  }
  size_=n_elements;

};


template<Integer Dim,Integer ManifoldDim> 
void ElemEntity<Simplex<Dim,ManifoldDim>,ManifoldDim>::init(const Mesh<Dim,ManifoldDim>& mesh)
{mesh_ptr_=std::make_shared<Mesh<Dim,ManifoldDim>>(mesh);} 







}
#endif //MARS_ENTITY_HPP