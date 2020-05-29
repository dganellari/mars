#ifndef MARS_CONTEXT_HPP
#define MARS_CONTEXT_HPP
#include "mars_base.hpp"
#include "mars_shape_functions_collection.hpp"
#include "mars_dirichlet_bc.hpp"
#include "mars_general_form.hpp"
#include "mars_evaluation.hpp"

namespace mars {

template<typename...Ts>
class Context;


template<typename BilinearForm, typename LinearForm, typename...DirichletBCs>
class Context<BilinearForm,LinearForm,DirichletBCs...>
{

 public:
    using Bilinear=GeneralForm<BilinearForm>;
    using Linear=GeneralForm<LinearForm>;
    using Coefficients=ShapeFunctionCoefficientsCollection<Bilinear,Linear>;
    using Maps=MapFromReferenceCollection<Bilinear,Linear>;
    using ShapeFunctions=ShapeFunctionsCollection<Bilinear,Linear>;
    using BCs=DirichletBoundaryConditionCollection<DirichletBCs...>;
    using MeshT=typename Bilinear::FunctionSpace::MeshT;
    using Elem= typename MeshT::Elem;
    static constexpr auto NLocalDofs=Maps::FunctionSpace::DofsDM::NLocalDofs;

    Context(const BilinearForm& bilinearform, const LinearForm&linearform, const DirichletBCs&...bcs):
    bilinear_form_(general_form(bilinearform)),
    linear_form_(general_form(linearform)),
    shape_coefficients_(shape_function_coefficients(bilinear_form_.spaces_ptr()->mesh(),bilinear_form_,linear_form_)),
    reference_maps_(reference_maps2(bilinear_form_,linear_form_)),
    shapefunctions_(shape_functions(shape_coefficients_,reference_maps_,bilinear_form_,linear_form_)),
    eval_bilinear_form_(Eval(bilinear_form_,shapefunctions_)),
    eval_linear_form_(Eval(linear_form_,shapefunctions_)),
    bcs_(bcs...)
    {
      // LinearForm ededee(5);

      // decltype(linearform.spaces_ptr()) dedeede4(33);
    }
    

    auto full_spaces_ptr(){return bilinear_form_.spaces_ptr();}
    auto full_spaces_ptr()const{return bilinear_form_.spaces_ptr();}

    
    template<typename SystemMat, typename Rhs>
    void assembly(SystemMat& A, Rhs& b,const Integer level=-1)
    {
      std::cout<<"assembly inside 1"<<std::endl;
     auto spaces_ptr=bilinear_form_.spaces_ptr()->spaces_ptr();
     auto& bisection=spaces_ptr->bisection();
     auto& tracker=bisection.tracker();
     // n_dofs_=spaces_ptr->n_dofs();
     std::cout<<"assembly inside"<<std::endl;
     if(level==-1)
      level_=tracker.current_iterate()-1;
    else
      level_=level;


     // //std::cout<<"level="<<level<<std::endl;
     // //std::cout<<"level="<<level_<<std::endl;
     auto& level_cumultive_n_dofs=full_spaces_ptr()->dofsdofmap().level_cumultive_n_dofs();
     std::cout<<"level_cumultive_n_dofs="<<std::endl;
     for(int i=0;i<level_cumultive_n_dofs.size();i++)
     {
      // for(int j=0; j<n_dofs_arr[i].size();j++)
        // //std::cout<<level_cumultive_n_dofs[i]<<" ";
      // //std::cout<<std::endl;
     }



     n_dofs_=0;
     // for(int i=0;i<n_dofs_arr.size();i++)
     // {
     //    n_dofs_+=n_dofs_arr[i][level_];
     // }
     // n_dofs_=n_dofs_arr[level_];
     // auto level_cumultive_n_dofs=spaces_ptr->dofsdofmap().level_cumultive_n_dofs();
     // //std::cout<<" level_cumultive_n_dofs =====>>>> "<<std::endl;
     // for(Integer i=0;i<level_cumultive_n_dofs.size();i++)
     //    //std::cout<<level_cumultive_n_dofs[i]<<" ";
     // //std::cout<<std::endl;

     n_dofs_=level_cumultive_n_dofs[level_];

     // n_dofs_=spaces_ptr->level_n_dofs_array(level_);
     n_dofs_=level_cumultive_n_dofs[level_];
     std::cout<<"n_dofs_="<<n_dofs_<<std::endl;
     auto mesh_ptr=full_spaces_ptr()->mesh_ptr();
     auto& mesh=full_spaces_ptr()->mesh();
     
     auto& n2e=full_spaces_ptr()->node2elem();

     // FiniteElem<Elem> FE(mesh_ptr);
     FiniteElem<Elem> FE(mesh);


     // auto& node_2_elem=full_spaces_ptr()->node2elem();
    // NodeToElem<MeshT> node_2_elem(mesh);
     // const auto& node2elem=node_2_elem.val();
     Integer max_cols=n2e.max_n_nodes();
     // //std::cout<<"------_______-----max_cols"<<max_cols<<std::endl;
     // // //std::cout<<node2elem.size()<<std::endl;
     // for(Integer i=0;i<node2elem.size();i++)
     // {
     //    const auto & n2e=node2elem[i];
     //  // //std::cout<<n2e.size()<<std::endl;
     //  // for(Integer j=0;j<n2e.size();j++)
     //    // //std::cout<<n2e[j]<<" "<<std::endl;
     //  // //std::cout<<std::endl;
     //  // if()
     //  if(max_cols<n2e.size())
     //  max_cols=n2e.size(); 
     // }
     // //std::cout<<" max_cols "<<max_cols<<std::endl;
     // //std::cout<<" max_n_nodes "<<n2e.max_n_nodes()<<std::endl;
     max_cols=min(NLocalDofs*max_cols,n_dofs_);
     


     // A.resize(n_dofs_,std::vector<Real>(n_dofs_));
     // //std::cout<<"------_______----- A init"<<std::endl;
     // //std::cout<<"n_dofs_=="<<n_dofs_<<std::endl;
     // //std::cout<<"NLocalDofs=="<<NLocalDofs<<std::endl;
     // //std::cout<<"max_cols=="<<max_cols<<std::endl;
     // NLocalDofs*max_number of_elements of a node
     A.init(n_dofs_,n_dofs_,max_cols);

     std::cout<<"------_______----- b init"<<std::endl;
     b.resize(n_dofs_);   
     constrained_dofs_.clear();
     constrained_mat_.clear();
     constrained_vec_.clear();

     constrained_dofs_.resize(n_dofs_,false);
     constrained_mat_.resize(n_dofs_,0);
     constrained_vec_.resize(n_dofs_,0);

     // shape_coefficients_.init(*mesh_ptr);
     // mesh.init_signed_normal();
     // auto& signed_normal=mesh.signed_normal();
     shape_coefficients_.init();
       // for(std::size_t el=0;el<mesh_ptr->n_elements();el++)
       for(std::size_t el=0;el<mesh.n_elements();el++)
       {
          // if(!mesh_ptr->is_active(el)) continue;
          // //std::cout<<"level=="<<level_<<std::endl;
          // if(!elem_belongs_to_level(mesh_ptr,el,level_,tracker)) continue;
          if(!elem_belongs_to_level(mesh,el,level_,tracker)) continue;

          // std::cout<<"------_______----- ELEMENT ID = "<<el<<". -----_______--------"<<std::endl;
          FE.init(el,level_);
          // //std::cout<<"fe jac"<<FE.jac()<<std::endl;
          // //std::cout<<"------_______----- qui1 -----_______--------"<<std::endl;
          shape_coefficients_.init(el);
          // //std::cout<<"------_______----- qui2 -----_______--------"<<std::endl;
          reference_maps_.init(FE);
          // //std::cout<<"------_______----- qui3 -----_______--------"<<std::endl;
          shapefunctions_.init(FE);
          // //std::cout<<"------_______----- qui4 -----_______--------"<<std::endl;
          eval_bilinear_form_.apply(A,FE);
          // //std::cout<<"------_______----- qui5 -----_______--------"<<std::endl;
          eval_linear_form_.apply(b,FE);




          // std::cout<<"------_______----- SURFACE-----_______--------"<<std::endl;

          if(FE.is_on_boundary())
          {
            for(std::size_t s=0;s<FE.n_side();s++)
              {
                // controlla, qui passiamo side_id, ma dovremmo avere label
                // dovresti costruire mesh coi label

                // std::cout<<"------_______----- BEGIN SIDE===="<<s<<" tag="<<FE.side_tags()[s]<<std::endl;
                FE.init_boundary(s);
                if(FE.is_side_on_boundary())
                {
                  // std::cout<<"------_______----- reference_maps_===="<<s<<std::endl;
                  reference_maps_.init_boundary(FE);
                  //std::cout<<"------_______----- shapefunctions_===="<<s<<std::endl;
                  shapefunctions_.init_boundary(FE);
                  // std::cout<<"------_______----- FE.elem_id="<<FE.elem_id()<<std::endl;                  
                  // std::cout<<"------_______----- BILINEAR BEGIN SIDE EVAL===="<<s<<" tag="<<FE.side_tags()[s]<<std::endl;
                  eval_bilinear_form_.apply_boundary(A,FE);
                  // A.print_val();
                  // std::cout<<"------_______----- LINEAR BEGIN SIDE EVAL===="<<s<<" tag="<<FE.side_tags()[s]<<std::endl;
                  // auto&elem =mesh.elem(el);
                  // auto&nodes=elem.nodes;
                  // std::cout<<"nodes="<<std::endl;
                  // for(Integer n=0;n<nodes.size();n++)
                  //   std::cout<<nodes[n]<<" ";
                  // std::cout<<std::endl;
                  // std::cout<<"J="<<std::endl;
                  //   std::cout<<FE.jac()<<" ";
                  // std::cout<<std::endl;       
                  // std::cout<<"side_J="<<std::endl;
                  //   std::cout<<FE.jac_side()<<" ";
                  // std::cout<<std::endl;                  
                  eval_linear_form_.apply_boundary(b,FE);
                }
                //std::cout<<"------_______----- END SIDE===="<<s<<std::endl;
              }


            for(std::size_t s=0;s<FE.n_side();s++)
              {
                // controlla, qui passiamo side_id, ma dovremmo avere label
                // dovresti costruire mesh coi label

                // std::cout<<"------_______----- ELEM===="<<el<<std::endl;
                // std::cout<<"------_______----- BEGIN SIDE===="<<s<<std::endl;

                
                if(FE.side_tags()[s]!=INVALID_INDEX)
                {
                  // std::cout<<"------_______----- ELEM===="<<el<<std::endl;
                  // std::cout<<"------_______----- INIT BOUNDARY===="<<FE.side_tags()[s]<<std::endl;
                  FE.init_boundary(s);
                  // std::cout<<"------_______----- END INIT BOUNDARY===="<<s<<std::endl;
                  bcs_.assembly(full_spaces_ptr(),constrained_dofs_,constrained_mat_,constrained_vec_,FE);
                  // reference_maps_.init_boundary(FE);
                  // shapefunctions_.init_boundary(FE);
                  // std::cout<<"------_______-----after bcs_.assembly"<<std::endl;
                  // eval_bilinear_form_.apply_boundary(A,FE);
                  // eval_linear_form_.apply_boundary(b,FE);
                  //std::cout<<"------_______----- END SIDE EVAL===="<<s<<std::endl;
                }
                // //std::cout<<"------_______----- END SIDE===="<<s<<std::endl;
              }
           // constrained_dofs_
         }
        
       }
      // A.print_val();

       // //std::cout<<"------RHS -------"<<std::endl;
       // for(std::size_t i=0;i<n_dofs_;i++)
       //    //std::cout<<b[i]<<std::endl;



       // //std::cout<<"------CONSTRAINED DOFS-------"<<std::endl;
       // for(std::size_t i=0;i<n_dofs_;i++)
       //    //std::cout<<i<<" "<<constrained_dofs_[i]<<std::endl;

       // //std::cout<<"------CONSTRAINED MAT-------"<<std::endl;
       // for(std::size_t i=0;i<n_dofs_;i++)
       //    //std::cout<<i<<" "<<constrained_mat_[i]<<std::endl;

       // //std::cout<<"------CONSTRAINED VEC-------"<<std::endl;
       // for(std::size_t i=0;i<n_dofs_;i++)
       //    //std::cout<<i<<" "<<constrained_vec_[i]<<std::endl;

       //std::cout<<"------END ASSEMBLY -------"<<std::endl;

    }
 


    
    void assembly(std::vector<Real>& b,const Integer level=-1)
    {


     auto spaces_ptr=bilinear_form_.spaces_ptr()->spaces_ptr();
     auto& bisection=spaces_ptr->bisection();
     auto& tracker=bisection.tracker();
     if(level==-1)
      level_=tracker.current_iterate()-1;
    else
      level_=level;

     auto& level_cumultive_n_dofs=full_spaces_ptr()->dofsdofmap().level_cumultive_n_dofs();




     n_dofs_=0;
     n_dofs_=level_cumultive_n_dofs[level_];
     n_dofs_=level_cumultive_n_dofs[level_];
     auto mesh_ptr=full_spaces_ptr()->mesh_ptr();
     auto& mesh=full_spaces_ptr()->mesh();
     
     FiniteElem<Elem> FE(mesh);


     std::cout<<"level = "<<level<<std::endl;
     std::cout<<"n_dofs_ = "<<n_dofs_<<std::endl;




     b.resize(n_dofs_);   
     constrained_dofs_.clear();
     constrained_mat_.clear();
     constrained_vec_.clear();

     constrained_dofs_.resize(n_dofs_,false);
     constrained_mat_.resize(n_dofs_,0);
     constrained_vec_.resize(n_dofs_,0);

     shape_coefficients_.init();
       for(std::size_t el=0;el<mesh.n_elements();el++)
       {
          if(!elem_belongs_to_level(mesh,el,level_,tracker)) continue;
          FE.init(el,level_);
          shape_coefficients_.init(el);
          reference_maps_.init(FE);
          shapefunctions_.init(FE);
          eval_linear_form_.apply(b,FE);

          std::cout<<"el=="<<el<<std::endl;

          if(FE.is_on_boundary())
          {
            for(std::size_t s=0;s<FE.n_side();s++)
              {
                FE.init_boundary(s);
                if(FE.is_side_on_boundary())
                {
                  reference_maps_.init_boundary(FE);
                  shapefunctions_.init_boundary(FE);
                  eval_linear_form_.apply_boundary(b,FE);
                }
              }


            for(std::size_t s=0;s<FE.n_side();s++)
              {              
                if(FE.side_tags()[s]!=INVALID_INDEX)
                {
                  FE.init_boundary(s);
                  bcs_.assembly(full_spaces_ptr(),constrained_dofs_,constrained_mat_,constrained_vec_,FE);
                }
              }
         }      
       }
    }



    template<typename SystemMat>
    void matrix_assembly(SystemMat& A,const Integer level=-1)
    {
     auto spaces_ptr=bilinear_form_.spaces_ptr()->spaces_ptr();
     auto& bisection=spaces_ptr->bisection();
     auto& tracker=bisection.tracker();
     // n_dofs_=spaces_ptr->n_dofs();
     //std::cout<<"assembly"<<std::endl;
     if(level==-1)
      level_=tracker.current_iterate()-1;
    else
      level_=level;


     //std::cout<<"level="<<level<<std::endl;
     //std::cout<<"level="<<level_<<std::endl;
     auto& level_cumultive_n_dofs=full_spaces_ptr()->dofsdofmap().level_cumultive_n_dofs();
     //std::cout<<"level_cumultive_n_dofs="<<std::endl;
     for(int i=0;i<level_cumultive_n_dofs.size();i++)
     {
      // for(int j=0; j<n_dofs_arr[i].size();j++)
        //std::cout<<level_cumultive_n_dofs[i]<<" ";
      //std::cout<<std::endl;
     }



     n_dofs_=0;
     // for(int i=0;i<n_dofs_arr.size();i++)
     // {
     //    n_dofs_+=n_dofs_arr[i][level_];
     // }
     // n_dofs_=n_dofs_arr[level_];
     // auto level_cumultive_n_dofs=spaces_ptr->dofsdofmap().level_cumultive_n_dofs();
     // //std::cout<<" level_cumultive_n_dofs =====>>>> "<<std::endl;
     // for(Integer i=0;i<level_cumultive_n_dofs.size();i++)
     //    //std::cout<<level_cumultive_n_dofs[i]<<" ";
     // //std::cout<<std::endl;

     n_dofs_=level_cumultive_n_dofs[level_];

     // n_dofs_=spaces_ptr->level_n_dofs_array(level_);
     n_dofs_=level_cumultive_n_dofs[level_];
     //std::cout<<"n_dofs_="<<n_dofs_<<std::endl;
     auto mesh_ptr=full_spaces_ptr()->mesh_ptr();
     auto& mesh=full_spaces_ptr()->mesh();
     
     auto& n2e=full_spaces_ptr()->node2elem();

     // FiniteElem<Elem> FE(mesh_ptr);
     FiniteElem<Elem> FE(mesh);


     // auto& node_2_elem=full_spaces_ptr()->node2elem();
    // NodeToElem<MeshT> node_2_elem(mesh);
     // const auto& node2elem=node_2_elem.val();
     Integer max_cols=n2e.max_n_nodes();
     // //std::cout<<"------_______-----llll"<<std::endl;
     // // //std::cout<<node2elem.size()<<std::endl;
     // for(Integer i=0;i<node2elem.size();i++)
     // {
     //    const auto & n2e=node2elem[i];
     //  // //std::cout<<n2e.size()<<std::endl;
     //  // for(Integer j=0;j<n2e.size();j++)
     //    // //std::cout<<n2e[j]<<" "<<std::endl;
     //  // //std::cout<<std::endl;
     //  // if()
     //  if(max_cols<n2e.size())
     //  max_cols=n2e.size(); 
     // }
     // //std::cout<<" max_cols "<<max_cols<<std::endl;
     // //std::cout<<" max_n_nodes "<<n2e.max_n_nodes()<<std::endl;
     max_cols=min(NLocalDofs*max_cols,n_dofs_);
     


     // A.resize(n_dofs_,std::vector<Real>(n_dofs_));
     //std::cout<<"------_______----- A init"<<std::endl;
     //std::cout<<"n_dofs_=="<<n_dofs_<<std::endl;
     //std::cout<<"NLocalDofs=="<<NLocalDofs<<std::endl;
     //std::cout<<"max_cols=="<<max_cols<<std::endl;
     A.init(n_dofs_,n_dofs_,max_cols);

     //std::cout<<"------_______----- b init"<<std::endl;
     // b.resize(n_dofs_);   
     constrained_dofs_.clear();
     constrained_mat_.clear();
     constrained_vec_.clear();

     constrained_dofs_.resize(n_dofs_,false);
     constrained_mat_.resize(n_dofs_,0);
     constrained_vec_.resize(n_dofs_,0);

     shape_coefficients_.init();
     for(std::size_t el=0;el<mesh.n_elements();el++)
       {
        if(!elem_belongs_to_level(mesh,el,level_,tracker)) continue;
        FE.init(el,level_);
        shape_coefficients_.init(el);
        reference_maps_.init(FE);
        shapefunctions_.init(FE);
        eval_bilinear_form_.apply(A,FE);



        if(FE.is_on_boundary())
        {
          for(std::size_t s=0;s<FE.n_side();s++)
            {
               FE.init_boundary(s);
              if(FE.is_side_on_boundary())
              {
                reference_maps_.init_boundary(FE);
                shapefunctions_.init_boundary(FE);
                eval_bilinear_form_.apply_boundary(A,FE);
              }
            }

          for(std::size_t s=0;s<FE.n_side();s++)
            {
            if(FE.side_tags()[s]!=INVALID_INDEX)
              {
                FE.init_boundary(s);
                bcs_.assembly(full_spaces_ptr(),constrained_dofs_,constrained_mat_,constrained_vec_,FE);

              }
            }
         }
        
       }

       //std::cout<<"end assembly"<<std::endl;

    }


    void build_boundary_info(std::vector<bool>& a, std::vector<Real>& b, std::vector<Real>& c,const Integer level=-1)
    {
     auto spaces_ptr=bilinear_form_.spaces_ptr()->spaces_ptr();
     auto& bisection=spaces_ptr->bisection();
     auto& tracker=bisection.tracker();

     if(level==-1)
      level_=tracker.current_iterate()-1;
     else
      level_=level;

     auto& level_cumultive_n_dofs=full_spaces_ptr()->dofsdofmap().level_cumultive_n_dofs();

     n_dofs_=level_cumultive_n_dofs[level_];

     auto& mesh=full_spaces_ptr()->mesh();

     FiniteElem<Elem> FE(mesh);
 
 
     a.clear();
     b.clear();
     c.clear();
     a.resize(n_dofs_,false);
     b.resize(n_dofs_,0);
     c.resize(n_dofs_,0);

     shape_coefficients_.init();
     //std::cout<<"BUILD BOUNDARY INFO"<<std::endl;
       for(std::size_t el=0;el<mesh.n_elements();el++)
       {
          if(!elem_belongs_to_level(mesh,el,level_,tracker)) continue;

          FE.init(el,level_);

          if(FE.is_on_boundary())
          {

            for(std::size_t s=0;s<FE.n_side();s++)
              {            
                if(FE.side_tags()[s]!=INVALID_INDEX)
                {
                  // //std::cout<<"el=="<<el<<" bound=="<<FE.side_tags()[s] <<std::endl;
                  FE.init_boundary(s);
                  bcs_.assembly(full_spaces_ptr(),a,b,c,FE);
                 }
              }
         }      
       }
    }

    void build_boundary_info(const std::vector<Integer>& levels)
    {


     Integer n_levels=levels.size();
     levels_.resize(n_levels);
     constrained_dofs_levels_.resize(n_levels);
     constrained_mat_levels_.resize(n_levels);
     constrained_vec_levels_.resize(n_levels);


     for(Integer i=0;i<n_levels;i++)
     {
      levels_[i]=levels[i];
      build_boundary_info(constrained_dofs_levels_[i],constrained_mat_levels_[i],constrained_vec_levels_[i],levels[i]);
     }

    }


    void build_boundary_info(const Integer level=-1)
    {
      build_boundary_info(constrained_dofs_,constrained_mat_,constrained_vec_,level);
    }


    // void build_boundary_info(const Integer level=-1)
    // {
    //  auto spaces_ptr=bilinear_form_.spaces_ptr()->spaces_ptr();
    //  auto& bisection=spaces_ptr->bisection();
    //  auto& tracker=bisection.tracker();

    //  if(level==-1)
    //   level_=tracker.current_iterate()-1;
    //  else
    //   level_=level;

    //  auto& level_cumultive_n_dofs=full_spaces_ptr()->dofsdofmap().level_cumultive_n_dofs();

    //  n_dofs_=level_cumultive_n_dofs[level_];

    //  auto& mesh=full_spaces_ptr()->mesh();

    //  FiniteElem<Elem> FE(mesh);
 
 
    //  constrained_dofs_.clear();
    //  constrained_mat_.clear();
    //  constrained_vec_.clear();
    //  constrained_dofs_.resize(n_dofs_,false);
    //  constrained_mat_.resize(n_dofs_,0);
    //  constrained_vec_.resize(n_dofs_,0);

    //  shape_coefficients_.init();
    //  //std::cout<<"BUILD BOUNDARY INFO"<<std::endl;
    //    for(std::size_t el=0;el<mesh.n_elements();el++)
    //    {
    //       if(!elem_belongs_to_level(mesh,el,level_,tracker)) continue;

    //       FE.init(el,level_);

    //       if(FE.is_on_boundary())
    //       {

    //         for(std::size_t s=0;s<FE.n_side();s++)
    //           {            
    //             if(FE.side_tags()[s]!=INVALID_INDEX)
    //             {
    //               // //std::cout<<"el=="<<el<<" bound=="<<FE.side_tags()[s] <<std::endl;
    //               FE.init_boundary(s);
    //               bcs_.assembly(full_spaces_ptr(),constrained_dofs_,constrained_mat_,constrained_vec_,FE);
    //              }
    //           }
    //      }      
    //    }
    // }


    auto& constrained_dofs_levels(){return constrained_dofs_levels_;}

    auto& constrained_dofs_levels()const {return constrained_dofs_levels_;}

    auto& constrained_mat_levels(){return constrained_mat_levels_;}

    auto& constrained_mat_levels()const {return constrained_mat_levels_;}

    auto& constrained_vec_levels(){return constrained_vec_levels_;}

    auto& constrained_vec_levels()const{return constrained_vec_levels_;}

    void constrained_quantities(std::vector<bool>& a, std::vector<Real>& b, std::vector<Real>& c)const 
    {
      a.resize(n_dofs_);
      b.resize(n_dofs_);
      c.resize(n_dofs_);

      for(Integer i=0;i<n_dofs_;i++)
      {
        a[i]=constrained_dofs_[i];
        b[i]=constrained_mat_[i];
        c[i]=constrained_vec_[i];
      }

    }




    template<typename SystemMat, typename Rhs>
    void apply_bc(SystemMat& A, Rhs& b)const
    {
      // //std::cout<<"------APPLY BC -------"<<std::endl;
      // //std::cout<<"------b.size="<<b.size()<<std::endl;
      // //std::cout<<"------n_dofs_="<<n_dofs_<<std::endl;
      // A.print_val();
      // //std::cout<<"constrained_dofs_ size="<<constrained_vec_.size()<<std::endl;

     // for(Integer i=0;i<n_dofs_;++i)
     // {
     //  //std::cout<<constrained_dofs_[i]<<std::endl;
     // }
     //  //std::cout<<"------constrained_vec_ BC -------"<<std::endl;
     // for(Integer i=0;i<n_dofs_;++i)
     // {
     //  //std::cout<< constrained_vec_[i]<<", "  <<std::endl;
     // }
     // //std::cout<<"------constrained_mat_ BC -------"<<std::endl;
     // for(Integer i=0;i<n_dofs_;++i)
     // {
     //  //std::cout<<constrained_mat_[i] << ", "   <<std::endl;
     // }
     for(Integer i=0;i<n_dofs_;++i)
     {
      // std::cout<<"i="<<i<<"/"<<n_dofs_<<std::endl;
      if(constrained_dofs_[i])
      {
        // std::cout<<"constrained = "<< constrained_vec_[i] << std::endl;
        
        b[i]=constrained_vec_[i];

        // for(Integer j=0;j<i;++j)
        //     { 
        //       A[i][j]=0;
        //     }
        // A[i][i]=constrained_mat_[i];
        // //std::cout<<"set_zero_row = "<< i << std::endl;
        A.set_zero_row(i);
        // //std::cout<<"equal = "<< i << std::endl;
        A.equal(constrained_mat_[i],i,i);
        // //std::cout<<"after equal = "<< i << std::endl;

        // for(Integer j=i+1;j<n_dofs_;++j)
        //     A[i][j]=0;
      
      }
      else
      {
        // //std::cout<<"not constrained = "<< constrained_vec_[i] << std::endl;
        A.row_static_condensation(i,constrained_vec_,constrained_dofs_,b[i]);
        // for(Integer j=0;j<n_dofs_;++j)
        //   if(constrained_dofs_[j])
        //   {
        //     b[i]-=A[i][j]*constrained_vec_[j];
        //     A[i][j]=0;
        //   }
      }


     }

     // A.print_val();
    }


    template<typename SystemMat, typename Rhs>
    void apply_zero_bc(SystemMat& A, Rhs& b,const Integer level)const
    {
      Real one=1.0;

     //  //std::cout<<"constrained_dofs_levels_[level][i]"<<std::endl;

     //  for(Integer i=0;i<constrained_dofs_levels_[level].size();++i)
     // {
     //  //std::cout<<constrained_dofs_levels_[level][i]<<std::endl;
     // }

      for(Integer i=0;i<constrained_dofs_levels_[level].size();++i)
     {
      if(constrained_dofs_levels_[level][i])
      {       
        A.set_zero_row(i);
        A.equal(one,i,i);
        b[i]=0;
      }
      else
      {
        A.row_static_condensation(i,constrained_dofs_levels_[level]);
      }
     }
    }



    template<typename SystemMat, typename Rhs>
    void apply_bc(SystemMat& A, Rhs& b,const Integer level)const
    {
     auto& level_cumultive_n_dofs=full_spaces_ptr()->dofsdofmap().level_cumultive_n_dofs();

     Integer n_dofs=level_cumultive_n_dofs[level_];

     for(Integer i=0;i<n_dofs;++i)
     {
      if(constrained_dofs_levels_[level][i])
      {        
        b[i]=constrained_vec_levels_[level][i];
        A.set_zero_row(i);
        A.equal(constrained_mat_levels_[level][i],i,i);
      }
      else
      {
        A.row_static_condensation(i,constrained_vec_levels_[level],constrained_dofs_levels_[level],b[i]);
      }
     }
    }



    template<typename SystemMat, typename Rhs>
    void apply_zero_bc_for_null_diagonal_element(SystemMat& A, Rhs& b)const
    {

     Real toll=1e-10;
     Real one=1.0;

     for(Integer i=0;i<b.size();++i)
     {       
        if(abs(A(i,i))<toll)
        {
           A.equal(one,i,i);
           b[i]=0.0;
        }
     }
    }

    template<typename SystemMat, typename Rhs>
    void apply_zero_bc_for_null_diagonal_element(SystemMat& A, Rhs& b, const std::vector<Integer>& w)const
    {

     Real toll=1e-10;
     Real big=1.0/toll;
     Real tmp;

     for(Integer i=0;i<b.size();++i)
     {  

      tmp=0.0;


        if(w[i]==0)
        {
           A.equal(big,i,i);
           b[i]=0.0;
        }
     }
    }

    template<typename SystemMat>
    void apply_zero_bc_for_null_diagonal_element(SystemMat& A)const
    {

     Real toll=1e-10;
     Real one=1.0;

     for(Integer i=0;i<A.rows();++i)
     {       
        if(abs(A(i,i))<toll)
        {
           A.equal(one,i,i);
        }
     }
    }


    template<typename SystemMat>
    void apply_zero_bc_to_matrix(SystemMat& A,const Integer level)
    {

      Real one=1.0;

      auto& level_cumultive_n_dofs=full_spaces_ptr()->dofsdofmap().level_cumultive_n_dofs();

      Integer n_dofs=level_cumultive_n_dofs[level];


      for(Integer i=0;i<n_dofs_;++i)
     {
      if(constrained_dofs_levels_[level][i])
      {       
        A.set_zero_row(i);
        A.equal(one,i,i);
      }
      else
      {
        A.row_static_condensation(i,constrained_dofs_levels_[level]);
      }
     }
    }



    template<typename SystemMat>
    void apply_zero_bc_to_matrix(SystemMat& A, const std::vector<bool>& constrained_dofs)
    {
      Real one=1.0;
      for(Integer i=0;i<n_dofs_;++i)
     {
      if(constrained_dofs[i])
      {       
        A.set_zero_row(i);
        A.equal(one,i,i);
      }
      else
      {
        A.row_static_condensation(i,constrained_dofs);
      }
     }
    }


    void apply_bc_to_vector(std::vector<Real>& b, Integer level=-1)const
    {
      Real one=1.0;

      std::cout<<"1 apply_bc_to_vector"<<std::endl;
      
      // auto& level_cumultive_n_dofs=full_spaces_ptr()->dofsdofmap().level_cumultive_n_dofs();
      if(level==-1)
        level=levels_[levels_.size()-1];
      std::cout<<"2 apply_bc_to_vector"<<std::endl;
      // //std::cout<<"level="<<level<<std::endl;
      std::cout<<"constrained_vec_levels_.size()="<<constrained_vec_levels_.size()<<std::endl;
      std::cout<<std::endl;
      std::cout<<"level="<<level<<std::endl;
      for(Integer i=0;i<levels_.size();i++)
        std::cout<<levels_[i]<<std::endl;
      std::cout<<std::endl;


      std::cout<<std::endl;
      
      auto& constrained_dofs=constrained_dofs_levels_[constrained_dofs_levels_.size()-1];
      auto& constrained_vec=constrained_vec_levels_[constrained_vec_levels_.size()-1];
      Integer ndofs=constrained_vec.size();

      if(b.size()<ndofs)
         {b.resize(ndofs,0);}

        for(Integer i=0;i<ndofs;++i)
       {
        if(constrained_dofs[i])
        {   

         b[i]=constrained_vec[i];
        std::cout<<i<<"/"<<ndofs<<"   "<<b[i]<<std::endl;

        }

       }
      
    }


    void apply_zero_bc_to_vector(std::vector<Real>& b, Integer level=-1)const
    {
      Real one=1.0;
      
      auto& level_cumultive_n_dofs=full_spaces_ptr()->dofsdofmap().level_cumultive_n_dofs();
      if(level==-1)
        level=levels_.size()-1;//[levels_.size()-1];
      // //std::cout<<"apply_zero_bc_to_vector"<<std::endl;
      // //std::cout<<"level="<<level<<std::endl;
      // //std::cout<<"constrained_vec_levels_.size()="<<constrained_vec_levels_.size()<<std::endl;
      Integer ndofs=constrained_vec_levels_[level].size();
      if(b.size()<ndofs)
         b.resize(ndofs,0);
      else
      {
        for(Integer i=0;i<ndofs;++i)
       {
        // //std::cout<<i<<"/"<<ndofs<<std::endl;
        if(constrained_vec_levels_[level][i])
        {       
         b[i]=0;//constrained_vec_levels_[level][i];
        }
       }
      }

    }

 private:
    SignedNormal<Elem> signed_normal_;
    Bilinear bilinear_form_;
    Linear linear_form_;
    Coefficients shape_coefficients_;
    Maps reference_maps_;
    ShapeFunctions shapefunctions_;
    Evaluation<Expression<Bilinear>,ShapeFunctions> eval_bilinear_form_;
    Evaluation<Expression<Linear>,ShapeFunctions> eval_linear_form_;
    BCs bcs_;
    Integer n_dofs_;
    std::vector<bool> constrained_dofs_;
    std::vector<Real> constrained_mat_;
    std::vector<Real> constrained_vec_;
    std::vector<Integer> levels_;
    std::vector<std::vector<bool>> constrained_dofs_levels_;
    std::vector<std::vector<Real>> constrained_mat_levels_;
    std::vector<std::vector<Real>> constrained_vec_levels_;
    Integer level_;
};




template<typename BilinearForm, typename LinearForm, typename...DirichletBcs>
constexpr auto create_context(const BilinearForm& bilinearform, const LinearForm&linearform, const DirichletBcs&...dirichlet_bcs)
{
  return Context<BilinearForm,LinearForm,DirichletBcs...>(bilinearform,linearform,dirichlet_bcs...);
}








template<typename ZeroForm>
class Context<ZeroForm>
{

 public:
    using Form=GeneralForm<ZeroForm>;
    using Coefficients=ShapeFunctionCoefficientsCollection<Form>;
    using Maps=MapFromReferenceCollection<Form>;
    using ShapeFunctions=ShapeFunctionsCollection<Form>;
    using Elem= typename Form::FunctionSpace::Elem;

    Context(const ZeroForm& zeroform)
    :
    zeroform_(general_form(zeroform))
    ,
    shape_coefficients_(shape_function_coefficients(zeroform_))
    ,
    reference_maps_(reference_maps2(zeroform_))
    ,
    shapefunctions_(shape_functions(shape_coefficients_,reference_maps_,zeroform_))
    ,
    eval_zero_form_(Eval(zeroform_,shapefunctions_)) 
    {

      // //std::cout<< "oooooooooooooooo "<<std::endl;
      // zeroform_.
      // typename Evaluation<Expression<Form>,ShapeFunctions>::template L2Products<0> dee(4,5);
      // typename Evaluation<Expression<Form>,ShapeFunctions>::template L2Products<1> de3e(4,5);
      // typename Form::TupleOfPairsNumbers defeeffe(5);
      // typename Form::TupleFunctionSpace eeee1(5,5,6,7,8,9,0);
      // typename Form::FunctionSpace eeee2(5,5,6,7,8,9,0);
      // typename Form::form eeee3(5,5,6,7,8,9,0);
      // typename Form::UniqueElementFunctionSpacesTupleType eeee4(5,5,6,7,8,9,0);
      // typename Form::TupleOfPairsNumbers eeee5(5,5,6,7,8,9,0);
    }



    void assembly(Real& A)
    {


     auto spaces_ptr=zeroform_.spaces_ptr()->spaces_ptr();
     n_dofs_=spaces_ptr->n_dofs();
     auto mesh_ptr=spaces_ptr->mesh_ptr();
     FiniteElem<Elem> FE(mesh_ptr);



     shape_coefficients_.init(*mesh_ptr);
     // //std::cout<<"------_______--6"<<std::endl;
       for(std::size_t el=0;el<mesh_ptr->n_elements();el++)
       {
          if(!mesh_ptr->is_active(el)) continue;
          
          // //std::cout<<"------_______----- ELEMENT ID = "<<el<<". -----_______--------"<<std::endl;
          FE.init(el);
          shape_coefficients_.init(el);
          reference_maps_.init(FE);
          shapefunctions_.init(FE);
          // //std::cout<<"eval_zero_form_ begin "<<std::endl;
          eval_zero_form_.apply(A,FE);
          // //std::cout<<"eval_zero_form_ end "<<std::endl;
          // //std::cout<<"------_______----- SURFACE-----_______--------"<<std::endl;

          if(FE.is_on_boundary())
          {
            for(std::size_t s=0;s<FE.n_side();s++)
              {
                // controlla, qui passiamo side_id, ma dovremmo avere label
                // dovresti costruire mesh coi label

                // //std::cout<<"------_______----- BEGIN SIDE===="<<s<<std::endl;
                FE.init_boundary(s);
                if(FE.is_side_on_boundary())
                {
                  reference_maps_.init_boundary(FE);
                  shapefunctions_.init_boundary(FE);
                  // //std::cout<<"------_______----- BEGIN SIDE EVAL===="<<s<<std::endl;


                  // eval_zero_form_.apply_boundary(A,FE);
                  // //std::cout<<"------_______----- END SIDE EVAL===="<<s<<std::endl;
                }
                // //std::cout<<"------_______----- END SIDE===="<<s<<std::endl;
              }


         }
        
       }
       }

        Real assembly()
        {
         Real A; 
         assembly(A);
         return A;}

 private:
    Form zeroform_;
    Coefficients shape_coefficients_;
    Maps reference_maps_;
    ShapeFunctions shapefunctions_;
    Evaluation<Expression<Form>,ShapeFunctions> eval_zero_form_;
    Integer n_dofs_;
};

template<typename ZeroForm>
constexpr auto create_context(const ZeroForm& zeroform)
{
  return Context<ZeroForm>(zeroform);
}


}
#endif