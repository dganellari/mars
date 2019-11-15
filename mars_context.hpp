#ifndef MARS_CONTEXT_HPP
#define MARS_CONTEXT_HPP
#include "mars_base.hpp"
#include "mars_shape_functions_collection.hpp"
#include "mars_dirichlet_bc.hpp"
#include "mars_general_form.hpp"
#include "mars_evaluation.hpp"

namespace mars {




template<typename BilinearForm, typename LinearForm, typename...DirichletBCs>
class Context
{

 public:
    using Bilinear=GeneralForm<BilinearForm>;
    using Linear=GeneralForm<LinearForm>;
    using Coefficients=ShapeFunctionCoefficientsCollection<Bilinear,Linear>;
    using Maps=MapFromReferenceCollection<Bilinear,Linear>;
    using ShapeFunctions=ShapeFunctionsCollection<Bilinear,Linear>;
    using BCs=DirichletBoundaryConditionCollection<DirichletBCs...>;
    using Elem= typename Bilinear::FunctionSpace::Elem;


    Context(const BilinearForm& bilinearform, const LinearForm&linearform, const DirichletBCs&...bcs):
    bilinear_form_(general_form(bilinearform)),
    linear_form_(general_form(linearform)),
    shape_coefficients_(shape_function_coefficients(bilinear_form_,linear_form_)),
    reference_maps_(reference_maps2(bilinear_form_,linear_form_)),
    shapefunctions_(shape_functions(shape_coefficients_,reference_maps_,bilinear_form_,linear_form_)),
    eval_bilinear_form_(Eval(bilinear_form_,shapefunctions_)),
    eval_linear_form_(Eval(linear_form_,shapefunctions_)),
    bcs_(bcs...)
    {}


    template<typename SystemMat, typename Rhs>
    void assembly(SystemMat& A, Rhs& b)
    {
     auto spaces_ptr=bilinear_form_.spaces_ptr()->spaces_ptr();
     n_dofs_=spaces_ptr->n_dofs();
     auto mesh_ptr=spaces_ptr->mesh_ptr();
     FiniteElem<Elem> FE(mesh_ptr);
     A.resize(n_dofs_,std::vector<Real>(n_dofs_));
     b.resize(n_dofs_);   
     constrained_dofs.resize(n_dofs_,false);
     constrained_mat.resize(n_dofs_,0);
     constrained_vec.resize(n_dofs_,0);


     shape_coefficients_.init(*mesh_ptr);
       for(std::size_t el=0;el<mesh_ptr->n_elements();el++)
       {

          FE.init(el);
          shape_coefficients_.init(el);
          reference_maps_.init(FE);
          shapefunctions_.init(FE);
          eval_bilinear_form_.apply(A,FE);
          eval_linear_form_.apply(b,FE);
          std::cout<<"------_______----- SURFACE-----_______--------"<<std::endl;

          if(FE.is_on_boundary())
          {
            for(std::size_t s=0;s<FE.n_side();s++)
              {
                // controlla, qui passiamo side_id, ma dovremmo avere label
                // dovresti costruire mesh coi label

                std::cout<<"------_______----- BEGIN SIDE===="<<s<<std::endl;
                FE.init_boundary(s);
                if(FE.is_side_on_boundary())
                {
                  reference_maps_.init_boundary(FE);
                  shapefunctions_.init_boundary(FE);
                  std::cout<<"------_______----- BEGIN SIDE EVAL===="<<s<<std::endl;
                  eval_bilinear_form_.apply_boundary(A,FE);
                  eval_linear_form_.apply_boundary(b,FE);
                  std::cout<<"------_______----- END SIDE EVAL===="<<s<<std::endl;
                }
                std::cout<<"------_______----- END SIDE===="<<s<<std::endl;
              }


            for(std::size_t s=0;s<FE.n_side();s++)
              {
                // controlla, qui passiamo side_id, ma dovremmo avere label
                // dovresti costruire mesh coi label

                std::cout<<"------_______----- ELEM===="<<el<<std::endl;
                std::cout<<"------_______----- BEGIN SIDE===="<<s<<std::endl;

                
                if(FE.side_tags()[s]!=INVALID_INDEX)
                {
                  std::cout<<"------_______----- INIT BOUNDARY===="<<s<<std::endl;
                  FE.init_boundary(s);
                  std::cout<<"------_______----- END INIT BOUNDARY===="<<s<<std::endl;
                  bcs_.assembly(constrained_dofs,constrained_mat,constrained_vec,FE);
                  // reference_maps_.init_boundary(FE);
                  // shapefunctions_.init_boundary(FE);
                  // std::cout<<"------_______----- BEGIN SIDE EVAL===="<<s<<std::endl;
                  // eval_bilinear_form_.apply_boundary(A,FE);
                  // eval_linear_form_.apply_boundary(b,FE);
                  // std::cout<<"------_______----- END SIDE EVAL===="<<s<<std::endl;
                }
                std::cout<<"------_______----- END SIDE===="<<s<<std::endl;
              }
           // constrained_dofs
        }
       }


       std::cout<<"------CONSTRAINED DOFS-------"<<std::endl;
       for(std::size_t i=0;i<n_dofs_;i++)
          std::cout<<i<<" "<<constrained_dofs[i]<<std::endl;

       std::cout<<"------CONSTRAINED MAT-------"<<std::endl;
       for(std::size_t i=0;i<n_dofs_;i++)
          std::cout<<i<<" "<<constrained_mat[i]<<std::endl;

       std::cout<<"------CONSTRAINED VEC-------"<<std::endl;
       for(std::size_t i=0;i<n_dofs_;i++)
          std::cout<<i<<" "<<constrained_vec[i]<<std::endl;

    }
 


    template<typename SystemMat, typename Rhs>
    void apply_bc(SystemMat& A, Rhs& b)
    {
      std::cout<<"------APPLY BC -------"<<std::endl;
     for(Integer i=0;i<n_dofs_;++i)
     {
      std::cout<<"i="<<i<<std::endl;
      if(constrained_dofs[i])
      {
        std::cout<<"constrained"<<std::endl;
        
        b[i]=constrained_vec[i];

        for(Integer j=0;j<i;++j)
            { 
              std::cout<<"??????????"<<A[i][j]<<std::endl;
              A[i][j]=0;
            }
        A[i][i]=constrained_mat[i];

        for(Integer j=i+1;j<n_dofs_;++j)
            A[i][j]=0;
      
      }
      else
      {
        std::cout<<"not constrained"<<std::endl;
        for(Integer j=0;j<n_dofs_;++j)
          if(constrained_dofs[j])
          {

            A[i][j]=0;
            b[i]-=A[i][j]*constrained_vec[j];
          }
      }


     }
    }


 private:
    Bilinear bilinear_form_;
    Linear linear_form_;
    Coefficients shape_coefficients_;
    Maps reference_maps_;
    ShapeFunctions shapefunctions_;
    Evaluation<Expression<Bilinear>,ShapeFunctions> eval_bilinear_form_;
    Evaluation<Expression<Linear>,ShapeFunctions> eval_linear_form_;
    BCs bcs_;
    Integer n_dofs_;
    std::vector<bool> constrained_dofs;
    std::vector<Real> constrained_mat;
    std::vector<Real> constrained_vec;
};




template<typename BilinearForm, typename LinearForm, typename...DirichletBcs>
constexpr auto create_context(const BilinearForm& bilinearform, const LinearForm&linearform, const DirichletBcs&...dirichlet_bcs)
{
  return Context<BilinearForm,LinearForm,DirichletBcs...>(bilinearform,linearform,dirichlet_bcs...);
}



}
#endif