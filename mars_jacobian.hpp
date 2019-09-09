#ifndef MARS_JACOBIAN_HPP
#define MARS_JACOBIAN_HPP
#include "mars_simplex.hpp"



namespace mars{




template<typename Elem>
class Jacobian;

template<Integer Dim, Integer ManifoldDim>
class Jacobian<Simplex<Dim, ManifoldDim>>
{
  public:

  Jacobian(const Mesh<Dim,ManifoldDim>&mesh):
  mesh_ptr_(std::make_shared<Mesh<Dim,ManifoldDim>>(mesh))
  {}
  
  constexpr void init_det(){detJ_=det(J_);}

  constexpr void init(const Integer id)
  {jacobian(mesh_ptr_->elem(id),mesh_ptr_->points(),J_);
   init_det();}

  constexpr auto & operator()()const {return J_;}
  constexpr auto & get_det()const {return detJ_;}

  private:  
  Simplex<Dim, ManifoldDim> simplex_;
  Matrix<Real, Dim, ManifoldDim> J_;
  Real detJ_;
  std::shared_ptr<Mesh<Dim,ManifoldDim>> mesh_ptr_;


};




    template<Integer Dim>
    inline constexpr Real dotofdots(const Vector<Real, Dim> &left, const Vector<Real,Dim> &right)
    {
        return dot(left,right);
    }

    template<Integer Rows,Integer Cols>
    inline constexpr Real dotofdots(const Matrix<Real, Rows,Cols> &left, const Matrix<Real,Rows,Cols> &right)
    {
        Real result = left(0,0)*right(0,0);
        
        for(Integer col = 1; col < Cols; ++col) 
            result += left(0,col)*right(0,col);
        for(Integer col = 0; col < Cols; ++col)
           for(Integer row = 1; row < Rows; ++row)  
            result += left(row,col)*right(row,col);  

        return result;
    }

    template<Integer Dim, typename T>
    inline constexpr Real dotofdots(const Vector<T, Dim> &left, const Vector<T, Dim> &right)
    {
        Real ret = dotofdots(left(0),right(0));
        for(Integer d = 1; d < Dim; ++d) {
            ret += dotofdots(left(d),right(d));
        }

        return ret;
    }

    template<typename T, Integer Rows,Integer Cols>
    inline constexpr Real dotofdots(const Matrix<T, Rows,Cols> &left, const Matrix<T,Rows,Cols> &right)
    {
        Real result = dotofdots(left(0,0),right(0,0));
        
        for(Integer col = 1; col < Cols; ++col) 
            result += dotofdots(left(0,col),right(0,col));
        for(Integer col = 0; col < Cols; ++col)
           for(Integer row = 1; row < Rows; ++row)  
            result += dotofdots(left(row,col),right(row,col));  

        return result;
    }



template<bool PositiveWeights,typename TestTrialSpaces,typename L2,typename Form>
class LocalMatrix;



template<typename TestTrialSpaces, typename Left,typename Right,Integer QR,typename Form>
class LocalMatrix<false,TestTrialSpaces,L2DotProductIntegral<Left,Right,QR>,Form>
// template<typename TestTrialSpaces, typename MeshT, typename Left,typename Right,Integer QR,typename Form>
// class LocalMatrix<false,TestTrialSpaces,L2DotProductIntegral<MeshT,Left,Right,QR>,Form>
{
    static_assert("negative weights not permitted");
};




// template<typename TestTrialSpaces, typename MeshT, typename Left,typename Right,Integer QR,typename Form>
// class LocalMatrix<true,TestTrialSpaces,L2DotProductIntegral<MeshT,Left,Right,QR>,Form>
template<typename TestTrialSpaces, typename Left,typename Right,Integer QR,typename Form>
class LocalMatrix<true,TestTrialSpaces,L2DotProductIntegral<Left,Right,QR>,Form>
{
 public:
 // using type= L2DotProductIntegral<MeshT,Left,Right,QR>;
 using type= L2DotProductIntegral<Left,Right,QR>;
 using QRule=typename type::QRule;
 using EvalLeft=Evaluation<Expression<Left>,QRule>;
 using EvalRight=Evaluation<Expression<Right>,QRule>;
 using subtype= OperatorType<type,ShapeFunctions2<Form>>;
 

 template<typename Elem,typename ShapeFunctions>
 void apply(subtype& mat, const Jacobian<Elem>& J, const ShapeFunctions& shape_functions)
 {
  const auto& detJ=J.get_det();
  eval_left_.apply(left_value_,shape_functions);
  eval_right_.apply(right_value_,shape_functions);

  // std::cout<<mat<<std::endl;
  // std::cout<<"left_value_"<<std::endl;
  // // std::cout<<left_value_<<std::endl;
  // std::cout<<"right_value_"<<std::endl;
  // std::cout<<right_value_<<std::endl;
  for(Integer ii=0;ii<mat.rows();ii++)
    for(Integer jj=0;jj<mat.cols();jj++)
    {   

        // std::cout<<"(ii,jj)=("<<ii<<","<<jj<<")"<<std::endl;
   


        // left value is trial, right value is test
        mat(ii,jj)=detJ*dotofdots(right_value_[ii],left_value_[jj]);
    }
    // std::cout<<detJ<<std::endl;
    // std::cout<<mat<<std::endl;


  mat=mat*detJ;
 }

  
 private:
 EvalLeft eval_left_;
 EvalRight eval_right_;
 OperatorType<Left,QRule> left_value_;
 OperatorType<Right,QRule> right_value_;

};





template<typename ...Ts>
class ConstructL2Inner;


// template<typename MeshT,typename Left1,typename Right1,Integer QR1>
// class ConstructL2Inner<L2DotProductIntegral<MeshT,Left1,Right1,QR1>>
                                   
// {
// public:
//    using T=L2DotProductIntegral<MeshT,Left1,Right1,QR1>;

//    static auto apply(const MeshT& mesh)
//    {
//     return L2Inner(mesh, typename T::Left(), typename T::Right());
//     }

// };
template<typename Left1,typename Right1,Integer QR1>
class ConstructL2Inner<L2DotProductIntegral<Left1,Right1,QR1>>
                                   
{
public:
   using T=L2DotProductIntegral<Left1,Right1,QR1>;

   static auto apply()
   {
    return L2Inner(typename T::Left(), typename T::Right());
    }

};

// template<typename MeshT,typename Left1,typename Right1,Integer QR1,
//                         typename Left2,typename Right2,Integer QR2>
// class ConstructL2Inner<Addition<Expression<L2DotProductIntegral<MeshT,Left1,Right1,QR1>>,
//                                 Expression<L2DotProductIntegral<MeshT,Left2,Right2,QR2>> > >
                                   
// {
// public:
//    using Left=L2DotProductIntegral<MeshT,Left1,Right1,QR1>;
//    using Right=L2DotProductIntegral<MeshT,Left2,Right2,QR2>;

//    static auto apply(const MeshT& mesh)
//    {
//     auto e1=L2Inner(mesh, typename Left::Left(), typename Left::Right());
//     auto e2=L2Inner(mesh, typename Right::Left(), typename Right::Right());
//    return e1+e2;
//     }

// };
template<typename Left1,typename Right1,Integer QR1,
         typename Left2,typename Right2,Integer QR2>
class ConstructL2Inner<Addition<Expression<L2DotProductIntegral<Left1,Right1,QR1>>,
                                Expression<L2DotProductIntegral<Left2,Right2,QR2>> > >
                                   
{
public:
   using Left=L2DotProductIntegral<Left1,Right1,QR1>;
   using Right=L2DotProductIntegral<Left2,Right2,QR2>;

   static auto apply()
   {
    auto e1=L2Inner(typename Left::Left(), typename Left::Right());
    auto e2=L2Inner(typename Right::Left(), typename Right::Right());
   return e1+e2;
    }

};
                         


// template<typename Left,typename Right>
// class ConstructL2Inner< Addition<Expression<Left>,Expression<Right>> >                             
// {
// public:
//    template<typename MeshT>
//    static auto apply(const MeshT& mesh)
//    {
//     auto e1=ConstructL2Inner<Left>::apply(mesh);
//     auto e2=ConstructL2Inner<Right>::apply(mesh);
//    return e1+e2;
//     }

// };
template<typename Left,typename Right>
class ConstructL2Inner< Addition<Expression<Left>,Expression<Right>> >                             
{
public:
   static auto apply()
   {
    auto e1=ConstructL2Inner<Left>::apply();
    auto e2=ConstructL2Inner<Right>::apply();
   return e1+e2;
    }

};


template<typename Tuple,Integer Nmax,Integer N, typename ShapeFunctions>
class EvalOfL2InnersAuxHelper;

template<typename Tuple,Integer Nmax, typename ShapeFunctions>
class EvalOfL2InnersAuxHelper<Tuple,Nmax,Nmax,ShapeFunctions>
{
 public:
    using single_type=Evaluation<Expression<GetType<Tuple,Nmax>>,ShapeFunctions>;
    using type=std::tuple<single_type>;
};

template<typename Tuple,Integer Nmax,Integer N, typename ShapeFunctions>
class EvalOfL2InnersAuxHelper
{
 public:
    using single_type=Evaluation<Expression<GetType<Tuple,N>>,ShapeFunctions>;
    using type=TupleCatType<std::tuple<single_type>,typename EvalOfL2InnersAuxHelper<Tuple,Nmax,N+1,ShapeFunctions>::type >;
};

template<typename Tuple,Integer N, typename ShapeFunctions>
using EvalOfL2InnersAux=typename EvalOfL2InnersAuxHelper<Tuple,TupleTypeSize<Tuple>::value-1,N,ShapeFunctions>::type;

template<typename Tuple, typename ShapeFunctions>
using EvalOfL2InnersType=EvalOfL2InnersAux<Tuple,0,ShapeFunctions>;

// template<typename Tuple,Integer N,typename MeshT, typename ShapeFunctions>
// constexpr std::enable_if_t<(N==TupleTypeSize<Tuple>::value-1), EvalOfL2InnersAux<Tuple,N,ShapeFunctions> > 
// EvalOfL2InnersHelper(const MeshT& mesh, ShapeFunctions& shape_functions)
// {   
//       using ens0=GetType<Tuple,N>;
//       using ens=Evaluation<Expression<ens0>,ShapeFunctions>;
//       return std::tuple<ens>(Eval(ConstructL2Inner<ens0>::apply(mesh),shape_functions));
// }
template<typename Tuple,Integer N,typename ShapeFunctions>
constexpr std::enable_if_t<(N==TupleTypeSize<Tuple>::value-1), EvalOfL2InnersAux<Tuple,N,ShapeFunctions> > 
EvalOfL2InnersHelper(ShapeFunctions& shape_functions)
{   
      using ens0=GetType<Tuple,N>;
      using ens=Evaluation<Expression<ens0>,ShapeFunctions>;
      return std::tuple<ens>(Eval(ConstructL2Inner<ens0>::apply(),shape_functions));
}
// template<typename Tuple,Integer N,typename MeshT, typename ShapeFunctions>
// constexpr std::enable_if_t< (N<TupleTypeSize<Tuple>::value-1), EvalOfL2InnersAux<Tuple,N,ShapeFunctions>>  
// EvalOfL2InnersHelper(const MeshT& mesh, ShapeFunctions& shape_functions)
// {
//       using ens0=GetType<Tuple,N>;
//       using ens=Evaluation<Expression<ens0>,ShapeFunctions>;
//       return std::tuple_cat(std::tuple<ens>(Eval(ConstructL2Inner<ens0>::apply(mesh),shape_functions)),
//                             EvalOfL2InnersHelper<Tuple,N+1>(mesh,shape_functions));
// }

template<typename Tuple,Integer N,typename ShapeFunctions>
constexpr std::enable_if_t< (N<TupleTypeSize<Tuple>::value-1), EvalOfL2InnersAux<Tuple,N,ShapeFunctions>>  
EvalOfL2InnersHelper(ShapeFunctions& shape_functions)
{
      using ens0=GetType<Tuple,N>;
      using ens=Evaluation<Expression<ens0>,ShapeFunctions>;
      return std::tuple_cat(std::tuple<ens>(Eval(ConstructL2Inner<ens0>::apply(),shape_functions)),
                            EvalOfL2InnersHelper<Tuple,N+1>(shape_functions));
}

// template<typename Tuple,typename MeshT, typename ShapeFunctions>
// constexpr auto EvalOfL2Inners(const MeshT& mesh, ShapeFunctions& shape_functions)
// {
//       return EvalOfL2InnersHelper<Tuple,0>(mesh,shape_functions);
// }
template<typename Tuple,typename ShapeFunctions>
constexpr auto EvalOfL2Inners(ShapeFunctions& shape_functions)
{
      return EvalOfL2InnersHelper<Tuple,0>(shape_functions);
}







template<typename GeneralForm >
class LocalMatrices
{
public:
    // using GeneralForm=typename EvaluationGeneralForm::type;
    using TupleOfPairsNumbers=typename GeneralForm::TupleOfPairsNumbers;
    static constexpr auto Nelem_dofs_array=GeneralForm::FunctionSpace::Nelem_dofs_array;

    template<Integer Nmax, Integer N>
    class AuxHelper;

    template<Integer Nmax>
    class AuxHelper<Nmax,Nmax>
    {
        public:
        using Numbers=GetType<TupleOfPairsNumbers,Nmax>;
        static constexpr Integer Rows=Nelem_dofs_array[GetType<Numbers,0>::value];
        static constexpr Integer Cols=Nelem_dofs_array[GetType<Numbers,1>::value];
        using type= std::tuple< Matrix<Real,Rows,Cols>>;
    };

    template<Integer Nmax, Integer N>
    class AuxHelper
    {
        public:
        using Numbers=GetType<TupleOfPairsNumbers,N>;
        static constexpr Integer Rows=Nelem_dofs_array[GetType<Numbers,0>::value];
        static constexpr Integer Cols=Nelem_dofs_array[GetType<Numbers,1>::value];
        using type=TupleCatType< std::tuple< Matrix<Real,Rows,Cols>>, typename AuxHelper<Nmax,N+1>::type >;
    };

    using type=typename AuxHelper<TupleTypeSize<TupleOfPairsNumbers>::value-1,0>::type;
    
};



template<typename EvaluationGeneralForm,typename ShapeFunctions>
class EvaluationOfL2Inners
{
public: 
    using GeneralForm=typename EvaluationGeneralForm::type;
    using LocalMatrices=typename LocalMatrices<GeneralForm>::type;    
    using L2Products=typename EvaluationGeneralForm::L2Products;
    using EvalOfL2InnersType=EvalOfL2InnersAux<L2Products,0,ShapeFunctions>;
    
    // template<typename MeshT>
    // EvaluationOfL2Inners(const MeshT&mesh, ShapeFunctions& shapefunctions):
    // eval_inners_(EvalOfL2Inners<L2Products>(mesh,shapefunctions))
    // {}

    EvaluationOfL2Inners(ShapeFunctions& shapefunctions):
    eval_inners_(EvalOfL2Inners<L2Products>(shapefunctions))
    {}

 template<Integer N,typename Output,typename Jacobian>
    constexpr void apply_aux_aux(Output& mat,Jacobian&J)
    {
     auto& local_mat=std::get<N>(mat_tuple_);
     local_mat.zero();
     auto & eval_N=std::get<N>(eval_inners_);
     eval_N.apply(local_mat,J);
     // std::cout<<local_mat<<std::endl;
    }   

 template<Integer Nmax,Integer N,typename Output,typename Jacobian>
    constexpr std::enable_if_t<(N==Nmax),void> apply_aux(Output& mat,Jacobian&J)
    {
     apply_aux_aux<N>(mat,J);
    }

 template<Integer Nmax,Integer N,typename Output,typename Jacobian>
    constexpr std::enable_if_t<(N<Nmax),void> apply_aux(Output& mat,Jacobian&J)
    {
     {apply_aux_aux<N>(mat,J);}
      apply_aux<Nmax,N+1>(mat,J);
     }

 template<typename Output,typename Jacobian>
    void apply(Output& mat,Jacobian&J)
    {apply_aux<TupleTypeSize<L2Products>::value-1,0>(mat,J);}

private:
    LocalMatrices mat_tuple_;
    EvalOfL2InnersType eval_inners_;    
};



}
#endif