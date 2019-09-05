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




template<typename TestTrialSpaces, typename MeshT, typename Left,typename Right,Integer QR,typename Form>
class LocalMatrix<false,TestTrialSpaces,L2DotProductIntegral<MeshT,Left,Right,QR>,Form>
{
    static_assert("negative weights not permitted");
};




template<typename TestTrialSpaces, typename MeshT, typename Left,typename Right,Integer QR,typename Form>
class LocalMatrix<true,TestTrialSpaces,L2DotProductIntegral<MeshT,Left,Right,QR>,Form>
{
 public:
 using type= L2DotProductIntegral<MeshT,Left,Right,QR>;
 using QRule=typename type::QRule;
 using EvalLeft=Evaluation<Expression<Left>,QRule>;
 using EvalRight=Evaluation<Expression<Right>,QRule>;
 using subtype= OperatorType<type,ShapeFunctions2<Form>>;
 

 template<typename Elem,typename ShapeFunctions>
 void apply(subtype& mat, const Jacobian<Elem>& J, const ShapeFunctions& shape_functions)
 {
  const auto& detJ=J.get_det();


  std::cout<<"L2DotProductIntegral"<<std::endl;
  eval_left_.apply(left_value_,shape_functions);
  eval_left_.apply(right_value_,shape_functions);

  // eval_right_.apply<1>(right_value_,shape_functions);

  for(Integer ii=0;ii<mat.rows();ii++)
    for(Integer jj=0;jj<mat.cols();jj++)
        mat(ii,jj)=detJ*dotofdots(left_value_[ii],right_value_[jj]);
    std::cout<<detJ<<std::endl;
    std::cout<<mat<<std::endl;


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


template<typename MeshT,typename Left1,typename Right1,Integer QR1>
class ConstructL2Inner<L2DotProductIntegral<MeshT,Left1,Right1,QR1>>
                                   
{
public:
   using T=L2DotProductIntegral<MeshT,Left1,Right1,QR1>;

   static auto apply(const MeshT& mesh)
   {
    return L2Inner(mesh, typename T::Left(), typename T::Right());
    }

};


template<typename MeshT,typename Left1,typename Right1,Integer QR1,
                        typename Left2,typename Right2,Integer QR2>
class ConstructL2Inner<Addition<Expression<L2DotProductIntegral<MeshT,Left1,Right1,QR1>>,
                                Expression<L2DotProductIntegral<MeshT,Left2,Right2,QR2>> > >
                                   
{
public:
   using Left=L2DotProductIntegral<MeshT,Left1,Right1,QR1>;
   using Right=L2DotProductIntegral<MeshT,Left2,Right2,QR2>;

   static auto apply(const MeshT& mesh)
   {
    auto e1=L2Inner(mesh, typename Left::Left(), typename Left::Right());
    auto e2=L2Inner(mesh, typename Right::Left(), typename Right::Right());
   return e1+e2;
    }

};



// template<typename MeshT,typename Left,
//                         typename Left2,typename Right2,Integer QR2>
// class ConstructL2Inner<Addition<Expression<Left>,
//                                 Expression<L2DotProductIntegral<MeshT,Left2,Right2,QR2>> > >
                                   
// {
// public:
//    using Right=L2DotProductIntegral<MeshT,Left2,Right2,QR2>;

//    static auto apply(const MeshT& mesh)
//    {
//     auto e1=ConstructL2Inner<Left>::apply(mesh);
//     auto e2=L2Inner(mesh, typename Right::Left(), typename Right::Right());
//    return e1+e2;
//     }

// };



// template<typename MeshT,typename Left1,typename Right1,Integer QR1,
//                         typename Right>
// class ConstructL2Inner<Addition<Expression<L2DotProductIntegral<MeshT,Left1,Right1,QR1>>,
//                                 Expression<Right> > >
                                   
// {
// public:
//    using Left=L2DotProductIntegral<MeshT,Left1,Right1,QR1>;

//    static auto apply(const MeshT& mesh)
//    {
//     auto e1=L2Inner(mesh, typename Left::Left(), typename Left::Right());
//     auto e2=ConstructL2Inner<Right>::apply(mesh);
//    return e1+e2;
//     }

// };

template<typename Left,typename Right>
class ConstructL2Inner< Addition<Expression<Left>,Expression<Right>> >                             
{
public:
   template<typename MeshT>
   static auto apply(const MeshT& mesh)
   {
    auto e1=ConstructL2Inner<Left>::apply(mesh);
    auto e2=ConstructL2Inner<Right>::apply(mesh);
   return e1+e2;
    }

};



template<typename Tuple,Integer Nmax,Integer N, typename ShapeFunctions>
class OKKAuxHelper;

template<typename Tuple,Integer Nmax, typename ShapeFunctions>
class OKKAuxHelper<Tuple,Nmax,Nmax,ShapeFunctions>
{
 public:
    using single_type=Evaluation<Expression<GetType<Tuple,Nmax>>,ShapeFunctions>;
    using type=std::tuple<single_type>;
};

template<typename Tuple,Integer Nmax,Integer N, typename ShapeFunctions>
class OKKAuxHelper
{
 public:
    using single_type=Evaluation<Expression<GetType<Tuple,N>>,ShapeFunctions>;
    using type=TupleCatType<std::tuple<single_type>,typename OKKAuxHelper<Tuple,Nmax,N+1,ShapeFunctions>::type >;
};

template<typename Tuple,Integer N, typename ShapeFunctions>
using OKKAux=typename OKKAuxHelper<Tuple,TupleTypeSize<Tuple>::value-1,N,ShapeFunctions>::type;


template<typename Tuple,Integer N,typename MeshT, typename ShapeFunctions>
constexpr std::enable_if_t<(N==TupleTypeSize<Tuple>::value-1), OKKAux<Tuple,N,ShapeFunctions> > 
OKK(const MeshT& mesh, ShapeFunctions& shape_functions)
{   
      using ens0=GetType<Tuple,N>;
      using ens=Evaluation<Expression<ens0>,ShapeFunctions>;
      return std::tuple<ens>(Eval(ConstructL2Inner<ens0>::apply(mesh),shape_functions));
}

template<typename Tuple,Integer N,typename MeshT, typename ShapeFunctions>
constexpr std::enable_if_t< (N<TupleTypeSize<Tuple>::value-1), OKKAux<Tuple,N,ShapeFunctions>>  
OKK(const MeshT& mesh, ShapeFunctions& shape_functions)
{
      using ens0=GetType<Tuple,N>;
      using ens=Evaluation<Expression<ens0>,ShapeFunctions>;
      return std::tuple_cat(std::tuple<ens>(Eval(ConstructL2Inner<ens0>::apply(mesh),shape_functions)),
                            OKK<Tuple,N+1>(mesh,shape_functions));
}


// template<typename Left,typename Right>
// class ConstructL2Inner< Subtraction<Expression<Left>,Expression<Right>> >                             
// {
// public:
//    template<typename MeshT>
//    static auto apply(const MeshT& mesh)
//    {
//     auto e1=ConstructL2Inner<Left>::apply(mesh);
//     auto e2=ConstructL2Inner<Right>::apply(mesh);
//    return e1-e2;
//     }

// };

}
#endif