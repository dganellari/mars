#ifndef MARS_QUADRATURE_RULES_HPP
#define MARS_QUADRATURE_RULES_HPP


namespace mars {

	
template<typename Elem_,Integer Order_>
class BaseQuadrature
{
 using Elem=Elem_;
 static constexpr Integer Order=Order;
};

template<typename Elem,Integer Order>
class GaussPoints;

template< Integer Dim >
class GaussPoints< Simplex<Dim,2> , 1>:BaseQuadrature<Simplex<Dim,2>,1>
{
private:
 Matrix<Real,1,2> qp_points_;
 Vector<Real,1> qp_weights_;
public:
  using qp_points_type=Matrix<Real,1,2>;
  static constexpr Integer NQPoints=1;
  const Matrix<Real,1,2> qp_points()const {return qp_points_;};
  const Vector<Real,1> qp_weights()const {return qp_weights_;};

  GaussPoints<Simplex<Dim,2>,1>():
  qp_points_({0.33333333333333, 0.33333333333333}),
  qp_weights_({1})
   {}; 

};



template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,2>:BaseQuadrature<Simplex<Dim,2>,2>
{
private:
  Matrix<Real,3,2>  qp_points_;
  Vector<Real,3> qp_weights_;
public:
  using qp_points_type=Matrix<Real,3,2>;
  static constexpr Integer NQPoints=3;
  const Matrix<Real,3,2> qp_points()const {return qp_points_;};
  const Vector<Real,3> qp_weights()const {return qp_weights_;};
  GaussPoints<Simplex<Dim,2>,2>():
  qp_points_({0.16666666666667, 0.16666666666667,
              0.16666666666667, 0.66666666666667,
              0.16666666666667, 0.16666666666667}),
  qp_weights_({0.33333333333333,0.33333333333333,0.33333333333333})
  {}
};




template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,3>:BaseQuadrature<Simplex<Dim,2>,3>
{
private:
  Matrix<Real,4,2>  qp_points_;
  Vector<Real,4> qp_weights_;
public:
  using qp_points_type=Matrix<Real,4,2>;
  static constexpr Integer NQPoints=4;
  const Matrix<Real,4,2> qp_points()const {return qp_points_;};
  const Vector<Real,4> qp_weights()const {return qp_weights_;};
  GaussPoints<Simplex<Dim,2>,3>():
  qp_points_({0.33333333333333, 0.33333333333333,
              0.20000000000000, 0.20000000000000,
              0.20000000000000, 0.60000000000000,
              0.60000000000000, 0.20000000000000}),
  qp_weights_({ -0.56250000000000,0.52083333333333,0.52083333333333,0.52083333333333})
  {}
};

template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,4>:BaseQuadrature<Simplex<Dim,2>,4>
{
private:
  Matrix<Real,6,2>  qp_points_;
  Vector<Real,6> qp_weights_;
public:
  using qp_points_type=Matrix<Real,6,2>;
  static constexpr Integer NQPoints=6;
  const Matrix<Real,6,2> qp_points()const {return qp_points_;};
  const Vector<Real,6> qp_weights()const {return qp_weights_;};
  GaussPoints<Simplex<Dim,2>,4>():
  qp_points_({0.44594849091597, 0.44594849091597,
              0.44594849091597, 0.10810301816807,
              0.10810301816807, 0.44594849091597, 
              0.09157621350977, 0.09157621350977, 
              0.09157621350977, 0.81684757298046,
              0.81684757298046, 0.09157621350977} ),
  qp_weights_({0.22338158967801, 
            0.22338158967801,
            0.22338158967801,
            0.10995174365532,
            0.10995174365532, 
            0.10995174365532})
  {}
};


template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,5>:BaseQuadrature<Simplex<Dim,2>,5>
{
private:
  Matrix<Real,7,2>  qp_points_;
  Vector<Real,7> qp_weights_;
public:
  using qp_points_type=decltype(qp_points_);
  static constexpr Integer NQPoints=7;
  const Matrix<Real,7,2> qp_points()const {return qp_points_;};
  const Vector<Real,7> qp_weights()const {return qp_weights_;};
  GaussPoints<Simplex<Dim,2>,5>():
  qp_points_(   {0.33333333333333, 0.33333333333333, 
                 0.47014206410511, 0.47014206410511,
                 0.47014206410511, 0.05971587178977,
                 0.05971587178977, 0.47014206410511, 
                 0.10128650732346, 0.10128650732346, 
                 0.10128650732346, 0.79742698535309, 
                 0.79742698535309, 0.10128650732346}  ),
  qp_weights_({0.22500000000000, 
            0.13239415278851, 
            0.13239415278851, 
            0.13239415278851, 
            0.12593918054483, 
            0.12593918054483, 
            0.12593918054483 })
  {}
};
}

#endif //MARS_QUADRATURE_RULES_HPP
