#ifndef MARS_QUADRATURE_RULES_HPP
#define MARS_QUADRATURE_RULES_HPP





namespace mars {









static constexpr Integer GaussianQuadrature = 1;

	
template<typename Elem_,Integer Order_>
class BaseQuadrature
{
public:
 using Elem=Elem_;
 static constexpr Integer Order=Order;
};


template<typename Elem,Integer QuadratureRuleType>
class MaxOrder;

template<Integer Dim>
class MaxOrder<Simplex<Dim,2>,GaussianQuadrature>
{
public:
  static constexpr Integer value=5;
};

template<typename Elem,Integer Order>
class GaussPoints;

static constexpr  Matrix<Real,1,2> QPSimplex1_qp_points(0.33333333333333, 0.33333333333333);
static constexpr  Matrix<Real,3,2> QPSimplex2_qp_points(0.16666666666667, 0.16666666666667,
                                                        0.16666666666667, 0.66666666666667,
                                                        0.16666666666667, 0.16666666666667);
static constexpr  Matrix<Real,4,2> QPSimplex3_qp_points(0.33333333333333, 0.33333333333333,
                                                        0.20000000000000, 0.20000000000000,
                                                        0.20000000000000, 0.60000000000000,
                                                        0.60000000000000, 0.20000000000000);
static constexpr  Matrix<Real,6,2> QPSimplex4_qp_points(0.44594849091597, 0.44594849091597,
                                                        0.44594849091597, 0.10810301816807,
                                                        0.10810301816807, 0.44594849091597, 
                                                        0.09157621350977, 0.09157621350977, 
                                                        0.09157621350977, 0.81684757298046,
                                                        0.81684757298046, 0.09157621350977);
static constexpr  Matrix<Real,7,2> QPSimplex5_qp_points(0.33333333333333, 0.33333333333333, 
                                                        0.47014206410511, 0.47014206410511,
                                                        0.47014206410511, 0.05971587178977,
                                                        0.05971587178977, 0.47014206410511, 
                                                        0.10128650732346, 0.10128650732346, 
                                                        0.10128650732346, 0.79742698535309, 
                                                        0.79742698535309, 0.10128650732346);


static constexpr  Vector<Real,1> QPSimplex1_qp_weights(1);                                                        
static constexpr  Vector<Real,3> QPSimplex2_qp_weights(0.33333333333333,
                                                       0.33333333333333,
                                                       0.33333333333333);
static constexpr  Vector<Real,4> QPSimplex3_qp_weights( -0.56250000000000,
                                                        0.52083333333333,
                                                        0.52083333333333,
                                                        0.52083333333333);
static constexpr  Vector<Real,6> QPSimplex4_qp_weights (0.22338158967801, 
                                                        0.22338158967801,
                                                        0.22338158967801,
                                                        0.10995174365532,
                                                        0.10995174365532, 
                                                        0.10995174365532);
static constexpr  Vector<Real,7> QPSimplex5_qp_weights (0.22500000000000, 
                                                        0.13239415278851, 
                                                        0.13239415278851, 
                                                        0.13239415278851, 
                                                        0.12593918054483, 
                                                        0.12593918054483, 
                                                        0.12593918054483 );









template< Integer Dim >
class GaussPoints< Simplex<Dim,2> , 1>: public BaseQuadrature<Simplex<Dim,2>,1>
{

public:
  using qp_points_type=Matrix<Real,1,2>;
  static constexpr Integer NQPoints=1;
  static constexpr  Matrix<Real,1,2> qp_points=QPSimplex1_qp_points;
  static constexpr  Vector<Real,1> qp_weights=QPSimplex1_qp_weights;
};


template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,2>: public BaseQuadrature<Simplex<Dim,2>,2>
{

public:
  using qp_points_type=Matrix<Real,3,2>;
  static constexpr Integer NQPoints=3;
  static constexpr  Matrix<Real,3,2> qp_points=QPSimplex2_qp_points;
  static constexpr  Vector<Real,3> qp_weights=QPSimplex2_qp_weights;
};





template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,3>: public BaseQuadrature<Simplex<Dim,2>,3>
{

public:
  using qp_points_type=Matrix<Real,4,2>;
  static constexpr Integer NQPoints=4;
  static constexpr  Matrix<Real,4,2> qp_points=QPSimplex3_qp_points;
  static constexpr  Vector<Real,4> qp_weights=QPSimplex3_qp_weights;

};

template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,4>: public BaseQuadrature<Simplex<Dim,2>,4>
{

public:
  using qp_points_type=Matrix<Real,6,2>;
  static constexpr Integer NQPoints=6;
  static constexpr  Matrix<Real,6,2> qp_points=QPSimplex4_qp_points;
  static constexpr  Vector<Real,6> qp_weights=QPSimplex4_qp_weights;

};


template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,5>: public BaseQuadrature<Simplex<Dim,2>,5>
{
public:
  using qp_points_type=Matrix<Real,7,2>;
  static constexpr Integer NQPoints=7;
  static constexpr  Matrix<Real,7,2> qp_points=QPSimplex5_qp_points;
  static constexpr  Vector<Real,7> qp_weights=QPSimplex5_qp_weights;

};


template< Integer Dim >
constexpr  Matrix<Real,1,2> GaussPoints< Simplex<Dim,2> , 1>::qp_points;
template< Integer Dim >
constexpr  Vector<Real,1> GaussPoints< Simplex<Dim,2> , 1>::qp_weights;

template< Integer Dim >
constexpr  Matrix<Real,3,2> GaussPoints< Simplex<Dim,2> , 2>::qp_points;
template< Integer Dim >
constexpr  Vector<Real,3> GaussPoints< Simplex<Dim,2> , 2>::qp_weights;

template< Integer Dim >
constexpr  Matrix<Real,4,2> GaussPoints< Simplex<Dim,2> , 3>::qp_points;
template< Integer Dim >
constexpr  Vector<Real,4> GaussPoints< Simplex<Dim,2> , 3>::qp_weights;

template< Integer Dim >
constexpr  Matrix<Real,6,2> GaussPoints< Simplex<Dim,2> , 4>::qp_points;
template< Integer Dim >
constexpr  Vector<Real,6> GaussPoints< Simplex<Dim,2> , 4>::qp_weights;

template< Integer Dim >
constexpr  Matrix<Real,7,2> GaussPoints< Simplex<Dim,2> , 5>::qp_points;
template< Integer Dim >
constexpr  Vector<Real,7> GaussPoints< Simplex<Dim,2> , 5>::qp_weights;







template<Integer QuadratureRuleType>
class QuadratureRule
{
public:
  template<typename Elem,Integer Order>
  using rule=GaussPoints< Elem , Order>;
};


template<>
class QuadratureRule<GaussianQuadrature>
{
public:
  template<typename Elem,Integer Order>
  using rule=GaussPoints< Elem , Order>;
};


template<Integer Order>
class GaussQP
{
public:
  template<typename Elem>
  using rule=GaussPoints< Elem , Order>;
};


}

#endif //MARS_QUADRATURE_RULES_HPP
