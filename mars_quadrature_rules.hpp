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
static constexpr  Matrix<Real,4,2> QPSimplex3_qp_points(0.15505102572168219018027159252941,0.17855872826361642311703513337422,
                                                        0.64494897427831780981972840747059,0.075031110222608118177475598324603,
                                                        0.15505102572168219018027159252941,0.66639024601470138670269327409637,
                                                        0.64494897427831780981972840747059,0.28001991549907407200279599420481);
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


static constexpr  Array<Real,1> QPSimplex1_qp_weights(1); 
static constexpr  Array<Real,1> QPSimplex1_qp_sqrt_abs_weights(1);                                                        
static constexpr  Array<Real,3> QPSimplex2_qp_weights(0.33333333333333,
                                                       0.33333333333333,
                                                       0.33333333333333);
static constexpr  Array<Real,3> QPSimplex2_qp_sqrt_abs_weights(
                                                       0.577350269189623,
                                                       0.577350269189623,
                                                       0.577350269189623);
static constexpr  Array<Real,4> QPSimplex3_qp_weights(0.15902069087198858469718450103758,
                                                      0.090979309128011415302815498962418,
                                                      0.15902069087198858469718450103758,
                                                      0.090979309128011415302815498962418);
static constexpr  Array<Real,4> QPSimplex3_qp_sqrt_abs_weights(
                                                     0.398773984698085,
                                                     0.301627765843948,
                                                     0.398773984698085,
                                                     0.301627765843948);
static constexpr  Array<Real,6> QPSimplex4_qp_weights (0.22338158967801, 
                                                        0.22338158967801,
                                                        0.22338158967801,
                                                        0.10995174365532,
                                                        0.10995174365532, 
                                                        0.10995174365532);
static constexpr  Array<Real,6> QPSimplex4_qp_sqrt_abs_weights(
                                                       0.472632615969328,
                                                       0.472632615969328,
                                                       0.472632615969328,
                                                       0.331589721878287,
                                                       0.331589721878287,
                                                       0.331589721878287);
static constexpr  Array<Real,7> QPSimplex5_qp_weights (0.22500000000000, 
                                                        0.13239415278851, 
                                                        0.13239415278851, 
                                                        0.13239415278851, 
                                                        0.12593918054483, 
                                                        0.12593918054483, 
                                                        0.12593918054483 );
static constexpr  Array<Real,7> QPSimplex5_qp_sqrt_abs_weights(
                                                       0.474341649025257, 
                                                       0.363860073089244, 
                                                       0.363860073089244, 
                                                       0.363860073089244, 
                                                       0.354879106943238, 
                                                       0.354879106943238, 
                                                       0.354879106943238);









template< Integer Dim >
class GaussPoints< Simplex<Dim,2> , 1>: public BaseQuadrature<Simplex<Dim,2>,1>
{

public:
  using qp_points_type=Matrix<Real,1,2>;
  static constexpr Integer NQPoints=1;
  static constexpr  Matrix<Real,1,2> qp_points=QPSimplex1_qp_points;
  static constexpr  Array<Real,1> qp_weights=QPSimplex1_qp_weights;
  static constexpr  Array<Real,1> qp_sqrt_abs_weights=QPSimplex1_qp_sqrt_abs_weights;
};


template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,2>: public BaseQuadrature<Simplex<Dim,2>,2>
{

public:
  using qp_points_type=Matrix<Real,3,2>;
  static constexpr Integer NQPoints=3;
  static constexpr  Matrix<Real,3,2> qp_points=QPSimplex2_qp_points;
  static constexpr  Array<Real,3> qp_weights=QPSimplex2_qp_weights;
  static constexpr  Array<Real,3> qp_sqrt_abs_weights=QPSimplex2_qp_sqrt_abs_weights;
};





template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,3>: public BaseQuadrature<Simplex<Dim,2>,3>
{

public:
  using qp_points_type=Matrix<Real,4,2>;
  static constexpr Integer NQPoints=4;
  static constexpr  Matrix<Real,4,2> qp_points=QPSimplex3_qp_points;
  static constexpr  Array<Real,4> qp_weights=QPSimplex3_qp_weights;
  static constexpr  Array<Real,4> qp_sqrt_abs_weights=QPSimplex3_qp_sqrt_abs_weights;
};

template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,4>: public BaseQuadrature<Simplex<Dim,2>,4>
{

public:
  using qp_points_type=Matrix<Real,6,2>;
  static constexpr Integer NQPoints=6;
  static constexpr  Matrix<Real,6,2> qp_points=QPSimplex4_qp_points;
  static constexpr  Array<Real,6> qp_weights=QPSimplex4_qp_weights;
  static constexpr  Array<Real,6> qp_sqrt_abs_weights=QPSimplex4_qp_sqrt_abs_weights;
};


template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,5>: public BaseQuadrature<Simplex<Dim,2>,5>
{
public:
  using qp_points_type=Matrix<Real,7,2>;
  static constexpr Integer NQPoints=7;
  static constexpr  Matrix<Real,7,2> qp_points=QPSimplex5_qp_points;
  static constexpr  Array<Real,7> qp_weights=QPSimplex5_qp_weights;
  static constexpr  Array<Real,7> qp_sqrt_abs_weights=QPSimplex5_qp_sqrt_abs_weights;
};


template< Integer Dim >
constexpr  Matrix<Real,1,2> GaussPoints< Simplex<Dim,2> , 1>::qp_points;
template< Integer Dim >
constexpr  Array<Real,1> GaussPoints< Simplex<Dim,2> , 1>::qp_weights;

template< Integer Dim >
constexpr  Matrix<Real,3,2> GaussPoints< Simplex<Dim,2> , 2>::qp_points;
template< Integer Dim >
constexpr  Array<Real,3> GaussPoints< Simplex<Dim,2> , 2>::qp_weights;

template< Integer Dim >
constexpr  Matrix<Real,4,2> GaussPoints< Simplex<Dim,2> , 3>::qp_points;
template< Integer Dim >
constexpr  Array<Real,4> GaussPoints< Simplex<Dim,2> , 3>::qp_weights;

template< Integer Dim >
constexpr  Matrix<Real,6,2> GaussPoints< Simplex<Dim,2> , 4>::qp_points;
template< Integer Dim >
constexpr  Array<Real,6> GaussPoints< Simplex<Dim,2> , 4>::qp_weights;

template< Integer Dim >
constexpr  Matrix<Real,7,2> GaussPoints< Simplex<Dim,2> , 5>::qp_points;
template< Integer Dim >
constexpr  Array<Real,7> GaussPoints< Simplex<Dim,2> , 5>::qp_weights;







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
