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

template<typename Elem,Integer Order>
class GaussPoints;

template<Integer Dim>
class MaxOrder<Simplex<Dim,1>,GaussianQuadrature>
{
public:
  static constexpr Integer value=10;
};




static constexpr Matrix<Real,1,1> QP1_Simplex1_qp_points(0.5000);
static constexpr Array<Real,1>    QP1_Simplex1_qp_weights(1.0);

static constexpr Matrix<Real,2,1> QP2_Simplex1_qp_points(7.886751345948129e-01,
                                                         2.113248654051871e-01);
static constexpr Array<Real,2>   QP2_Simplex1_qp_weights(0.5, 
                                                         0.5);

static constexpr Matrix<Real,3,1> QP3_Simplex1_qp_points(8.872983346207417e-01,
                                                         5.000000000000000e-01,
                                                         1.127016653792582e-01);
static constexpr Array<Real,3>   QP3_Simplex1_qp_weights(2.777777777777772e-01,
                                                         4.444444444444444e-01,  
                                                         2.777777777777772e-01);


static constexpr Matrix<Real,4,1> QP4_Simplex1_qp_points(9.305681557970262e-01,
                                                         6.699905217924281e-01,
                                                         3.300094782075719e-01,
                                                         6.943184420297371e-02);
static constexpr Array<Real,4>   QP4_Simplex1_qp_weights(1.739274225687272e-01,
                                                         3.260725774312730e-01,
                                                         3.260725774312730e-01,
                                                         1.739274225687272e-01);



static constexpr Matrix<Real,5,1> QP5_Simplex1_qp_points(9.530899229693319e-01,
                                                         7.692346550528415e-01,
                                                         5.000000000000000e-01,
                                                         2.307653449471584e-01,
                                                         4.691007703066802e-02);
static constexpr Array<Real,5>   QP5_Simplex1_qp_weights(1.184634425280946e-01,
                                                         2.393143352496831e-01,
                                                         2.844444444444444e-01,
                                                         2.393143352496831e-01,
                                                         1.184634425280946e-01);


static constexpr Matrix<Real,6,1> QP6_Simplex1_qp_points(9.662347571015760e-01,
                                                         8.306046932331322e-01,
                                                         6.193095930415985e-01,
                                                         3.806904069584016e-01,
                                                         1.693953067668678e-01,
                                                         3.376524289842397e-02);

static constexpr Array<Real,6>   QP6_Simplex1_qp_weights(8.566224618958525e-02,
                                                         1.803807865240693e-01,
                                                         2.339569672863455e-01,
                                                         2.339569672863455e-01,
                                                         1.803807865240693e-01,
                                                         8.566224618958525e-02);


static constexpr  Matrix<Real,7,1> QP7_Simplex1_qp_points(9.745539561713792e-01,
                                                         8.707655927996972e-01,
                                                         7.029225756886985e-01,
                                                         5.000000000000000e-01,
                                                         2.970774243113014e-01,
                                                         1.292344072003028e-01,
                                                         2.544604382862076e-02);

static constexpr Array<Real,7>   QP7_Simplex1_qp_weights(6.474248308443482e-02,
                                                         1.398526957446384e-01,
                                                         1.909150252525595e-01,
                                                         2.089795918367347e-01,
                                                         1.909150252525595e-01,
                                                         1.398526957446384e-01,
                                                         6.474248308443482e-02);



static constexpr  Matrix<Real,8,1> QP8_Simplex1_qp_points(9.801449282487682e-01,
                                                         8.983332387068134e-01,
                                                         7.627662049581645e-01,
                                                         5.917173212478249e-01,
                                                         4.082826787521751e-01,
                                                         2.372337950418355e-01,
                                                         1.016667612931866e-01,
                                                         1.985507175123186e-02);

static constexpr Array<Real,8>   QP8_Simplex1_qp_weights(5.061426814518841e-02,
                                                         1.111905172266872e-01,
                                                         1.568533229389437e-01,
                                                         1.813418916891811e-01,
                                                         1.813418916891811e-01,
                                                         1.568533229389437e-01,
                                                         1.111905172266872e-01,
                                                         5.061426814518841e-02);


static constexpr  Matrix<Real,9,1> QP9_Simplex1_qp_points(9.840801197538130e-01,
                                                          9.180155536633179e-01,
                                                          8.066857163502952e-01,
                                                          6.621267117019045e-01,
                                                          5.000000000000000e-01,
                                                          3.378732882980955e-01,
                                                          1.933142836497048e-01,
                                                          8.198444633668212e-02,
                                                          1.591988024618696e-02);

static constexpr  Array<Real,9>   QP9_Simplex1_qp_weights(4.063719418078732e-02,
                                                          9.032408034742870e-02,
                                                          1.303053482014678e-01,
                                                          1.561735385200015e-01,
                                                          1.651196775006299e-01,
                                                          1.561735385200015e-01,
                                                          1.303053482014678e-01,
                                                          9.032408034742870e-02,
                                                          4.063719418078732e-02);


static constexpr  Matrix<Real,10,1> QP10_Simplex1_qp_points(9.869532642585859e-01,
                                                            9.325316833444923e-01,
                                                            8.397047841495122e-01,
                                                            7.166976970646236e-01,
                                                            5.744371694908156e-01,
                                                            4.255628305091844e-01,
                                                            2.833023029353764e-01,
                                                            1.602952158504878e-01,
                                                            6.746831665550773e-02,
                                                            1.304673574141413e-02);

static constexpr   Array<Real,10>  QP10_Simplex1_qp_weights(3.333567215434401e-02,
                                                            7.472567457529028e-02,
                                                            1.095431812579910e-01,
                                                            1.346333596549981e-01,
                                                            1.477621123573765e-01,
                                                            1.477621123573765e-01,
                                                            1.346333596549981e-01,
                                                            1.095431812579910e-01,
                                                            7.472567457529028e-02,
                                                            3.333567215434401e-02);



template< Integer Dim >
class GaussPoints< Simplex<Dim,1> , 1>: public BaseQuadrature<Simplex<Dim,1>,1>
{
public:
  using qp_points_type=Matrix<Real,1,1>;
  static constexpr Integer NQPoints=1;
  static constexpr  Matrix<Real,1,1> qp_points=QP1_Simplex1_qp_points;
  static constexpr  Array<Real,1> qp_weights=QP1_Simplex1_qp_weights;
};

template< Integer Dim >
class GaussPoints< Simplex<Dim,1> , 2>: public BaseQuadrature<Simplex<Dim,1>,2>
{
public:
  using qp_points_type=Matrix<Real,2,1>;
  static constexpr Integer NQPoints=2;
  static constexpr  Matrix<Real,2,1> qp_points=QP2_Simplex1_qp_points;
  static constexpr  Array<Real,2> qp_weights=QP2_Simplex1_qp_weights;
};

template< Integer Dim >
class GaussPoints< Simplex<Dim,1> , 3>: public BaseQuadrature<Simplex<Dim,1>,3>
{
public:
  using qp_points_type=Matrix<Real,3,1>;
  static constexpr Integer NQPoints=3;
  static constexpr  Matrix<Real,3,1> qp_points=QP3_Simplex1_qp_points;
  static constexpr  Array<Real,3> qp_weights=QP3_Simplex1_qp_weights;
};

template< Integer Dim >
class GaussPoints< Simplex<Dim,1> , 4>: public BaseQuadrature<Simplex<Dim,1>,4>
{
public:
  using qp_points_type=Matrix<Real,4,1>;
  static constexpr Integer NQPoints=4;
  static constexpr  Matrix<Real,4,1> qp_points=QP4_Simplex1_qp_points;
  static constexpr  Array<Real,4> qp_weights=QP4_Simplex1_qp_weights;
};

template< Integer Dim >
class GaussPoints< Simplex<Dim,1> , 5>: public BaseQuadrature<Simplex<Dim,1>,5>
{
public:
  using qp_points_type=Matrix<Real,5,1>;
  static constexpr Integer NQPoints=5;
  static constexpr  Matrix<Real,5,1> qp_points=QP5_Simplex1_qp_points;
  static constexpr  Array<Real,5> qp_weights=QP5_Simplex1_qp_weights;
};

template< Integer Dim >
class GaussPoints< Simplex<Dim,1> , 6>: public BaseQuadrature<Simplex<Dim,1>,6>
{
public:
  using qp_points_type=Matrix<Real,6,1>;
  static constexpr Integer NQPoints=6;
  static constexpr  Matrix<Real,6,1> qp_points=QP6_Simplex1_qp_points;
  static constexpr  Array<Real,6> qp_weights=QP6_Simplex1_qp_weights;
};

template< Integer Dim >
class GaussPoints< Simplex<Dim,1> , 7>: public BaseQuadrature<Simplex<Dim,1>,7>
{
public:
  using qp_points_type=Matrix<Real,7,1>;
  static constexpr Integer NQPoints=7;
  static constexpr  Matrix<Real,7,1> qp_points=QP7_Simplex1_qp_points;
  static constexpr  Array<Real,7> qp_weights=QP7_Simplex1_qp_weights;
};

template< Integer Dim >
class GaussPoints< Simplex<Dim,1> , 8>: public BaseQuadrature<Simplex<Dim,1>,8>
{
public:
  using qp_points_type=Matrix<Real,8,1>;
  static constexpr Integer NQPoints=8;
  static constexpr  Matrix<Real,8,1> qp_points=QP8_Simplex1_qp_points;
  static constexpr  Array<Real,8> qp_weights=QP8_Simplex1_qp_weights;
};

template< Integer Dim >
class GaussPoints< Simplex<Dim,1> , 9>: public BaseQuadrature<Simplex<Dim,1>,9>
{
public:
  using qp_points_type=Matrix<Real,9,1>;
  static constexpr Integer NQPoints=9;
  static constexpr  Matrix<Real,9,1> qp_points=QP9_Simplex1_qp_points;
  static constexpr  Array<Real,9> qp_weights=QP9_Simplex1_qp_weights;
};

template< Integer Dim >
class GaussPoints< Simplex<Dim,1> , 10>: public BaseQuadrature<Simplex<Dim,1>,10>
{
public:
  using qp_points_type=Matrix<Real,10,1>;
  static constexpr Integer NQPoints=10;
  static constexpr  Matrix<Real,10,1> qp_points=QP10_Simplex1_qp_points;
  static constexpr  Array<Real,10> qp_weights=QP10_Simplex1_qp_weights;
};



template< Integer Dim >
constexpr  Matrix<Real,1,1> GaussPoints< Simplex<Dim,1> , 1>::qp_points;
template< Integer Dim >
constexpr  Array<Real,1> GaussPoints< Simplex<Dim,1> , 1>::qp_weights;


template< Integer Dim >
constexpr  Matrix<Real,2,1> GaussPoints< Simplex<Dim,1> , 2>::qp_points;
template< Integer Dim >
constexpr  Array<Real,2> GaussPoints< Simplex<Dim,1> , 2>::qp_weights;


template< Integer Dim >
constexpr  Matrix<Real,3,1> GaussPoints< Simplex<Dim,1> , 3>::qp_points;
template< Integer Dim >
constexpr  Array<Real,3> GaussPoints< Simplex<Dim,1> , 3>::qp_weights;


template< Integer Dim >
constexpr  Matrix<Real,4,1> GaussPoints< Simplex<Dim,1> , 4>::qp_points;
template< Integer Dim >
constexpr  Array<Real,4> GaussPoints< Simplex<Dim,1> , 4>::qp_weights;


template< Integer Dim >
constexpr  Matrix<Real,5,1> GaussPoints< Simplex<Dim,1> , 5>::qp_points;
template< Integer Dim >
constexpr  Array<Real,5> GaussPoints< Simplex<Dim,1> , 5>::qp_weights;


template< Integer Dim >
constexpr  Matrix<Real,6,1> GaussPoints< Simplex<Dim,1> , 6>::qp_points;
template< Integer Dim >
constexpr  Array<Real,6> GaussPoints< Simplex<Dim,1> , 6>::qp_weights;


template< Integer Dim >
constexpr  Matrix<Real,7,1> GaussPoints< Simplex<Dim,1> , 7>::qp_points;
template< Integer Dim >
constexpr  Array<Real,7> GaussPoints< Simplex<Dim,1> , 7>::qp_weights;


template< Integer Dim >
constexpr  Matrix<Real,8,1> GaussPoints< Simplex<Dim,1> , 8>::qp_points;
template< Integer Dim >
constexpr  Array<Real,8> GaussPoints< Simplex<Dim,1> , 8>::qp_weights;


template< Integer Dim >
constexpr  Matrix<Real,9,1> GaussPoints< Simplex<Dim,1> , 9>::qp_points;
template< Integer Dim >
constexpr  Array<Real,9> GaussPoints< Simplex<Dim,1> , 9>::qp_weights;


template< Integer Dim >
constexpr  Matrix<Real,10,1> GaussPoints< Simplex<Dim,1> , 10>::qp_points;
template< Integer Dim >
constexpr  Array<Real,10> GaussPoints< Simplex<Dim,1> , 10>::qp_weights;





template<Integer Dim>
class MaxOrder<Simplex<Dim,2>,GaussianQuadrature>
{
public:
  static constexpr Integer value=5;
};



static constexpr  Matrix<Real,1,2> QP1_Simplex2_qp_points(0.33333333333333, 0.33333333333333);
static constexpr  Matrix<Real,3,2> QP2_Simplex2_qp_points(0.16666666666667, 0.16666666666667,
                                                        0.16666666666667, 0.66666666666667,
                                                        0.16666666666667, 0.16666666666667);
static constexpr  Matrix<Real,4,2> QP3_Simplex2_qp_points(0.15505102572168219018027159252941,0.17855872826361642311703513337422,
                                                        0.64494897427831780981972840747059,0.075031110222608118177475598324603,
                                                        0.15505102572168219018027159252941,0.66639024601470138670269327409637,
                                                        0.64494897427831780981972840747059,0.28001991549907407200279599420481);
static constexpr  Matrix<Real,6,2> QP4_Simplex2_qp_points(0.44594849091597, 0.44594849091597,
                                                        0.44594849091597, 0.10810301816807,
                                                        0.10810301816807, 0.44594849091597, 
                                                        0.09157621350977, 0.09157621350977, 
                                                        0.09157621350977, 0.81684757298046,
                                                        0.81684757298046, 0.09157621350977);
static constexpr  Matrix<Real,7,2> QP5_Simplex2_qp_points(0.33333333333333, 0.33333333333333, 
                                                        0.47014206410511, 0.47014206410511,
                                                        0.47014206410511, 0.05971587178977,
                                                        0.05971587178977, 0.47014206410511, 
                                                        0.10128650732346, 0.10128650732346, 
                                                        0.10128650732346, 0.79742698535309, 
                                                        0.79742698535309, 0.10128650732346);


static constexpr  Array<Real,1> QP1_Simplex2_qp_weights(1); 
static constexpr  Array<Real,1> QP1_Simplex2_qp_sqrt_abs_weights(1);                                                        
static constexpr  Array<Real,3> QP2_Simplex2_qp_weights(0.33333333333333,
                                                       0.33333333333333,
                                                       0.33333333333333);
static constexpr  Array<Real,3> QP2_Simplex2_qp_sqrt_abs_weights(
                                                       0.577350269189623,
                                                       0.577350269189623,
                                                       0.577350269189623);
static constexpr  Array<Real,4> QP3_Simplex2_qp_weights(0.15902069087198858469718450103758,
                                                      0.090979309128011415302815498962418,
                                                      0.15902069087198858469718450103758,
                                                      0.090979309128011415302815498962418);
static constexpr  Array<Real,4> QP3_Simplex2_qp_sqrt_abs_weights(
                                                     0.398773984698085,
                                                     0.301627765843948,
                                                     0.398773984698085,
                                                     0.301627765843948);
static constexpr  Array<Real,6> QP4_Simplex2_qp_weights (0.22338158967801, 
                                                        0.22338158967801,
                                                        0.22338158967801,
                                                        0.10995174365532,
                                                        0.10995174365532, 
                                                        0.10995174365532);
static constexpr  Array<Real,6> QP4_Simplex2_qp_sqrt_abs_weights(
                                                       0.472632615969328,
                                                       0.472632615969328,
                                                       0.472632615969328,
                                                       0.331589721878287,
                                                       0.331589721878287,
                                                       0.331589721878287);
static constexpr  Array<Real,7> QP5_Simplex2_qp_weights (0.22500000000000, 
                                                        0.13239415278851, 
                                                        0.13239415278851, 
                                                        0.13239415278851, 
                                                        0.12593918054483, 
                                                        0.12593918054483, 
                                                        0.12593918054483 );
static constexpr  Array<Real,7> QP5_Simplex2_qp_sqrt_abs_weights(
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
  static constexpr  Matrix<Real,1,2> qp_points=QP1_Simplex2_qp_points;
  static constexpr  Array<Real,1> qp_weights=QP1_Simplex2_qp_weights;
  static constexpr  Array<Real,1> qp_sqrt_abs_weights=QP1_Simplex2_qp_sqrt_abs_weights;
};


template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,2>: public BaseQuadrature<Simplex<Dim,2>,2>
{

public:
  using qp_points_type=Matrix<Real,3,2>;
  static constexpr Integer NQPoints=3;
  static constexpr  Matrix<Real,3,2> qp_points=QP2_Simplex2_qp_points;
  static constexpr  Array<Real,3> qp_weights=QP2_Simplex2_qp_weights;
  static constexpr  Array<Real,3> qp_sqrt_abs_weights=QP2_Simplex2_qp_sqrt_abs_weights;
};





template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,3>: public BaseQuadrature<Simplex<Dim,2>,3>
{

public:
  using qp_points_type=Matrix<Real,4,2>;
  static constexpr Integer NQPoints=4;
  static constexpr  Matrix<Real,4,2> qp_points=QP3_Simplex2_qp_points;
  static constexpr  Array<Real,4> qp_weights=QP3_Simplex2_qp_weights;
  static constexpr  Array<Real,4> qp_sqrt_abs_weights=QP3_Simplex2_qp_sqrt_abs_weights;
};

template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,4>: public BaseQuadrature<Simplex<Dim,2>,4>
{

public:
  using qp_points_type=Matrix<Real,6,2>;
  static constexpr Integer NQPoints=6;
  static constexpr  Matrix<Real,6,2> qp_points=QP4_Simplex2_qp_points;
  static constexpr  Array<Real,6> qp_weights=QP4_Simplex2_qp_weights;
  static constexpr  Array<Real,6> qp_sqrt_abs_weights=QP4_Simplex2_qp_sqrt_abs_weights;
};


template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,5>: public BaseQuadrature<Simplex<Dim,2>,5>
{
public:
  using qp_points_type=Matrix<Real,7,2>;
  static constexpr Integer NQPoints=7;
  static constexpr  Matrix<Real,7,2> qp_points=QP5_Simplex2_qp_points;
  static constexpr  Array<Real,7> qp_weights=QP5_Simplex2_qp_weights;
  static constexpr  Array<Real,7> qp_sqrt_abs_weights=QP5_Simplex2_qp_sqrt_abs_weights;
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
