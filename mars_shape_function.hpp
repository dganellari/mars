// Ho una classe ElemIntegral< Elem, QuadratureRule, Order> che ha un metodo integral che prende una funzione f,
// e calcola l'integrale di f sull'elemento (di riferimento?), quindi:
// ElemIntegral< Elem, QuadratureRule, Order>::integral(f)
// oppure f
// ElemIntegral< Elem, QuadratureRule, Order>::integral(f,actual_elem)
// per calcolare l'integrale sull'elemento attuale.

// Quindi, integral prende i pesi e li moltiplica per f valutata nei punti di quadratura.

// Devo dare f come input generico. Quindi devo poter creare f(p)=(A+1)(f_1(p))*f_2(p)*f_3(p).
// Oppure f=(A+1)(f_1) * f_2 * f_3



// class Function;

// class FunctionToIntegrate
// {

//   private:
//          Function f_;

//     FunctionToIntegrate operator * (const FunctionToIntegrate& rhs) const
//     {
//         FunctionToIntegrate p;
//         p.f_=this->f_ * rhs.f_;
//         p.a = this;
//         p.b = &rhs;
//         return p;
//     }

// }

// function(Point p)
// {return {this->f_(p)* rhs.f_(p);}


// template<typename T>
// {


  
// }

//   class gradient() prende una funzione di base generica  ene ritorna il gradiente

//   grad(phi)(p)
//   cosi posso fare:
//   gradient(shape_function)            
//   oppure
//   A= (I+div)(shape_function)
//   dove definisco I(shape_function)=shape_function e div(shape_function)=div_shape_function







#ifndef MARS_SHAPE_FUNCTION_HPP
#define MARS_SHAPE_FUNCTION_HPP

#include "mars_simplex.hpp"
#include "mars_base_elementfunctionspace.hpp"
#include "mars_vector.hpp"

namespace mars{




template<typename TensorType>
class DotProduct;


template<>
class DotProduct<Real>
{

};

template<Integer Dim>
class DotProduct<Vector<Real,Dim>>
{

 public:
 static const Real compute(const Vector<Real,Dim>& v1, const Vector<Real,Dim>& v2, const double& alpha=1)
 {
  Real result=0.0;
  for(Integer nn=0;nn<Dim;nn++)
  {
   result+=v1[nn]*v2[nn];
  }
  return result*alpha;
 };
};

template<Integer Rows, Integer Cols>
class DotProduct<Matrix<Real,Rows,Cols>>
{

 public:


 static const Real compute(const Matrix<Real,Rows,Cols>& m1, const Matrix<Real,Rows,Cols>& m2, const Real& alpha=1)
 {
  Real result=0.0;
  for(Integer nn=0;nn<Rows;nn++)
    {
      for(Integer mm=0;mm<Cols;mm++)
      {
       result+=m1(nn,mm)*m2(nn,mm);
      }
    }
  return result*alpha;
 }

  static const Real compute(const Matrix<Real,Rows,Cols>& m1, const Matrix<Real,Rows,Cols>& m2, const Matrix<Real,Rows,Cols>& m3)
 {
    Real result=0.0;
    for(Integer nn=0;nn<Rows;nn++)
      {
        for(Integer mm=0;mm<Cols;mm++)
        {
         result+=m1(nn,mm)*m2(nn,mm)*m3(nn,mm);
        }
      } 
    return result;
   }
};




template<typename TensorType, Integer NComponents=1>
class QPVector;


template<Integer Rows>
class QPVector<Vector<Real,Rows>,1>
{
 public:
 // compute (U,V)
 static const Vector<Real,Rows> compute(const Vector<Real,Rows>& U, 
                                        const Vector<Real,Rows>& V)
 {
  Vector<Real,Rows> result;
  for(Integer mm=0;mm<Rows;mm++)
       result[mm]=U[mm]*V[mm];
  return result;
 }
};


template<Integer Rows, Integer Cols, Integer NComponents>
class QPVector<Matrix<Real,Rows,Cols>,NComponents>
{

 public:

 // compute (U,V)
 static const Vector<Real,Rows> compute(const Matrix<Real,Rows,Cols>& U, 
                                        const Matrix<Real,Rows,Cols>& V)
 {
  Vector<Real,Rows> result;

  for(Integer mm=0;mm<Rows;mm++)
    {
      result[mm]=0;
      for(Integer nn=0;nn<Cols;nn++)
      {
       result[mm]+=U(mm,nn)*V(mm,nn);
      }
    }
  return result;
 }


 // compute (U,V)
 static const Vector<Real,Rows> compute(const Real& alpha, 
                                        const Matrix<Real,Rows,Cols>& U,         
                                        const Matrix<Real,Rows,Cols>& V)
 {
  Vector<Real,Rows> result;

  for(Integer mm=0;mm<Rows;mm++)
    {
      result[mm]=0;
      for(Integer nn=0;nn<Cols;nn++)
      {
       result[mm]+=U(mm,nn)*V(mm,nn);
      }
      result[mm]*=alpha;
    }
  return result;
 }

 // compute (A U, V),A=const
 static const Vector<Real,Rows> compute(const Matrix<Real,Cols,Cols>& A, const Matrix<Real,Rows,Cols>& U, 
                                                                         const Matrix<Real,Rows,Cols>& V)
 {
  Vector<Real,Rows> result;
  Vector<Real,Cols> rowU;
  Vector<Real,Cols> rowV;
  for(Integer mm=0;mm<Rows;mm++)
    {
      result[mm]=0;
      U.get_row(mm,rowU);
      V.get_row(mm,rowV);
      rowU=A*rowU;
      for(Integer nn=0;nn<Cols;nn++)
      {
       result[mm]+=rowU(nn)*rowV(nn);
      }
    }
  return result;
 }

 // compute (A U, V),A=const
 static const Vector<Real,Rows> compute(                                 const Matrix<Real,Rows,Cols>& U, 
                                        const Matrix<Real,Cols,Cols>& B, const Matrix<Real,Rows,Cols>& V)
 {
  Vector<Real,Rows> result;
  Vector<Real,Cols> rowU;
  Vector<Real,Cols> rowV;
  for(Integer mm=0;mm<Rows;mm++)
    {
      result[mm]=0;
      U.get_row(mm,rowU);
      V.get_row(mm,rowV);
      rowV=B*rowV;
      for(Integer nn=0;nn<Cols;nn++)
      {
       result[mm]+=rowU(nn)*rowV(nn);
      }
    }
  return result;
 }

 // compute (A U,B V),A=const, B =const
 static const Vector<Real,Rows> compute(const Matrix<Real,Cols,Cols>& A, const Matrix<Real,Rows,Cols>& U, 
                                        const Matrix<Real,Cols,Cols>& B, const Matrix<Real,Rows,Cols>& V)
 {
  Vector<Real,Rows> result;
  Vector<Real,Cols> rowU;
  Vector<Real,Cols> rowV;
  for(Integer mm=0;mm<Rows;mm++)
    {
      result[mm]=0;
      U.get_row(mm,rowU);
      V.get_row(mm,rowV);
      rowU=A*rowU;
      rowV=B*rowV;
      for(Integer nn=0;nn<Cols;nn++)
      {
       result[mm]+=rowU(nn)*rowV(nn);
      }
    }
  return result;
 }

 // compute (A U,B V),A=const, B =const

 static const Matrix< Vector<Real,Rows>, NComponents,NComponents> compute(
                                        const Vector< Matrix<Real,NComponents,NComponents>,Rows>& A,
                                        const Vector< Matrix<Real,NComponents,NComponents>,Rows>& B)
 {
  Matrix< Vector<Real,Rows>, NComponents,NComponents> mat;
  for(Integer qp=0;qp<Rows;qp++)
      {
       auto mat_qp=transpose(A[qp])*B[qp];
       for(Integer ii=0; ii< NComponents; ii++)
           for(Integer jj=0; jj< NComponents; jj++)
              mat(ii,jj)[qp]=mat_qp(ii,jj); 
      }
  return mat;
 }

};


template<Integer Rows>
using QPVecV=QPVector<Vector<Real,Rows>,1>;

template<Integer Rows,Integer Cols,Integer NComponents=1>
using QPVecM=QPVector<Matrix<Real,Rows,Cols>,NComponents>;

















template<typename Elem,Integer Order>
class GaussPoints;

template< Integer Dim >
class GaussPoints< Simplex<Dim,2> , 1>
{
private:
Matrix<Real,1,2> qp_points_;
Vector<Real,1> weights_;

public:
  static constexpr Integer n_qp_points=1;
  const Matrix<Real,1,2> qp_points()const {return qp_points_;};
  const Vector<Real,1> weights()const {return weights_;};

  GaussPoints<Simplex<Dim,2>,1>():
  qp_points_({0.33333333333333, 0.33333333333333}),
  weights_({1})
   {}; 

};



template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,2>
{
private:
  Matrix<Real,3,2>  qp_points_;
  Vector<Real,3> weights_;

public:
  static constexpr Integer n_qp_points=3;
  const Matrix<Real,3,2> qp_points()const {return qp_points_;};
  const Vector<Real,3> weights()const {return weights_;};
  GaussPoints<Simplex<Dim,2>,2>():
  qp_points_({0.16666666666667, 0.16666666666667,
              0.16666666666667, 0.66666666666667,
              0.16666666666667, 0.16666666666667}),
  weights_({0.33333333333333,0.33333333333333,0.33333333333333})
  {}
};




template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,3>
{
private:
  Matrix<Real,4,2>  qp_points_;
  Vector<Real,4> weights_;

public:
  static constexpr Integer n_qp_points=4;
  const Matrix<Real,4,2> qp_points()const {return qp_points_;};
  const Vector<Real,4> weights()const {return weights_;};
  GaussPoints<Simplex<Dim,2>,3>():
  qp_points_({0.33333333333333, 0.33333333333333,
              0.20000000000000, 0.20000000000000,
              0.20000000000000, 0.60000000000000,
              0.60000000000000, 0.20000000000000}),
  weights_({ -0.56250000000000,0.52083333333333,0.52083333333333,0.52083333333333})
  {}
};

template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,4>
{
private:
  Matrix<Real,6,2>  qp_points_;
  Vector<Real,6> weights_;

public:
  static constexpr Integer n_qp_points=6;
  const Matrix<Real,6,2> qp_points()const {return qp_points_;};
  const Vector<Real,6> weights()const {return weights_;};
  GaussPoints<Simplex<Dim,2>,4>():
  qp_points_({0.44594849091597, 0.44594849091597,
              0.44594849091597, 0.10810301816807,
              0.10810301816807, 0.44594849091597, 
              0.09157621350977, 0.09157621350977, 
              0.09157621350977, 0.81684757298046,
              0.81684757298046, 0.09157621350977} ),
  weights_({0.22338158967801, 
            0.22338158967801,
            0.22338158967801,
            0.10995174365532,
            0.10995174365532, 
            0.10995174365532})
  {}
};


template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,5>
{
private:
  Matrix<Real,7,2>  qp_points_;
  Vector<Real,7> weights_;

public:
  static constexpr Integer n_qp_points=7;
  const Matrix<Real,7,2> qp_points()const {return qp_points_;};
  const Vector<Real,7> weights()const {return weights_;};
  GaussPoints<Simplex<Dim,2>,5>():
  qp_points_(   {0.33333333333333, 0.33333333333333, 
                 0.47014206410511, 0.47014206410511,
                 0.47014206410511, 0.05971587178977,
                 0.05971587178977, 0.47014206410511, 
                 0.10128650732346, 0.10128650732346, 
                 0.10128650732346, 0.79742698535309, 
                 0.79742698535309, 0.10128650732346}  ),
  weights_({0.22500000000000, 
            0.13239415278851, 
            0.13239415278851, 
            0.13239415278851, 
            0.12593918054483, 
            0.12593918054483, 
            0.12593918054483 })
  {}
};










template<Integer Ndofs_, Integer ManifoldDim, Integer ShapeFunctionDim>
class BaseShapeFunction
{public:
  static constexpr Integer Ndofs=Ndofs_;
  virtual ~BaseShapeFunction(){};
  virtual const Matrix<Real,Ndofs,ShapeFunctionDim>
  phi(const Vector<Real,ManifoldDim>&)=0;



  template<Integer NQPpoints>
  const Vector<Matrix<Real,NQPpoints,ShapeFunctionDim>,Ndofs>
  shape_function(const Matrix<Real,NQPpoints,2>& qp_points)
  {
        Vector<Matrix<Real,NQPpoints,ShapeFunctionDim>,Ndofs> shapefunction;
        Vector<Real,ManifoldDim> point;
        for(Integer nn=0;nn<Ndofs;nn++)
          {
            for(Integer qp=0;qp<NQPpoints;qp++)
            {
              qp_points.get_row(qp,point);
              auto& sf=this->phi(point);
              for(Integer dim=0;dim<ShapeFunctionDim;dim++)
              shapefunction[nn](qp,dim)=sf(nn,dim);  
            }            
          }

        return shapefunction;};
};




template<typename Elem,typename BaseFunctionSpace>
class ShapeFunction;








template<Integer NComponents_>
class ShapeFunction<Simplex<2,2>, Lagrange1<NComponents_> > : 
public BaseShapeFunction<3,2,1>
{
public:
  static constexpr Integer ShapeFunctionDim=1;
  static constexpr Integer NComponents=NComponents_;
  static constexpr Integer Ndofs=3;

  virtual const Matrix<Real,3,1>
  phi(const Vector<Real,2>& point)
       {const auto& xi=point[0];
        const auto& eta=point[1];
        const Real zeta = 1. - xi - eta;
        const Matrix<Real,3,1> shapefunction{zeta,xi,eta};
        return shapefunction;};


 const Matrix<Real,3,2> 
 grad_phi(const Vector<Real,2>& point)
       {const Matrix<Real,3,2> grad{-1,-1,1,0,0,1} ;
        return grad;};
};


template<Integer NComponents_>
class ShapeFunction<Simplex<2,2>, Lagrange2<NComponents_> > : 
public BaseShapeFunction<6,2,1>
{
public:
  static constexpr Integer ShapeFunctionDim=1;
  static constexpr Integer NComponents=NComponents_;
  static constexpr Integer Ndofs=6;

      virtual const Matrix<Real,6,1>
      phi(const Vector<Real,2>& point) 
      {
          const auto& xi=point[0];
          const auto& eta=point[1];
          const Real zeta = 1. - xi - eta;
          Matrix<Real,6,1> shape_function{2.*zeta*(zeta-0.5),
                                        2.*xi*(xi-0.5),
                                        2.*eta*(eta-0.5),
                                        4.*zeta*xi,
                                        4.*xi*eta, 
                                        4.*eta*zeta };
          return shape_function;
      };


 const Matrix<Real,3,2> 
 grad_phi(const Vector<Real,2>& point)
       {const Matrix<Real,3,2> grad{-1,-1,1,0,0,1} ;
        return grad;};
};





template<Integer NComponents_>
class ShapeFunction<Simplex<2,2>, RT0<NComponents_> > : 
public BaseShapeFunction<3,2,2>
{
public:
  static constexpr Integer ShapeFunctionDim=2;
  static constexpr Integer NComponents=NComponents_;
  static constexpr Integer Ndofs=3;

  virtual const  Matrix<Real,3,2>
  phi(const Vector<Real,2>& point)
       {const auto& xi=point[0];
        const auto& eta=point[1];
        const Matrix<Real,3,2> shapefunction{xi,eta-1,xi,eta,1-xi,-eta};
        return shapefunction;};

 const Real
 grad_phi(const Vector<Real,2>& point)
       {assert("gradient is not defined for RT elements");};
};













template<typename T, Integer Rows1,Integer Cols1,Integer Rows2,Integer Cols2>
Matrix<Real, Rows1*Rows2, Cols1*Cols2 > tensorproduct(const Matrix<T, Rows1, Cols1>& mat1,const Matrix<T, Rows2, Cols2>& mat2 )
{
 Matrix<Real, Rows1*Rows2, Cols1*Cols2 > mat;

  for(Integer i1=0;i1<Rows1;i1++)
  {
    for(Integer j1=0;j1<Cols1;j1++)
    {
      for(Integer i2=0;i2<Rows2;i2++)
      {
       for(Integer j2=0;j2<Cols2;j2++)
        {
          mat(i1*Rows2+i2 ,j1*Cols2+j2)= DotProduct<T>::compute(mat1(i1,j1),mat2(i2,j2));
        }
      }
    }
  }
  return mat;
}



template <typename T>
T sort_indexes(const T &v) {
  // initialize original index locations
  T idx;
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

template <typename T>
T T_sort(const T &v, const T& index) {
  // initialize original index locations
  T sortv;
  for(Integer ii=0;ii<v.size();ii++)
      sortv[ii]=v[index[ii]];
  return sortv;
}

template <typename T>
T T_sort(const T &v, Integer* index) {
  // initialize original index locations
  T sortv;
  for(Integer ii=0;ii<v.size();ii++)
      sortv[ii]=v[index[ii]];
  return sortv;
}


// template <typename T>
// T sort_indexes(const T &v,const Integer& size) {
//   // initialize original index locations
//   T idx(size);
//   std::iota(idx.begin(), idx.begin()+size, 0);

//   // sort indexes based on comparing values in v
//   std::sort(idx.begin(), idx.end(),
//        [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

//   return idx;
// }

// template <typename T>
// std::vector<T> sub_vector(const std::vector<T> &v,const std::vector<Integer>& index) {


//   Integer size=index.size();
//   std::vector<Integer> subvec(size);
//   for(Integer ii=0;ii<size;ii++)
//       subvec[ii]=v[index[ii]];

//   return subvec;
// }


template <std::size_t SDim, std::size_t TDim, typename S>
std::array<S,TDim> sub_vector(const std::array<S,SDim> &v,const std::array<S,TDim>& index) {

  //static_assert(IndexDim<=Dim,"in sub_vector the length of the index vector must be smaller than the one of the vector");
  std::array<S,TDim> subvec;
  for(Integer ii=0;ii<TDim;ii++)
      subvec[ii]=v[index[ii]];

  return subvec;
}




template <typename TrialFunction,typename TestFunction, typename MeshT >
void MassIntegrator(const MeshT& mesh)
{




std::cout<<"massintegrator1"<<std::endl;
using Elem=typename TrialFunction::Elem;
const Integer Dim=Elem::Dim;
const Integer ManifoldDim=Elem::ManifoldDim;    
using Point    = mars::Vector<Real, Dim>;

constexpr Integer n_trial_dofs=FunctionSpaceDofsPerElem<TrialFunction>::value;
constexpr Integer n_test_dofs=FunctionSpaceDofsPerElem<TestFunction>::value;


std::array<std::array<Real,n_trial_dofs>,n_test_dofs> matrix;

auto& trial_entities_nums= TrialFunction::entities_nums;  
auto& trial_dofs_per_entity= TrialFunction::dofs_per_entity;
auto& test_entities_nums= TestFunction::entities_nums;  
auto& test_dofs_per_entity= TestFunction::dofs_per_entity;


Integer row=0;

for(Integer nn1=0;nn1<test_entities_nums;nn1++)
{
   // constexpr Integer test_entity_dim=test_dofs_per_entity[nn1];
   // const auto& n_test_entity=ElemEntityCombinations<Elem,test_entity_dim>::value;
   // for(Integer mm1=0;mm1<n_test_entity;mm1++)
   // {
   //    Integer col=0;
   //      for(Integer nn2=0;nn2<trial_entities_nums;nn2++)
   //      {
   //          const auto& trial_entity_dim=trial_dofs_per_entity[nn2];
   //          const auto& n_trial_entity=ElemEntityCombinations<Elem,trial_entity_dim>::value;
   //          for(Integer mm2=0;mm2<n_trial_entity;mm2++)
   //          {
   //          // for(Integer pp=0;pp<n_qp_points;pp++)
   //          //      mass(matrix[dato][col],test,trial, point p);
   //          col++;
   //          }
   //      }  
   // }
}

using trialfun=Elem2FunctionSpace<TrialFunction>;
constexpr Integer QPOrder=4;
GaussPoints<Elem,QPOrder> gauss;

//Matrix<Real,Ndofs*NComponents,Ndofs*NComponents> mass;

ShapeFunction<Elem, trialfun> shapef;
auto sfall2=shapef.shape_function(gauss.qp_points());
std::cout<<" inheritance"<<std::endl;
sfall2[0].describe(std::cout);
sfall2[1].describe(std::cout);
sfall2[2].describe(std::cout);
std::cout<<"massintegrator2"<<std::endl;

constexpr Integer ShapeFunctionDim=ShapeFunction<Elem, trialfun>::ShapeFunctionDim;
constexpr Integer Ndofs=ShapeFunction<Elem, trialfun>::Ndofs;
constexpr Integer n_qp_points=GaussPoints<Elem,QPOrder>::n_qp_points; 
constexpr Integer NComponents=ShapeFunction<Elem, trialfun>::NComponents;
Matrix<Real,NComponents,NComponents> matrix_a{1.0};
Matrix<Real,NComponents,NComponents> matrix_b{1.0};


//Matrix<Real,NComponents,NComponents> matrix_a{1.0,0.0,0.0,1.0};
//Matrix<Real,NComponents,NComponents> matrix_b{1.0,2.0,2.0,3.0};
Matrix<Real,Ndofs,Ndofs> mass;

std::cout<<"massintegrator3"<<std::endl;

Vector<decltype(matrix_a),n_qp_points> mata(matrix_a);
Vector<decltype(matrix_b),n_qp_points> matb(matrix_b);
std::cout<<" mata"<<std::endl;

for(Integer qp=0;qp<n_qp_points;qp++)
    mata[qp].describe(std::cout);

Matrix< Vector<Real,n_qp_points>, NComponents,NComponents> qpmat=QPVecM<n_qp_points,ShapeFunctionDim,NComponents>::compute(mata,mata);
Matrix<Vector<Real,n_qp_points> ,Ndofs,Ndofs> mat_mass;


// vector, n_qp_points long, whose components are the matrix A evaluated in different qp points
//Vector< Matrix<Real,NComponents,NComponents> , n_qp_points> matA;
for(Integer nn=0;nn<Ndofs;nn++)
   for(Integer mm=0;mm<Ndofs;mm++)
   {

    // Quello che voglio fare alla fine e' una combinazione di matrici. Per H1 ad esempio:
    // (A gradu, gradv) + (B u,v) = (weights,  A GradUGradV+BUV)
    // (A divS,B divT)
    // dove BUV e' il vettore 
    // un elemento e' sempre il prodotto scalare tra wieghts e un'altro vettore lungo n_qp_points
    // quindi se ho il prodotto tra due Matrix(n_qppoints,Dim) mat, devo restituire
    // Vector(n_qppoints ) vec, dove vec(qp)=sum( mat1(qp,ii)*mat2(qp,ii))

    
    auto& qpvec=QPVecM<n_qp_points,ShapeFunctionDim>::compute( sfall2[nn], sfall2[mm]);
    mat_mass(nn,mm)=QPVecV<n_qp_points>::compute(qpvec,gauss.weights());
    mass(nn,mm)=DotProduct<Vector<Real,n_qp_points>>::compute(qpvec,gauss.weights());
   }
std::cout<<"mass"<<std::endl;
mass.describe(std::cout);
std::cout<<"mat_mass"<<std::endl;
mat_mass.describe(std::cout);
std::cout<<"qpmat"<<std::endl;
qpmat.describe(std::cout);
std::cout<<"ECCO2"<<std::endl;
auto mattensor=tensorproduct(mat_mass,qpmat);
mattensor.describe(std::cout);


for(Integer ee=0;ee<mesh.n_elements();ee++)
{
Elem elem=mesh.elem(ee);
std::array<Integer,ManifoldDim+1> nodes=elem.nodes;
std::vector<Vector<Real,Dim>> points;
std::vector<Vector<Real,Dim>> facepoints(ManifoldDim);
std::vector<Vector<Real,Dim>> elempoints(ManifoldDim+1);

Integer facenodes_tmp[ManifoldDim];
std::array< Integer, ManifoldDim> facenodes;

Integer elemnodes[ManifoldDim+1];
std::vector<Integer> reordernodes(ManifoldDim+1);



for(Integer nn=0;nn<ManifoldDim+1 ;nn++)
  {
    Integer tmp[1];
    Combinations<ManifoldDim + 1,1>::generate(nn,tmp);
    //std::cout<<tmp[0]<<" ";
    elemnodes[nn]=tmp[0];
  }
  std::cout<<std::endl;
std::reverse(std::begin(elemnodes), std::end(elemnodes));
for(Integer nn=0;nn<ManifoldDim+1 ;nn++)
  {
    //std::cout<<elemnodes[nn]<<" ";
  }
 std::cout<<"elem=="<<ee<<std::endl;

for(Integer mm=0;mm<nodes.size();mm++)
 {points.push_back(mesh.point(nodes[mm]));
  std::cout<<points[mm];
 }
 //std::cout<<std::endl;


 

for(Integer mm=0;mm<nodes.size();mm++)
   {
    //std::cout<<nodes[mm]<<" ";
  }
    

for(Integer mm=0;mm<ManifoldDim+1;mm++)
   {
    Combinations<ManifoldDim + 1,ManifoldDim>::generate(mm,facenodes_tmp);

    std::copy(std::begin(facenodes_tmp), std::end(facenodes_tmp), std::begin(facenodes));
    // for(Integer nn=0;nn<nodes.size();nn++)
    //     std::cout<<"nodes="<<nodes[nn]<<std::endl;
    // for(Integer nn=0;nn<facenodes.size();nn++)
    //     std::cout<<"facenodes="<<facenodes[nn]<<std::endl;      

    std::array<Integer,ManifoldDim> subvec=sub_vector(nodes,facenodes);
    auto sorted=sort_indexes(subvec);
    auto tmp_sort=T_sort(facenodes,&sorted[0]);


    // std::cout<<"elemnodes[mm]=="<<elemnodes[mm]<<std::endl;
    // std::cout<<"subvec: "<<std::endl;
    // for(Integer nn=0;nn<subvec.size();nn++)
    //    std::cout<<subvec[nn]<<" ";
    // std::cout<<"sorted: "<<std::endl;
    // for(Integer nn=0;nn<sorted.size();nn++)
    //    std::cout<<sorted[nn]<<" ";
    // std::cout<<std::endl;


    for(Integer nn=0;nn<ManifoldDim;nn++)
        reordernodes[nn]=tmp_sort[nn];
    reordernodes[ManifoldDim]=elemnodes[mm];
    std::cout<<"reordernodes: "<<std::endl;
    for(Integer nn=0;nn<reordernodes.size();nn++)
       std::cout<<reordernodes[nn]<<", ";

    //std::cout<<std::endl;
    //std::cout<<"facepoints: "<<std::endl;
    for(Integer nn=0;nn<ManifoldDim;nn++)
    {
      facepoints[nn]=points[facenodes[nn]];
      //std::cout<<facepoints[nn]<<std::endl;
    }
    for(Integer nn=0;nn<ManifoldDim+1;nn++)
    {
      elempoints[nn]=points[reordernodes[nn]];
      }
    Triangle2 tri2;
    tri2.nodes[0] = 0;
    tri2.nodes[1] = 1;
    tri2.nodes[2] = 2;
    std::cout<<std::endl<<"normal "<<std::endl;
    auto n2 = normal(tri2, elempoints);
    for(Integer nn=0;nn<ManifoldDim;nn++)
        std::cout<<n2[nn]<<" ";
    std::cout<<std::endl;
          }
std::cout<<std::endl;

}


    Triangle2 tri2;
    tri2.nodes[0] = 0;
    tri2.nodes[1] = 1;
    tri2.nodes[2] = 2;


    Tetrahedron4 tet; 
    tet.nodes[0] = 0;
    tet.nodes[1] = 1;
    tet.nodes[2] = 2;
    tet.nodes[3] = 3;

    std::vector<Vector2r> points2(
    {
      
      
       { 1., 0. },
       { 0., 0. },
       { 0., 1. }
    });

    std::vector<Vector4r> points4(
    {
      { 2., 0., 0., 0. },
      { 0., 2., 0., 0. },
      { 0., 0., 2., 0. },
      { 0., 0., 0., 2. },
      { 0., 0., 0., 0. }
      
      
      
      
    });
    auto n2 = normal(tri2, points2);
    auto n4 = normal(tet, points4);

    n2.describe(std::cout);
    n4.describe(std::cout);
    std::cout << n4 << std::endl;

};


// template< Integer Dim>
// shape_function<Simplex<Dim,2>, RaviartThomasFE, 0 >(const std::array<double,2>& point, std::array<double,3>)
// {
// const auto& xi=point[0];
// const auto& eta=point[1];
// const double zeta = 1. - xi - eta;

// shape_function[0]={xi,eta-1};
// shape_function[0]={xi-1,eta};
// shape_function[0]={xi,eta};
// }





}



#endif