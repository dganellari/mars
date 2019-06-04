#ifndef MARS_EXAMPLES_HPP
#define MARS_EXAMPLES_HPP

#include "mars_simplex.hpp"
#include "mars_connectivity.hpp"
#include "mars_functionspace_dofmap.hpp"
#include "mars_functionspace.hpp"
#include "mars_shape_function.hpp"
#include "mars_matrix.hpp"



namespace mars{


using std::cout;
using std::endl;


































// Base Function: f(x) = x
template<typename T,typename...Parameters>
class Expression
{
    protected:
    
    struct Implementation
    {
    virtual ~Implementation() {};
    virtual T evaluate(const Parameters &... params)
    {// not used
    T t; return t;};
    };

    public:
    Expression()
    :   self_(std::make_shared<Implementation>())
    {}

    Expression(const std::shared_ptr<Implementation>& self)
    :   self_(self)
    {}

    T operator () (const Parameters &... params) const { return self_->evaluate(params...); }

    private:
    std::shared_ptr<Implementation> self_;
};




// Unary Function: u(-f(x))
template<typename T,typename...Parameters>
class ExpressionUnaryMinus : public Expression<T,Parameters...>
{
public:


    protected:
    struct Implementation : Expression<T,Parameters...>::Implementation
    {
        Expression<T,Parameters...> f;


        Implementation(const Expression<T,Parameters...>& f1):   
        f(f1)
        {};

        virtual T evaluate(const Parameters&... params) override
        {return (-f(params...)); }
    };

    public:
    ExpressionUnaryMinus(const Expression<T,Parameters...>& f)
    :   Expression<T,Parameters...>(std::make_shared<Implementation>(f))
    {}
};





// Binary Function: u(f(x) * g(x))
template<typename T,typename...Parameters>
class ExpressionBinaryMultiply :  public Expression<T,Parameters...>
{
    protected:
    struct Implementation : Expression<T,Parameters...>::Implementation
    {
        Expression<T,Parameters...> f;
        Expression<T,Parameters...> g;
        Implementation(const Expression<T,Parameters...>& f1,const Expression<T,Parameters...>& g1)
        :   f(f1), g(g1)
        {};

        virtual T evaluate(const Parameters &...params) override { return f(params...) * g(params...); }
    };

    public:
    ExpressionBinaryMultiply<T,Parameters...>(const Expression<T,Parameters...>& f,const Expression<T,Parameters...>& g)
    :   Expression<T,Parameters...>(std::make_shared<Implementation>(f, g))
    {}
};


// Binary Function: u(f(x) + g(x))
template<typename T,typename...Parameters>
class ExpressionBinaryAdd : public Expression<T,Parameters...>
{
    protected:
    struct Implementation : Expression<T,Parameters...>::Implementation
    {
        Expression<T,Parameters...> f;
        Expression<T,Parameters...> g;
        Implementation(const Expression<T,Parameters...>& f,const Expression<T,Parameters...>& g)
        :   f(f), g(g)
        {};

        virtual T evaluate(const Parameters &...params) override { return f(params...) + g(params...); }
    };

    public:
    ExpressionBinaryAdd(const Expression<T,Parameters...>& f, const Expression<T,Parameters...>& g)
    :   Expression<T,Parameters...>(std::make_shared<Implementation>(f, g))
    {}
};

// Binary Function: u(f(x) + g(x))
template<typename T,typename...Parameters>
class ExpressionBinarySubtract : public Expression<T,Parameters...>
{
    protected:
    struct Implementation : Expression<T,Parameters...>::Implementation
    {
        Expression<T,Parameters...> f;
        Expression<T,Parameters...> g;
        Implementation(Expression<T,Parameters...> f1, Expression<T,Parameters...> g1)
        :   f(f1), g(g1)
        {};

        virtual T evaluate(const Parameters &...params) override { return f(params...) - g(params...); }
    };

    public:
    ExpressionBinarySubtract(const Expression<T,Parameters...>& f,const Expression<T,Parameters...>& g)
    :   Expression<T,Parameters...>(std::make_shared<Implementation>(f, g))
    {}
};

// Binary Function: u(f(x) * g(x))
template<typename T,typename...Parameters>
class ExpressionRightScalarMultiply : public Expression<T,Parameters...>
{
    protected:
    struct Implementation : Expression<T,Parameters...>::Implementation
    {
        Expression<T,Parameters...> f;
        Real alpha_;
        Implementation(const Expression<T,Parameters...>& f, const Real& alpha):   
        f(f), 
        alpha_(alpha)
        {};


        virtual T evaluate(const Parameters &...params) override { return f(params...) * alpha_; }
    };

    public:
    ExpressionRightScalarMultiply(const Expression<T,Parameters...>& f, const Real& alpha)
    :   Expression<T,Parameters...>(std::make_shared<Implementation>(f, alpha))
    {}

};

// Binary Function: u(f(x) * g(x))
template<typename T,typename...Parameters>
class ExpressionLeftScalarMultiply : public Expression<T,Parameters...>
{
    protected:
    struct Implementation : Expression<T,Parameters...>::Implementation
    {
        Expression<T,Parameters...> f;
        Real alpha_;
        Implementation( const Real& alpha, const Expression<T,Parameters...>& f):   
        f(f), 
        alpha_(alpha)
        {};


        virtual T evaluate(const Parameters &...params) override { return f(params...) * alpha_; }
    };

    public:
    ExpressionLeftScalarMultiply(const Real& alpha, const Expression<T,Parameters...>& f)
    :   Expression<T,Parameters...>(std::make_shared<Implementation>(alpha,f))
    {}

};










// Binary Function: u(f(x) + g(x))
// template<typename Point,Integer Rows>
// class ExpressionBinaryMultiply<Vector<Real,Rows>,Point>
// {
//     protected:
//     template<Integer Cols>
//     struct Implementation : public Expression<Matrix<Real,Rows,Cols>,Point>::Implementation, public Expression<Vector<Real,Cols>,Point>::Implementation
//     {
//         Expression<Matrix<Real,Rows,Cols>,Point> mat;
//         Expression<Vector<Real,Cols>,Point> vec;
//         Implementation(const Expression<Matrix<Real,Rows,Cols>,Point>& mat1,const Expression<Vector<Real,Cols>,Point>& vec1)
//         :   mat(mat1), vec(vec1)
//         {};

//         virtual Vector<Real,Rows> evaluate(const Point &point) override { return mat(point) * vec(point); }
//     };

//     public:
//     template<Integer Cols>
//     ExpressionBinaryMultiply(const Expression<Matrix<Real,Rows,Cols>,Point>& mat,const Expression<Vector<Real,Cols>,Point>& vec)
//     :   Expression<Vector<Real,Rows>,Point>(std::make_shared<Implementation>(mat,vec))
//     {}
// };

template<typename Point,Integer Rows>
class ExpressionBinaryMultiply<Vector<Real,Rows>,Point>: public Expression<Vector<Real,Rows>,Point>
{
    protected:
    template<Integer Cols>
    struct Implementation: Expression<Vector<Real,Rows>,Point>::Implementation
    {   Expression<Matrix<Real,Rows,Cols>,Point> mat;
        Expression<Vector<Real,Cols>,Point> vec;
        Implementation(const Expression<Matrix<Real,Rows,Cols>,Point>& mat1,const Expression<Vector<Real,Cols>,Point>& vec1)
        :   mat(mat1),vec(vec1)
        {};

        virtual Vector<Real,Rows> evaluate(const Point &point) override { 
            return mat(point) * vec(point); }
    };

    public:
    template<Integer Cols>
    ExpressionBinaryMultiply(const Expression<Matrix<Real,Rows,Cols>,Point>& mat,const Expression<Vector<Real,Cols>,Point>& vec)
    :   Expression<Vector<Real,Rows>,Point>(std::make_shared<Implementation<Cols>>(mat,vec))
    {}
};

template<typename Point, Integer Rows,Integer Cols>
inline ExpressionBinaryMultiply<Vector<Real,Rows>,Point> operator *
(const Expression<Matrix<Real,Rows,Cols>,Point>& mat,const Expression<Vector<Real,Cols>,Point>& vec) 
{ return ExpressionBinaryMultiply<Vector<Real,Rows>,Point>(mat,vec);};




template<typename T,typename...Parameters>
inline ExpressionUnaryMinus<T,Parameters...> operator - 
(const Expression<T,Parameters...>& f) { return ExpressionUnaryMinus<T,Parameters...>(f); }

template<typename T,typename...Parameters>
inline ExpressionBinarySubtract<T,Parameters...> operator - 
(const Expression<T,Parameters...>& f,const Expression<T,Parameters...>& g) { return ExpressionBinarySubtract<T,Parameters...>(f,g); }


template<typename T,typename...Parameters>
inline ExpressionBinaryMultiply<T,Parameters...> operator * 
(const Expression<T,Parameters...>& f,const Expression<T,Parameters...>& g) { return ExpressionBinaryMultiply<T,Parameters...>(f,g); }

template<typename T,typename...Parameters>
inline ExpressionRightScalarMultiply<T,Parameters...> operator *
(const Expression<T,Parameters...>& f,const Real& alpha) { return ExpressionRightScalarMultiply<T,Parameters...>(f,alpha);}

template<typename T,typename...Parameters>
inline ExpressionLeftScalarMultiply<T,Parameters...> operator *
(const Real& alpha, const Expression<T,Parameters...>& f) { return ExpressionLeftScalarMultiply<T,Parameters...>(alpha,f);}



template<typename T,typename...Parameters>
inline ExpressionBinaryAdd<T,Parameters...> operator + 
(const Expression<T,Parameters...>& f,const Expression<T,Parameters...>& g) { return ExpressionBinaryAdd<T,Parameters...>(f,g); }







template<typename Point>
class ExpressionMatrix: public Expression<Matrix<Real,3,3>,Point>
{ 

public:
    static constexpr Integer Rows=3;
    static constexpr Integer Cols=3;

     protected:
    struct Implementation: Expression<Matrix<Real,Rows,Cols>,Point>::Implementation
    {
        Expression<Matrix<Real,Rows,Cols>,Point> f;
        Implementation(Expression<Matrix<Real,Rows,Cols>,Point> f1):   
        f(f1)
        {};

        virtual Matrix<Real,3,3> evaluate(const Point& point) override
        { std::cout<<point<<std::endl;
          mat_(0,0)=point[0];            mat_(0,1)=point[0]+point[1];        mat_(0,2)=point[0]+point[2];
          mat_(1,0)=point[1]+point[0];   mat_(1,1)=point[1]+point[1];        mat_(1,2)=point[1]+point[2];
          mat_(2,0)=point[2]+point[0];   mat_(2,1)=point[2]+point[1];        mat_(2,2)=point[2]+point[2];
         return mat_; };
    protected:
    Matrix<Real,Rows,Cols> mat_;    
    };  
public: 




    ExpressionMatrix():Expression<Matrix<Real,Rows,Cols>,Point>(std::make_shared<Implementation>(Expression<Matrix<Real,Rows,Cols>,Point>())){};

    ExpressionMatrix(Expression<Matrix<Real,Rows,Cols>,Point> f)
    :   Expression<Matrix<Real,Rows,Cols>,Point>(std::make_shared<Implementation>(f))
    {}

};

template<typename Point>
class ExpressionVector: public Expression<Vector<Real,3>,Point>
{ 

public:
    static constexpr Integer Rows=3;

     protected:
    struct Implementation: Expression<Vector<Real,3>,Point>::Implementation
    {
        Expression<Vector<Real,3>,Point> f;
        Implementation(Expression<Vector<Real,3>,Point> f1):   
        f(f1)
        {};

        virtual Vector<Real,Rows>  evaluate(const Point& point) override
        { 
          vec_[0]=point[0];          
          vec_[1]=point[1];
          vec_[2]=point[2];
         return vec_; };
    protected:
    Vector<Real,Rows> vec_;    
    };  
public: 




    ExpressionVector():Expression<Vector<Real,3>,Point>(std::make_shared<Implementation>(Expression<Vector<Real,3>,Point>())){};

    ExpressionVector(Expression<Vector<Real,3>,Point> f)
    :   Expression<Vector<Real,3>,Point>(std::make_shared<Implementation>(f))
    {}

};

template<typename QuadratureRule, typename Elem,typename BaseFunctionSpace, typename QP>
class ExpressionShapeFunction: public Expression<ShapeFunctionOperator4< QuadratureRule, Elem, BaseFunctionSpace>,QP>
{ 

public:
    static constexpr Integer Rows=3;

     protected:
    struct Implementation: Expression<ShapeFunctionOperator4< QuadratureRule, Elem, BaseFunctionSpace>,QP>::Implementation
    {
        Expression<ShapeFunctionOperator4< QuadratureRule, Elem, BaseFunctionSpace>,QP> f;
        Implementation(Expression<ShapeFunctionOperator4< QuadratureRule, Elem, BaseFunctionSpace>,QP> f1):   
        f(f1)
        {};

        virtual Vector<Real,Rows>  evaluate(const QP& qp_points, const GradientOperator& grad) override
        { 
          vec_[0]=qp_points[0];          
          vec_[1]=qp_points[1];
          vec_[2]=qp_points[2];
         return vec_; };
        virtual Vector<Real,Rows>  evaluate(const QP& qp_points, const IdentityOperator& identity) override
        { 
          vec_[0]=qp_points[0];          
          vec_[1]=qp_points[1];
          vec_[2]=qp_points[2];
         return vec_; };

    protected:
    Vector<Real,Rows> vec_;    
    };  
public: 




    ExpressionShapeFunction():Expression<ShapeFunctionOperator4< QuadratureRule, Elem, BaseFunctionSpace>,QP>
    (std::make_shared<Implementation>(Expression<ShapeFunctionOperator4< QuadratureRule, Elem, BaseFunctionSpace>,QP>())){};

    ExpressionShapeFunction(Expression<ShapeFunctionOperator4< QuadratureRule, Elem, BaseFunctionSpace>,QP> f)
    :   Expression<ShapeFunctionOperator4< QuadratureRule, Elem, BaseFunctionSpace>,QP>(std::make_shared<Implementation>(f))
    {}

};




















void cents_example()
{
    using Point=Vector<Real,3>;
    Point point{0,1,2};


    ExpressionMatrix<Point> em;
    ExpressionVector<Point> vec;
    auto emm=-em;
    auto prod=emm*em;
    auto summa=emm+em;
    auto right=emm*3.0;
    auto half=0.5*emm*3.0+right*prod;
    auto again=emm-em;
    auto again2=em*vec;
    auto again3=3*(again2+vec);
    std::cout<<"eee---"<<emm(point)<<std::endl;
    std::cout<<"prod---"<<prod(point)<<std::endl;
    std::cout<<"sum---"<<summa(point)<<std::endl;
    std::cout<<"right---"<<right(point)<<std::endl;
    std::cout<<"half---"<<half(point)<<std::endl;
    std::cout<<"again---"<<again(point)<<std::endl;
    std::cout<<"vec---"<<vec(point)<<std::endl;
    std::cout<<"again2---"<<again2(point)<<std::endl;
    std::cout<<"again3---"<<again3(point)<<std::endl;

    constexpr Integer ManifoldDim=2;
    constexpr Integer Dim=2;
    using MeshT=mars::Mesh<Dim, ManifoldDim>;
    MeshT mesh;
    using Elem = typename MeshT::Elem; 
    read_mesh("../data/beam-tri.MFEM", mesh);
    constexpr Integer QPOrder=4;
    constexpr Integer NQPoints=GaussPoints<Elem,QPOrder>::NQPoints; 
    using QP=typename GaussPoints<Elem,QPOrder>::qp_points_type;
    GaussPoints<Elem,QPOrder> gauss;
    const auto& qp_points=gauss.qp_points();

    IdentityOperator identity;

    // ExpressionShape<GaussPoints<Elem,QPOrder>, Elem, Lagrange1<2>,QP> a;
    // auto a0=a(identity,qp_points);


    // //Operator::identity;
    // std::cout<<a0<<std::endl;
    // auto bb=-a; 

};

// template<typename TrialFunction, typename TestFunction>
// void Assembling()
// {
// auto trial_fe=TrialFunction::FEFamily;
// auto test_fe=TestFunction::FEFamily;
// auto trial_order=TrialFunction::Order;
// auto test_order=TestFunction::Order;
// };















void functionspaces_example2D()
{

    constexpr Integer ManifoldDim=2;
    constexpr Integer Dim=2;
    using MeshT=mars::Mesh<Dim, ManifoldDim>;
    MeshT mesh;
    using Elem = typename MeshT::Elem; 
    read_mesh("../data/beam-tri.MFEM", mesh);
 
    NodeToElem<Elem> node2elem3(mesh);
    auto node2elem=node2elem3.val();

    for(Integer nn=0; nn<node2elem.size();nn++)
        {std::cout<<"node=="<<nn<<"     ";
         for(Integer mm=0; mm<node2elem[nn].size();mm++)
          std::cout<<node2elem[nn][mm]<<" ";
         std::cout<<std::endl;
      }
    ElemEntity<Elem,0> nodes(mesh,node2elem);

    std::cout<<"node 2 elem=="<<std::endl;
    for(Integer nn=0; nn<nodes.entity_2_elem().size();nn++)
         std::cout<<nodes.entity_2_elem(nn)[0]<<" ";
    std::cout<<std::endl;
    for(Integer nn=0; nn<nodes.elem_2_entity().size();nn++)
        {
         std::cout<<"elem="<<nn<< std::endl;
         for(Integer mm=0; mm<nodes.elem_2_entity(nn).size();mm++)
            std::cout<<nodes.elem_2_entity(nn)[mm]<<" ";
        } 


    using FSspace= FunctionSpace< MeshT, Lagrange2<1>,RT0<1>>;
    FSspace FEspace(mesh);

   std::cout<<"n_dofs="<<FEspace.n_dofs()<< std::endl;
    std::cout<<std::endl;
    for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
    {
     auto &elem_id=elem_iter;
     std::cout<<std::endl<<" elem =="<<elem_iter<<std::endl;
     for(Integer nn=0;nn<FEspace.dofmap(elem_id).size();nn++)
     {
        std::cout<<FEspace.dofmap(elem_id)[nn]<< "  ";
     }
    } 
std::cout<<std::endl;
for(Integer ss=0;ss<FEspace.n_subspaces();ss++)
       std::cout<<"components of space["<<ss<<"]=="<<FEspace.components(ss)<<std::endl;
std::cout<<std::endl;

for(Integer ss=0;ss<FEspace.n_subspaces();ss++)
{
   std::cout<<"dofs of space ="<<ss<<std::endl;
   for(Integer cc=0;cc<FEspace.components(ss);cc++)
    {std::cout<<"component ="<<cc<<"   "<<std::endl;
      auto& vec=FEspace.space_dofs(ss,cc);
      for(Integer mm=0;mm<FEspace.n_dofs(ss,cc);mm++)
          std::cout<<vec[mm]<<" ";
   }
   std::cout<<std::endl;

}
 

    for(Integer mm=0;mm<FEspace.offset().size();mm++)
     {
      std::cout<<"offset space="<<mm<<std::endl;
      for(Integer nn=0;nn<FEspace.offset()[mm].size();nn++)
         {
            std::cout<< FEspace.offset()[mm][nn]<<" ";
         }
     }
     std::cout<<std::endl;


 auto spaceinf=FEspace.space_info();
 for(Integer mm=0;mm<spaceinf.size();mm++)
     {
      std::cout<<"Space=="<<mm<<std::endl;
      for(Integer nn=0;nn<spaceinf[mm].size();nn++)
        std::cout<<spaceinf[mm][nn]<<" ";
      std::cout<<std::endl;
     }
 using fespace1= ElemFunctionSpace<Simplex<2,2>, Lagrange1<1>>;
 using fespace2= ElemFunctionSpace<Simplex<2,2>, Lagrange2<2>>;
 using fespace3= ElemFunctionSpace<Simplex<2,2>, RT0<1>>;
 std::cout<<"fespace"<<fespace2::FEFamily<<std::endl;
 std::cout<<"fespace"<<fespace2::Order<<std::endl;
 BilinearFormIntegrator<4,FSspace,FSspace,GradientOperator,IdentityOperator>(FEspace,FEspace);

 }


void normals_example3D()
{
    constexpr Integer ManifoldDim=3;
    constexpr Integer Dim=3;   
    using MeshT=mars::Mesh<Dim, ManifoldDim>;
    MeshT mesh;
    using Elem = typename MeshT::Elem;   
    read_mesh("../data/beam-tet.MFEM", mesh); 
    auto sn=SignedNormal<Simplex<Dim,ManifoldDim>>(mesh);
    auto n=sn();
    auto alpha=sn.sign();
    sn.print(mesh);
}

void normals_example4D()
{
    constexpr Integer ManifoldDim=4;
    constexpr Integer Dim=4;   
    using MeshT=mars::Mesh<Dim, ManifoldDim>;
    MeshT mesh;
    using Elem = typename MeshT::Elem;   
    read_mesh("../data/pentatope_2.MFEM", mesh);
    auto sn=SignedNormal<Simplex<Dim,ManifoldDim>>(mesh);
    auto n=sn();
    auto alpha=sn.sign();
    for(Integer e=0;e<mesh.n_elements();e++)
        {
           std::cout<<"elem=="<<e<<std::endl; 
           std::cout<<"normals:"<<std::endl; 
           for(Integer f=0;f<n[e].size();f++)
               n[e][f].describe(std::cout);
           std::cout<<"signs:"<<std::endl;
           std::cout<<alpha[e]<<std::endl; 
        }    
    sn.print(mesh);

}



void functionspaces_example4D()
{

    constexpr Integer ManifoldDim=4;
    constexpr Integer Dim=4;
    constexpr Integer FEFamily=LagrangeFE;
    constexpr Integer Order=3;
    using MeshT=mars::Mesh<Dim, ManifoldDim>;
    MeshT mesh;
    using Elem = typename MeshT::Elem; 
    read_mesh("../data/pentatope_2.MFEM", mesh);
    //read_mesh("../data/cube_6.MFEM", mesh);
    //read_mesh("../data/square_2_def.MFEM", mesh);
   

    NodeToElem<Elem> node2elem3(mesh);
    auto node2elem=node2elem3.val();
    //auto vec=dofmap<ElemLagrange1<Elem>   , ElemLagrange3<Elem,1,1> >(mesh);

    // Primitive minusone(-1);
    // Primitive zero(0);
    // Composite first(1); 
    // Composite second(2); 
    // Composite third(3);
    // Composite fifth(5);  

    // first.add(zero);//std::make_shared<decltype(zero)>(zero));
    // third.add(minusone); 
    // third.add(zero); 
    // first.add(third);//std::make_shared<decltype(third)>(third));
    // second.add(first);
    // second.traverse(); 
    // fifth.add(minusone,zero,second,third);
    // std::cout<<std::endl<<"-----------------------"<<std::endl;
    // fifth.traverse();
    // std::cout<<std::endl<<"-----------------------"<<std::endl;
    // fifth[0]->traverse();
    // std::cout<<std::endl<<"-----------------------"<<std::endl;
    // fifth[1]->traverse();
    // std::cout<<std::endl<<"-----------------------"<<std::endl;
    // fifth[2]->traverse();
    // std::cout<<std::endl<<"-----------------------"<<std::endl;

    // auto fourth=second[0];
    // fourth.traverse(); 
    // cout<<endl<<" nnn = "<< *nnn<<endl;
     Leaf a0(3);
     Leaf a1(1);
     Leaf a2(2);
     Tree t1(-10);
     Tree t2(-20);
     a0.print();
     a1.print();
     a2.print();
     a2[0].print();
     t2.add(a0);
     t2.add(a1);

     t1.add(t2);
     t1.add(a0);
     t1.add(a1);
     t1.add(a2);
     
     auto ecco=t1[0][0];
     ecco.print();
    // ecco.print();
    // tree t1;
    
    // t1.add(a1);
    // auto ecc=t1[0];
    // std::cout<<"-----------------------"<<a1.val()<<"   "<<ecc->val()<<std::endl;

    // a1.print();
    // auto c1=t1[0];
    // std::cout<<"MA CI ARRIVO?"<<std::endl;
    // std::cout<<a1.val()<<" "<<a2.val()<<" "<<t3.val_vec().size()<<" "<<t3.val_vec()[0]->val()<<std::endl;
    // // constexpr auto n_spaces=2;
    // constexpr auto dofs_per_elem=DofsPerElemNums1<Elem,RT0<1>,Lagrange3<1>>::value;
    // std::vector<std::array<Integer,dofs_per_elem>> dofmap_vec;
    // std::array<std::vector<Integer>,n_spaces> offset;
    // dofmap1<RT0<1>,Lagrange3<1>>(mesh,dofmap_vec,offset);
    

    FunctionSpace< MeshT, Lagrange1<2>, Lagrange2<1> > FEspace(mesh);

    // auto eeh=FunctionSpaceSystem(FEspace);


   std::cout<<"n_dofs="<<FEspace.n_dofs()<< std::endl;
    std::cout<<std::endl;
    for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
    {
     auto &elem_id=elem_iter;
     std::cout<<"elem_id="<<elem_id<<", number of dofs=s"<<FEspace.dofmap(elem_id).size()<<std::endl;
     for(Integer nn=0;nn<FEspace.dofmap(elem_id).size();nn++)
     {
        std::cout<<FEspace.dofmap(elem_id)[nn]<<" ";
     }
     std::cout<<std::endl;
    } 
std::cout<<std::endl;
for(Integer ss=0;ss<FEspace.n_subspaces();ss++)
       std::cout<<"components of space["<<ss<<"]=="<<FEspace.components(ss)<<std::endl;
std::cout<<std::endl;

for(Integer ss=0;ss<FEspace.n_subspaces();ss++)
{
   std::cout<<"dofs of space ="<<ss<<std::endl;
   for(Integer cc=0;cc<FEspace.components(ss);cc++)
    {std::cout<<"component ="<<cc<<"   "<<std::endl;
      auto& vec=FEspace.space_dofs(ss,cc);
      for(Integer mm=0;mm<FEspace.n_dofs(ss,cc);mm++)
          std::cout<<vec[mm]<<" ";
      std::cout<<std::endl;
   }
   std::cout<<std::endl;

}
 

    for(Integer mm=0;mm<FEspace.offset().size();mm++)
     {
      std::cout<<"offset space="<<mm<<std::endl;
      for(Integer nn=0;nn<FEspace.offset()[mm].size();nn++)
         {
            std::cout<< FEspace.offset()[mm][nn]<<" ";
         }
     }
     std::cout<<std::endl;


   for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
    {std::cout<<std::endl;
      const auto size=FEspace.dofmap(1,elem_iter).size();
      std::cout<<"elem_iter="<<elem_iter<<std::endl;
      for(Integer nn=0;nn<size;nn++)
          std::cout<<FEspace.dofmap(1,elem_iter)[nn]<<" ";
      std::cout<<std::endl;
     }

 std::cout<<std::endl;
 }





void connectivity_example5D()
{

    constexpr Integer ManifoldDim=4;
    constexpr Integer Dim=4;
    using MeshT=mars::Mesh<Dim, ManifoldDim>;
    MeshT mesh;
    using Elem = typename MeshT::Elem; 
    read_mesh("../data/pentatope_2.MFEM", mesh);
    //read_mesh("../data/cube_6.MFEM", mesh);
    //read_mesh("../data/square_2_def.MFEM", mesh)    
    

    NodeToElem<Elem> node2elem3(mesh);
    const auto node2elem=node2elem3.val();
    const auto const_entities_tuple=EntitiesOfFunctionSpace<Elem,GeneralSpace,0>(mesh,node2elem);
    const auto node=std::get<0>(const_entities_tuple);
    const auto edge=std::get<1>(const_entities_tuple);
    const auto triangle=std::get<2>(const_entities_tuple);
    constexpr Integer entity_from_index=7;
    constexpr Integer entitydim_from=2;
    constexpr Integer subentitydim_from=1;
    constexpr Integer entitydim_to=1;

    Integer entity_e[3];
    Integer entity_t[3];
     ElemConnectivity<Elem,entitydim_from,subentitydim_from,entitydim_to> conn_e2t(mesh,node2elem);
     
        
        
    const auto &entity_2_elem_e=edge.entity_2_elem();
    const auto &elem2entity_e=edge.elem_2_entity();
        
    const auto &entity_2_elem_t=triangle.entity_2_elem();
    const auto &elem2entity_t=triangle.elem_2_entity();
    
    cout<<"ENTITY FROM #elems="<<triangle.size()<<endl;
    cout<<"ENTITY FROM dimension="<<entitydim_from<<endl;

        const auto & elemii=mesh.elem(entity_2_elem_t[entity_from_index][0]);
        Combinations<ManifoldDim + 1, triangle.num_of_points()>::generate(entity_2_elem_t[entity_from_index][1],entity_t);
        for(int jj=0;jj<triangle.num_of_points();jj++)
           cout<<elemii.nodes[entity_t[jj]]<<" ";    
        cout<<endl;


    cout<<"SUBENTITY FROM dimension="<<subentitydim_from<<endl;
    cout<<"ENTITY TO dimension="<<entitydim_to<<endl;
    cout<<"ENTITY FROM of interest family="<<entity_from_index<<endl;
    
 

 
    cout<<"ENTITY TO #elems="<<edge.size()<<endl;
    for(int ii=0;ii<edge.size();ii++)
    {
    
        cout<<endl;
        for(int jj=0;jj<1;jj++)
           cout<<entity_2_elem_e[ii][jj]<<" ";
        cout<<"--------ENTITY TO elem id="<<ii<<"   ";
        const auto & elemii=mesh.elem(entity_2_elem_e[ii][0]);
        Combinations<ManifoldDim + 1, edge.num_of_points()>::generate(entity_2_elem_e[ii][1],entity_e);
        for(int jj=0;jj<edge.num_of_points();jj++)
           cout<<elemii.nodes[entity_e[jj]]<<" ";    
    }   

    const auto & connection=conn_e2t.compute(triangle,entity_from_index,edge);

    cout<<endl;
    for(Integer ii=0;ii<connection.size();ii++)
       cout<<" CONNECTIONS: ="<<connection[ii]<<endl;

 

     }


}



#endif




// Base Function: f(x) = x
// class Function
// {
//     protected:
//     struct Implementation
//     {
//         virtual ~Implementation() {}
//         virtual double evaluate(double x) const { return x; }
//     };

//     public:
//     Function()
//     :   self_(std::make_shared<Implementation>())
//     {}

//     double operator () (double x) const { return self_->evaluate(x); }

//     protected:
//     Function(std::shared_ptr<Implementation> self)
//     :   self_(self)
//     {}

//     private:
//     std::shared_ptr<Implementation> self_;
// };



// class Function2: public Function
// {
//      protected:
//     struct Implementation: Function::Implementation
//     {
//         Function f;
//         Implementation(Function f):   
//         f(f)
//         {};
//         virtual double evaluate(double x) const override{return 3.333*x; }
//     };  
// public: 
//     // double operator () (double x) const  { return 2*x; };
//     // Function2(Function f)
//     // :   Function(std::make_shared<Implementation>(f))
//     // {}
//     Function2()  
//     :   Function(std::make_shared<Implementation>(Function()))
//     {};
// };

// // Unary Function: u(-f(x))
// class UnaryMinus : public Function
// {
//     protected:
//     struct Implementation : Function::Implementation
//     {
//         Function f;


//         Implementation(Function f):   
//         f(f)
//         {};

//         virtual double evaluate(double x) const override { return -f(x); }
//     };

//     public:
//     UnaryMinus(Function f)
//     :   Function(std::make_shared<Implementation>(f))
//     {}
// };

// // Binary Function: u(f(x) + g(x))
// class BinaryAdd : public Function
// {
//     protected:
//     struct Implementation : Function::Implementation
//     {
//         Function f;
//         Function g;
//         Implementation(Function f, Function g)
//         :   f(f), g(g)
//         {};

//         virtual double evaluate(double x) const override { return f(x) + g(x); }
//     };

//     public:
//     BinaryAdd(Function f, Function g)
//     :   Function(std::make_shared<Implementation>(f, g))
//     {}
// };

// // Binary Function: u(f(x) * g(x))
// class BinaryMultiply : public Function
// {
//     protected:
//     struct Implementation : Function::Implementation
//     {
//         Function f;
//         Function g;
//         Implementation(Function f, Function g)
//         :   f(f), g(g)
//         {};

//         virtual double evaluate(double x) const override { return f(x) * g(x); }
//     };

//     public:
//     BinaryMultiply(Function f, Function g)
//     :   Function(std::make_shared<Implementation>(f, g))
//     {}
// };

// // Binary Function: u(f(x) * g(x))
// class scalarMultiply : public Function
// {
//     protected:
//     struct Implementation : Function::Implementation
//     {
//         Function f;
//         double alpha_;
//         Implementation(Function f, double alpha):   
//         f(f), 
//         alpha_(alpha)
//         {};


//         virtual double evaluate(double x) const override { return f(x) * alpha_; }
//     };

//     public:
//     scalarMultiply(Function f, double alpha)
//     :   Function(std::make_shared<Implementation>(f, alpha))
//     {}

// };

// inline scalarMultiply operator * (Function f,double alpha) { return scalarMultiply(f,alpha); }
// inline scalarMultiply operator * (double alpha,Function f) { return scalarMultiply(f,alpha); }
// inline UnaryMinus operator - (Function f) { return UnaryMinus(f); }
// inline BinaryAdd operator + (Function f, Function g) { return BinaryAdd(f, g); }
// inline BinaryMultiply operator * (Function f, Function g) { return BinaryMultiply(f, g); }
