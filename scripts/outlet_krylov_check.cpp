// GPT/Codex, 2026-09-10. Runs the production FGMRES controller with independent dense operations.
#include "../backend/distributed/unstructured/fem/mars_outlet_fgmres.hpp"
#include <iostream>
#include <numeric>
#include <stdexcept>

struct HostOps
{
    int n, depth, calls = 0, reports = 0;
    bool varying = false, reject = false, corrupt_check = false;
    std::vector<double> matrix, storage;
    HostOps(int size, int restart) : n(size), depth(restart), matrix(n*n), storage(n*(2*depth+5)) {}
    double* solution() { return storage.data(); }
    double* rhs() { return solution()+n; }
    double* residual() { return solution()+2*n; }
    double* work() { return solution()+3*n; }
    double* basis(int i) { return solution()+(4+i)*n; }
    double* direction(int i) { return solution()+(5+depth+i)*n; }
    void zero(double* v) { std::fill(v,v+n,0); }
    void copy(const double* a, double* b) { std::copy(a,a+n,b); }
    void scale(double a, double* v) { for (int i=0;i<n;++i) v[i]*=a; }
    void axpy(double a, const double* x, double* y) { for (int i=0;i<n;++i) y[i]+=a*x[i]; }
    double norm(const double* v) { return std::sqrt(std::inner_product(v,v+n,v,0.0)); }
    void apply(const double* x, double* y)
    {
        for (int i=0;i<n;++i)
        {
            y[i]=0;
            for (int j=0;j<n;++j) y[i]+=matrix[i*n+j]*x[j];
        }
        if (corrupt_check && x==solution()) y[0]+=1;
    }
    bool precondition(const double* v, double* z)
    {
        ++calls;
        for (int i=0;i<n;++i) z[i]=v[i]*(varying ? 0.5+((calls+2*i)%7) : 1.0);
        return !reject;
    }
    void orthogonalize(double* w, int count, double* coefficients)
    {
        std::fill(coefficients,coefficients+count,0.0);
        for (int pass=0;pass<2;++pass)
        {
            std::vector<double> partial(count);
            for (int j=0;j<count;++j) partial[j]=std::inner_product(basis(j),basis(j)+n,w,0.0);
            for (int j=0;j<count;++j) { coefficients[j]+=partial[j]; axpy(-partial[j],basis(j),w); }
        }
    }
    void report(int, double) { ++reports; }
};

int main()
{
    int checks=0;
    auto require=[&](bool ok, const char* message) {
        ++checks;
        if (!ok) throw std::runtime_error(message);
    };
    auto solve=[](HostOps& ops, int cap=100, double tol=1e-12) {
        return outlet_fgmres(ops,ops.depth,cap,tol,std::numeric_limits<double>::epsilon());
    };
    for (bool varying : {false,true})
        for (double units : {1e-100,1.0,1e100})
        {
            HostOps ops(4,4);
            ops.varying=varying;
            // Nonsymmetric, with unequal row scales (the device path uses volume row scaling).
            ops.matrix={4,3,-1,0, 0,2,1,1, 1,0,3,-2, 0,0,1,2};
            const double exact[4]={units,-2*units,.3*units,4*units};
            ops.apply(exact,ops.rhs());
            auto result=solve(ops);
            require(result.converged && result.iterations<=4,"dense variable-preconditioner solve failed");
            for (int i=0;i<4;++i) require(std::abs((ops.solution()[i]-exact[i])/units)<1e-11,"wrong Z-basis reconstruction");
            ops.apply(ops.solution(),ops.work()); ops.axpy(-1,ops.rhs(),ops.work());
            require(ops.norm(ops.work())/ops.norm(ops.rhs())<=1e-12,"reported convergence without true residual");
        }
    {
        HostOps ops(2,4); ops.matrix={1,0,0,1}; ops.rhs()[0]=1; ops.rhs()[1]=-3;
        auto result=solve(ops,1);
        require(result.converged && result.iterations==1,"happy breakdown or cap shorter than restart failed");
    }
    {
        HostOps ops(2,2); ops.matrix={1,0,0,1};
        auto result=solve(ops);
        require(result.converged && result.iterations==0 && ops.calls==0,"zero RHS invokes preconditioner");
        ops.rhs()[0]=1; ops.reject=true;
        require(!solve(ops).converged,"failed preconditioner accepted");
    }
    {
        HostOps ops(2,2); ops.rhs()[0]=1;
        require(!solve(ops).converged,"singular zero operator accepted");
        ops.matrix={1,0,0,1}; ops.corrupt_check=true;
        require(!solve(ops).converged,"Hessenberg residual bypassed true residual check");
        ops.corrupt_check=false; ops.matrix[0]=INFINITY;
        require(!solve(ops).converged,"nonfinite action accepted");
    }
    {
        HostOps ops(5,2);
        for (int i=0;i<5;++i) { ops.matrix[i*5+i]=2+.1*i; ops.rhs()[i]=1+i; }
        auto limited=solve(ops,3,1e-15);
        require(!limited.converged && limited.iterations==3,"partial final restart exceeds or truncates iteration cap");
        auto result=solve(ops,100);
        require(result.converged && result.iterations>2 && ops.reports>1,"restart failed");
    }
    {
        // Positive eigenvalues do not make a single preconditioned direction a reliable descent method.
        HostOps ops(2,2); ops.matrix={1,10,0,2}; ops.rhs()[1]=1;
        double r[2]={0,1}, ar[2];
        for (int i=0;i<20;++i)
        {
            ops.apply(r,ar);
            const double numerator=r[0]*ar[0]+r[1]*ar[1];
            const double denominator=ar[0]*ar[0]+ar[1]*ar[1];
            const double omega=std::max(0.0,numerator/denominator);
            ops.axpy(-omega,ar,r);
        }
        require(ops.norm(r)>.9,"stagnation fixture does not stagnate");
        require(solve(ops).converged,"FGMRES failed stagnation fixture");
    }
    std::cout << "PASS: " << checks << " production FGMRES controller checks\n";
}
