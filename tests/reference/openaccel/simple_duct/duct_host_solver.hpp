#pragma once
// Test-only linear solve policy for the host DistributedSimpleRunner on duct meshes. The gates'
// dense LU does not scale past a few thousand rows, so this gathers the owned rows by solver DOF
// onto rank 0 and runs right-preconditioned restarted GMRES with ILU(0), then scatters. It
// stops at a relative true residual of 1e-12 (the runner re-checks 1e-13 + 1e-10 |b| itself)
// and returns the same verdict on every rank. Production uses Hypre on the device instead.
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>
#include <vector>

namespace duct {

struct HostMatrix {
    std::vector<int> offsets, columns; std::vector<double> values;
    void allocate(int rows,int,int nnz) { offsets.assign(std::size_t(rows)+1,0); columns.assign(std::size_t(nnz),0); values.assign(std::size_t(nnz),0); }
    int* rowOffsetsPtr() { return offsets.data(); } const int* rowOffsetsPtr() const { return offsets.data(); }
    int* colIndicesPtr() { return columns.data(); } const int* colIndicesPtr() const { return columns.data(); }
    double* valuesPtr() { return values.data(); } const double* valuesPtr() const { return values.data(); }
};

// Scalar CSR with sorted columns and a located diagonal.
struct Csr {
    int n=0; std::vector<int> start{0}, column, diagonal; std::vector<double> value;
    void multiply(const double* x,double* y) const {
        for (int i=0;i<n;++i) { double s=0; for (int k=start[i];k<start[i+1];++k) s+=value[k]*x[column[k]]; y[i]=s; }
    }
};

// ILU(0) in place: strictly lower part holds L (unit diagonal), the rest U.
inline bool ilu0(Csr& a) {
    std::vector<int> where(a.n>0?std::size_t(a.n):0,-1);
    for (int i=0;i<a.n;++i) {
        for (int k=a.start[i];k<a.start[i+1];++k) where[a.column[k]]=k;
        for (int k=a.start[i];k<a.diagonal[i];++k) {
            const int j=a.column[k];
            const double factor=a.value[k]/=a.value[a.diagonal[j]];
            for (int m=a.diagonal[j]+1;m<a.start[j+1];++m) { const int w=where[a.column[m]]; if (w>=0) a.value[w]-=factor*a.value[m]; }
        }
        for (int k=a.start[i];k<a.start[i+1];++k) where[a.column[k]]=-1;
        const double d=a.value[a.diagonal[i]];
        if (!std::isfinite(d) || d==0) return false;
    }
    return true;
}
inline void ilu_apply(const Csr& f,const double* r,double* z) {
    for (int i=0;i<f.n;++i) { double s=r[i]; for (int k=f.start[i];k<f.diagonal[i];++k) s-=f.value[k]*z[f.column[k]]; z[i]=s; }
    for (int i=f.n-1;i>=0;--i) {
        double s=z[i]; for (int k=f.diagonal[i]+1;k<f.start[i+1];++k) s-=f.value[k]*z[f.column[k]];
        z[i]=s/f.value[f.diagonal[i]];
    }
}
inline double norm(const std::vector<double>& v) { double s=0; for (double a:v) s+=a*a; return std::sqrt(s); }

// GMRES(m), right ILU(0) preconditioning, x0=0. Returns the iteration count, or -1 when the
// true residual has not reached rtol*|b| within max_iterations.
inline int gmres(const Csr& a,const Csr& m,const std::vector<double>& b,std::vector<double>& x,double rtol,int restart,int max_iterations) {
    const int n=a.n; x.assign(std::size_t(n),0.0);
    const double target=rtol*norm(b);
    if (target==0) return 0;
    const std::size_t rows=std::size_t(n), m1=std::size_t(restart);
    std::vector<double> r(b), w(rows), z(rows);
    std::vector<std::vector<double>> v(m1+1,std::vector<double>(rows));
    std::vector<double> h((m1+1)*m1), cs(m1), sn(m1), g(m1+1), y(m1);
    int iterations=0;
    for (;;) {
        const double beta=norm(r);
        if (!std::isfinite(beta)) return -1;
        if (beta<=target) return iterations;
        if (iterations>=max_iterations) return -1;
        for (int i=0;i<n;++i) v[0][i]=r[i]/beta;
        std::fill(g.begin(),g.end(),0.0); g[0]=beta;
        int j=0;
        for (;j<restart && iterations<max_iterations;++j,++iterations) {
            ilu_apply(m,v[j].data(),z.data()); a.multiply(z.data(),w.data());
            for (int i=0;i<=j;++i) {   // modified Gram-Schmidt
                double d=0; for (int k=0;k<n;++k) d+=w[k]*v[i][k];
                h[std::size_t(i)*restart+j]=d; for (int k=0;k<n;++k) w[k]-=d*v[i][k];
            }
            const double hn=norm(w); h[std::size_t(j+1)*restart+j]=hn;
            if (hn>0) for (int k=0;k<n;++k) v[j+1][k]=w[k]/hn;
            for (int i=0;i<j;++i) {
                const double a0=h[std::size_t(i)*restart+j], a1=h[std::size_t(i+1)*restart+j];
                h[std::size_t(i)*restart+j]=cs[i]*a0+sn[i]*a1; h[std::size_t(i+1)*restart+j]=-sn[i]*a0+cs[i]*a1;
            }
            const double a0=h[std::size_t(j)*restart+j], r0=std::hypot(a0,hn);
            if (!(r0>0)) return -1;
            cs[j]=a0/r0; sn[j]=hn/r0; h[std::size_t(j)*restart+j]=r0;
            g[j+1]=-sn[j]*g[j]; g[j]*=cs[j];
            if (std::abs(g[j+1])<=.5*target || hn==0) { ++j; ++iterations; break; }
        }
        for (int i=j-1;i>=0;--i) {
            double s=g[i]; for (int k=i+1;k<j;++k) s-=h[std::size_t(i)*restart+k]*y[k];
            y[i]=s/h[std::size_t(i)*restart+i];
        }
        std::fill(w.begin(),w.end(),0.0);
        for (int i=0;i<j;++i) for (int k=0;k<n;++k) w[k]+=y[i]*v[i][k];
        ilu_apply(m,w.data(),z.data());
        for (int k=0;k<n;++k) x[k]+=z[k];
        a.multiply(x.data(),w.data());
        for (int k=0;k<n;++k) r[k]=b[k]-w[k];   // restart from the true residual
    }
}

struct SolveStats { long long solves=0, iterations=0; int last=0; };

template<int C> struct GatheredSolve {
    MPI_Comm comm;
    std::vector<double> b, x;
    SolveStats stats;
    double rtol=1e-12;
    explicit GatheredSolve(MPI_Comm c):comm(c) {}
    double* rhs(std::size_t rows) { b.assign(rows,0); return b.data(); }
    const double* rhs() const { return b.data(); }
    const double* solution() const { return x.data(); }
    std::size_t size() const { return x.size(); }
    template<class System> bool operator()(const System& s) {
        int rank=0, ranks=1; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks);
        const int rows=s.rows(); const long long total=C*(long long)s.solver_nodes();
        const int* o=s.matrix().rowOffsetsPtr(); const int* c=s.matrix().colIndicesPtr(); const double* v=s.matrix().valuesPtr();
        const auto& map=s.solver_dof_map();
        std::vector<int> length(rows>0?std::size_t(rows):0); std::vector<long long> column; std::vector<double> value;
        for (int i=0;i<rows;++i) {
            length[i]=o[i+1]-o[i];
            for (int k=o[i];k<o[i+1];++k) { column.push_back((long long)map[c[k]]); value.push_back(v[k]); }
        }
        int count=int(value.size());
        std::vector<int> counts(ranks), displs(ranks), row_counts(ranks), row_displs(ranks);
        MPI_Gather(&count,1,MPI_INT,counts.data(),1,MPI_INT,0,comm);
        MPI_Gather(&rows,1,MPI_INT,row_counts.data(),1,MPI_INT,0,comm);
        long long all=0, all_rows=0;
        for (int q=0;q<ranks;++q) { displs[q]=int(all); all+=counts[q]; row_displs[q]=int(all_rows); all_rows+=row_counts[q]; }
        const bool root=rank==0;
        std::vector<int> glength(std::size_t(root?all_rows:0)); std::vector<long long> gcolumn(std::size_t(root?all:0));
        std::vector<double> gvalue(gcolumn.size()), gb(glength.size()), gx;
        // Solver rows are contiguous per rank in rank order, so rank order is global row order.
        MPI_Gatherv(length.data(),rows,MPI_INT,glength.data(),row_counts.data(),row_displs.data(),MPI_INT,0,comm);
        MPI_Gatherv(column.data(),count,MPI_LONG_LONG,gcolumn.data(),counts.data(),displs.data(),MPI_LONG_LONG,0,comm);
        MPI_Gatherv(value.data(),count,MPI_DOUBLE,gvalue.data(),counts.data(),displs.data(),MPI_DOUBLE,0,comm);
        MPI_Gatherv(b.data(),rows,MPI_DOUBLE,gb.data(),row_counts.data(),row_displs.data(),MPI_DOUBLE,0,comm);
        int ok=1, iterations=0;
        if (root) {
            ok=all_rows==total;
            Csr a; a.n=int(all_rows); std::size_t k=0;
            std::vector<std::pair<int,double>> row;
            for (int i=0;i<a.n && ok;++i) {
                row.clear();
                for (int m=0;m<glength[i];++m,++k) {
                    ok=ok && gcolumn[k]>=0 && gcolumn[k]<total && std::isfinite(gvalue[k]);
                    row.push_back({int(gcolumn[k]),gvalue[k]});
                }
                std::sort(row.begin(),row.end());
                int diagonal=-1;
                for (std::size_t m=0;m<row.size();++m) {
                    if (m && row[m].first==row[m-1].first) { a.value.back()+=row[m].second; continue; }
                    if (row[m].first==i) diagonal=int(a.column.size());
                    a.column.push_back(row[m].first); a.value.push_back(row[m].second);
                }
                ok=ok && diagonal>=0; a.diagonal.push_back(diagonal); a.start.push_back(int(a.column.size()));
            }
            Csr m=a;
            ok=ok && ilu0(m);
            if (ok) { iterations=gmres(a,m,gb,gx,rtol,std::min(a.n,60),20000); ok=iterations>=0; }
        }
        MPI_Bcast(&ok,1,MPI_INT,0,comm); MPI_Bcast(&iterations,1,MPI_INT,0,comm);
        x.assign(std::size_t(rows),0);
        if (ok) MPI_Scatterv(gx.data(),row_counts.data(),row_displs.data(),MPI_DOUBLE,x.data(),rows,MPI_DOUBLE,0,comm);
        ++stats.solves; stats.iterations+=iterations; stats.last=iterations;
        return ok!=0;
    }
};

} // namespace duct
