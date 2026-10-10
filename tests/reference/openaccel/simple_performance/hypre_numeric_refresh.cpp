#include "hypre_host_shim.hpp"
#include "../../../../backend/distributed/unstructured/solvers/mars_hypre_pressure_recovery.hpp"
// Inspect resource identity without exposing production test hooks.
#define private public
#include "host_gmres.hpp"
#undef private
#include "pressure_profile_fixture.hpp"
#include "../distributed_matrix/gate_problem.hpp"
// Sequential Hypre headers use MPI aliases; the wrapper's host stubs remain active.
#undef MPI_Comm
#undef MPI_COMM_WORLD
#undef MPI_INT
#undef MPI_DOUBLE
#undef MPI_MAX
#undef MPI_SUM
#undef MPI_Comm_rank
#undef MPI_Allreduce
#undef MPI_Barrier
#undef MPI_Abort
#undef MPI_Wtime
#include "../../../../backend/distributed/unstructured/fem/segregated/mars_segregated_compensated_dot.hpp"
#include "../../../../backend/distributed/unstructured/fem/segregated/mars_segregated_pressure_refinement.hpp"

using Solver = mars::fem::HypreGMRESSolver<double, int, mars::HostTestTag>;
using Matrix = Solver::Matrix;

struct RefinementNorms { double residual2=0; bool finite=true,passed=false; };

void check(bool okay, const char* message) {
    if (!okay) throw std::runtime_error(message);
}

template<class Function> void expect_failure(Function function) {
    bool rejected = false;
    try { function(); }
    catch (const std::runtime_error&) { rejected = true; }
    check(rejected, "invalid input was accepted");
}

void make_graph(Matrix& matrix, int n) {
    matrix.column_count = n + 1;
    matrix.offsets = {0};
    matrix.columns.clear();
    for (int row = 0; row < n; ++row) {
        // Reversed local numbering tests the supplied map. Invalid entries must
        // stay excluded regardless of their large numeric values.
        matrix.columns.push_back(n - 1 - row);
        if (row > 0) matrix.columns.push_back(n - row);
        if (row + 1 < n) matrix.columns.push_back(n - 2 - row);
        matrix.columns.push_back(-1);
        matrix.columns.push_back(n);
        matrix.columns.push_back(n + 7);
        matrix.offsets.push_back(matrix.columns.size());
    }
    matrix.values.resize(matrix.columns.size());
}

void fill_system(Matrix& matrix, const std::vector<HYPRE_BigInt>& map, int epoch,
                 std::vector<double>& rhs, std::vector<double>& truth) {
    const int n = matrix.numRows();
    rhs.assign(n, 0);
    truth.resize(n);
    for (int i = 0; i < n; ++i) truth[i] = 2 + std::sin(0.11 * i + 0.3 * epoch);
    for (int row = 0; row < n; ++row) {
        for (int slot = matrix.offsets[row]; slot < matrix.offsets[row + 1]; ++slot) {
            const int local = matrix.columns[slot];
            const int col = local >= 0 && local < int(map.size()) ? map[local] : -1;
            double value = 1e12;
            if (col >= 0) {
                value = col == row ? 4.0 + 0.01 * row + 0.8 * epoch
                    : (epoch == 1 || (epoch == 3 && row % 3 == 0)) ? 0.0
                    : col < row ? -0.5 - 0.1 * epoch : -1.0;
                rhs[row] += value * truth[col];
            }
            matrix.values[slot] = value;
        }
    }
}

void verify(const Matrix& matrix, const std::vector<HYPRE_BigInt>& map,
            const std::vector<double>& b, const std::vector<double>& x,
            const std::vector<double>& truth) {
    double residual2 = 0, rhs2 = 0, error = 0;
    for (int row = 0; row < matrix.numRows(); ++row) {
        double residual = b[row];
        for (int slot = matrix.offsets[row]; slot < matrix.offsets[row + 1]; ++slot) {
            const int local = matrix.columns[slot];
            if (local >= 0 && local < int(map.size()) && map[local] >= 0)
                residual -= matrix.values[slot] * x[map[local]];
        }
        residual2 += residual * residual;
        rhs2 += b[row] * b[row];
        error = std::max(error, std::abs(x[row] - truth[row]));
    }
    check(std::sqrt(residual2 / rhs2) < 2e-10, "independent true residual failed");
    check(error < 2e-9, "manufactured solution error failed");
}

void scaled_systems() {
    Matrix matrix;
    make_graph(matrix, 32);
    std::vector<HYPRE_BigInt> map(33);
    for (int i=0;i<32;++i) map[i]=31-i;
    map.back()=-1;
    for (bool cached : {false,true}) {
        Solver solver(0,300,1e-10);
        solver.setVerbose(false);
        solver.setAMGCoarseRelaxType(18);
        solver.enable_true_residual_check();
        if (cached) solver.enable_fixed_graph_updates();
        for (double scale : {1e-9,1.0,1e15}) {
            std::vector<double> b,truth,x(32,0);
            fill_system(matrix,map,0,b,truth);
            for (double& value:matrix.values) value*=scale;
            for (double& value:b) value*=scale;
            check(solver.solve(matrix,b,x,0,32,0,32,map), "scaled system rejected");
            verify(matrix,map,b,x,truth);
            check(solver.getLastFinalResidual()<1e-10, "explicit residual not reported");
            HYPRE_Real target=0,reported=0;
            HYPRE_Int restart=0,minimum=0,limit=0,iterations=0,print=0;
            const bool flex=solver.useFlexGmres_;
            (flex?HYPRE_FlexGMRESGetTol:HYPRE_GMRESGetTol)(solver.solver_,&target);
            (flex?HYPRE_FlexGMRESGetKDim:HYPRE_GMRESGetKDim)(solver.solver_,&restart);
            (flex?HYPRE_FlexGMRESGetMinIter:HYPRE_GMRESGetMinIter)(solver.solver_,&minimum);
            (flex?HYPRE_FlexGMRESGetMaxIter:HYPRE_GMRESGetMaxIter)(solver.solver_,&limit);
            (flex?HYPRE_FlexGMRESGetPrintLevel:HYPRE_GMRESGetPrintLevel)(solver.solver_,&print);
            (flex?HYPRE_FlexGMRESGetNumIterations:HYPRE_GMRESGetNumIterations)(solver.solver_,&iterations);
            (flex?HYPRE_FlexGMRESGetFinalRelativeResidualNorm:HYPRE_GMRESGetFinalRelativeResidualNorm)(solver.solver_,&reported);
            check(target==solver.tolerance_ && restart==solver.kDim_ && minimum==3 && limit==solver.maxIter_ && print==0,
                  "wrong Krylov configuration API");
            check(iterations==solver.getLastIterations() && iterations>0 && reported<1e-10,
                  "wrong Krylov result API");
        }
        std::vector<double> b(32,0),truth,x(32,0);
        check(solver.solve(matrix,b,x,0,32,0,32,map), "explicit zero RHS rejected");
        check(solver.getLastFinalResidual()==0, "explicit zero residual was stale");
        fill_system(matrix,map,0,b,truth);
        Solver stalled(0,1,1e-14,Solver::JACOBI,1);
        stalled.setVerbose(false);
        stalled.enable_true_residual_check();
        if (cached) stalled.enable_fixed_graph_updates();
        check(!stalled.solve(matrix,b,x,0,32,0,32,map), "unconverged solve accepted");
        check(stalled.getLastFinalResidual()>1e-14, "failed solve hid its true residual");
    }
}

void spmv_backend_selection() {
    Matrix matrix;
    make_graph(matrix,32);
    std::vector<HYPRE_BigInt> map(33);
    for (int i=0;i<32;++i) map[i]=31-i;
    map.back()=-1;
    for (const char* choice : {static_cast<const char*>(nullptr),"0","1"}) {
        if (choice) setenv("MARS_HYPRE_SPMV_VENDOR",choice,1);
        else unsetenv("MARS_HYPRE_SPMV_VENDOR");
        Solver solver(0,300,1e-12);
        solver.setVerbose(false);
        solver.setAMGCoarseRelaxType(18);
        solver.enable_fixed_graph_updates();
        solver.enable_true_residual_check(1e-13,1e-10);
        const int before=host_spmv_set_calls;
        for (int epoch=0;epoch<2;++epoch) {
            std::vector<double> b,truth,x(32,0);
            fill_system(matrix,map,epoch,b,truth);
            check(solver.solve(matrix,b,x,0,32,0,32,map),"SpMV selection solve rejected");
            verify(matrix,map,b,x,truth);
            check(host_spmv_set_calls==before+1,"SpMV selection missing or repeated on a cached solve");
            check(host_spmv_last_request==(choice?choice[0]-'0':0),"wrong default or explicit SpMV policy");
        }
    }
    for (const char* bad : {"","-1","2","native"}) {
        setenv("MARS_HYPRE_SPMV_VENDOR",bad,1);
        Solver solver;
        const int before=host_spmv_set_calls;
        expect_failure([&] { solver.configure_spmv(); });
        check(host_spmv_set_calls==before,"invalid SpMV option reached Hypre");
    }
    setenv("MARS_HYPRE_SPMV_VENDOR","0",1);
    Solver solver;
    const int before=host_spmv_set_calls;
    host_spmv_rank_disagreement=true;
    expect_failure([&] { solver.configure_spmv(); });
    check(!host_spmv_rank_disagreement && host_spmv_set_calls==before,
          "inconsistent rank selection reached Hypre");
    unsetenv("MARS_HYPRE_SPMV_VENDOR");
    std::cout<<"PASS: native SpMV default, vendor opt-in, cached solves and collective option checks\n";
}

void residual_workspace() {
    const int n=257;
    Matrix matrix;
    make_graph(matrix,n);
    std::vector<HYPRE_BigInt> map(n+1);
    for (int i=0;i<n;++i) map[i]=n-1-i;
    map.back()=-1;
    Solver solver(0,300,1e-12);
    solver.setVerbose(false);
    solver.setAMGCoarseRelaxType(18);
    solver.enable_fixed_graph_updates();
    solver.enable_true_residual_check(1e-13,1e-10);
    for (int epoch=0;epoch<3;++epoch) {
        std::vector<double> b,truth,x(n,0),got_b(n),got_x(n);
        fill_system(matrix,map,epoch,b,truth);
        check(solver.solve(matrix,b,x,0,n,0,n,map),"residual fixture solve failed");
        for (double rhs_scale : {0.,1e-9,1.,1e9}) {
            for (double x_scale : {0.,1.,1.-1e-6}) {
                auto rhs=b,guess=truth;
                for (double& value:rhs) value*=rhs_scale;
                for (double& value:guess) value*=rhs_scale==0?x_scale:rhs_scale*x_scale;
                HYPRE_IJVectorSetValues(solver.b_hypre_,n,solver.d_row_global_.data(),rhs.data());
                HYPRE_IJVectorSetValues(solver.x_hypre_,n,solver.d_row_global_.data(),guess.data());
                HYPRE_IJVectorAssemble(solver.b_hypre_);
                HYPRE_IJVectorAssemble(solver.x_hypre_);
                long double residual2=0,rhs2=0,scale2=0;
                for (int row=0;row<n;++row) {
                    long double residual=rhs[row],scale=std::abs(residual);
                    for (int slot=matrix.offsets[row];slot<matrix.offsets[row+1];++slot) {
                        const int local=matrix.columns[slot];
                        if (local<0 || local>=int(map.size()) || map[local]<0) continue;
                        const long double term=static_cast<long double>(matrix.values[slot])*guess[map[local]];
                        residual-=term; scale+=std::abs(term);
                    }
                    residual2+=residual*residual; scale2+=scale*scale;
                    rhs2+=static_cast<long double>(rhs[row])*rhs[row];
                }
                const double expected=std::sqrt(residual2),bnorm=std::sqrt(rhs2);
                const double roundoff=64*std::numeric_limits<double>::epsilon()*std::sqrt(scale2);
                for (double poison : {1e100,std::numeric_limits<double>::quiet_NaN()}) {
                    // Cached scratch must be overwritten, including after zero-RHS solves.
                    HYPRE_ParVectorSetConstantValues(solver.par_r_,poison);
                    const double relative=solver.true_relative_residual();
                    check(std::isfinite(relative),"poisoned residual workspace escaped");
                    check(std::abs(solver.last_absolute_residual_-expected)<=roundoff,
                          "residual disagrees with independent original CSR");
                    check(std::abs(solver.last_rhs_norm_-bnorm)<=roundoff,"RHS norm changed");
                    check(std::abs(relative-(bnorm>0?expected/bnorm:expected))
                          <=4*roundoff/(bnorm>0?bnorm:1),"relative/absolute residual convention changed");
                    HYPRE_IJVectorGetValues(solver.b_hypre_,n,solver.d_row_global_.data(),got_b.data());
                    HYPRE_IJVectorGetValues(solver.x_hypre_,n,solver.d_row_global_.data(),got_x.data());
                    check(got_b==rhs && got_x==guess,"residual calculation modified b or x");
                    check(HYPRE_GetError()==0,"residual fixture API failure");
                }
                if (epoch==0 && rhs_scale==1. && x_scale==1.-1e-6) {
                    const auto old_ij=solver.r_hypre_;
                    const auto old_par=solver.par_r_;
                    const double old_absolute=solver.last_absolute_residual_;
                    const double old_rhs=solver.last_rhs_norm_;
                    auto other_b=rhs,other_x=guess;
                    other_b[0]+=1.; other_x[1]-=.5;
                    HYPRE_ParVectorSetConstantValues(solver.par_r_,2.);
                    const auto audit=solver.audit_true_residual(other_b,other_x);
                    check(std::abs(audit.stored_norm-2*std::sqrt(double(n)))<1e-12,
                          "audit missed overwritten residual storage");
                    check(std::abs(audit.rhs_difference-1.)<1e-12
                          && std::abs(audit.solution_difference-.5)<1e-12,
                          "audit missed changed input copies");
                    for (double norm : {audit.synchronized_norm,audit.copy_matvec_norm,audit.fresh_workspace_norm})
                        check(std::abs(norm-expected)<=roundoff,"audit recomputation disagrees with CSR");
                    check(solver.r_hypre_==old_ij && solver.par_r_==old_par
                          && solver.last_absolute_residual_==old_absolute && solver.last_rhs_norm_==old_rhs,
                          "audit changed cached handles or acceptance evidence");
                    HYPRE_IJVectorGetValues(solver.b_hypre_,n,solver.d_row_global_.data(),got_b.data());
                    HYPRE_IJVectorGetValues(solver.x_hypre_,n,solver.d_row_global_.data(),got_x.data());
                    check(got_b==rhs && got_x==guess,"audit modified b or x");
                }
            }
        }
    }
    std::cout<<"PASS: residual workspace overwrite, immutable inputs and independent CSR norms\n";
}

void mixed_residual_acceptance() {
    Matrix matrix;
    make_graph(matrix,32);
    std::vector<HYPRE_BigInt> map(33);
    for (int i=0;i<32;++i) map[i]=31-i;
    map.back()=-1;
    auto norms=[&](const std::vector<double>& b,const std::vector<double>& x) {
        double r2=0,b2=0;
        for (int row=0;row<32;++row) {
            double r=b[row];
            for (int slot=matrix.offsets[row];slot<matrix.offsets[row+1];++slot) {
                const int local=matrix.columns[slot];
                if (local>=0 && local<int(map.size()) && map[local]>=0)
                    r-=matrix.values[slot]*x[map[local]];
            }
            r2+=r*r; b2+=b[row]*b[row];
        }
        return std::make_pair(std::sqrt(r2),std::sqrt(b2));
    };
    // Force a successful Hypre early stop at the supplied guess. The wrapper
    // must decide from b-Ax, not the success code or Hypre's stopping criterion.
    const char* saved_env=std::getenv("MARS_HYPRE_ABSTOL");
    const bool had_env=saved_env!=nullptr;
    const std::string saved=had_env?saved_env:"";
    setenv("MARS_HYPRE_ABSTOL","1e100",1);
    for (bool cached : {false,true}) {
        Solver solver(0,300,1e-12);
        solver.setVerbose(false);
        solver.setAMGCoarseRelaxType(18);
        if (cached) solver.enable_fixed_graph_updates();
        for (double scale : {1e-6,1.,1e6}) {
            std::vector<double> b,truth,x;
            fill_system(matrix,map,0,b,truth);
            for (double& value:matrix.values) value*=scale;
            for (double& value:b) value*=scale;
            for (double error : {1e-11,1e-9,1e-5}) {
                auto guess=[&] { x=truth; for (double& value:x) value*=1-error; };
                const bool expected=error==1e-11 || (scale==1e-6 && error==1e-9);
                guess();
                solver.enable_true_residual_check(1e-13,1e-10);
                check(solver.solve(matrix,b,x,0,32,0,32,map)==expected,"mixed acceptance mismatch");
                const auto [r,bnorm]=norms(b,x);
                check((r<=1e-13+1e-10*bnorm)==expected,"independent mixed residual mismatch");
                check(solver.getLastIterations()==0,"early-stop fixture performed iterations");
                check(solver.getLastFinalResidual()>1e-12,"fixture did not miss the Krylov target");
                check(solver.tolerance_==1e-12,"acceptance changed the Krylov target");
                guess();
                solver.enable_true_residual_check();
                check(!solver.solve(matrix,b,x,0,32,0,32,map),"strict residual mode changed");
            }
        }
        std::vector<double> b,truth,x(32,1);
        fill_system(matrix,map,0,b,truth);
        std::fill(b.begin(),b.end(),0);
        const double unit_residual=norms(b,x).first;
        solver.enable_true_residual_check(1e-13,1e-10);
        setenv("MARS_HYPRE_RESIDUAL_AUDIT","1",1);
        for (double factor : {0.,0.5,2.}) {
            std::fill(x.begin(),x.end(),factor*1e-13/unit_residual);
            const int synchronizations=host_device_synchronizations;
            check(solver.solve(matrix,b,x,0,32,0,32,map)==(factor<=1),"zero RHS mixed acceptance mismatch");
            check((norms(b,x).first<=1e-13)==(factor<=1),"independent zero RHS check failed");
            check((host_device_synchronizations>synchronizations)==(factor>1),
                  "audit synchronized an accepted solve or skipped a rejection");
        }
        unsetenv("MARS_HYPRE_RESIDUAL_AUDIT");
        expect_failure([&] { solver.enable_true_residual_check(-1.,1e-10); });
        expect_failure([&] { solver.enable_true_residual_check(1e-13,-1.); });
        expect_failure([&] { solver.enable_true_residual_check(0.,0.); });
        expect_failure([&] { solver.enable_true_residual_check(std::numeric_limits<double>::infinity(),1e-10); });
        expect_failure([&] { solver.enable_true_residual_check(1e-13,std::numeric_limits<double>::quiet_NaN()); });
    }
    if (had_env) setenv("MARS_HYPRE_ABSTOL",saved.c_str(),1);
    else unsetenv("MARS_HYPRE_ABSTOL");
    std::cout<<"PASS: mixed residual acceptance, strict mode and zero RHS (fresh/cached)\n";
}

void explicit_stopping_tolerances() {
    Matrix matrix;
    make_graph(matrix,32);
    std::vector<HYPRE_BigInt> map(33);
    for (int i=0;i<32;++i) map[i]=31-i;
    map.back()=-1;
    const char* previous=std::getenv("MARS_HYPRE_ABSTOL");
    const bool had_env=previous!=nullptr;
    const std::string saved=had_env?previous:"";
    for (bool cached : {false,true}) for (const char* environment : {"0","1e100"}) {
        setenv("MARS_HYPRE_ABSTOL",environment,1);
        for (double absolute : {0.,0.125}) {
            Solver solver(0,300,1e-6,Solver::JACOBI);
            solver.setVerbose(false);
            if (cached) solver.enable_fixed_graph_updates();
            solver.set_stopping_tolerances(1e-10,absolute);
            solver.enable_true_residual_check(absolute,1e-10,true);
            int calls=0;
            for (int epoch=0;epoch<2;++epoch) {
                std::vector<double> b,truth;
                fill_system(matrix,map,epoch,b,truth);
                const double bnorm=std::sqrt(std::inner_product(b.begin(),b.end(),b.begin(),0.));
                for (double factor : {0.5,2.}) {
                    const double initial=absolute>0?factor*absolute:1e-3*bnorm;
                    auto x=truth;
                    for (double& value:x) value*=1-initial/bnorm;
                    check(solver.solve(matrix,b,x,0,32,0,32,map),"explicit stopping target rejected");
                    ++calls;
                    double residual2=0;
                    for (int row=0;row<32;++row) {
                        double residual=b[row];
                        for (int slot=matrix.offsets[row];slot<matrix.offsets[row+1];++slot) {
                            const int local=matrix.columns[slot];
                            if (local>=0 && local<int(map.size()) && map[local]>=0)
                                residual-=matrix.values[slot]*x[map[local]];
                        }
                        residual2+=residual*residual;
                    }
                    check(std::sqrt(residual2)<=std::max(absolute,1e-10*bnorm),
                          "explicit stopping target failed original CSR residual");
                    check((solver.getLastIterations()==0)==(absolute>0 && factor<1),
                          "explicit atol did not override the environment");
                    HYPRE_Real target=0;
                    (solver.useFlexGmres_?HYPRE_FlexGMRESGetTol:HYPRE_GMRESGetTol)(solver.solver_,&target);
                    check(target==1e-10 && HYPRE_GetError()==0,"explicit relative target was not forwarded");
                    expect_failure([&] { solver.set_stopping_tolerances(1e-8,0.5); });
                    check(solver.tolerance_==1e-10 && solver.stopping_absolute_tolerance_==absolute,
                          "late stopping configuration changed the active target");
                }
            }
            check(solver.get_graph_build_count()==(cached?1:calls),"explicit target changed graph reuse");
        }
    }
    Solver invalid;
    invalid.set_stopping_tolerances(1e-8,0.125);
    const double infinity=std::numeric_limits<double>::infinity(), nan=std::numeric_limits<double>::quiet_NaN();
    for (double relative : {-1.,0.,1.,infinity,nan})
        expect_failure([&] { invalid.set_stopping_tolerances(relative,0.); });
    for (double absolute : {-1.,infinity,nan})
        expect_failure([&] { invalid.set_stopping_tolerances(1e-10,absolute); });
    check(invalid.tolerance_==1e-8 && invalid.stopping_absolute_tolerance_==0.125,
          "invalid stopping configuration changed the target");
    if (had_env) setenv("MARS_HYPRE_ABSTOL",saved.c_str(),1);
    else unsetenv("MARS_HYPRE_ABSTOL");
    std::cout<<"PASS: explicit stopping targets, environment overrides, cached refresh and invalid/late setters\n";
}

void maximum_residual_acceptance() {
    Matrix matrix;
    matrix.column_count=4; matrix.offsets={0,1,2,3,4}; matrix.columns={0,1,2,3}; matrix.values={1,1,1,1};
    const std::vector<HYPRE_BigInt> map={0,1,2,3};
    for (bool cached : {false,true}) {
        Solver solver(0,300,1e-12,Solver::JACOBI);
        solver.setVerbose(false);
        if (cached) solver.enable_fixed_graph_updates();
        // Preserve exact dyadic candidates so acceptance alone determines the verdict.
        solver.set_stopping_tolerances(1e-12,1e100);
        auto candidate=[&](double rhs,double guess,double absolute,double relative,bool maximum,bool expected) {
            std::vector<double> b(4,rhs),x(4,guess);
            if (maximum) solver.enable_true_residual_check(absolute,relative,true);
            else solver.enable_true_residual_check(absolute,relative);
            check(solver.solve(matrix,b,x,0,4,0,4,map)==expected,"maximum/additive acceptance mismatch");
            check(solver.getLastIterations()==0 && x==std::vector<double>(4,guess),
                  "acceptance fixture did not preserve its initial candidate");
            const double residual=2*std::abs(rhs-guess), bnorm=2*std::abs(rhs);
            const double limit=maximum?std::max(absolute,relative*bnorm):absolute+relative*bnorm;
            check((residual<=limit)==expected && solver.last_absolute_residual_==residual,
                  "maximum/additive verdict disagrees with exact residual");
        };
        for (bool maximum : {true,false}) {
            candidate(1.,31./32,0.125,0.0625,maximum,true);
            candidate(1.,29./32,0.125,0.0625,maximum,!maximum);
            candidate(1.,26./32,0.125,0.0625,maximum,false);
        }
        candidate(1.,29./32,0.25,0.03125,true,true);  // absolute term controls
        candidate(1.,29./32,0.03125,0.125,true,true); // relative term controls
        candidate(1.,0.75,0.25,0.03125,true,false);
        candidate(1.,0.75,0.03125,0.125,true,false);
        candidate(0.,0.,0.125,0.0625,true,true);
        candidate(0.,0.03125,0.125,0.0625,true,true);
        candidate(0.,0.125,0.125,0.0625,true,false);
    }
    std::cout<<"PASS: exact maximum versus additive acceptance, absolute/relative targets and zero RHS\n";
}

void prepared_corrections() {
    for (bool flex:{false,true}) for (bool cached:{false,true}) {
        setenv("MARS_HYPRE_FLEXGMRES",flex?"1":"0",1);
        setenv("MARS_HYPRE_ABSTOL","37",1);
        setenv("MARS_HYPRE_MINITER","0",1);
        Solver solver(0,40,1e-12,Solver::JACOBI,10);
        solver.setVerbose(false); solver.enable_true_residual_check(0.,1e-10,true);
        if (cached) solver.enable_fixed_graph_updates();
        Matrix a; a.column_count=8; a.offsets={0};
        std::vector<HYPRE_BigInt> map(8); std::iota(map.begin(),map.end(),0);
        for (int i=0;i<8;++i) { a.columns.push_back(i); a.values.push_back(2.); a.offsets.push_back(i+1); }
        std::vector<double> b(8,1),x(8,0),delta(8,0);
        check(!solver.solve(a,b,x,0,8,0,8,map),"loose backend stop escaped original residual check");
        auto controls=solver.prepared_controls(); controls.minimum=7; solver.set_prepared_controls(controls);
        check(controls.minimum==7 && controls.absolute==37,"test did not capture environment controls");
        const auto matrix=solver.parcsr_A_; const auto precond=solver.precond_;
        const int setups=solver.get_setup_count(),graphs=solver.get_graph_build_count();
        check(solver.solve_prepared_correction(b,delta,10),"prepared correction failed");
        check(solver.getLastIterations()>0 && solver.getLastIterations()<=10,"correction budget failed");
        check(solver.parcsr_A_==matrix && solver.precond_==precond
            && solver.get_setup_count()==setups && solver.get_graph_build_count()==graphs,"correction rebuilt matrix or preconditioner");
        for (double d:delta) check(std::abs(d-.5)<4*std::numeric_limits<double>::epsilon(),"wrong diagonal correction");
        check(solver.check_prepared_solution(b,delta),"corrected original equation rejected");
        const auto after=solver.prepared_controls();
        check(after.minimum==controls.minimum && after.maximum==controls.maximum
            && after.relative==controls.relative && after.absolute==controls.absolute
            && solver.residual_relative_tolerance_==1e-10 && solver.residual_absolute_tolerance_==0,
            "correction did not restore stopping and acceptance controls");
        solver.set_prepared_controls({controls.relative,controls.absolute,0,controls.maximum});
        x.assign(8,0);
        check(!solver.solve(a,b,x,0,8,0,8,map),"correction controls leaked to the next original solve");
        solver.set_prepared_controls(controls);
        expect_failure([&] { solver.solve_prepared_correction(b,delta,0); });
        expect_failure([&] { solver.solve_prepared_correction(b,delta,41); });
        auto bad=b; bad[0]=std::numeric_limits<double>::quiet_NaN();
        expect_failure([&] { solver.solve_prepared_correction(bad,delta,10); });
        const auto failed=solver.prepared_controls();
        check(failed.minimum==controls.minimum && failed.maximum==controls.maximum
            && failed.relative==controls.relative && failed.absolute==controls.absolute,"failed correction leaked controls");
    }
    unsetenv("MARS_HYPRE_ABSTOL"); unsetenv("MARS_HYPRE_MINITER"); unsetenv("MARS_HYPRE_FLEXGMRES");
    std::cout<<"PASS: prepared corrections reuse matrix/setup, enforce budget and restore GMRES/FlexGMRES controls\n";
}

void refinement_policy() {
    struct Operations {
        int mode=0,attempts=0,kept=0,restored=0,checked=0;
        double best=1;
        RefinementNorms defect() { return {best,mode!=6,false}; }
        mars::segregated::PressureCorrectionResult correct(int remaining) {
            ++attempts;
            return {mode!=3,mode==4?0:mode==5?remaining+1:1};
        }
        RefinementNorms candidate() {
            const double residual=mode==1?best:mode==2?2*best:best*.01;
            return {residual,mode!=7,residual<.1};
        }
        bool verify() { ++checked; return mode!=8; }
        void keep() { ++kept; best*=.01; }
        void restore() { ++restored; }
    };
    for (int mode=0;mode<9;++mode) {
        Operations op; op.mode=mode;
        const auto result=mars::segregated::refine_pressure(op,2,8);
        check(result.accepted==(mode==0),"refinement falsely accepted a failed/stalled/uncertain candidate");
        check(result.rounds<=3 && op.attempts<=3 && result.iterations<=8,"refinement exceeded a bound");
        check(op.restored==int(!result.accepted),"best candidate was not restored on failure");
        if (mode==1 || mode==2 || mode==7) check(op.kept==0,"nonimproving or nonfinite candidate replaced best");
        if (mode==8) check(op.checked==3 && op.kept==3,"original-equation checks were bypassed");
    }
    Operations capped;
    check(!mars::segregated::refine_pressure(capped,8,8).accepted && capped.attempts==0,"exhausted budget launched correction");
    Operations exhausted; exhausted.mode=8;
    const auto result=mars::segregated::refine_pressure(exhausted,7,8);
    check(!result.accepted && result.iterations==8 && exhausted.attempts==1,"remaining budget not enforced");
    std::cout<<"PASS: bounded correction policy rejects stalls, nonfinite values, failed verification and exhausted budgets\n";
}

void stagnation_correction() {
    setenv("MARS_HYPRE_MINITER","0",1);
    Matrix a; a.column_count=8; a.offsets={0};
    std::vector<double> b(8,0),x(8,0),defect(8),delta(8),trial(8);
    std::vector<HYPRE_BigInt> map(8); std::iota(map.begin(),map.end(),0);
    for (int i=0;i<8;++i) {
        const int k=i/2; const double scale=1.+k/16.;
        a.columns.push_back(2*k); a.columns.push_back(2*k+1);
        a.values.push_back(scale*(i%2?-1.:1.));
        a.values.push_back(scale*(i%2?1.+0x1p-24:-1.));
        a.offsets.push_back(a.columns.size()); b[i]=i%2?0.:scale;
    }
    Solver solver(0,200,1e-10,Solver::JACOBI,100);
    solver.setVerbose(false); solver.enable_fixed_graph_updates();
    solver.enable_true_residual_check(0.,1e-10,true); solver.set_stopping_tolerances(1e-10,0.);
    check(!solver.solve(a,b,x,0,8,0,8,map),"public early-stop fixture no longer reproduces rejection");
    const auto saved_a=a.values,saved_b=b;
    struct Operations {
        Solver& solver; Matrix& a; std::vector<double> &b,&x,&r,&delta,&trial;
        RefinementNorms measure(const std::vector<double>& values,bool bounded=false) {
            double r2=0,b2=0;
            for (int i=0;i<8;++i) {
                mars::segregated::CompensatedDot dot;
                for (int k=a.offsets[i];k<a.offsets[i+1];++k) dot.product(a.values[k],values[a.columns[k]]);
                dot.product(-1.,b[i]); r[i]=-dot.value();
                const double absolute=std::abs(r[i])+(bounded?dot.error_bound():0);
                r2+=(bounded?2:1)*absolute*absolute; b2+=b[i]*b[i];
            }
            return {r2,std::isfinite(r2),r2<=1e-20*b2};
        }
        auto defect() { return measure(x); }
        mars::segregated::PressureCorrectionResult correct(int remaining) {
            std::fill(delta.begin(),delta.end(),0);
            const int before=solver.getLastIterations();
            const bool accepted=solver.solve_prepared_correction(r,delta,remaining);
            return {accepted,solver.getLastIterations()-before};
        }
        auto candidate() { for (int i=0;i<8;++i) trial[i]=x[i]+delta[i]; return measure(trial); }
        bool verify() { return solver.check_prepared_solution(b,trial) && measure(trial,true).passed; }
        void keep() { x.swap(trial); }
        void restore() { solver.check_prepared_solution(b,x); }
    } op{solver,a,b,x,defect,delta,trial};
    const int before=solver.getLastIterations(),setups=solver.get_setup_count();
    const auto result=mars::segregated::refine_pressure(op,before,solver.get_max_iterations());
    check(result.accepted && result.rounds>0 && result.iterations<=200,"bounded correction did not recover public stagnation case");
    check(solver.get_setup_count()==setups && a.values==saved_a && b==saved_b,"refinement modified original equation or setup");
    // These dyadic blocks have exactly representable solutions; this oracle
    // does not rely on either residual implementation to detect a false pass.
    for (int i=0;i<8;++i) check(x[i]==(i%2?0x1p24:1.+0x1p24),"recovered solution differs from exact dyadic oracle");
    unsetenv("MARS_HYPRE_MINITER");
    std::cout<<"PASS: public Hypre early stop recovered against exact solution without tolerance or setup changes\n";
}

void fill_profile_system(const dmatrix_gate::Problem& problem,Matrix& matrix,unsigned seed,
                         std::vector<double>& rhs,std::vector<double>& truth) {
    matrix.column_count=problem.nodes+1;
    matrix.offsets={0}; matrix.columns.clear(); matrix.values.clear();
    rhs.resize(problem.nodes); truth.resize(problem.nodes);
    int all_weak=0;
    for (int row=0;row<problem.nodes;++row) {
        const int g=problem.solver_to_global[row];
        double row_sum=0,diag=0;
        for (int u:problem.neighbors[g]) {
            const double value=problem.value(g,u,0,0,seed);
            matrix.columns.push_back(problem.nodes-1-problem.solver_node[u]);
            matrix.values.push_back(value);
            row_sum+=value;
            if (u==g) diag=value;
            else if (problem.options.diffusion) {
                check(value<0 && value==problem.value(u,g,0,0,seed),"diffusion coupling is not symmetric negative");
            }
        }
        matrix.offsets.push_back(matrix.columns.size());
        if (problem.options.diffusion) check(diag>0 && row_sum>0,"diffusion matrix lost strict diagonal dominance");
        all_weak+=std::abs(row_sum)>.9*std::abs(diag);
        rhs[row]=problem.product(g,0,seed,10+seed);
        truth[row]=problem.solution(g,0,10+seed);
    }
    check(all_weak==(problem.options.diffusion?0:problem.nodes),"fixture strength-filter expectation changed");
}

void explicit_pressure_profile() {
    namespace settings=mars::fem::pressure_settings;
    // Also retain the old all-weak fixture: maxlevels alone cannot force coarsening.
    for (bool diffusion:{false,true}) for (bool flex:{false,true}) for (bool one_level:{false,true}) {
        dmatrix_gate::Options options; options.diffusion=diffusion;
        dmatrix_gate::Problem problem(1,1,options);
        const int n=problem.nodes;
        Matrix matrix;
        std::vector<HYPRE_BigInt> map(n+1);
        for (int i=0;i<n;++i) map[i]=n-1-i;
        map.back()=-1;
        // Deliberately oppose the profile; pressure must override only this instance.
        setenv("MARS_HYPRE_FLEXGMRES",flex?"0":"1",1);
        setenv("MARS_HYPRE_MINITER","3",1);
        Solver pressure(0,300,1e-12),momentum(0,300,1e-12);
        pressure.setVerbose(false); momentum.setVerbose(false);
        const auto profile=pressure_profile_fixture(flex,one_level);
        auto bad=profile; bad["coarserelax"]=9;
        expect_failure([&] { settings::configure_gpu_profile(pressure,bad); });
        bad=profile; bad["miniter"]=201;
        expect_failure([&] { settings::configure_gpu_profile(pressure,bad); });
        settings::configure_gpu_profile(pressure,profile);
        pressure.enable_fixed_graph_updates();
        for (unsigned seed:{1u,2u}) {
            std::vector<double> b,truth,x(n,0),other(n,0);
            fill_profile_system(problem,matrix,seed,b,truth);
            check(pressure.solve(matrix,b,x,0,n,0,n,map),"profile solve rejected");
            verify(matrix,map,b,x,truth);
            pressure.inspect_prepared([&](auto solver,auto amg,bool flexible) {
                const auto actual=settings::snapshot(solver,amg,flexible);
                std::ostringstream detail;
                check(flexible==flex && pressure_profile_controls_match(profile,actual,detail),"pressure profile not applied after setup");
                check(pressure_profile_levels_match(actual,one_level || !diffusion),"pressure profile hierarchy depth incorrect");
                std::cout<<"profile diffusion="<<diffusion<<" flexible="<<flex<<" maxlevels="<<profile.at("maxlevels")
                         <<" round="<<seed<<" actual_levels="<<actual.at("effective_levels")<<'\n';
                if (diffusion && !flex && !one_level && seed==1) {
                    for (const auto& [key,value]:profile) {
                        const auto name=key=="relax_down"?"effective_relax_1":key=="relax_up"?"effective_relax_2":key;
                        auto changed=actual; changed[name]=value+1;
                        std::ostringstream mismatch;
                        check(!pressure_profile_controls_match(profile,changed,mismatch)
                              && mismatch.str().find(name)!=std::string::npos,"changed setting escaped the profile check");
                    }
                    for (const auto* key:{"effective_relax_3","kdim"}) {
                        auto missing=actual; missing.erase(key);
                        std::ostringstream mismatch;
                        check(!pressure_profile_controls_match(profile,missing,mismatch),"missing setting escaped the profile check");
                    }
                    auto collapsed=actual; collapsed["effective_levels"]=1;
                    check(pressure_profile_controls_match(profile,collapsed,detail)
                          && !pressure_profile_levels_match(collapsed,false),"collapsed hierarchy was confused with controls");
                    collapsed.erase("effective_levels");
                    check(!pressure_profile_levels_match(collapsed,false),"missing hierarchy depth passed");
                }
            });
            check(momentum.solve(matrix,b,other,0,n,0,n,map),"unmodified momentum solve rejected");
            momentum.inspect_prepared([&](auto solver,auto amg,bool flexible) {
                const auto actual=settings::snapshot(solver,amg,flexible);
                check(flexible!=flex && actual.at("miniter")==3 && actual.at("maxiter")==300
                    && actual.at("kdim")==30,"pressure controls leaked into momentum");
            });
        }
        check(pressure.get_graph_build_count()==1 && pressure.get_numeric_update_count()==2,"profile cache was rebuilt");
        expect_failure([&] { settings::configure_gpu_profile(pressure,profile); });
    }
    unsetenv("MARS_HYPRE_FLEXGMRES"); setenv("MARS_HYPRE_MINITER","0",1);
    std::cout<<"PASS: pressure profiles, one/multiple levels, GMRES/FlexGMRES and cached updates\n";
}

void retained_pressure_checks() {
    namespace recovery=mars::fem::pressure_recovery;
    for(bool cached:{false,true}) {
        Matrix matrix; matrix.column_count=1; matrix.offsets={0,1}; matrix.columns={0}; matrix.values={3.};
        std::vector<HYPRE_BigInt> map{0};
        std::vector<double> b{std::nextafter(3.,INFINITY)},x(1,0),low;
        Solver solver(0,20,1e-20,Solver::BOOMERAMG,5);
        solver.setVerbose(false); solver.set_stopping_tolerances(1e-20,1e6);
        solver.enable_true_residual_check(1e-17,1e-20,true);
        if(cached) solver.enable_fixed_graph_updates();
        check(!solver.solve(matrix,b,x,0,1,0,1,map),"loose initial stop did not fail true residual");
        solver.set_prepared_controls({1e-20,1e-17,0,20});
        const int setups=solver.get_setup_count(),graphs=solver.get_graph_build_count();
        const auto result=solver.recover_prepared_expansion(x,low);
        check(result.accepted && result.rounds>0 && result.rounds<=4 && low.size()==1 && low[0]!=0,
              "wrapper did not retain a passing nonzero remainder");
        const double rounded=std::fma(3.,x[0],-b[0]);
        check(std::abs(rounded)>1e-17 && std::abs(rounded+3*low[0])<1e-17,
              "retained scalar candidate did not remove update rounding error");
        const auto restored=solver.prepared_controls();
        check(restored.relative==1e-20 && restored.absolute==1e-17 && restored.minimum==0 && restored.maximum==20
            && solver.get_setup_count()==setups && solver.get_graph_build_count()==graphs,
              "expansion changed original controls or rebuilt AMG");
        auto* a=reinterpret_cast<hypre_ParCSRMatrix*>(solver.parcsr_A_);
        auto* px=reinterpret_cast<hypre_ParVector*>(solver.par_x_);
        recovery::Vector values(px,HYPRE_MEMORY_HOST);
        recovery::Residual evaluator(a,px,HYPRE_MEMORY_HOST);
        for(double value:{0.,1e-200,1e-160,1.}) {
            recovery::checked(hypre_ParVectorSetConstantValues(values.p,value));
            const auto norm=evaluator.norm(values.p);
            check(norm.finite() && norm.lower<=value && norm.upper>=value,
                  "shared norm failed to enclose a scalar with underflowing square");
        }
    }
}

int main() {
    setenv("MARS_HYPRE_MINITER", "0", 1);
    unsetenv("MARS_HYPRE_FLEXGMRES");
    unsetenv("MARS_HYPRE_RESIDUAL_AUDIT");
    try {
        retained_pressure_checks();
        explicit_pressure_profile();
        for (const char* flexible : {"0","1"}) {
            setenv("MARS_HYPRE_FLEXGMRES",flexible,1);
            setenv("MARS_HYPRE_MINITER","3",1);
            scaled_systems();
            spmv_backend_selection();
            setenv("MARS_HYPRE_MINITER","0",1);
            mixed_residual_acceptance();
            explicit_stopping_tolerances();
            maximum_residual_acceptance();
            residual_workspace();
            std::cout<<"PASS: Krylov API and residual checks, MARS_HYPRE_FLEXGMRES="<<flexible<<'\n';
        }
        unsetenv("MARS_HYPRE_FLEXGMRES");
        prepared_corrections();
        refinement_policy();
        stagnation_correction();
        Matrix matrix;
        make_graph(matrix, 160);
        std::vector<HYPRE_BigInt> map(161);
        for (int i = 0; i < 160; ++i) map[i] = 159 - i;
        map.back() = -1;
        std::vector<double> b, truth, x(160, 0), baseline_x(160, 0);
        Solver persistent(0, 300, 1e-10), baseline(0, 300, 1e-10);
        expect_failure([&] { persistent.inspect_prepared([](auto,auto,bool){}); });
        persistent.setVerbose(false);
        baseline.setVerbose(false);
        baseline.enable_timing();
        persistent.enable_fixed_graph_updates();
        persistent.enable_timing();
        void *ij = nullptr, *parcsr = nullptr, *solver = nullptr, *amg = nullptr;
        void *rhs = nullptr, *solution = nullptr;
        const void *packed = nullptr, *slots = nullptr;
        for (int epoch = 0; epoch < 8; ++epoch) {
            fill_system(matrix, map, epoch, b, truth);
            std::fill(x.begin(), x.end(), 0.1 * epoch);
            std::fill(baseline_x.begin(), baseline_x.end(), 0.1 * epoch);
            const int barriers = host_barriers;
            const int reductions = host_reductions;
            check(persistent.solve(matrix, b, x, 0, 160, 0, 160, map), "persistent solve failed");
            bool inspected=false;
            persistent.inspect_prepared([&](auto krylov,auto preconditioner,bool flexible) {
                inspected=krylov==persistent.solver_ && preconditioner==persistent.precond_ && !flexible;
            });
            check(inspected,"pressure capture did not inspect the prepared solver");
            check(host_barriers == barriers, "fixed graph path added a barrier");
            if (epoch > 0) check(host_reductions - reductions == 9, "steady update collectives changed");
            verify(matrix, map, b, x, truth);
            check(baseline.solve(matrix, b, baseline_x, 0, 160, 0, 160, map), "baseline solve failed");
            check(host_barriers == barriers + 2, "legacy barrier behavior changed");
            verify(matrix, map, b, baseline_x, truth);
            check(baseline.get_last_timing().packing_seconds > 0, "baseline packing time absent");
            for (int i = 0; i < 160; ++i)
                check(std::abs(x[i] - baseline_x[i]) < 2e-9, "fresh/persistent field mismatch");
            if (epoch == 0) {
                ij = persistent.A_hypre_; parcsr = persistent.parcsr_A_;
                solver = persistent.solver_; amg = persistent.precond_;
                rhs = persistent.b_hypre_; solution = persistent.x_hypre_;
                packed = persistent.d_packed_values_.data(); slots = persistent.d_source_slots_.data();
            }
            check(ij == persistent.A_hypre_ && parcsr == persistent.parcsr_A_
                && solver == persistent.solver_ && amg == persistent.precond_
                && rhs == persistent.b_hypre_ && solution == persistent.x_hypre_
                && packed == persistent.d_packed_values_.data() && slots == persistent.d_source_slots_.data(),
                "fixed graph resources changed identity");
            check(persistent.get_graph_build_count() == 1, "graph unexpectedly rebuilt");
            check(persistent.get_numeric_update_count() == epoch + 1, "numeric refresh skipped");
            check(persistent.get_setup_count() == epoch + 1, "numeric AMG setup skipped");
            const auto& timing = persistent.get_last_timing();
            check(timing.prepare_seconds >= timing.packing_seconds && timing.packing_seconds > 0
                && timing.setup_seconds > 0 && timing.solve_seconds > 0 && timing.finish_seconds > 0,
                "timing scopes invalid");
        }
        check(baseline.get_graph_build_count() == 8, "baseline did not rebuild each time");

        // Value-buffer relocation alone does not invalidate a fixed graph.
        auto relocated_values = matrix.values;
        matrix.values.swap(relocated_values);
        fill_system(matrix, map, 8, b, truth);
        check(persistent.solve(matrix, b, x, 0, 160, 0, 160, map), "relocated values failed");
        verify(matrix, map, b, x, truth);
        check(persistent.get_graph_build_count() == 1, "value relocation rebuilt graph");

        // Moving graph storage is detected even when the dimensions stay fixed.
        auto relocated_columns = matrix.columns;
        matrix.columns.swap(relocated_columns);
        check(persistent.solve(matrix, b, x, 0, 160, 0, 160, map), "graph relocation failed");
        check(persistent.get_graph_build_count() == 2, "graph relocation was not detected");
        verify(matrix, map, b, x, truth);

        // In-place map edits require explicit invalidation under the API contract.
        std::swap(map[0], map[1]);
        persistent.invalidate_setup();
        fill_system(matrix, map, 9, b, truth);
        check(persistent.solve(matrix, b, x, 0, 160, 0, 160, map), "invalidated map failed");
        verify(matrix, map, b, x, truth);
        check(persistent.get_graph_build_count() == 3, "explicit invalidation did not rebuild");
        persistent.setPointBlock(2);
        check(persistent.solve(matrix, b, x, 0, 160, 0, 160, map), "configuration invalidation failed");
        check(persistent.get_graph_build_count() == 4, "configuration did not rebuild");
        verify(matrix, map, b, x, truth);

        // Frozen reuse keeps its old contract and setup count.
        persistent.setPointBlock(1);
        persistent.enable_reuse();
        check(persistent.solve(matrix, b, x, 0, 160, 0, 160, map), "frozen setup failed");
        const int frozen_setups = persistent.get_setup_count();
        for (double& value : b) value *= 2;
        for (double& value : truth) value *= 2;
        check(persistent.solve(matrix, b, x, 0, 160, 0, 160, map), "frozen solve failed");
        verify(matrix, map, b, x, truth);
        check(persistent.get_setup_count() == frozen_setups, "frozen mode rebuilt setup");
        persistent.enable_fixed_graph_updates();
        fill_system(matrix, map, 10, b, truth);
        check(persistent.solve(matrix, b, x, 0, 160, 0, 160, map), "mode transition failed");
        verify(matrix, map, b, x, truth);

        // New partitions and maps rebuild the graph and vector ranges together.
        const int graphs_before_resize = persistent.get_graph_build_count();
        make_graph(matrix, 192);
        map.resize(193);
        for (int i = 0; i < 192; ++i) map[i] = 191 - i;
        map.back() = -1;
        fill_system(matrix, map, 11, b, truth);
        x.assign(192, 0);
        check(persistent.solve(matrix, b, x, 0, 192, 0, 192, map), "partition resize failed");
        verify(matrix, map, b, x, truth);
        check(persistent.get_graph_build_count() == graphs_before_resize + 1,
              "partition resize did not rebuild");

        const double saved = matrix.values[0];
        matrix.values[0] = std::numeric_limits<double>::quiet_NaN();
        expect_failure([&] { persistent.solve(matrix, b, x, 0, 192, 0, 192, map); });
        matrix.values[0] = saved;
        const double saved_b = b[0];
        b[0] = std::numeric_limits<double>::infinity();
        expect_failure([&] { persistent.solve(matrix, b, x, 0, 192, 0, 192, map); });
        b[0] = saved_b;
        auto short_rhs = b;
        short_rhs.resize(2);
        expect_failure([&] { persistent.solve(matrix, short_rhs, x, 0, 192, 0, 192, map); });
        expect_failure([&] { persistent.solve<HYPRE_BigInt>(matrix, b, x, 0, 192, 0, 192, map); });
        persistent.setPrecondMatrix(&matrix);
        expect_failure([&] { persistent.solve(matrix, b, x, 0, 192, 0, 192, map); });
        persistent.setPrecondMatrix(nullptr);
        check(persistent.solve(matrix, b, x, 0, 192, 0, 192, map), "post-rejection rebuild failed");
        verify(matrix, map, b, x, truth);

        // The same fixed-graph path also supports the diagonal preconditioner.
        Solver jacobi(0, 300, 1e-10, Solver::JACOBI);
        jacobi.setVerbose(false);
        jacobi.enable_fixed_graph_updates();
        for (int epoch = 0; epoch < 3; ++epoch) {
            fill_system(matrix, map, epoch, b, truth);
            std::fill(x.begin(), x.end(), 0);
            check(jacobi.solve(matrix, b, x, 0, 192, 0, 192, map), "Jacobi refresh failed");
            verify(matrix, map, b, x, truth);
        }
        check(jacobi.get_setup_count() == 3 && jacobi.get_graph_build_count() == 1,
              "Jacobi lifecycle failed");
        b.assign(192, 0);
        truth.assign(192, 0);
        x.assign(192, 1);
        std::vector<double> fresh_zero_x(192, 1);
        Solver fresh_jacobi(0, 300, 1e-10, Solver::JACOBI);
        fresh_jacobi.setVerbose(false);
        const bool refreshed_zero = jacobi.solve(matrix, b, x, 0, 192, 0, 192, map);
        const bool fresh_zero = fresh_jacobi.solve(matrix, b, fresh_zero_x, 0, 192, 0, 192, map);
        check(refreshed_zero == fresh_zero, "zero RHS warm-guess acceptance changed");
        for (int i = 0; i < 192; ++i)
            check(std::abs(x[i] - fresh_zero_x[i]) < 2e-9, "zero RHS warm-guess parity failed");
        x.assign(192, 0);
        check(jacobi.solve(matrix, b, x, 0, 192, 0, 192, map), "zero RHS refresh failed");
        for (double value : x) check(std::abs(value) < 1e-9, "zero RHS retained stale solution");
        check(jacobi.getLastIterations() == 0 && jacobi.getLastFinalResidual() == 0,
              "zero-iteration result retained previous residual");
        const auto residual_scratch = jacobi.r_hypre_;
        fill_system(matrix, map, 1, b, truth);
        x = truth;
        check(jacobi.solve(matrix, b, x, 0, 192, 0, 192, map), "exact nonzero initial guess failed");
        verify(matrix, map, b, x, truth);
        check(jacobi.getLastIterations() == 0 && jacobi.r_hypre_ == residual_scratch,
              "zero-iteration residual scratch was not reused");
        Matrix empty;
        make_graph(empty, 0);
        const std::vector<HYPRE_BigInt> empty_map;
        std::vector<double> empty_rhs, empty_x;
        expect_failure([&] { jacobi.solve(empty, empty_rhs, empty_x, 0, 0, 0, 0, empty_map); });
        std::cout << "PASS: scale-independent acceptance and rejected stalled solves, refreshed values/zeros, true residuals, fresh-solve parity, stable resources, "
                     "relocation, resize, invalidation, frozen/Jacobi modes, rejected bad inputs, "
                     "timing and collective counts\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
