#pragma once
#include "mfem.hpp"
#include "../bp3_phases/phases.hpp"
#include <memory>

// Preserve the memory preference of the driver's true-DOF operators.
namespace bp3_measure
{
using namespace mfem;
inline void sync_device() { MFEM_DEVICE_SYNC; }

class TimedAmg : public HypreBoomerAMG
{
public:
    using HypreBoomerAMG::HypreBoomerAMG;
    using HypreBoomerAMG::Mult;
    void Mult(const Vector& b, Vector& x) const override
    {
        bp3_phase::entry_wait();
        bp3_phase::Scope coarse(bp3_phase::Coarse, false, true);
        bp3_phase::Scope hypre(bp3_phase::Hypre, false, true);
        HypreBoomerAMG::Mult(b, x);
    }
};

class TimedOperator : public Operator
{
    std::unique_ptr<Operator> op_;

public:
    explicit TimedOperator(Operator* op)
        : Operator(op->Height(), op->Width())
        , op_(op)
    {
    }
    MemoryClass GetMemoryClass() const override { return op_->GetMemoryClass(); }
    void MultTranspose(const Vector& b, Vector& x) const override
    {
        bp3_phase::Scope phase(bp3_phase::Operator, true, true);
        op_->MultTranspose(b, x);
    }
    // Multigrid uses additive actions; retain their scratch reuse and batched implementations.
    void AddMult(const Vector& b, Vector& x, const real_t a = 1.0) const override
    {
        bp3_phase::Scope phase(bp3_phase::Operator, true, true);
        op_->AddMult(b, x, a);
    }
    void AddMultTranspose(const Vector& b, Vector& x, const real_t a = 1.0) const override
    {
        bp3_phase::Scope phase(bp3_phase::Operator, true, true);
        op_->AddMultTranspose(b, x, a);
    }
    void ArrayMult(const Array<const Vector*>& b, Array<Vector*>& x) const override
    {
        bp3_phase::Scope phase(bp3_phase::Operator, true, true);
        op_->ArrayMult(b, x);
    }
    void ArrayMultTranspose(const Array<const Vector*>& b, Array<Vector*>& x) const override
    {
        bp3_phase::Scope phase(bp3_phase::Operator, true, true);
        op_->ArrayMultTranspose(b, x);
    }
    void ArrayAddMult(const Array<const Vector*>& b, Array<Vector*>& x, const real_t a = 1.0) const override
    {
        bp3_phase::Scope phase(bp3_phase::Operator, true, true);
        op_->ArrayAddMult(b, x, a);
    }
    void ArrayAddMultTranspose(const Array<const Vector*>& b, Array<Vector*>& x, const real_t a = 1.0) const override
    {
        bp3_phase::Scope phase(bp3_phase::Operator, true, true);
        op_->ArrayAddMultTranspose(b, x, a);
    }
    void AssembleDiagonal(Vector& diagonal) const override { op_->AssembleDiagonal(diagonal); }
    void Mult(const Vector& b, Vector& x) const override
    {
        bp3_phase::Scope phase(bp3_phase::Operator, true, true);
        op_->Mult(b, x);
    }
};

class TimedSmoother : public Solver
{
    std::unique_ptr<Solver> solver_;

public:
    explicit TimedSmoother(Solver* solver)
        : Solver(solver->Height())
        , solver_(solver)
    {
    }
    MemoryClass GetMemoryClass() const override { return solver_->GetMemoryClass(); }
    void ArrayMult(const Array<const Vector*>& b, Array<Vector*>& x) const override
    {
        bp3_phase::Scope phase(bp3_phase::Smoother, true, true);
        solver_->iterative_mode = iterative_mode;
        solver_->ArrayMult(b, x);
    }
    void ArrayMultTranspose(const Array<const Vector*>& b, Array<Vector*>& x) const override
    {
        bp3_phase::Scope phase(bp3_phase::Smoother, true, true);
        solver_->iterative_mode = iterative_mode;
        solver_->ArrayMultTranspose(b, x);
    }
    void SetOperator(const Operator& op) override
    {
        solver_->SetOperator(op);
        height = solver_->Height();
        width  = solver_->Width();
    }
    void Mult(const Vector& b, Vector& x) const override
    {
        bp3_phase::Scope phase(bp3_phase::Smoother, true, true);
        solver_->iterative_mode = iterative_mode;
        solver_->Mult(b, x);
    }
    void MultTranspose(const Vector& b, Vector& x) const override
    {
        bp3_phase::Scope phase(bp3_phase::Smoother, true, true);
        solver_->iterative_mode = iterative_mode;
        solver_->MultTranspose(b, x);
    }
};

} // namespace bp3_measure
