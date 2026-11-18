#include "phase_operators.hpp"
#include <iostream>

using namespace mfem;
void require(bool ok)
{
    if (!ok)
    {
        std::cerr << "Wrapper contract failed\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}
class TestOperator : public Operator
{
public:
    mutable int additive = 0, additive_transpose = 0, batches = 0;
    TestOperator()
        : Operator(2)
    {
    }
    MemoryClass GetMemoryClass() const override { return MemoryClass::DEVICE; }
    void Mult(const Vector& b, Vector& x) const override
    {
        x = b;
        x *= 2;
    }
    void MultTranspose(const Vector& b, Vector& x) const override
    {
        x = b;
        x *= 3;
    }
    void AddMult(const Vector& b, Vector& x, const real_t a = 1.0) const override
    {
        ++additive;
        x.Add(2 * a, b);
    }
    void AddMultTranspose(const Vector& b, Vector& x, const real_t a = 1.0) const override
    {
        ++additive_transpose;
        x.Add(3 * a, b);
    }
    void ArrayMult(const Array<const Vector*>& b, Array<Vector*>& x) const override
    {
        ++batches;
        Operator::ArrayMult(b, x);
    }
    void ArrayMultTranspose(const Array<const Vector*>& b, Array<Vector*>& x) const override
    {
        ++batches;
        Operator::ArrayMultTranspose(b, x);
    }
    void ArrayAddMult(const Array<const Vector*>& b, Array<Vector*>& x, const real_t a = 1.0) const override
    {
        ++batches;
        Operator::ArrayAddMult(b, x, a);
    }
    void ArrayAddMultTranspose(const Array<const Vector*>& b, Array<Vector*>& x, const real_t a = 1.0) const override
    {
        ++batches;
        Operator::ArrayAddMultTranspose(b, x, a);
    }
    void AssembleDiagonal(Vector& d) const override { d = 2; }
};
class TestSmoother : public Solver
{
public:
    TestSmoother()
        : Solver(2)
    {
    }
    MemoryClass GetMemoryClass() const override { return MemoryClass::DEVICE; }
    void SetOperator(const Operator& op) override
    {
        height = op.Height();
        width  = op.Width();
    }
    void Mult(const Vector& b, Vector& x) const override
    {
        if (iterative_mode)
            x += b;
        else
            x = b;
    }
    void MultTranspose(const Vector& b, Vector& x) const override
    {
        Mult(b, x);
        x *= 3;
    }
};
int main(int argc, char** argv)
{
    Mpi::Init(argc, argv);
    auto* inner = new TestOperator;
    bp3_measure::TimedOperator op(inner);
    bp3_measure::TimedSmoother smoother(new TestSmoother);
    require(op.GetMemoryClass() == MemoryClass::DEVICE && smoother.GetMemoryClass() == MemoryClass::DEVICE);
    Vector b(2), x(2);
    b = 1.;
    x = 0.;
    op.Mult(b, x);
    require(x(0) == 2.);
    op.MultTranspose(b, x);
    require(x(1) == 3.);
    op.AssembleDiagonal(x);
    require(x(0) == 2.);
    x = 4.;
    op.AddMult(b, x, -.5);
    require(x(0) == 3. && inner->additive == 1);
    op.AddMultTranspose(b, x, 2.);
    require(x(0) == 9. && inner->additive_transpose == 1);
    Vector b2(2), x2(2);
    b2 = 2.;
    Array<const Vector*> inputs(2);
    Array<Vector*> outputs(2);
    inputs[0] = &b; inputs[1] = &b2;
    outputs[0] = &x; outputs[1] = &x2;
    op.ArrayMult(inputs, outputs);
    require(x(0) == 2. && x2(0) == 4. && inner->batches == 1);
    op.ArrayMultTranspose(inputs, outputs);
    require(x(0) == 3. && x2(0) == 6. && inner->batches == 2);
    op.ArrayAddMult(inputs, outputs, -.5);
    require(x(0) == 2. && x2(0) == 4. && inner->batches == 3 && inner->additive == 3);
    op.ArrayAddMultTranspose(inputs, outputs, 2.);
    require(x(0) == 8. && x2(0) == 16. && inner->batches == 4 && inner->additive_transpose == 3);
    smoother.iterative_mode = false;
    smoother.Mult(b, x);
    require(x(0) == 1.);
    smoother.iterative_mode = true;
    smoother.Mult(b, x);
    require(x(1) == 2.);
    smoother.MultTranspose(b, x);
    require(x(0) == 9.);
    smoother.ArrayMult(inputs, outputs);
    require(x(0) == 10. && x2(0) == 18.);
    smoother.ArrayMultTranspose(inputs, outputs);
    require(x(0) == 33. && x2(0) == 60.);
    smoother.SetOperator(op);
    require(smoother.Height() == 2 && smoother.Width() == 2);
    std::cout << "PASS memory class, additive and batched dispatch, transpose and iterative-mode contracts\n";
}
