#include "node_export.hpp"
#include <vector>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    using mars_reference::NodeExport;
    const std::vector<double> diagonal = {16, 9, -7, 3, 32, 4, -8, 5, 48};
    {
        NodeExport writer(MPI_COMM_WORLD, "momentum.node");
        int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) writer.node(42, {
            {"density", std::vector<double>{2}}, {"volume", std::vector<double>{3}},
            {"pseudo_dt", std::vector<double>{.5}}, {"mass_divergence", std::vector<double>{-4}},
            {"velocity", std::vector<double>{1,-2,3}}, {"pressure_gradient", std::vector<double>{2,4,-1}},
            {"force", std::vector<double>{1,2,3}}, {"source", std::vector<double>{-1,1,2}},
            {"coriolis", std::vector<double>(9,0)}
        }, {{"lhs", std::vector<double>{16,0,0,0,16,0,0,0,16}}, {"rhs", std::vector<double>{-10,5,6}}});
    }
    {
        NodeExport writer(MPI_COMM_WORLD, "momentum.relaxation");
        int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) writer.node(42, {{"lhs", std::vector<double>{4,9,-7,3,8,4,-8,5,12}},
                                      {"alpha", std::vector<double>{.25}}}, {{"lhs", diagonal}});
    }
    {
        NodeExport writer(MPI_COMM_WORLD, "momentum.influence");
        int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::vector<double> row = {-4,1e8,-1e8,1e8,-8,1e8,1e8,1e8,-12};
        row.insert(row.end(), diagonal.begin(), diagonal.end());
        if (rank == 0) writer.node(42, {{"volume", std::vector<double>{8}}, {"row_blocks", row},
            {"diagonal_block", std::vector<double>{1}}, {"consistent", std::vector<double>{1}},
            {"fractional_step", std::vector<double>{0}}, {"transient", std::vector<double>{0}},
            {"small", std::vector<double>{std::numeric_limits<double>::epsilon()}}},
            {{"d", std::vector<double>{.5,.25,1.0/6}}, {"d_tilde", std::vector<double>{2.0/3,1.0/3,2.0/9}}});
    }
    {
        NodeExport writer(MPI_COMM_WORLD, "momentum.boundary_relaxation");
        int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) writer.node(42, {{"rhs", std::vector<double>{8,-4,0}}, {"factor", std::vector<double>{.75}}},
                                     {{"rhs", std::vector<double>{6,-3,0}}});
    }
    MPI_Finalize();
}
