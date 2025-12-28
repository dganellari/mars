#ifndef MY_NODE_GRAPH_H
#define MY_NODE_GRAPH_H

#include <CRSNodeGraph.h>
#include <fstream>

namespace accel
{

class myNodeGraph : public ::linearSolver::CRSNodeGraph
{
public:
    myNodeGraph(const std::string& filename,
                const MPI_Comm comm,
                const ::linearSolver::GraphLayout layout,
                const int verbose = 0)
        : ::linearSolver::CRSNodeGraph(comm, layout), filename_(filename),
          verbose_(verbose)
    {
        this->buildGraph();
        this->modifyHost();
        this->syncToDevice();
    }

protected:
    std::string filename_;
    const int verbose_;

    void buildGraph_() override
    {
        std::ifstream f(filename_, std::ios::binary);
        if (!f)
        {
            throw std::runtime_error("ERROR, could not read graph file");
        }

        if (verbose_ > 0)
        {
            std::cout << "deserializing" << std::endl;
        }
        deserialize(f);
        if (verbose_ > 0)
        {
            std::cout << "graph file read successfully" << std::endl;
        }

        assert(row_ptr_.size() > 1);
        assert(static_cast<Index>(primary_indices_.size()) == row_ptr_.back());
        n_owned_nodes_ = static_cast<Index>(row_ptr_.size()) - 1;
        is_built_ = true;

        computeAuxiliaryData_();
        computePackInfos_();

#ifndef NDEBUG
        // sanity checks:
        assert(is_built_);
        assert(n_owned_nodes_ > 0);
        assert(n_ghost_nodes_ != ~0);
        assert(global_row_offset_ != ~0);
        assert(global_number_nodes_ != ~0);
        assert(global_number_indices_ != ~0ull);

        // a.)
        Index sum_nnz = 0;
        assert(n_owned_nodes_ == static_cast<Index>(row_nnz_owned_.size()));
        assert(n_owned_nodes_ == static_cast<Index>(row_nnz_ghost_.size()));
        for (Index i = 0; i < n_owned_nodes_; i++)
        {
            const Index row_sum_nnz = nnzOwned(i) + nnzGhost(i);
            sum_nnz += row_sum_nnz;
            assert(row_sum_nnz == row_ptr_[i + 1] - row_ptr_[i]);
        }
        assert(sum_nnz == this->nIndices());
        // b.)
        Index sum_received = 0;
        // for (const PackInfo& p : pack_infos_)
        for (size_t k = 0; k < pack_infos_.size(); ++k)
        {
            const PackInfo& p = pack_infos_[k];
            sum_received += static_cast<Index>(p.recv_idx.size());
        }
        assert(sum_received <= n_ghost_nodes_);
        if (comm_ != MPI_COMM_NULL)
        {
            for (int r = 0; r < size_; r++)
            {
                if (r == rank_ && sum_received < n_ghost_nodes_)
                {
                    std::cout
                        << "WARNING (rank=" << rank_
                        << "): receiving fewer ghost nodes than number of "
                           "ghosts in graph (receiving="
                        << sum_received << "; n_ghosts=" << n_ghost_nodes_
                        << ")" << std::endl;
                }
                MPI_Barrier(comm_);
            }
        }
        // c.)
        assert(n_owned_nodes_ ==
               static_cast<Index>(diagonal_row_offset_.size()));
        for (size_t k = 0; k < diagonal_row_offset_.size(); ++k)
        // for (const Index i : diagonal_row_offset_)
        {
            const Index i = diagonal_row_offset_[k];
            assert(i >= 0);
        }
#endif /* NDEBUG */
    }
};

} // namespace accel
#endif
