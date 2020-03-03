/* Copyright (c) 2016, Eidgenössische Technische Hochschule Zürich and
Forschungszentrum Jülich GmbH.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */


#include <algorithm>
#include <string>
#include <vector>

#include <mars_distributed_context.hpp>
//#include <mars_threading.hpp>

namespace mars {

struct dry_run_context_impl {

    explicit dry_run_context_impl(unsigned num_ranks, unsigned num_cells_per_tile):
        num_ranks_(num_ranks), num_cells_per_tile_(num_cells_per_tile) {};

    gathered_vector<unsigned int>
    gather_gids(const std::vector<unsigned int> &local_gids) const
    {
        using count_type = typename gathered_vector<unsigned int>::count_type;

        count_type local_size = local_gids.size();

        std::vector<unsigned int> gathered_gids;
        gathered_gids.reserve(local_size * num_ranks_);

        for (count_type i = 0; i < num_ranks_; i++)
        {
            gathered_gids.insert(gathered_gids.end(), local_gids.begin(), local_gids.end());
        }

        for (count_type i = 0; i < num_ranks_; i++)
        {
            for (count_type j = i * local_size; j < (i + 1) * local_size; j++)
            {
                gathered_gids[j] += num_cells_per_tile_ * i;
            }
        }

        std::vector<count_type> partition;
        for (count_type i = 0; i <= num_ranks_; i++)
        {
            partition.push_back(static_cast<count_type>(i * local_size));
        }

        return gathered_vector<unsigned int>(std::move(gathered_gids), std::move(partition));
    }

    ViewVectorType<unsigned int>
    scatter_gids(const ViewVectorType<unsigned int> global, const ViewVectorType<unsigned int> local) const
    {
        return ViewVectorType<unsigned int>("dry_run", 0);
    }

    void
    scatterv_gids(const ViewVectorType<unsigned int> global, const ViewVectorType<unsigned int> local, 
                const std::vector<int>& counts) const
    {
    }

    int id() const { return 0; }

    int size() const { return num_ranks_; }

    template <typename T>
    T min(T value) const { return value; }

    template <typename T>
    T max(T value) const { return value; }

    template <typename T>
    T sum(T value) const { return value * num_ranks_; }

    template <typename T>
    std::vector<T> gather(T value, int) const {
        return std::vector<T>(num_ranks_, value);
    }

    void barrier() const {}

    std::string name() const { return "dryrun"; }

    unsigned num_ranks_;
    unsigned num_cells_per_tile_;
};

std::shared_ptr<distributed_context> make_dry_run_context(unsigned num_ranks, unsigned num_cells_per_tile) {
    return std::make_shared<distributed_context>(dry_run_context_impl(num_ranks, num_cells_per_tile));
}

} // namespace mars
