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

#pragma once

#include <assert.h>
#include <algorithm>
#include <cstddef>
#include <vector>
#include "mars_globals.hpp"

namespace mars {

    template <typename T>
    class gathered_vector {
    public:
        using value_type = T;
        using count_type = Integer;

        gathered_vector(std::vector<value_type> &&v, std::vector<count_type> &&p)
            : values_(std::move(v)), partition_(std::move(p)) {
            assert(std::is_sorted(partition_.begin(), partition_.end()));
            assert(partition_.back() == values_.size());
        }

        /// the partition of distribution
        const std::vector<count_type> &partition() const { return partition_; }

        /// the number of entries in the gathered vector in partition i
        count_type count(std::size_t i) const { return partition_[i + 1] - partition_[i]; }

        /// the values in the gathered vector
        const std::vector<value_type> &values() const { return values_; }

        /// the size of the gathered vector
        std::size_t size() const { return values_.size(); }

    private:
        std::vector<value_type> values_;
        std::vector<count_type> partition_;
    };

}  // namespace mars
