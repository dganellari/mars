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

#include <cmath>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "mars_base.hpp"

#include <mars_pp_util.hpp>

#include "mars_gathered_vector.hpp"
#ifdef MARS_ENABLE_KOKKOS
#include "mars_utils_kokkos.hpp"
#endif

namespace mars {

#define MARS_PUBLIC_COLLECTIVES_(T)                    \
    T min(T value) const { return impl_->min(value); } \
    T max(T value) const { return impl_->max(value); } \
    T sum(T value) const { return impl_->sum(value); } \
    std::vector<T> gather(T value, int root) const { return impl_->gather(value, root); }

#define MARS_INTERFACE_COLLECTIVES_(T) \
    virtual T min(T value) const = 0;  \
    virtual T max(T value) const = 0;  \
    virtual T sum(T value) const = 0;  \
    virtual std::vector<T> gather(T value, int root) const = 0;

#define MARS_WRAP_COLLECTIVES_(T)                                \
    T min(T value) const override { return wrapped.min(value); } \
    T max(T value) const override { return wrapped.max(value); } \
    T sum(T value) const override { return wrapped.sum(value); } \
    std::vector<T> gather(T value, int root) const override { return wrapped.gather(value, root); }

#ifdef _WIN32
    // Type-list does not compile on MSVC

#define MAKE_COLLECTIVE_FUN(F_) \
    F_(float)                   \
    F_(double)                  \
    F_(int)                     \
    F_(unsigned)                \
    F_(long)                    \
    F_(unsigned long)
    //\
    // F_(long long)               \
    // F_(unsigned long long)

#define GENERATE_COLLECTIVE_PUBLIC MAKE_COLLECTIVE_FUN(MARS_PUBLIC_COLLECTIVES_)
#define GENERATE_COLLECTIVE_INTERFACE MAKE_COLLECTIVE_FUN(MARS_INTERFACE_COLLECTIVES_)
#define GENERATE_COLLECTIVE_WRAP MAKE_COLLECTIVE_FUN(MARS_WRAP_COLLECTIVES_)

#else  //_WIN32

#define MARS_COLLECTIVE_TYPES_ float, double, int, unsigned, long, unsigned long, long long, unsigned long long

#define GENERATE_COLLECTIVE_PUBLIC MARS_PP_FOREACH(MARS_PUBLIC_COLLECTIVES_, MARS_COLLECTIVE_TYPES_);
#define GENERATE_COLLECTIVE_INTERFACE MARS_PP_FOREACH(MARS_INTERFACE_COLLECTIVES_, MARS_COLLECTIVE_TYPES_)
#define GENERATE_COLLECTIVE_WRAP MARS_PP_FOREACH(MARS_WRAP_COLLECTIVES_, MARS_COLLECTIVE_TYPES_)

#endif  // _WIN32

#ifdef MARS_ENABLE_KOKKOS
#define MARS_PUBLIC_PTOP_(T)                                                                                  \
    ViewObject<T> min(ViewObject<T> value) const { return impl_->min(value); }                                \
    ViewObject<T> max(ViewObject<T> value) const { return impl_->max(value); }                                \
    void i_send_recv_all_to_all(const std::vector<T> &send_count, std::vector<T> &receive_count) {            \
        impl_->i_send_recv_all_to_all(send_count, receive_count);                                             \
    }                                                                                                         \
    void i_send_recv_vec(const std::vector<T> &send_count, std::vector<T> &receive_count) {                   \
        impl_->i_send_recv_vec(send_count, receive_count);                                                    \
    }                                                                                                         \
    void i_send_recv_view_to_all(const ViewVectorType<T> &dest,                                               \
                                 const Integer *dest_displ,                                                   \
                                 const ViewVectorType<T> &src,                                                \
                                 const Integer *src_displ) {                                                  \
        impl_->i_send_recv_view_to_all(dest, dest_displ, src, src_displ);                                     \
    }                                                                                                         \
    void i_send_recv_view(const ViewVectorType<T> &dest,                                                      \
                          const Integer *dest_displ,                                                          \
                          const ViewVectorType<T> &src,                                                       \
                          const Integer *src_displ) {                                                         \
        impl_->i_send_recv_view(dest, dest_displ, src, src_displ);                                            \
    }                                                                                                         \
    void gather_all_view(T value, const ViewVectorType<T> &buffer) { impl_->gather_all_view(value, buffer); } \
    void gather_all_view(const ViewVectorType<T> &values, ViewVectorType<T> &buffer) {                  \
        impl_->gather_all_view(values, buffer);                                                               \
    }

#define MARS_INTERFACE_PTOP_(T)                                                                                     \
    virtual ViewObject<T> min(ViewObject<T> value) const = 0;                                                       \
    virtual ViewObject<T> max(ViewObject<T> value) const = 0;                                                       \
    virtual void i_send_recv_all_to_all(const std::vector<T> &send_count, std::vector<T> &receive_count) const = 0; \
    virtual void i_send_recv_vec(const std::vector<T> &send_count, std::vector<T> &receive_count) const = 0;        \
    virtual void i_send_recv_view_to_all(const ViewVectorType<T> &dest,                                             \
                                         const Integer *dest_displ,                                                 \
                                         const ViewVectorType<T> &src,                                              \
                                         const Integer *src_displ) const = 0;                                       \
    virtual void i_send_recv_view(const ViewVectorType<T> &dest,                                                    \
                                  const Integer *dest_displ,                                                        \
                                  const ViewVectorType<T> &src,                                                     \
                                  const Integer *src_displ) const = 0;                                              \
    virtual void gather_all_view(T value, const ViewVectorType<T> &buffer) const = 0;                               \
    virtual void gather_all_view(const ViewVectorType<T> &values, ViewVectorType<T> &buffer) const = 0;

#define MARS_WRAP_PTOP_(T)                                                                                        \
    ViewObject<T> min(ViewObject<T> value) const override { return wrapped.min(value); }                          \
    ViewObject<T> max(ViewObject<T> value) const override { return wrapped.max(value); }                          \
    void i_send_recv_all_to_all(const std::vector<T> &send_count, std::vector<T> &receive_count) const override { \
        wrapped.i_send_recv_all_to_all(send_count, receive_count);                                                \
    }                                                                                                             \
    void i_send_recv_vec(const std::vector<T> &send_count, std::vector<T> &receive_count) const override {        \
        wrapped.i_send_recv_vec(send_count, receive_count);                                                       \
    }                                                                                                             \
    void i_send_recv_view_to_all(const ViewVectorType<T> &dest,                                                   \
                                 const Integer *dest_displ,                                                       \
                                 const ViewVectorType<T> &src,                                                    \
                                 const Integer *src_displ) const override {                                       \
        wrapped.i_send_recv_view_to_all(dest, dest_displ, src, src_displ);                                        \
    }                                                                                                             \
    void i_send_recv_view(const ViewVectorType<T> &dest,                                                          \
                          const Integer *dest_displ,                                                              \
                          const ViewVectorType<T> &src,                                                           \
                          const Integer *src_displ) const override {                                              \
        wrapped.i_send_recv_view(dest, dest_displ, src, src_displ);                                               \
    }                                                                                                             \
    void gather_all_view(T value, const ViewVectorType<T> &buffer) const override {                               \
        wrapped.gather_all_view(value, buffer);                                                                   \
    }                                                                                                             \
    void gather_all_view(const ViewVectorType<T> &values, ViewVectorType<T> &buffer) const override {             \
        wrapped.gather_all_view(values, buffer);                                                                  \
    }

#ifdef _WIN32
    // Type-list does not compile on MSVC

#define MAKE_PTOP_FUN(F_) \
    F_(double)            \
    F_(float)             \
    F_(int)               \
    F_(unsigned)          \
    F_(short)

#define GENERATE_PTOP_PUBLIC MAKE_PTOP_FUN(MARS_PUBLIC_PTOP_);
#define GENERATE_PTOP_INTERFACE MAKE_PTOP_FUN(MARS_INTERFACE_PTOP_)
#define GENERATE_PTOP_WRAP MAKE_PTOP_FUN(MARS_WRAP_PTOP_)

#else  //_WIN32

#define MARS_PTOP_TYPES_ double, Integer, float, int, unsigned, short, uint64_t
#define GENERATE_PTOP_PUBLIC MARS_PP_FOREACH(MARS_PUBLIC_PTOP_, MARS_PTOP_TYPES_);
#define GENERATE_PTOP_INTERFACE MARS_PP_FOREACH(MARS_INTERFACE_PTOP_, MARS_PTOP_TYPES_)
#define GENERATE_PTOP_WRAP MARS_PP_FOREACH(MARS_WRAP_PTOP_, MARS_PTOP_TYPES_)

#endif  //_WIN32

#endif

    // Defines the concept/interface for a distributed communication context.
    //
    // Uses value-semantic type erasure to define the interface, so that
    // types that implement the interface can use duck-typing, without having
    // to inherit from distributed_context.
    //
    // For the simplest example of a distributed_context implementation,
    // see local_context, which is the default context.

    class distributed_context {
    public:
        using gid_vector = std::vector<Integer>;

#ifdef MARS_ENABLE_KOKKOS
        using local_sfc = ViewVectorType<Integer>;
#endif
        // default constructor uses a local context: see below.
        distributed_context();

        template <typename Impl>
        distributed_context(Impl &&impl) : impl_(new wrap<Impl>(std::forward<Impl>(impl))) {}

        distributed_context(distributed_context &&other) = default;
        distributed_context &operator=(distributed_context &&other) = default;

        gathered_vector<Integer> gather_gids(const gid_vector &local_gids) const {
            return impl_->gather_gids(local_gids);
        }

#ifdef MARS_ENABLE_KOKKOS
        local_sfc scatter_gids(const local_sfc global, const local_sfc local) const {
            return impl_->scatter_gids(global, local);
        }

        void scatterv_gids(const local_sfc global, const local_sfc local, const std::vector<int> &counts) const {
            impl_->scatterv_gids(global, local, counts);
        }

        void broadcast(const ViewVectorType<Integer> global) const { impl_->broadcast(global); }
#endif

        int id() const { return impl_->id(); }

        int size() const { return impl_->size(); }

        void barrier() const { impl_->barrier(); }

        std::string name() const { return impl_->name(); }

        // MARS_PP_FOREACH(MARS_PUBLIC_COLLECTIVES_, MARS_COLLECTIVE_TYPES_);
        GENERATE_COLLECTIVE_PUBLIC

#ifdef MARS_ENABLE_KOKKOS
        // MARS_PP_FOREACH(MARS_PUBLIC_PTOP_, MARS_PTOP_TYPES_);
        GENERATE_PTOP_PUBLIC
#endif

        std::vector<std::string> gather(std::string value, int root) const { return impl_->gather(value, root); }

    private:
        struct interface {
            virtual int id() const = 0;
            virtual int size() const = 0;
            virtual void barrier() const = 0;
            virtual std::string name() const = 0;

            virtual gathered_vector<Integer> gather_gids(const gid_vector &local_gids) const = 0;

            // MARS_PP_FOREACH(MARS_INTERFACE_COLLECTIVES_, MARS_COLLECTIVE_TYPES_)
            GENERATE_COLLECTIVE_INTERFACE

#ifdef MARS_ENABLE_KOKKOS
            virtual local_sfc scatter_gids(const local_sfc global, const local_sfc local) const = 0;
            virtual void scatterv_gids(const local_sfc global,
                                       const local_sfc local,
                                       const std::vector<int> &counts) const = 0;
            virtual void broadcast(const ViewVectorType<Integer> global) const = 0;
            // MARS_PP_FOREACH(MARS_INTERFACE_PTOP_, MARS_PTOP_TYPES_)
            GENERATE_PTOP_INTERFACE

#endif
            virtual std::vector<std::string> gather(std::string value, int root) const = 0;

            virtual ~interface() {}
        };

        template <typename Impl>
        struct wrap : interface {
            explicit wrap(const Impl &impl) : wrapped(impl) {}
            explicit wrap(Impl &&impl) : wrapped(std::move(impl)) {}

#ifdef MARS_ENABLE_KOKKOS
            virtual local_sfc scatter_gids(const local_sfc global, const local_sfc local) const override {
                return wrapped.scatter_gids(global, local);
            }
            virtual void scatterv_gids(const local_sfc global,
                                       const local_sfc local,
                                       const std::vector<int> &counts) const override {
                wrapped.scatterv_gids(global, local, counts);
            }
            virtual void broadcast(const ViewVectorType<Integer> global) const override { wrapped.broadcast(global); }

            // MARS_PP_FOREACH(MARS_WRAP_PTOP_, MARS_PTOP_TYPES_)
            GENERATE_PTOP_WRAP

#endif
            virtual gathered_vector<Integer> gather_gids(const gid_vector &local_gids) const override {
                return wrapped.gather_gids(local_gids);
            }

            int id() const override { return wrapped.id(); }
            int size() const override { return wrapped.size(); }
            void barrier() const override { wrapped.barrier(); }
            std::string name() const override { return wrapped.name(); }

            // MARS_PP_FOREACH(MARS_WRAP_COLLECTIVES_, MARS_COLLECTIVE_TYPES_)
            GENERATE_COLLECTIVE_WRAP

            std::vector<std::string> gather(std::string value, int root) const override {
                return wrapped.gather(value, root);
            }

            Impl wrapped;
        };

        std::unique_ptr<interface> impl_;
    };

    struct local_context {
        using gid_vector = std::vector<Integer>;
#ifdef MARS_ENABLE_KOKKOS
        using local_sfc = ViewVectorType<Integer>;
        local_sfc scatter_gids(const local_sfc global, const local_sfc local) const {
            return ViewVectorType<Integer>("local_context_view", 0);
        }

        void scatterv_gids(const local_sfc global, const local_sfc local, const std::vector<int> &counts) const {}

        template <typename T>
        void i_send_recv_all_to_all(const std::vector<T> &send_count, std::vector<T> &receive_count) const {}

        template <typename T>
        void i_send_recv_vec(const std::vector<T> &send_count, std::vector<T> &receive_count) const {}

        template <typename T>
        void i_send_recv_view_to_all(const ViewVectorType<T> &dest,
                                     const Integer *dest_displ,
                                     const ViewVectorType<T> &src,
                                     const Integer *src_displ) const {}
        template <typename T>
        void i_send_recv_view(const ViewVectorType<T> &dest,
                              const Integer *dest_displ,
                              const ViewVectorType<T> &src,
                              const Integer *src_displ) const {}

        template <typename T>
        void gather_all_view(T value, const ViewVectorType<T> &buffer) const {}

        template <typename T>
        void gather_all_view(const ViewVectorType<T> &values, ViewVectorType<T> &buffer) const {}

        void broadcast(const ViewVectorType<Integer> global) const {}

        template <typename T>
        ViewObject<T> min(ViewObject<T> value) const {
            return value;
        }

        template <typename T>
        ViewObject<T> max(ViewObject<T> value) const {
            return value;
        }

#endif
        gathered_vector<Integer> gather_gids(const std::vector<Integer> &local_gids) const {
            using count_type = typename gathered_vector<Integer>::count_type;
            return gathered_vector<Integer>(std::vector<Integer>(local_gids),
                                            {0u, static_cast<count_type>(local_gids.size())});
        }

        int id() const { return 0; }

        int size() const { return 1; }

        template <typename T>
        T min(T value) const {
            return value;
        }

        template <typename T>
        T max(T value) const {
            return value;
        }

        template <typename T>
        T sum(T value) const {
            return value;
        }

        template <typename T>
        std::vector<T> gather(T value, int) const {
            return {std::move(value)};
        }

        void barrier() const {}

        std::string name() const { return "local"; }
    };

    inline distributed_context::distributed_context() : distributed_context(local_context()) {}

    using distributed_context_handle = std::shared_ptr<distributed_context>;

    inline distributed_context_handle make_local_context() { return std::make_shared<distributed_context>(); }

    distributed_context_handle make_dry_run_context(unsigned num_ranks, unsigned num_cells_per_rank);

    // MPI context creation functions only provided if built with MPI support.
    template <typename MPICommType>
    distributed_context_handle make_mpi_context(MPICommType);

}  // namespace mars
