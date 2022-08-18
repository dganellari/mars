#ifndef GENERATION_MARS_UTILS_KOKKOS_HPP_
#define GENERATION_MARS_UTILS_KOKKOS_HPP_

#include "mars_config.hpp"

#ifdef WITH_KOKKOS_KERNELS
#include "KokkosKernels_default_types.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#endif

#ifdef WITH_KOKKOS
#include "Kokkos_Layout.hpp"
#include "Kokkos_UnorderedMap.hpp"
#endif
#include "mars_err.hpp"
#include "mars_globals.hpp"

namespace mars {

    // return size of an array as a compile-time constant.
    template <typename T, std::size_t N, std::size_t M>
    constexpr std::size_t arraySize(T (&)[N][M]) noexcept {
        return N;
    }

// Defined for host operations that need to be performed always on the host
#ifdef KOKKOS_ENABLE_OPENMP
    using KokkosHostSpace = Kokkos::HostSpace;
    using KokkosHostExecSpace = Kokkos::OpenMP;
#else  // Serial
    using KokkosHostSpace = Kokkos::HostSpace;
    using KokkosHostExecSpace = Kokkos::Serial;
#endif

// Default Kokkos Exec space.
#if defined(MARS_USE_CUDA)
#if defined(MARS_USE_CUDAUVM)
    using KokkosSpace = Kokkos::CudaUVMSpace;
#else
    using KokkosSpace = Kokkos::CudaSpace;
#endif  // MARS_USE_CUDAUVM
    using KokkosLayout = Kokkos::LayoutLeft;
#define MARS_LAMBDA_REF [&] __device__;
#else  // MARS_USE_CUDA
#ifdef KOKKOS_ENABLE_OPENMP
    using KokkosSpace = Kokkos::HostSpace;
    using KokkosLayout = Kokkos::LayoutRight;
#else   // Serial
    using KokkosSpace = Kokkos::HostSpace;
    using KokkosLayout = Kokkos::LayoutRight;
#endif  // KOKKOS_ENABLE_OPENMP
#define MARS_LAMBDA_REF [&]
#endif  // MARS_USE_CUDA

    template <typename T>
    using ViewVectorTypeStride = Kokkos::View<T*, Kokkos::LayoutStride, KokkosSpace>;

    template <typename T>
    using ViewVectorHost = Kokkos::View<T*, Kokkos::LayoutRight, KokkosHostSpace>;

    template <typename T>
    using ViewVectorType = Kokkos::View<T*, KokkosLayout, KokkosSpace>;

    template <typename T>
    using ViewMatrixTypeLeft = Kokkos::View<T**, Kokkos::LayoutLeft, KokkosSpace>;

    template <typename T>
    using ViewMatrixType = Kokkos::View<T**, KokkosLayout, KokkosSpace>;

    template <typename T, Integer yDim_>
    using ViewMatrixTypeRC = Kokkos::View<T* [yDim_], KokkosLayout, KokkosSpace>;

    template <typename T, Integer YDim_>
    using ViewMatrixTexture =
        Kokkos::View<T* [YDim_], KokkosLayout, KokkosSpace, Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

    template <typename T, Integer XDim_, Integer YDim_>
    using ViewMatrixTextureC =
        Kokkos::View<T[XDim_][YDim_], KokkosLayout, KokkosSpace, Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

    template <typename T, Integer XDim_>
    using ViewVectorTextureC =
        Kokkos::View<T[XDim_], KokkosLayout, KokkosSpace, Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

    template <typename T>
    using ViewVectorTexture = Kokkos::View<T*, KokkosLayout, KokkosSpace, Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

    template <typename T, Integer XDim_>
    using ViewVectorTypeC = Kokkos::View<T[XDim_], KokkosLayout, KokkosSpace>;

    template <typename T>
    using ViewVectorTypeU = Kokkos::View<T*, KokkosSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    template <typename T>
    using ViewObject = Kokkos::View<T[1], KokkosSpace>;

    template <typename T, class space>
    using ViewObjectU = Kokkos::View<T[1], space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    template <typename... T>
    using ViewsTuple = std::tuple<ViewVectorType<T>...>;
    /*

    template<typename T>
    using ViewObjectUH = Kokkos::View<T[1], KokkosHostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    */

    template <typename Key, typename Value>
    using UnorderedMap = Kokkos::UnorderedMap<Key, Value, KokkosSpace>;

    template <typename T, Integer YDim_>
    struct IndexView {
        ViewMatrixTexture<T, YDim_> view;
        int index;

        KOKKOS_INLINE_FUNCTION
        IndexView(ViewMatrixTexture<T, YDim_> v, int idx) : view(v), index(idx) {}

        KOKKOS_INLINE_FUNCTION
        T& operator[](int i) { return view(index, i); }
    };

    class KokkosImplementation {
        std::string name = "kokkos";
    };

    class DistributedImplementation {
        std::string name = "distributed";
    };

    // copy matrix from host data to the host mirror view and then deep copy to the device texture view.
    template <typename T, Integer xDim_, Integer yDim_>
    void copy_matrix_from_host(std::vector<std::vector<T>> hostData,
                               ViewMatrixTextureC<T, xDim_, yDim_> map_side_to_nodes,
                               const int xDim,
                               const int yDim) {
        using namespace Kokkos;

        typename ViewMatrixTextureC<T, xDim_, yDim_>::HostMirror h_view = create_mirror_view(map_side_to_nodes);

        parallel_for(
            MDRangePolicy<Rank<2>, KokkosHostExecSpace>({0, 0}, {xDim, yDim}),
            KOKKOS_LAMBDA(int i, int j) { h_view(i, j) = hostData[i][j]; });

        Kokkos::deep_copy(map_side_to_nodes, h_view);
    }

    /*template<Integer Dim, Integer ManifoldDim, class Point_>
    void remove_extra_nodes(Mesh<Dim, ManifoldDim, Point_>& mesh,
                    std::vector<Vector<Real, Dim> >& np, const std::vector<bool>& active) {

            int count = 0;
            for (unsigned int i = 0; i < active.size(); ++i) {
                    if (active[i]) {
                            np[count] = mesh.point(i);
                            ++count;
                    }

            }

            mesh.setPoints(move(np));

    }*/

    template <typename T>
    void compact_scan(const ViewVectorType<bool>& pred, const ViewVectorType<T>& scan, const ViewVectorType<T>& out) {
        using namespace Kokkos;
        const Integer size = pred.extent(0);

        parallel_for(
            "compact scan", size, KOKKOS_LAMBDA(const Integer i) {
                if (pred(i) == 1) {
                    out(scan(i)) = i;
                }
            });
    }

    template <typename T>
    void inclusive_scan(const Integer start, const Integer end, ViewVectorType<T> in_) {
        using namespace Kokkos;

        parallel_scan(
            RangePolicy<>(start, end), KOKKOS_LAMBDA(const int& i, T& upd, const bool& final) {
                // Load old value in case we update it before accumulating
                const T val_i = in_(i);

                upd += val_i;

                if (final) {
                    in_(i) = upd;  // only update array on final pass
                }
            });
    }

    template <typename T>
    void inclusive_scan(const Integer start, const Integer end, const ViewVectorType<T> in_, ViewVectorType<T> out_) {
        using namespace Kokkos;

        parallel_scan(
            RangePolicy<>(start, end), KOKKOS_LAMBDA(const int& i, T& upd, const bool& final) {
                // Load old value in case we update it before accumulating
                const T val_i = in_(i);

                upd += val_i;

                if (final) {
                    out_(i) = upd;  // only update array on final pass
                }
            });
    }

    /* template<typename T, typename U>
    void incl_excl_scan_strided(const Integer start, const Integer end,
                            const U in_, T out_)
    {
            using namespace Kokkos;

            parallel_scan (RangePolicy<>(start , end ),	KOKKOS_LAMBDA (const int& i,
                                    Integer& upd, const bool& final)
            {
                    // Load old value in case we update it before accumulating
                    const Integer val_i = in_(i);

                    upd += val_i;

                    if (final)
                    {
                            out_(i+1) = upd; // To have both ex and inclusive in the same output.
                    }
            });
    }
     */
    // works for all cases including those with strided access views coming from auto = subview...
    template <typename H, typename U>
    void incl_excl_scan(const Integer start, const Integer end, const U in_, H out_) {
        using namespace Kokkos;

        using T = typename H::value_type;

        parallel_scan(
            RangePolicy<>(start, end), KOKKOS_LAMBDA(const int& i, T& upd, const bool& final) {
                // Load old value in case we update it before accumulating
                const T val_i = in_(i);

                upd += val_i;

                if (final) {
                    out_(i + 1) = upd;  // To have both ex and inclusive in the same output.
                }
            });
    }

    // not a very performant scan since its access is not coalesced. However OK for small arrays.
    template <typename T, typename U>
    void column_scan(const Integer end, const Integer col_idx, const ViewMatrixType<U>& in_, ViewVectorType<T>& out_) {
        using namespace Kokkos;

        parallel_scan(
            RangePolicy<>(0, end), KOKKOS_LAMBDA(const int& i, T& upd, const bool& final) {
                // Load old value in case we update it before accumulating
                const T val_i = in_(i, col_idx);

                upd += val_i;

                if (final) {
                    out_(i + 1) = upd;  // To have both ex and inclusive in the same output.
                }
            });
    }

    // not a very performant scan since its access is not coalesced. However OK for small arrays.
    template <typename T, typename U>
    void row_scan(const Integer end, const Integer row_idx, const ViewMatrixTypeLeft<U>& in_, ViewVectorType<T>& out_) {
        using namespace Kokkos;

        parallel_scan(
            RangePolicy<>(0, end), KOKKOS_LAMBDA(const int& i, T& upd, const bool& final) {
                // Load old value in case we update it before accumulating
                const T val_i = in_(row_idx, i);

                upd += val_i;

                if (final) {
                    out_(i + 1) = upd;  // To have both ex and inclusive in the same output.
                }
            });
    }

    template <typename T>
    void complex_inclusive_scan(const Integer start,
                                const Integer end,
                                ViewVectorType<T> index_count_,
                                ViewVectorType<T> pt_count_) {
        using namespace Kokkos;

        // paralel scan on the index_count view for both columns at the same time misusing the complex instead.
        /*	parallel_scan (RangePolicy<>(1 , nr_elements +1 ),	KOKKOS_LAMBDA (const int& i,
                                        complex<Integer>& upd, const bool& final)
                {*/
        parallel_scan(
            RangePolicy<>(start, end), KOKKOS_LAMBDA(const int& i, complex<T>& upd, const bool& final) {
                // Load old value in case we update it before accumulating
                const T val_i = index_count_(i);
                const T val_ip = pt_count_(i);

                /*upd += val_i;
                upd_ip += val_ip;*/
                upd.real() += val_i;
                upd.imag() += val_ip;

                if (final) {
                    index_count_(i) = upd.real();  // only update array on final pass
                    pt_count_(i) = upd.imag();
                    /*index_count_(i,0) = upd;
                    index_count_(i,1) = upd_ip;*/
                }
                // For exclusive scan, change the update value after
                // updating array, like we do here. For inclusive scan,
                // change the update value before updating array.
                /*upd.real() += val_i0;
                upd.imag() += val_i1;*/
            });
    }

    template <typename T>
    void exclusive_scan(const Integer start, const Integer end, ViewVectorType<T> in_) {
        using namespace Kokkos;

        parallel_scan(
            RangePolicy<>(start, end), KOKKOS_LAMBDA(const Integer& i, T& upd, const bool& final) {
                // Load old value in case we update it before accumulating
                const T val_i = in_(i);

                if (final) {
                    in_(i) = upd;  // only update array on final pass
                }

                upd += val_i;
            });
    }

    template <typename T>
    void exclusive_bool_scan(const Integer start, const Integer end, ViewVectorType<T> out_, ViewVectorType<bool> in_) {
        using namespace Kokkos;

        parallel_scan(
            RangePolicy<>(start, end), KOKKOS_LAMBDA(const unsigned int& i, T& upd, const bool& final) {
                // Load old value in case we update it before accumulating
                const T val_i = in_(i);

                if (final) {
                    out_(i) = upd;  // only update array on final pass
                }

                upd += val_i;
            });
    }

    // returns the prefix sum of C into a mirror view
    template <typename C>
    void make_scan_index_mirror(const ViewVectorType<Integer>::HostMirror& out, C const& c) {
        static_assert(std::is_integral<typename C::value_type>::value, "make_index only applies to integral types");

        out(0) = 0;
        std::partial_sum(c.begin(), c.end(), out.data() + 1);
    }

    // Segmented scan on a bool view using hierachical Parallelism. Cub lib has the best impl.
    template <typename F>
    void segmented_scan(const Integer teams, ViewVectorType<bool> in_, F f) {
        Kokkos::parallel_for(
            "naive scan",
            Kokkos::TeamPolicy<>(teams, Kokkos::AUTO),
            MARS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& teamMember) {
                Integer i = teamMember.league_rank();
                if (in_(i) == 1) {
                    Integer segmentSum = 0;
                    Kokkos::parallel_reduce(
                        Kokkos::TeamVectorRange(teamMember, i),
                        [=](const Integer j, Integer& innerUpdate) { innerUpdate += in_(j); },
                        segmentSum);

                    f(segmentSum, i);
                }
            });
    }

}  // namespace mars
#endif /* GENERATION_MARS_UTILS_KOKKOS_HPP_ */
