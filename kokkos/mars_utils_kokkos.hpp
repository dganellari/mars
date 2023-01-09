#ifndef GENERATION_MARS_UTILS_KOKKOS_HPP_
#define GENERATION_MARS_UTILS_KOKKOS_HPP_

#include "mars_config.hpp"

#ifdef MARS_ENABLE_KOKKOS_KERNELS
#include "KokkosKernels_default_types.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#endif

#ifdef MARS_ENABLE_KOKKOS
#include "Kokkos_Core.hpp"
#include "Kokkos_UnorderedMap.hpp"
#ifdef MARS_ENABLE_CUDA
#include <cub/cub.cuh>  // or equivalently <cub/device/device_radix_sort.cuh>
#endif
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
#if defined(MARS_ENABLE_CUDA)
#if defined(MARS_ENABLE_CUDAUVM)
    using KokkosSpace = Kokkos::CudaUVMSpace;
#else
    using KokkosSpace = Kokkos::CudaSpace;
#endif  // MARS_ENABLE_CUDAUVM
    using KokkosLayout = Kokkos::LayoutLeft;
#define MARS_LAMBDA_REF [&] __device__;
#elif defined(MARS_ENABLE_HIP)
    using KokkosSpace = Kokkos::Experimental::HIPSpace;
    using KokkosLayout = Kokkos::LayoutLeft;
#define MARS_LAMBDA_REF [&] __device__;
#else
#ifdef KOKKOS_ENABLE_OPENMP
    using KokkosSpace = Kokkos::HostSpace;
    using KokkosLayout = Kokkos::LayoutRight;
#else   // Serial
    using KokkosSpace = Kokkos::HostSpace;
    using KokkosLayout = Kokkos::LayoutRight;
#endif  // KOKKOS_ENABLE_OPENMP
#define MARS_LAMBDA_REF [&]
#endif  // MARS_ENABLE_CUDA

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

        //Parallel for generates warning for [] operator as called from a device function.
        /* parallel_for(
            MDRangePolicy<Rank<2>, KokkosHostExecSpace>({0, 0}, {xDim, yDim}),
            MARS_LAMBDA(int i, int j) { h_view(i, j) = hostData[i][j]; }); */
        for (int i = 0; i < xDim; ++i) {
            for (int j = 0; j < yDim; ++j) {
                h_view(i, j) = hostData[i][j];
            }
        }

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

#ifdef MARS_ENABLE_CUDA

    // Sorts keys into ascending order. (~2N auxiliary storage required)
    template <typename V>
    V cub_radix_sort(V in) {
        // Declare, allocate, and initialize device-accessible pointers for sorting data
        auto size = in.extent(0);
        V out("out radix sort data", size);

        // Determine temporary device storage requirements
        void* d_temp_storage = NULL;
        size_t temp_storage_bytes = 0;
        cub::DeviceRadixSort::SortKeys(d_temp_storage, temp_storage_bytes, in.data(), out.data(), size);
        // Allocate temporary storage
        cudaMalloc(&d_temp_storage, temp_storage_bytes);
        // Run sorting operation
        cub::DeviceRadixSort::SortKeys(d_temp_storage, temp_storage_bytes, in.data(), out.data(), size);
        return out;
    }

    // Sorts keys into ascending order. (~N auxiliary storage required).
    template <typename T>
    ViewVectorType<T> buffer_cub_radix_sort(ViewVectorType<T> in) {
        // Declare, allocate, and initialize device-accessible pointers for sorting data
        auto size = in.extent(0);
        ViewVectorType<T> out("out radix sort data", size);

        // Create a DoubleBuffer to wrap the pair of device pointers
        cub::DoubleBuffer<T> d_keys(in.data(), out.data());

        // Determine temporary device storage requirements
        void* d_temp_storage = NULL;
        size_t temp_storage_bytes = 0;
        cub::DeviceRadixSort::SortKeys(d_temp_storage, temp_storage_bytes, d_keys, size);
        // Allocate temporary storage
        cudaMalloc(&d_temp_storage, temp_storage_bytes);
        // Run sorting operation
        cub::DeviceRadixSort::SortKeys(d_temp_storage, temp_storage_bytes, d_keys, size);

        //selector 0,1 based on the active buffer.
        if (d_keys.selector)
            return out;
        else
            return in;
    }

    // Block-sorting CUDA kernel
    template <typename T>
    __global__ void BlockSortKernel(T* d_in, T* d_out) {
        using namespace cub;

        // Specialize BlockRadixSort, BlockLoad, and BlockStore for 128 threads
        // owning 16 T items each
        typedef BlockRadixSort<T, 128, 16> BlockRadixSort;
        typedef BlockLoad<T, 128, 16, BLOCK_LOAD_TRANSPOSE> BlockLoad;
        typedef BlockStore<T, 128, 16, BLOCK_STORE_TRANSPOSE> BlockStore;

        // Allocate shared memory
        __shared__ union {
            typename BlockRadixSort::TempStorage sort;
            typename BlockLoad::TempStorage load;
            typename BlockStore::TempStorage store;
        } temp_storage;

        int block_offset = blockIdx.x * (128 * 16);  // OffsetT for this block's ment

        // Obtain a segment of 2048 consecutive keys that are blocked across threads
        T thread_keys[16];
        BlockLoad(temp_storage.load).Load(d_in + block_offset, thread_keys);
        __syncthreads();

        // Collectively sort the keys
        BlockRadixSort(temp_storage.sort).Sort(thread_keys);
        __syncthreads();

        // Store the sorted segment
        BlockStore(temp_storage.store).Store(d_out + block_offset, thread_keys);
    }
#endif


    /* ***************************distributed utils************************************************** */


    template <typename T, Integer Dim>
    void fill_view_vector(ViewVectorTextureC<T, Dim> v, const T value[Dim]) {
        using namespace Kokkos;

        typename ViewVectorTextureC<T, Dim>::HostMirror h_pt = create_mirror_view(v);

        for (Integer i = 0; i < Dim; ++i) h_pt(i) = value[i];

        deep_copy(v, h_pt);
    }

    template <typename T, Integer XDim, Integer YDim>
    void fill_view_matrix(ViewMatrixTextureC<T, XDim, YDim> m, const T value[XDim][YDim]) {
        using namespace Kokkos;

        typename ViewMatrixTextureC<T, XDim, YDim>::HostMirror h_pt = create_mirror_view(m);

        for (Integer i = 0; i < XDim; ++i) {
            for (Integer j = 0; j < YDim; ++j) {
                h_pt(i, j) = value[i][j];
            }
        }

        deep_copy(m, h_pt);
    }

    template <typename T>
    void print_view(const ViewVectorType<T>& view) {
        using namespace Kokkos;

        parallel_for(
            "print view", view.extent(0), KOKKOS_LAMBDA(const Integer i) {
                std::cout << "i: " << i << " value: " << view(i) << std::endl;
            });
    }

    struct resize_view_functor {
        resize_view_functor(std::string d, size_t s) : _desc(d), _size(s) {}
        template <typename ElementType>
        void operator()(ElementType& el, std::size_t I) const {
            el = ElementType(_desc + std::to_string(I), _size);
        }

        std::string _desc;
        size_t _size;
    };

    struct resize_functor {
        resize_functor(size_t s) : _size(s) {}
        template <typename ElementType>
        void operator()(ElementType& el, std::size_t I) const {
            el = ElementType(_size);
        }
        size_t _size;
    };

    struct print_functor {
        template <typename ElementType>
        void operator()(ElementType& el) const {
            std::cout << el << std::endl;
        }
    };

    /* special implementation of the binary search considering found
    if an element is between current and next proc value. */
    template <typename T>
    MARS_INLINE_FUNCTION Integer
    find_owner_processor(const ViewVectorType<T> view, const T enc_oc, const int offset, Integer guess) {
        const int last_index = view.extent(0) / offset - 1;
        int first_proc = 0;
        int last_proc = last_index;

        while (first_proc <= last_proc && first_proc != last_index) {
            T current = view(offset * guess);
            T next = view(offset * (guess + 1));

            if (enc_oc >= current && enc_oc < next) {
                return guess;
            } else if (enc_oc < current) {
                last_proc = guess - 1;
                guess = (first_proc + last_proc + 1) / 2;
            } else if (enc_oc >= next) {
                first_proc = guess + 1;
                guess = (first_proc + last_proc) / 2;
            }
        }

        return -1;
    }

    // standard binary search
    template <typename T>
    MARS_INLINE_FUNCTION Integer binary_search(const T* view, Integer f, Integer l, const T enc_oc) {
        while (f <= l) {
            Integer guess = (l + f) / 2;

            T current = *(view + guess);

            if (enc_oc == current) {
                return guess;
            } else if (enc_oc < current) {
                l = guess - 1;
            } else {
                f = guess + 1;
            }
        }

        return -1;
    }

    // Trilinos way of doing abs max and min
    template <class T, class H>
    struct AbsMinOp {
        MARS_INLINE_FUNCTION
        static T apply(const T& val1, const H& val2) {
            const auto abs1 = Kokkos::ArithTraits<T>::abs(val1);
            const auto abs2 = Kokkos::ArithTraits<H>::abs(val2);
            return abs1 < abs2 ? T(abs1) : H(abs2);
        }
    };

    template <typename SC>
    struct atomic_abs_min {
        MARS_INLINE_FUNCTION
        void operator()(SC& dest, const SC& src) const {
            Kokkos::Impl::atomic_fetch_oper(AbsMinOp<SC, SC>(), &dest, src);
        }
    };

    template <class T, class H>
    struct AbsMaxOp {
        MARS_INLINE_FUNCTION
        static T apply(const T& val1, const H& val2) {
            const auto abs1 = Kokkos::ArithTraits<T>::abs(val1);
            const auto abs2 = Kokkos::ArithTraits<H>::abs(val2);
            return abs1 > abs2 ? T(abs1) : H(abs2);
        }
    };

    template <typename SC>
    struct atomic_abs_max {
        KOKKOS_INLINE_FUNCTION
        void operator()(SC& dest, const SC& src) const {
            Kokkos::Impl::atomic_fetch_oper(AbsMaxOp<SC, SC>(), &dest, src);
        }
    };

    template <typename H>
    struct AtomicOp {
        AtomicOp(H f) : func(f) {}

        MARS_INLINE_FUNCTION
        void operator()(double& dest, const double& src) const { Kokkos::Impl::atomic_fetch_oper(func, &dest, src); }

        H func;
    };

    template <typename H, typename S>
    MARS_INLINE_FUNCTION void atomic_op(H f, S& dest, const S& src) {
        Kokkos::Impl::atomic_fetch_oper(f, &dest, src);
    }

    // max plus functor
    template <typename T>
    class MaxPlus {
    public:
        // Kokkos reduction functors need the value_type typedef.
        // This is the type of the result of the reduction.
        typedef T value_type;

        MaxPlus(const ViewVectorType<T> x) : x_(x) {}

        // This is helpful for determining the right index type,
        // especially if you expect to need a 64-bit index.
        typedef typename ViewVectorType<T>::size_type size_type;

        KOKKOS_INLINE_FUNCTION void operator()(const size_type i,
                                               value_type& update) const {  // max-plus semiring equivalent of "plus"
            if (update < x_(i)) {
                update = x_(i);
            }
        }

        // "Join" intermediate results from different threads.
        // This should normally implement the same reduction
        // operation as operator() above. Note that both input
        // arguments MUST be declared volatile.
        KOKKOS_INLINE_FUNCTION void join(
            volatile value_type& dst,
            const volatile value_type& src) const {  // max-plus semiring equivalent of "plus"
            if (dst < src) {
                dst = src;
            }
        }

        // Tell each thread how to initialize its reduction result.
        KOKKOS_INLINE_FUNCTION void init(value_type& dst) const {  // The identity under max is -Inf.
            dst = Kokkos::reduction_identity<value_type>::max();
        }

    private:
        ViewVectorType<T> x_;
    };

    //**********************************dof handler utils******************************************

    // The staggered grid implementation forbids unifying these two methods.
    // They work for both general and staggered dof handlers (unified for handlers).
    // However not possible to unify for owned and local. Either one or the other.
    template <Integer Label, typename H>
    ViewVectorType<Integer> compact_owned_dofs(const H& dof_handler, ViewVectorType<Integer>& locally_owned_dofs) {
        using namespace Kokkos;

        const Integer local_size = dof_handler.template get_owned_dof_size<1>();

        ViewVectorType<bool> dof_predicate("label_dof_predicate", local_size);
        Kokkos::parallel_for(
            "separateowneddoflabels", local_size, KOKKOS_LAMBDA(const Integer i) {
                const Integer local = dof_handler.template get_owned_dof<1>(i);
                if (dof_handler.template get_label<1>(local) & Label) {
                    dof_predicate(i) = 1;
                }
            });

        /* perform a scan on the dof predicate*/
        ViewVectorType<Integer> owned_dof_map("owned_dof_scan", local_size + 1);
        incl_excl_scan(0, local_size, dof_predicate, owned_dof_map);

        auto vol_subview = subview(owned_dof_map, local_size);
        auto h_vs = create_mirror_view(vol_subview);
        // Deep copy device view to host view.
        deep_copy(h_vs, vol_subview);

        locally_owned_dofs = ViewVectorType<Integer>("locally_owned_dofs", h_vs());

        parallel_for(
            local_size, KOKKOS_LAMBDA(const Integer i) {
                if (dof_predicate(i) == 1) {
                    Integer vindex = owned_dof_map(i);
                    const Integer local = dof_handler.template get_owned_dof<1>(i);
                    locally_owned_dofs(vindex) = local;
                }
            });

        return owned_dof_map;
    }

    template <Integer Label, typename H>
    ViewVectorType<Integer> compact_local_dofs(const H& dof_handler, ViewVectorType<Integer>& local_dofs) {
        using namespace Kokkos;

        const Integer local_size = dof_handler.get_base_dof_size();

        ViewVectorType<bool> dof_predicate("label_dof_predicate", local_size);
        Kokkos::parallel_for(
            "separatelocaldoflabels", local_size, KOKKOS_LAMBDA(const Integer i) {
                const Integer local = dof_handler.get_local_dof(i);
                if (dof_handler.template get_label<1>(local) & Label) {
                    dof_predicate(i) = 1;
                }
            });

        /* perform a scan on the dof predicate*/
        ViewVectorType<Integer> local_dof_map("local_dof_scan", local_size + 1);
        incl_excl_scan(0, local_size, dof_predicate, local_dof_map);

        auto vol_subview = subview(local_dof_map, local_size);
        auto h_vs = create_mirror_view(vol_subview);
        // Deep copy device view to host view.
        deep_copy(h_vs, vol_subview);

        local_dofs = ViewVectorType<Integer>("local_dofs", h_vs());

        /* Compact the predicate into the volume and face dofs views */
        parallel_for(
            local_size, KOKKOS_LAMBDA(const Integer i) {
                if (dof_predicate(i) == 1) {
                    Integer vindex = local_dof_map(i);
                    const Integer local = dof_handler.get_local_dof(i);
                    local_dofs(vindex) = local;
                }
            });

        return local_dof_map;
    }

    template <typename DH, typename F>
    ViewVectorType<bool> build_sfc_to_local_predicate(DH dofhandler,
                                                      F f,
                                                      const Integer local_size,
                                                      const ViewVectorType<Integer> in) {
        ViewVectorType<bool> predicate("separated_predicate", local_size);

        Kokkos::parallel_for(
            "separatedoflabelss", local_size, KOKKOS_LAMBDA(const Integer i) {
                const Integer sfc = in(i);
                /* if (f(dofhandler.sfc_to_local(sfc))) { */
                if (f(sfc)) {
                    predicate(i) = 1;
                }
            });

        return predicate;
    }

    template <typename DH, typename F>
    inline ViewVectorType<bool> compact_sfc_to_local(DH dofhandler,
                                                     F f,
                                                     const ViewVectorType<Integer> in,
                                                     ViewVectorType<Integer>& out) {
        using namespace Kokkos;

        auto local_size = in.extent(0);
        auto in_predicate = build_sfc_to_local_predicate<DH, F>(dofhandler, f, local_size, in);

        /* perform a scan on the separated dof predicate*/
        auto bscan = ViewVectorType<Integer>("boudnary_separated_scan", local_size + 1);
        incl_excl_scan(0, local_size, in_predicate, bscan);

        auto vol_subview = subview(bscan, local_size);
        auto h_vs = create_mirror_view(vol_subview);
        // Deep copy device view to host view.
        deep_copy(h_vs, vol_subview);

        out = ViewVectorType<Integer>("local_separated_dofs", h_vs());
        /* ViewVectorType<Integer> bvds = out; */

        /* Compact the predicate into the separated and face dofs views */
        parallel_for(
            local_size, KOKKOS_LAMBDA(const Integer i) {
                if (in_predicate(i) == 1) {
                    Integer vindex = bscan(i);
                    out(vindex) = dofhandler.sfc_to_local(in(i));
                }
            });
        return in_predicate;
    }

    inline ViewVectorType<Integer> count_sfc_to_local(const ViewVectorType<Integer> in_scan,
                                                      const ViewVectorType<bool> in_predicate) {
        using namespace Kokkos;

        const Integer count_size = in_scan.extent(0) - 1;

        ViewVectorType<Integer> count("count_to_local", count_size);
        ViewVectorType<Integer> scan("scan_to_local", count_size + 1);
        /* Compact the predicate into the separated and face dofs views */
        parallel_for(
            in_predicate.extent(0), KOKKOS_LAMBDA(const Integer i) {
                if (in_predicate(i) == 1) {
                    Integer proc = find_owner_processor(in_scan, i, 1, 0);
                    atomic_increment(&count(proc));
                }
            });

        incl_excl_scan(0, count_size, count, scan);

        return scan;
    }

}  // namespace mars
#endif /* GENERATION_MARS_UTILS_KOKKOS_HPP_ */
