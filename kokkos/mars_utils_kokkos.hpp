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
#include <Kokkos_Sort.hpp>
#ifdef MARS_ENABLE_CUDA
#include <cub/cub.cuh>  // or equivalently <cub/device/device_radix_sort.cuh>
#include <thrust/unique.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
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


    //generate a kokkos function that takes a view as input and returns the last value of the view by coping it to the host.
    template <typename T>
    T get_last_value(ViewVectorType<T> view) {
        using namespace Kokkos;

        auto sub_view = subview(view, view.extent(0) - 1);
        auto h_view = create_mirror_view(sub_view);
        deep_copy(h_view, sub_view);

        return h_view();
    }

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
    void complex_inclusive_scan(const Integer start,
                                const Integer end,
                                ViewVectorType<T> index_count_,
                                ViewVectorType<T> pt_count_) {
        using namespace Kokkos;
        // paralel scan on the index_count view for both columns at the same time misusing the complex instead.
        // Warning: Kokkos complex only works on floating point types
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


    //write a function that prints a view in parallel with kokkos
    //
    template <typename T>
    void print_view(ViewVectorType<T> view, const std::string& name) {
        using namespace Kokkos;

        parallel_for(
            "print view", view.extent(0), KOKKOS_LAMBDA(const int i) { printf("%s[%d] = %d\n", name.c_str(), i, view(i)); });
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

    // works for all cases including those with strided access views coming from auto = subview...
    template <typename H, typename U, typename S, typename J>
    void incl_excl_scan(const S start, const J end, const U in_, H out_) {
        using namespace Kokkos;

        using T = typename H::value_type;

        parallel_scan(
            RangePolicy<IndexType<J>>(start, end), KOKKOS_LAMBDA(const J i, T& upd, const bool final) {
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

    //generate a kokkos function that takes a view as input  then creates a mirror host view, fills it up and copies it to the gpu memory.
    template <typename V, typename H, typename U>
    void host_scan_to_device(const V& view, H& h_view, U const& data) {
        h_view = Kokkos::create_mirror_view(view);
        make_scan_index_mirror(h_view, data);
        Kokkos::deep_copy(view, h_view);
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

    template <typename T>
    void kokkos_sort(ViewVectorType<T> d_in) {
        Kokkos::sort(d_in);
    }

    template <typename T>
    ViewVectorType<T> kokkos_unique(const ViewVectorType<T>& d_in) {
        // Create a view to hold the unique flags
        ViewVectorType<bool> flags("flags", d_in.extent(0));

        // Mark the unique elements
        Kokkos::parallel_for(
            "mark_unique", Kokkos::RangePolicy<>(0, d_in.extent(0)), KOKKOS_LAMBDA(int i) {
                if (i == 0) {
                    flags(i) = true;  // The first element is always unique
                } else {
                    flags(i) = (d_in(i) != d_in(i - 1));
                }
            });

        // Count the number of unique elements
        int unique_count = 0;
        Kokkos::parallel_reduce(
            "count_unique",
            d_in.extent(0),
            KOKKOS_LAMBDA(int i, int& count) {
                if (flags(i)) count++;
            },
            unique_count);

        // Create a view to hold the unique elements
        ViewVectorType<T> d_out("d_out", unique_count);

        // Copy the unique elements
        Kokkos::parallel_scan(
            "copy_unique", d_in.extent(0), KOKKOS_LAMBDA(int i, int& index, bool final) {
                if (flags(i)) {
                    if (final) d_out(index) = d_in(i);
                    index++;
                }
            });

        return d_out;
    }

#ifdef MARS_ENABLE_CUDA

    template <typename T>
    void cub_inclusive_scan(ViewVectorType<T> data_) {
        // Determine temporary device storage requirements
        void* d_temp_storage = nullptr;
        size_t temp_storage_bytes = 0;
        cub::DeviceScan::InclusiveSum(d_temp_storage, temp_storage_bytes, data_.data(), data_.data(), data_.extent(0));
        // Allocate temporary storage
        cudaMalloc(&d_temp_storage, temp_storage_bytes);
        // Run inclusive scan operation
        cub::DeviceScan::InclusiveSum(d_temp_storage, temp_storage_bytes, data_.data(), data_.data(), data_.extent(0));
    }

    template <typename T>
    void cub_unique(ViewVectorType<T> d_in, ViewVectorType<T> d_out, ViewVectorType<T> d_num_selected_out) {
        auto num_items = d_in.extent(0);

        void* d_temp_storage = nullptr;
        size_t temp_storage_bytes = 0;
        cub::DeviceSelect::Unique(
            d_temp_storage, temp_storage_bytes, d_in.data(), d_out.data(), d_num_selected_out.data(), num_items);
        // Allocate temporary storage
        cudaMalloc(&d_temp_storage, temp_storage_bytes);
        // Run selection
        cub::DeviceSelect::Unique(
            d_temp_storage, temp_storage_bytes, d_in.data(), d_out.data(), d_num_selected_out.data(), num_items);
    }

    template <typename T>
    ViewVectorType<T> thrust_unique(const ViewVectorType<T>& d_in) {
        thrust::device_ptr<T> d_ptr_in(d_in.data());
        //do unique and count the unique elements
        auto end_unique = thrust::unique(d_ptr_in, d_ptr_in + d_in.extent(0));
        int unique_count = thrust::distance(d_ptr_in, end_unique);
        //allocate the output size based on the unique count.
        ViewVectorType<T> d_out(Kokkos::view_alloc("d_out no init", Kokkos::WithoutInitializing), unique_count);
        //copy from thrust to d_out
        thrust::copy(d_ptr_in, d_ptr_in + unique_count, d_out.data());
        return d_out;
    }

    //using thrust unique_count for newer thrust versions
    /* template <typename T>
    void thrust_unique(ViewVectorType<T> d_in) {
        thrust::device_ptr<T> d_ptr_in(d_in.data());
        auto count = thrust::unique_count(d_ptr_in, d_ptr_in + d_in.extent(0), thrust::equal_to<T>());

        ViewVectorType<T> d_out("unique data", count);
        thrust::device_ptr<T> d_ptr_out(d_out.data());

        thrust::unique_copy(d_ptr_in, d_ptr_in + d_in.extent(0), d_ptr_out);
        return d_out;
    } */

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

    // Sorts keys into ascending order. (~N auxiliary storage required).
    template <typename T>
    void buffer_cub_radix_sort_pairs(ViewVectorType<T>& keys, ViewVectorType<T>& values) {
        // Declare, allocate, and initialize device-accessible pointers for sorting data
        auto size = keys.extent(0);
        ViewVectorType<T> out_keys("out radix sort data", size);
        ViewVectorType<T> out_values("out radix sort values", size);

        // Create a DoubleBuffer to wrap the pair of device pointers
        cub::DoubleBuffer<T> d_keys(keys.data(), out_keys.data());
        cub::DoubleBuffer<T> d_values(values.data(), out_values.data());

        // Determine temporary device storage requirements
        void* d_temp_storage = NULL;
        size_t temp_storage_bytes = 0;
        cub::DeviceRadixSort::SortKeys(d_temp_storage, temp_storage_bytes, d_keys, d_values, size);
        // Allocate temporary storage
        cudaMalloc(&d_temp_storage, temp_storage_bytes);
        // Run sorting operation
        cub::DeviceRadixSort::SortKeys(d_temp_storage, temp_storage_bytes, d_keys, d_values, size);

        cudaFree(d_temp_storage);

        keys = d_keys.selector ? out_keys : keys;
        values = d_values.selector ? out_values : values;
    }

    // Sorts keys into ascending order. (~N auxiliary storage required).
    template <typename T>
    ViewVectorType<T> buffer_cub_segmented_radix_sort(ViewVectorType<T> in, ViewVectorType<T> d_scan) {
        // Declare, allocate, and initialize device-accessible pointers for sorting data
        auto size = in.extent(0);
        ViewVectorType<T> out("out radix sort data", size);

        auto num_segments = d_scan.extent(0) - 1;

        // Create a DoubleBuffer to wrap the pair of device pointers
        cub::DoubleBuffer<T> d_keys(in.data(), out.data());

        // Determine temporary device storage requirements
        void* d_temp_storage = NULL;
        size_t temp_storage_bytes = 0;
        cub::DeviceSegmentedRadixSort::SortKeys(d_temp_storage, temp_storage_bytes, d_keys, size, num_segments, d_scan.data(), d_scan.data() + 1);
        // Allocate temporary storage
        cudaMalloc(&d_temp_storage, temp_storage_bytes);
        // Run sorting operation
        cub::DeviceSegmentedRadixSort::SortKeys(d_temp_storage, temp_storage_bytes, d_keys, size, num_segments, d_scan.data(), d_scan.data() + 1);

        cudaFree(d_temp_storage);
        //selector 0,1 based on the active buffer.
        if (d_keys.selector)
            return out;
        else
            return in;
    }

    // Function to perform segmented unique operation using Thrust with transforming the original values to get to different segment values and work with the whole array instead of segmentes.
    template <typename T>
    ViewVectorType<T> thrust_segmented_unique(ViewVectorType<T>& d_in,
                                              ViewVectorType<T>& d_count,
                                              ViewVectorType<T> d_scan,
                                              int max_element) {
        auto size = d_in.extent(0);
        auto num_segments = d_scan.extent(0) - 1;

        // Create flags to denote membership of each element to the respective segment
        thrust::device_vector<int> flags(size);
        for (int i = 0; i < num_segments; ++i) {
            thrust::fill(flags.begin() + d_scan(i), flags.begin() + d_scan(i + 1), i);
        }

        // Transform data
        thrust::device_ptr<T> d_ptr_in(d_in.data());
        thrust::transform(d_ptr_in, d_ptr_in + size, flags.begin(), d_ptr_in, [max_element] __device__(T x, int flag) {
            return x + flag * 2 * max_element;
        });

        // Sort the transformed data
        /* buffer_cub_radix_sort_pairs(d_in, flags); */
        thrust::sort_by_key(d_ptr_in, d_ptr_in + size, flags.begin());

        // Apply thrust::unique_by_key
        thrust::device_vector<int> flags_out(size);
        auto end = thrust::unique_by_key(d_ptr_in, d_ptr_in + size, flags.begin(), flags_out.begin());
        int unique_count = thrust::distance(d_ptr_in, end.first);

        // Transform data back
        thrust::transform(
            d_ptr_in, d_ptr_in + unique_count, flags_out.begin(), d_ptr_in, [max_element] __device__(T x, int flag) {
                return x - flag * 2 * max_element;
            });

        // Allocate the output size based on the unique count
        ViewVectorType<T> d_out(Kokkos::view_alloc("d_out no init", Kokkos::WithoutInitializing), unique_count);
        thrust::copy(d_ptr_in, d_ptr_in + unique_count, d_out.data());

        // Use reduce_by_key to get unique keys and their counts
        thrust::device_vector<int> counts(unique_count, 1);
        thrust::device_vector<int> unique_keys(unique_count);
        thrust::device_vector<int> unique_counts(unique_count);
        auto reduce_end = thrust::reduce_by_key(flags_out.begin(),
                                                flags_out.begin() + unique_count,
                                                counts.begin(),
                                                unique_keys.begin(),
                                                unique_counts.begin());

        // Use the unique keys and counts to fill out d_count
        thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(unique_keys.begin(), unique_counts.begin())),
                         thrust::make_zip_iterator(thrust::make_tuple(reduce_end.first, reduce_end.second)),
                         [d_count_data = d_count.data()] __device__(const thrust::tuple<int, int>& t) {
                             d_count_data[thrust::get<0>(t)] = thrust::get<1>(t);
                         });

        return d_out;
    }

// Function to perform sort and unique operation
template <typename T>
ViewVectorType<T> segmented_sort_unique(ViewVectorType<T>& d_in, ViewVectorType<T>& d_count, ViewVectorType<int> d_scan) {
    // Determine the maximum element in the input data
    thrust::device_ptr<T> d_ptr_in(d_in.data());
    T max_element = *thrust::max_element(d_ptr_in, d_ptr_in + d_in.extent(0));

    // Perform unique and count the unique elements
    auto unique = thrust_segmented_unique(d_in, d_count, d_scan, max_element);
    return unique;
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

}  // namespace mars
#endif /* GENERATION_MARS_UTILS_KOKKOS_HPP_ */
