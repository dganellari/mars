#ifndef MARS_STATIC_MATH_KOKKOS_HPP
#define MARS_STATIC_MATH_KOKKOS_HPP

#include <array>
#include <iostream>
#include <vector>

#include "mars_base.hpp"

#include "mars_static_math.hpp"
#include "mars_sub_view.hpp"
#include "mars_utils_kokkos.hpp"

namespace mars {

    template <Integer N, Integer K>
    class CombinationsAux<N, K, KokkosImplementation> {
    public:
        static constexpr Integer NChooseK = Factorial<N>::value / (Factorial<K>::value * Factorial<N - K>::value);

        static MARS_INLINE_FUNCTION void apply(ViewMatrixTextureC<Integer, NChooseK, K> combs) {
            TempArray<Integer, K> data;
            // printf("Applying to combs\n");
            Integer comb_index = 0;
            apply(data, 0, 0, combs, comb_index);
        }

    private:
        static MARS_INLINE_FUNCTION void apply(
            TempArray<Integer, K>& data,  // call this with the shared memory view shared pointer
            const Integer index,
            const Integer i,
            ViewMatrixTextureC<Integer, NChooseK, K> combs,
            Integer& comb_index) {
            if (index == K) {
                // std::copy(data, data+K, combs[comb_index++]);
                for (Integer j = 0; j < K; ++j) {
                    combs(comb_index, j) = data[j];
                }
                // printf("Applying to combs with idx=%li\n", comb_index);
                ++comb_index;
                return;
            }

            if (i >= N) {
                return;
            }

            data[index] = i;

            apply(data, index + 1, i + 1, combs, comb_index);

            // current is excluded, replace it with next (Note that
            // i+1 is passed, but index is not changed)
            apply(data, index, i + 1, combs, comb_index);
        }
    };

    template <Integer N, Integer ChooseM>
    class Combinations<N, ChooseM, KokkosImplementation> {
    public:
        static constexpr Integer value =
            Factorial<N>::value / (Factorial<ChooseM>::value * Factorial<N - ChooseM>::value);
        ViewMatrixTextureC<Integer, value, ChooseM> combs;

        struct GenerateComb {
            ViewMatrixTextureC<Integer, value, ChooseM> combs;

            GenerateComb(ViewMatrixTextureC<Integer, value, ChooseM> cmb) : combs(cmb) {}

            MARS_INLINE_FUNCTION
            void operator()(int i) const { CombinationsAux<N, ChooseM, KokkosImplementation>::apply(combs); }
        };

        // reinitilize=false should only be used in the destructor of the bisection for deallocation of combs.
        static Combinations& instance(bool reinitilize = true) {
            static Combinations instance_;
            if (instance_.combs.data() == nullptr && reinitilize) {
                instance_.combs = ViewMatrixTextureC<Integer, value, ChooseM>("combs");
                Kokkos::parallel_for(1, GenerateComb(instance_.combs));
            }
            return instance_;
        }

        template <typename T>
        MARS_INLINE_FUNCTION static void choose(const Integer k, const SubView<Integer, N>& in, T* out) {
            assert(k < value);

            for (Integer i = 0; i < ChooseM; ++i) {
                assert(instance().combs(k, i) < N);
                assert(instance().combs(k, i) >= 0);

                out[i] = in[instance().combs(k, i)];
            }
        }

    private:
        Combinations() : combs("combinations") { Kokkos::parallel_for(1, GenerateComb(combs)); }
    };

}  // namespace mars

#endif  // MARS_STATIC_MATH_HPP
