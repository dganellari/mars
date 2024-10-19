/*
 * MIT License
 *
 * Copyright (c) 2021 CSCS, ETH Zurich
 *               2021 University of Basel
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <gtest/gtest.h>
#include <random>
#include "mars_distributed_octant.hpp"
#include "mars_sfc_code.hpp"

using namespace mars;

//! @brief tests numKeys random 3D points for encoding/decoding consistency
template<class KeyType>
void inversionTest2D()
{
    int numKeys  = 1000;
    int maxCoord = (1 << maxTreeLevel<KeyType>{}) - 1;

    std::mt19937 gen;
    std::uniform_int_distribution<unsigned> distribution(0, maxCoord);

    auto getRand = [&distribution, &gen]() { return distribution(gen); };

    std::vector<unsigned> x(numKeys);
    std::vector<unsigned> y(numKeys);

    std::generate(begin(x), end(x), getRand);
    std::generate(begin(y), end(y), getRand);

    for (int i = 0; i < numKeys; ++i)
    {
        KeyType hilbertKey = encode_hilbert_2D<KeyType>(x[i], y[i]);
        auto octant = decode_hilbert_2D<KeyType>(hilbertKey);
        EXPECT_EQ(x[i], octant.x);
        EXPECT_EQ(y[i], octant.y);
    }
}

TEST(HilbertCode, inversion2D)
{
    inversionTest2D<unsigned>();
    inversionTest2D<uint64_t>();
}

//! @brief 2D test the 2D Hilbert curve of order 1 and 2
template<class KeyType>
void Hilbert2D(int order)
{
    for (unsigned xi = 0; xi < pow(2, order); ++xi)
    {
        for (unsigned yi = 0; yi < pow(2, order); ++yi)
        {
            KeyType key = encode_hilbert_2D<KeyType>(xi, yi);
            auto t      = decode_hilbert_2D<KeyType>(key);
            EXPECT_EQ(xi, t.x);
            EXPECT_EQ(yi, t.y);
        }
    }
}

TEST(HilbertCode, Hilbert2D)
{
    Hilbert2D<unsigned>(2);
    Hilbert2D<uint64_t>(1);
}

//! @brief test the curve on the first eight octants
template<class KeyType>
void firstOrderCurve()
{
    constexpr unsigned hilbertToMorton[8] = {0, 1, 3, 2, 6, 7, 5, 4};

    for (unsigned xi = 0; xi < 2; ++xi)
    {
        for (unsigned yi = 0; yi < 2; ++yi)
        {
            for (unsigned zi = 0; zi < 2; ++zi)
            {
                unsigned L1Range      = (1 << maxTreeLevel<KeyType>{}) / 2;
                unsigned mortonOctant = 4 * xi + 2 * yi + zi;

                {
                    KeyType hilbertKey     = encode_hilbert_3D<KeyType>(L1Range * xi, L1Range * yi, L1Range * zi);
                    unsigned hilbertOctant = octalDigit(hilbertKey, 1);
                    EXPECT_EQ(mortonOctant, hilbertToMorton[hilbertOctant]);
                }
                {
                    KeyType hilbertKey     = encode_hilbert_3D<KeyType>(L1Range * xi + L1Range - 1, L1Range * yi + L1Range - 1,
                                                           L1Range * zi + L1Range - 1);
                    unsigned hilbertOctant = octalDigit(hilbertKey, 1);
                    EXPECT_EQ(mortonOctant, hilbertToMorton[hilbertOctant]);
                }
            }
        }
    }
}

TEST(HilbertCode, firstOrderCurve)
{
    firstOrderCurve<unsigned>();
    firstOrderCurve<uint64_t>();
}

//! @brief verifies continuity properties across consecutive octants at all levels
template<class KeyType>
void continuityTest()
{
    for (unsigned level = 1; level < maxTreeLevel<KeyType>{}; ++level)
    {
        // on the highest level, we can only check 7 octant crossings
        int maxOctant = (level > 1) ? 8 : 7;

        for (int octant = 0; octant < maxOctant; ++octant)
        {
            KeyType lastKey      = (octant + 1) * nodeRange<KeyType>(level) - 1;
            KeyType firstNextKey = lastKey + 1;

            auto lastOctant = decode_hilbert_3D<KeyType>(lastKey);

            auto nextOctant = decode_hilbert_3D<KeyType>(firstNextKey);

            // the points in 3D space should be right next to each other, i.e. delta == 1
            // this is a property that the Z-curve does not have
            int delta = std::abs(int(lastOctant.x) - int(nextOctant.x)) + std::abs(int(lastOctant.y) - int(nextOctant.y)) +
                        std::abs(int(lastOctant.z) - int(nextOctant.z));

            EXPECT_EQ(delta, 1);
        }
    }
}

TEST(HilbertCode, continuity)
{
    continuityTest<unsigned>();
    continuityTest<uint64_t>();
}

//! @brief tests numKeys random 3D points for encoding/decoding consistency
template<class KeyType>
void inversionTest()
{
    int numKeys  = 1000;
    int maxCoord = (1 << maxTreeLevel<KeyType>{}) - 1;

    std::mt19937 gen;
    std::uniform_int_distribution<unsigned> distribution(0, maxCoord);

    auto getRand = [&distribution, &gen]() { return distribution(gen); };

    std::vector<unsigned> x(numKeys);
    std::vector<unsigned> y(numKeys);
    std::vector<unsigned> z(numKeys);

    std::generate(begin(x), end(x), getRand);
    std::generate(begin(y), end(y), getRand);
    std::generate(begin(z), end(z), getRand);

    for (int i = 0; i < numKeys; ++i)
    {
        KeyType hilbertKey = encode_hilbert_3D<KeyType>(x[i], y[i], z[i]);

        auto octant = decode_hilbert_3D<KeyType>(hilbertKey);
        EXPECT_EQ(x[i], octant.x);
        EXPECT_EQ(y[i], octant.y);
        EXPECT_EQ(z[i], octant.z);
    }
}

TEST(HilbertCode, inversion)
{
    inversionTest<unsigned>();
    inversionTest<uint64_t>();
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
