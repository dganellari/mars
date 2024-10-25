#include <gtest/gtest.h>
#include "mars_sfc_code.hpp"
#include <iostream>

using namespace mars;



TEST(MortonCode, EncodeDecode2D) {
    unsigned x = 5, y = 10;
    auto code = encode_morton_2D(x, y);
    auto decoded = decode_morton_2D<unsigned>(code);
    EXPECT_EQ(decoded.x, x);
    EXPECT_EQ(decoded.y, y);
}

template<class KeyType>
void encode_morton(){
    KeyType x = 5, y = 10, z = 15;
    auto code = encode_morton_3D<KeyType>(x, y, z);
    auto decoded = decode_morton_3D<KeyType>(code);
    EXPECT_EQ(decoded.x, x);
    EXPECT_EQ(decoded.y, y);
    EXPECT_EQ(decoded.z, z);
}

TEST(MortonCode, EncodeDecode3D) {
    encode_morton<unsigned>();
    encode_morton<uint64_t>();
}

template<class KeyType>
void decode_morton(){
    KeyType x = 4, y = 2, z = 5;

    KeyType code = 340;
    auto octant = decode_morton_3D<KeyType>(code);
    EXPECT_EQ(x, octant.x);
    EXPECT_EQ(y, octant.y);
    EXPECT_EQ(z, octant.z);
}

TEST(MortonCode, DecodeMorton3D) {
    decode_morton<unsigned>();
    decode_morton<uint64_t>();
}

TEST(MortonCode, decodeMorton64)
{
    uint64_t code = 0x7FFFFFFFFFFFFFFFlu;
    auto octant = decode_morton_3D<uint64_t>(code);
    EXPECT_EQ((1u << 21u) - 1u, octant.z);
    EXPECT_EQ((1u << 21u) - 1u, octant.y);
    EXPECT_EQ((1u << 21u) - 1u, octant.x);

    code = 0x1249249241249249;
    EXPECT_EQ((1u << 21u) - 512u - 1u, decode_morton_3D<uint64_t>(code).x);

    code = 0b0111lu << (20u * 3);
    auto octant2 = decode_morton_3D<uint64_t>(code);
    EXPECT_EQ(1u << 20u, octant2.z);
    EXPECT_EQ(1u << 20u, octant2.y);
    EXPECT_EQ(1u << 20u, octant2.x);

    code = 0b0011lu << (20u * 3);
    auto octant3 = decode_morton_3D<uint64_t>(code);
    EXPECT_EQ(0, octant3.z);
    EXPECT_EQ(1u << 20u, octant3.y);
    EXPECT_EQ(1u << 20u, octant3.x);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
