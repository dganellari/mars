#include <gtest/gtest.h>
#include "mars_sfc_code.hpp"

using namespace mars;

template <class KeyType>
void encode_morton3D() {
    constexpr unsigned tree_level = 3;
    std::array<unsigned, 3> key = {5, 3, 6};
    KeyType morton_key = encode_morton_3D<KeyType>(key[0], key[1], key[2]);
    EXPECT_EQ(morton_key, pad(KeyType(0b101011110), 9));
}

TEST(MortonCode, imorton3D) {
    encode_morton3D<unsigned>();
    encode_morton3D<uint64_t>();
}

TEST(MortonCode, idecodeMorton32) {
    unsigned x = 5;
    unsigned y = 2;
    unsigned z = 4;

    unsigned code = 340;
    auto octant = decode_morton_3D(code);
    EXPECT_EQ(x, octant.x);
    EXPECT_EQ(y, octant.y);
    EXPECT_EQ(z, octant.z);
}

TEST(MortonCode, idecodeMorton64)
{
    std::size_t code = 0x7FFFFFFFFFFFFFFFlu;
    auto octant = decode_morton_3D(code);
    EXPECT_EQ((1u << 21u) - 1u, octant.x);
    EXPECT_EQ((1u << 21u) - 1u, octant.y);
    EXPECT_EQ((1u << 21u) - 1u, octant.z);

    code = 0x1249249241249249;
    EXPECT_EQ((1u << 21u) - 512u - 1u, decode_morton_3D(code).z);

    code = 0b0111lu << (20u * 3);
    auto octant2 = decode_morton_3D(code);
    EXPECT_EQ(1u << 20u, octant2.x);
    EXPECT_EQ(1u << 20u, octant2.y);
    EXPECT_EQ(1u << 20u, octant2.z);

    code = 0b0011lu << (20u * 3);
    auto octant3 = decode_morton_3D(code);
    EXPECT_EQ(0, octant3.x);
    EXPECT_EQ(1u << 20u, octant3.y);
    EXPECT_EQ(1u << 20u, octant3.z);
}

template<class KeyType>
void mortonNeighbors()
{
    //                      input    ref. out  treeLevel dx   dy   dz
    std::vector<std::tuple<KeyType, KeyType, unsigned, int, int, int>> codes{
        {pad(KeyType(0b000111111), 9), pad(KeyType(0b000111011), 9), 3, -1, 0, 0},
        {pad(KeyType(0b000111111), 9), pad(KeyType(0b100011011), 9), 3, 1, 0, 0},
        {pad(KeyType(0b000111111), 9), pad(KeyType(0b000111101), 9), 3, 0, -1, 0},
        {pad(KeyType(0b000111111), 9), pad(KeyType(0b010101101), 9), 3, 0, 1, 0},
        {pad(KeyType(0b000111111), 9), pad(KeyType(0b000111110), 9), 3, 0, 0, -1},
        {pad(KeyType(0b000111111), 9), pad(KeyType(0b001110110), 9), 3, 0, 0, 1},
        // over/underflow tests
        {pad(KeyType(0b100111111), 9), pad(KeyType(0b000011011), 9), 3, 1, 0, 0},  // PBC sensitive
        {pad(KeyType(0b000011011), 9), pad(KeyType(0b100111111), 9), 3, -1, 0, 0}, // PBC sensitive
        {pad(KeyType(0b011), 3), pad(KeyType(0b111), 3), 1, 1, 0, 0},
        {pad(KeyType(0b111), 3), pad(KeyType(0b011), 3), 1, 1, 0, 0},  // PBC sensitive
        {pad(KeyType(0b011), 3), pad(KeyType(0b111), 3), 1, -1, 0, 0}, // PBC sensitive
        // diagonal offset
        {pad(KeyType(0b000111111), 9), pad(KeyType(0b000111000), 9), 3, -1, -1, -1},
        {pad(KeyType(0b000111000), 9), pad(KeyType(0b000111111), 9), 3, 1, 1, 1},
        {pad(KeyType(0b000111111), 9), pad(KeyType(0b111111111), 9), 3, -4, -4, -4}, // PBC sensitive
        {pad(KeyType(0b000111000), 9), pad(KeyType(0b000000000), 9), 3, 6, 6, 6},    // PBC sensitive
    };

    auto computeCode = [](auto t)
    {
        auto [input, reference, level, dx, dy, dz] = t;
        IBox ibox                                  = mortonIBox(input, level);
        return sfcNeighbor<MortonKey<KeyType>>(ibox, level, dx, dy, dz);
    };

    std::vector<KeyType> probes(codes.size());
    std::transform(begin(codes), end(codes), begin(probes), computeCode);

    for (std::size_t i = 0; i < codes.size(); ++i)
    {
        EXPECT_EQ(std::get<1>(codes[i]), probes[i]);
    }
}

TEST(MortonCode, mortonNeighbor)
{
    mortonNeighbors<unsigned>();
    mortonNeighbors<uint64_t>();
}

template<class KeyType>
void mortonIBox()
{
    constexpr unsigned maxCoord = 1u << maxTreeLevel<KeyType>{};
    {
        KeyType nodeStart = nodeRange<KeyType>(0) - 1;
        KeyType nodeEnd   = nodeRange<KeyType>(0);
        IBox ibox         = mortonIBoxKeys(nodeStart, nodeEnd);
        IBox refBox{maxCoord - 1, maxCoord};
        EXPECT_EQ(ibox, refBox);
    }
}

TEST(MortonCode, makeIBox)
{
    mortonIBox<unsigned>();
    mortonIBox<uint64_t>();
}
