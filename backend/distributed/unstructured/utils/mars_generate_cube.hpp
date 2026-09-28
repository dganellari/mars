#pragma once

// Procedural per-rank structured-cube generation, drop-in for
// readMeshWithElementPartitioning (mars_read_mesh_binary.hpp): each rank computes
// ONLY its element slice of an Ncells^3 hex cube on [0,1]^3 -- no mesh file. This
// is the prerequisite for trillion-DOF runs, where an O(10 TB) mesh file is
// infeasible. The output format is identical to the file reader:
//   - per-rank coords for the nodes this rank's elements reference (local-indexed)
//   - 8 connectivity arrays remapped to those local node indices
// cstone's sync() then merges nodes shared across ranks by SFC key from coords,
// exactly as for a file mesh.
//
// ELEMENT->RANK SPLIT is a 3D CUBIC BLOCK (brick), not a flat-index slab. The rank
// count is factored Px*Py*Pz ~ numRanks and rank r owns a contiguous (ex,ey,ez)
// brick. WHY: a flat contiguous split (rank r owns elements [r*per,(r+1)*per)) is a
// thin SLAB at high rank count -- once per-rank elems < Ncells^2 the slab spans less
// than one ex-layer but still references TWO FULL transverse node planes (2*(N+1)^2),
// a per-rank NODE FLOOR that does NOT shrink as you add ranks and blows the per-GPU
// domain-build ceiling at trillion scale. A near-cubic brick's node count is
// surface-dominated and DOES drop as ranks grow. cstone re-partitions by SFC anyway;
// the initial split only needs to bound the per-rank host+device footprint.
//
// NODE NUMBERING is CLOSED-FORM, no host sort/unique/lower_bound. A brick's referenced
// nodes are exactly the dense lattice block [bx0..bx1]x[by0..by1]x[bz0..bz1], so the
// local node id of a corner is direct arithmetic (decode global -> local) -- O(1) per
// connectivity entry, no temporary multi-GB buffers, no ~100s/rank host sort. This is
// what makes the generator tractable at hundreds of millions of elements per rank.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <tuple>
#include <type_traits>
#include <vector>

namespace mars {

// Opt-in irregularity for the procedural cube. Default OFF -> byte-identical to a
// perfect grid so existing --ncells runs are unaffected. When ON, node coords are
// warped/deformed by a PURE FUNCTION of the global integer index (cx,cy,cz), so a
// node shared across ranks gets identical coords with NO communication -> the mesh
// stays conforming. Two levers:
//   warp:   non-uniform monotone per-axis spacing (spreads element VOLUMES; gives
//           volume-COV irregularity for the FEM/operator path). Monotone -> hexes
//           stay convex, positive Jacobian. NOTE: warp is essentially halo-NEUTRAL --
//           it does not move the SFC partition surface; keep it only if you want the
//           non-uniform volume distribution, not for halo realism.
//   deform: large-wavelength sinusoidal bend of the lattice. THIS is the halo lever.
//           Amplitude is a fraction of the DOMAIN (not of local spacing h) and the
//           wavelength is fixed in DOMAIN units, so the geometric distortion -- and
//           thus the SFC partition/halo bump -- is INDEPENDENT of Ncells. A
//           perturbation tied to local h is sub-cell and the SFC re-finds compact
//           locality at scale (cosmetic at trillion N); a domain-scale bend does not.
//
// VALIDITY (positive Jacobian): inversion is governed by the displacement GRADIENT,
// = deformAmp * 2*pi * deformWaves, NOT by deformAmp alone. Empirically the cube
// stays non-inverted while deformAmp*2*pi*deformWaves <= ~1.0 (clean at 1.0, hexes
// start inverting around 1.25). The default (deformAmp=0.04, deformWaves=3 ->
// product 0.75) is comfortably valid. generateCubeElementPartition asserts this
// bound so a user who raises either field past the cliff fails LOUDLY instead of
// silently producing inverted hexes (the in-memory ctor does not check Jacobian sign).
struct CubeIrregularity {
    bool   warp   = false;  // non-uniform monotone tensor spacing (volume COV; halo-neutral)
    bool   deform = false;  // large-wavelength domain-scale sinusoidal bend (the halo lever)
    double warpStrength   = 0.6;   // in [0,1): 0 = uniform, ->1 = strongly graded
    double deformAmp      = 0.04;  // fraction of the DOMAIN extent (absolute, NOT *h)
    // wavelengths across the domain (domain-scale). Distortion is amp*2*pi*waves and
    // is N-independent, so the halo bump holds at trillion scale. Validity bound:
    // deformAmp * 2*pi * deformWaves <= ~1.0.
    int    deformWaves    = 3;
};

// Monotone per-axis warp on the normalized coordinate t in [0,1]. Smooth-sinh-like
// grading: clusters nodes toward the ends. Strictly increasing for s in [0,1) ->
// preserves node order per axis (conforming, no inversion). s=0 returns t (uniform).
template<typename RealType>
inline RealType warpAxis(RealType t, double s)
{
    if (s <= 0.0) return t;
    const double a = 2.0 * t - 1.0;                 // [-1,1]
    const double w = std::sin(0.5 * M_PI * a);      // monotone S-curve on [-1,1]
    const double blended = (1.0 - s) * a + s * w;   // monotone for s in [0,1)
    return static_cast<RealType>(0.5 * (blended + 1.0));
}

// Factor numRanks into Px*Py*Pz == numRanks, picking the triple of divisors with the
// smallest max/min aspect ratio so each rank's element brick is as cubic as possible
// -> surface-dominated halo. Always exact (Px*Py*Pz == numRanks) since we only search
// over actual divisors. A prime numRanks degenerates to 1x1xnumRanks (a slab), but
// real runs use composite counts (the trillion target 20000 = 25x25x32 is near-cubic);
// ranks are never left idle.
inline std::array<size_t, 3> factorRankGrid(int numRanks)
{
    size_t best[3] = {1, 1, (size_t)numRanks};
    double bestScore = 1e300;
    for (size_t px = 1; px <= (size_t)numRanks; ++px) {
        if (numRanks % (int)px != 0) continue;
        size_t rem = (size_t)numRanks / px;
        for (size_t py = 1; py <= rem; ++py) {
            if (rem % py != 0) continue;
            size_t pz = rem / py;
            // score = aspect spread; lower is more cubic
            size_t mx = std::max({px, py, pz}), mn = std::min({px, py, pz});
            double score = (double)mx / (double)mn;
            if (score < bestScore) { bestScore = score; best[0] = px; best[1] = py; best[2] = pz; }
        }
    }
    return {best[0], best[1], best[2]};
}

// Rank grid for a box of cells[0] x cells[1] x cells[2] hexes: the divisor triple
// Px*Py*Pz == numRanks whose brick (cells per rank along each axis) is the most
// cubic, with no axis split into more pieces than it has cells. A one-cell-thick
// box therefore splits in x and y only.
inline std::array<size_t, 3> factorRankGrid(int numRanks, const std::array<size_t, 3>& cells)
{
    std::array<size_t, 3> best{0, 0, 0};
    double bestScore = 1e300;
    for (size_t px = 1; px <= (size_t)numRanks; ++px) {
        if (numRanks % (int)px != 0) continue;
        const size_t rem = (size_t)numRanks / px;
        for (size_t py = 1; py <= rem; ++py) {
            if (rem % py != 0) continue;
            const size_t p[3] = {px, py, rem / py};
            double lo = 1e300, hi = 0;
            bool fits = true;
            for (int d = 0; d < 3; ++d) {
                if (p[d] > cells[d]) { fits = false; break; }
                const double w = (double)cells[d] / (double)p[d];
                lo = std::min(lo, w);
                hi = std::max(hi, w);
            }
            if (fits && hi / lo < bestScore) { bestScore = hi / lo; best = {p[0], p[1], p[2]}; }
        }
    }
    if (best[0] == 0) {
        std::fprintf(stderr, "factorRankGrid: %d ranks cannot split %zu x %zu x %zu cells\n",
                     numRanks, cells[0], cells[1], cells[2]);
        std::abort();
    }
    return best;
}

// 1D block range [lo,hi) of cell indices for rank coordinate rc of P along an axis of
// N cells. Even split with the remainder spread over the first ranks.
inline void axisBlock(size_t N, size_t P, size_t rc, size_t& lo, size_t& hi)
{
    const size_t base = N / P, extra = N % P;
    lo = rc * base + std::min(rc, extra);
    hi = lo + base + (rc < extra ? 1 : 0);
}

// Returns (nodeCount, elementCount, x, y, z, localConn[8]) for this rank's brick of a
// cells[0] x cells[1] x cells[2] hex lattice split on rank grid P. coord(cx, cy, cz, x,
// y, z) maps a GLOBAL lattice index to its position; it must be a pure function of the
// index so nodes shared by two ranks get identical coords. coords are local-node-indexed;
// localConn holds local node indices in hex corner order matching generate_hex_cube.py:
//   0:(i,j,k) 1:(i+1,j,k) 2:(i+1,j+1,k) 3:(i,j+1,k)
//   4:(i,j,k+1) 5:(i+1,j,k+1) 6:(i+1,j+1,k+1) 7:(i,j+1,k+1)
template<typename RealType, typename KeyType, typename Coord>
inline auto generateBrickPartition(const std::array<size_t, 3>& cells,
                                   const std::array<size_t, 3>& P, int rank, Coord coord)
{
    // Key-range guard: a 32-bit KeyType would silently overflow gid past 4.29e9.
    {
        const long double maxGid = (long double)(cells[0] + 1) * (long double)(cells[1] + 1)
                                 * (long double)(cells[2] + 1);
        const long double keyMax = (long double)std::numeric_limits<
            typename std::make_unsigned<KeyType>::type>::max();
        if (maxGid - 1.0L > keyMax) {
            std::fprintf(stderr,
                "generateBrickPartition: %zu x %zu x %zu cells need node ids up to %.0Lf "
                "but KeyType holds only %.0Lf -- use a wider KeyType.\n",
                cells[0], cells[1], cells[2], maxGid - 1.0L, keyMax);
            std::abort();
        }
    }

    // rank -> (rcx,rcy,rcz). Ranks >= Px*Py*Pz own an empty brick.
    const size_t usedRanks = P[0] * P[1] * P[2];
    size_t bx0 = 0, bx1 = 0, by0 = 0, by1 = 0, bz0 = 0, bz1 = 0;
    if ((size_t)rank < usedRanks) {
        const size_t rcx = (size_t)rank / (P[1] * P[2]);
        const size_t rem = (size_t)rank % (P[1] * P[2]);
        const size_t rcy = rem / P[2];
        const size_t rcz = rem % P[2];
        axisBlock(cells[0], P[0], rcx, bx0, bx1);
        axisBlock(cells[1], P[1], rcy, by0, by1);
        axisBlock(cells[2], P[2], rcz, bz0, bz1);
    }
    const size_t nex = (bx1 > bx0) ? bx1 - bx0 : 0;
    const size_t ney = (by1 > by0) ? by1 - by0 : 0;
    const size_t nez = (bz1 > bz0) ? bz1 - bz0 : 0;
    const size_t elementCount = nex * ney * nez;

    // Node lattice block this brick references: [bx0..bx1] x [by0..by1] x [bz0..bz1]
    // (inclusive on the high end -> nlx = nex+1 node planes per axis).
    const size_t nlx = (nex > 0) ? nex + 1 : 0;
    const size_t nly = (ney > 0) ? ney + 1 : 0;
    const size_t nlz = (nez > 0) ? nez + 1 : 0;
    const size_t nodeCount = nlx * nly * nlz;

    // local node id from a GLOBAL lattice index (cx,cy,cz) inside the brick block.
    auto localNode = [&](size_t cx, size_t cy, size_t cz) -> KeyType {
        return static_cast<KeyType>(((cx - bx0) * nly + (cy - by0)) * nlz + (cz - bz0));
    };

    // --- connectivity: closed-form local ids, no sort/unique/lower_bound ---
    std::array<std::vector<KeyType>, 8> lconn;
    for (auto& v : lconn) v.resize(elementCount);
    {
        size_t li = 0;
        for (size_t ex = bx0; ex < bx1; ++ex)
            for (size_t ey = by0; ey < by1; ++ey)
                for (size_t ez = bz0; ez < bz1; ++ez, ++li) {
                    lconn[0][li] = localNode(ex,     ey,     ez);
                    lconn[1][li] = localNode(ex + 1, ey,     ez);
                    lconn[2][li] = localNode(ex + 1, ey + 1, ez);
                    lconn[3][li] = localNode(ex,     ey + 1, ez);
                    lconn[4][li] = localNode(ex,     ey,     ez + 1);
                    lconn[5][li] = localNode(ex + 1, ey,     ez + 1);
                    lconn[6][li] = localNode(ex + 1, ey + 1, ez + 1);
                    lconn[7][li] = localNode(ex,     ey + 1, ez + 1);
                }
    }

    std::vector<RealType> x(nodeCount), y(nodeCount), z(nodeCount);
    {
        size_t i = 0;
        for (size_t cx = bx0; cx < bx0 + nlx; ++cx)
            for (size_t cy = by0; cy < by0 + nly; ++cy)
                for (size_t cz = bz0; cz < bz0 + nlz; ++cz, ++i)
                    coord(cx, cy, cz, x[i], y[i], z[i]);
    }

    return std::make_tuple(nodeCount, elementCount, std::move(x), std::move(y), std::move(z), std::move(lconn));
}

// This rank's CUBIC BRICK of an Ncells^3 hex cube on [0,1]^3 (rank grid from
// factorRankGrid(numRanks)); see generateBrickPartition for the returned tuple.
template<typename RealType, typename KeyType>
inline auto generateCubeElementPartition(size_t Ncells, int rank, int numRanks,
                                         CubeIrregularity irr = {})
{
    // Validity guard: deform distortion is amp*2*pi*waves; inverts past ~1.0.
    if (irr.deform) {
        const double distortion = irr.deformAmp * 2.0 * M_PI * (double)irr.deformWaves;
        if (distortion > 1.0) {
            std::fprintf(stderr,
                "generateCubeElementPartition: deform distortion deformAmp*2pi*deformWaves "
                "= %.3f > 1.0 -> inverted hexes. Reduce deformAmp (%.3f) or deformWaves (%d).\n",
                distortion, irr.deformAmp, irr.deformWaves);
            std::abort();
        }
    }

    const RealType h = RealType(1) / static_cast<RealType>(Ncells);
    const bool irregular = irr.warp || irr.deform;
    const double twoPi = 2.0 * M_PI;
    auto coord = [&](size_t cx, size_t cy, size_t cz, RealType& x, RealType& y, RealType& z) {
        if (!irregular) {
            x = static_cast<RealType>(cx) * h;
            y = static_cast<RealType>(cy) * h;
            z = static_cast<RealType>(cz) * h;
            return;
        }
        RealType tx = static_cast<RealType>(cx) / static_cast<RealType>(Ncells);
        RealType ty = static_cast<RealType>(cy) / static_cast<RealType>(Ncells);
        RealType tz = static_cast<RealType>(cz) / static_cast<RealType>(Ncells);
        if (irr.warp) {
            tx = warpAxis<RealType>(tx, irr.warpStrength);
            ty = warpAxis<RealType>(ty, irr.warpStrength);
            tz = warpAxis<RealType>(tz, irr.warpStrength);
        }
        RealType px = tx, py = ty, pz = tz;
        if (irr.deform) {
            // domain-scale bend: amplitude is an absolute fraction of the
            // domain, wavelength fixed in domain units -> distortion (and
            // halo bump) independent of Ncells. Each axis displaced by the
            // OTHER two axes' phase so the deformation is genuinely 3D.
            const double A  = irr.deformAmp;
            const double kw = static_cast<double>(irr.deformWaves);
            px += static_cast<RealType>(A * std::sin(twoPi * kw * ty) * std::cos(twoPi * kw * tz));
            py += static_cast<RealType>(A * std::sin(twoPi * kw * tz) * std::cos(twoPi * kw * tx));
            pz += static_cast<RealType>(A * std::sin(twoPi * kw * tx) * std::cos(twoPi * kw * ty));
        }
        x = px; y = py; z = pz;
    };
    return generateBrickPartition<RealType, KeyType>({Ncells, Ncells, Ncells},
                                                     factorRankGrid(numRanks), rank, coord);
}

// This rank's brick of a box [lo, hi] meshed with cells[0] x cells[1] x cells[2] hexes,
// uniform in x and z. yGrading in [0,1) clusters the y spacing toward both y faces
// (warpAxis); 0 is uniform.
template<typename RealType, typename KeyType>
inline auto generateBoxElementPartition(const std::array<size_t, 3>& cells,
                                        const std::array<double, 3>& lo,
                                        const std::array<double, 3>& hi,
                                        double yGrading, int rank, int numRanks)
{
    auto coord = [&](size_t cx, size_t cy, size_t cz, RealType& x, RealType& y, RealType& z) {
        const double ty = warpAxis<double>((double)cy / (double)cells[1], yGrading);
        x = static_cast<RealType>(lo[0] + (hi[0] - lo[0]) * ((double)cx / (double)cells[0]));
        y = static_cast<RealType>(lo[1] + (hi[1] - lo[1]) * ty);
        z = static_cast<RealType>(lo[2] + (hi[2] - lo[2]) * ((double)cz / (double)cells[2]));
    };
    return generateBrickPartition<RealType, KeyType>(cells, factorRankGrid(numRanks, cells), rank, coord);
}

} // namespace mars
