// Realistic micro-kernels for silent-OT generation (sender, regular noise,
// GF(2) correlations): GGM PPRF expansion and the correlation-robust hash.
// The GGM and hash formulas match libOTe (RegularPprf.h, SilentOtExtSender.cpp);
// the kernels differ only in batching and traversal, which never change output
// except where a variant is explicitly marked UNPROVEN.
//
// AES runs on VAES ymm registers: each ymm holds two 128-bit blocks, so a node
// of W interleaved trees occupies W/2 ymm. AES-128 costs 10 AES-unit
// instructions per block (9 aesenc + aesenclast); the round-key whitening is a
// vpxor. Zen 4 issues 2 VAES ymm per cycle, i.e. 4 block-rounds per cycle.
#pragma once
#include <immintrin.h>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <new>
#include <utility>

namespace otgen {
using u64 = std::uint64_t;
#define OTGEN_INLINE inline __attribute__((always_inline))
#define OTGEN_NOINLINE __attribute__((noinline))

// ---------------------------------------------------------------- storage
struct BlockBuf {
    __m128i* p = nullptr;
    std::size_t n = 0;
    BlockBuf() = default;
    explicit BlockBuf(std::size_t blocks) : n(blocks) {
        p = static_cast<__m128i*>(::operator new[](blocks * 16, std::align_val_t{64}));
        std::memset(static_cast<void*>(p), 0, blocks * 16);
    }
    BlockBuf(BlockBuf&& o) noexcept : p(std::exchange(o.p, nullptr)), n(o.n) {}
    BlockBuf& operator=(BlockBuf&& o) noexcept { std::swap(p, o.p); std::swap(n, o.n); return *this; }
    BlockBuf(const BlockBuf&) = delete;
    ~BlockBuf() { if (p) ::operator delete[](p, std::align_val_t{64}); }
    __m256i* y() const { return reinterpret_cast<__m256i*>(p); }
};

// ---------------------------------------------------------------- AES-128
inline __m128i expandStep(__m128i k, __m128i g) {
    g = _mm_shuffle_epi32(g, 0xff);
    k = _mm_xor_si128(k, _mm_slli_si128(k, 4));
    k = _mm_xor_si128(k, _mm_slli_si128(k, 4));
    k = _mm_xor_si128(k, _mm_slli_si128(k, 4));
    return _mm_xor_si128(k, g);
}
struct Aes {
    __m128i k[11];
    __m256i k2[11];
    explicit Aes(__m128i key) {
        k[0] = key;
#define OTGEN_EXP(i, rc) k[i] = expandStep(k[i - 1], _mm_aeskeygenassist_si128(k[i - 1], rc))
        OTGEN_EXP(1, 0x01); OTGEN_EXP(2, 0x02); OTGEN_EXP(3, 0x04); OTGEN_EXP(4, 0x08);
        OTGEN_EXP(5, 0x10); OTGEN_EXP(6, 0x20); OTGEN_EXP(7, 0x40); OTGEN_EXP(8, 0x80);
        OTGEN_EXP(9, 0x1b); OTGEN_EXP(10, 0x36);
#undef OTGEN_EXP
        for (int i = 0; i < 11; ++i) k2[i] = _mm256_broadcastsi128_si256(k[i]);
    }
    __m128i enc1(__m128i x) const { // scalar reference
        x = _mm_xor_si128(x, k[0]);
        for (int r = 1; r < 10; ++r) x = _mm_aesenc_si128(x, k[r]);
        return _mm_aesenclast_si128(x, k[10]);
    }
};
template<int N>
OTGEN_INLINE void enc(const Aes& a, __m256i* x) {
#pragma GCC unroll 64
    for (int i = 0; i < N; ++i) x[i] = _mm256_xor_si256(x[i], a.k2[0]);
#pragma GCC unroll 9
    for (int r = 1; r < 10; ++r) {
#pragma GCC unroll 64
        for (int i = 0; i < N; ++i) x[i] = _mm256_aesenc_epi128(x[i], a.k2[r]);
    }
#pragma GCC unroll 64
    for (int i = 0; i < N; ++i) x[i] = _mm256_aesenclast_epi128(x[i], a.k2[10]);
}
// Two keys interleaved for ILP: x under a, y under b.
template<int N>
OTGEN_INLINE void enc2(const Aes& a, const Aes& b, __m256i* x, __m256i* y) {
#pragma GCC unroll 64
    for (int i = 0; i < N; ++i) { x[i] = _mm256_xor_si256(x[i], a.k2[0]); y[i] = _mm256_xor_si256(y[i], b.k2[0]); }
#pragma GCC unroll 9
    for (int r = 1; r < 10; ++r) {
#pragma GCC unroll 64
        for (int i = 0; i < N; ++i) {
            x[i] = _mm256_aesenc_epi128(x[i], a.k2[r]);
            y[i] = _mm256_aesenc_epi128(y[i], b.k2[r]);
        }
    }
#pragma GCC unroll 64
    for (int i = 0; i < N; ++i) {
        x[i] = _mm256_aesenclast_epi128(x[i], a.k2[10]);
        y[i] = _mm256_aesenclast_epi128(y[i], b.k2[10]);
    }
}

// ---------------------------------------------------------------- GGM
// libOTe RegularPprfSender::expandOne, GF(2) blocks, interleaved output.
// Internal level, per parent x (fixed key K):
//     a = AES_K(x);  child0 = aesenc(a, x);  child1 = a +_64 x
// Leaf level, per parent x (fixed keys K0, K1):
//     TwoAes:   leaf_keep = AES_Kkeep(x) ^ x                  (libOTe)
//     AesRound: leaf0 = aesenc(AES_K0(x), x); leaf1 = AES_K0(x) +_64 x
//               UNPROVEN: the internal-level PRG reused at the leaves. Its
//               leaves become OT correlations; needs a security argument
//               before any deployment. Measured only to bound the speed-up.
// Every level accumulates left/right sums per tree (sent to the receiver).
// Output: tree t (of 8 in a sub-group) leaf c is at subgroupBase + 8c + t.
enum class Leaf { TwoAes, AesRound };

// Rounds on the AES units per noise element (leaf), for the floor.
constexpr double aesOpsPerLeaf(Leaf l, u64 depth) {
    // internal parents per tree: 2^(depth-1)-1 at 11 ops; leaf parents 2^(depth-1).
    const double half = double(u64(1) << (depth - 1));
    const double leafOps = l == Leaf::TwoAes ? 20.0 : 11.0;
    return ((half - 1) * 11.0 + half * leafOps) / (2.0 * half);
}

template<int V>
OTGEN_INLINE void internalRange(const Aes& key, const __m256i* __restrict in, __m256i* __restrict ch,
                                u64 n, __m256i* s0m, __m256i* s1m) {
    __m256i s0[V], s1[V];
#pragma GCC unroll 64
    for (int v = 0; v < V; ++v) { s0[v] = s0m[v]; s1[v] = s1m[v]; }
    for (u64 p = 0; p < n; ++p) {
        __m256i x[V], a[V];
#pragma GCC unroll 64
        for (int v = 0; v < V; ++v) { x[v] = _mm256_load_si256(in + p * V + v); a[v] = x[v]; }
        enc<V>(key, a);
#pragma GCC unroll 64
        for (int v = 0; v < V; ++v) {
            const __m256i c0 = _mm256_aesenc_epi128(a[v], x[v]);
            const __m256i c1 = _mm256_add_epi64(a[v], x[v]);
            _mm256_store_si256(ch + (2 * p) * V + v, c0);
            _mm256_store_si256(ch + (2 * p + 1) * V + v, c1);
            s0[v] = _mm256_xor_si256(s0[v], c0);
            s1[v] = _mm256_xor_si256(s1[v], c1);
        }
    }
#pragma GCC unroll 64
    for (int v = 0; v < V; ++v) { s0m[v] = s0[v]; s1m[v] = s1[v]; }
}

// in: n leaf parents, global leaf-parent index p0.. ; group: the W-tree group's
// first output block (sub-group s starts at group + 8*s*domain).
template<int V, Leaf L, bool Stream = false>
OTGEN_INLINE void leafRange(const Aes& k0, const Aes& k1, const __m256i* __restrict in, u64 n, u64 p0,
                            __m128i* __restrict group, u64 domain, __m256i* s0m, __m256i* s1m) {
    __m256i s0[V], s1[V];
#pragma GCC unroll 64
    for (int v = 0; v < V; ++v) { s0[v] = s0m[v]; s1[v] = s1m[v]; }
    for (u64 p = 0; p < n; ++p) {
        __m256i x[V], a[V], b[V];
#pragma GCC unroll 64
        for (int v = 0; v < V; ++v) { x[v] = _mm256_load_si256(in + p * V + v); a[v] = x[v]; b[v] = x[v]; }
        if constexpr (L == Leaf::TwoAes) {
            if constexpr (V <= 8) enc2<V>(k0, k1, a, b); // both keys interleaved: 2V states
            else { enc<V>(k0, a); enc<V>(k1, b); }       // 32 states would not fit in registers
        } else enc<V>(k0, a);
        const u64 c = 2 * (p0 + p);
#pragma GCC unroll 64
        for (int v = 0; v < V; ++v) {
            __m256i l0, l1;
            if constexpr (L == Leaf::TwoAes) {
                l0 = _mm256_xor_si256(a[v], x[v]);
                l1 = _mm256_xor_si256(b[v], x[v]);
            } else {
                l0 = _mm256_aesenc_epi128(a[v], x[v]);
                l1 = _mm256_add_epi64(a[v], x[v]);
            }
            __m128i* sub = group + 8 * (v / 4) * domain;
            auto* d0 = reinterpret_cast<__m256i*>(sub + 8 * c) + (v % 4);
            auto* d1 = reinterpret_cast<__m256i*>(sub + 8 * (c + 1)) + (v % 4);
            if constexpr (Stream) { _mm256_stream_si256(d0, l0); _mm256_stream_si256(d1, l1); }
            else { _mm256_store_si256(d0, l0); _mm256_store_si256(d1, l1); }
            s0[v] = _mm256_xor_si256(s0[v], l0);
            s1[v] = _mm256_xor_si256(s1[v], l1);
        }
    }
#pragma GCC unroll 64
    for (int v = 0; v < V; ++v) { s0m[v] = s0[v]; s1m[v] = s1[v]; }
}

// Resumable expander: step(n) processes the next n parents (internal or leaf)
// in libOTe's breadth-first order, group after group, wrapping around. Used
// whole for the standalone kernel and in chunks for the fused experiment.
template<int W, Leaf L>
struct Ggm {
    static_assert(W % 8 == 0 && W <= 32);
    static constexpr int V = W / 2;
    u64 depth, trees, domain, groups;
    Aes key, leaf0, leaf1;
    __m128i seed;
    BlockBuf levels[24];
    __m128i* out;           // trees*domain blocks
    BlockBuf encSums;       // per group, per level: 2*W blocks (sent values)
    alignas(32) __m256i sums[2][V];
    u64 group = 0, d = 0, p = 0, parents = 0;

    Ggm(u64 depth_, u64 trees_, __m128i* out_)
        : depth(depth_), trees(trees_), domain(u64(1) << depth_), groups(trees_ / W),
          key(_mm_set_epi64x(0x1234, 0x5678)), leaf0(_mm_set_epi64x(0x9abc, 0xdef0)),
          leaf1(_mm_set_epi64x(0x1111, 0x2222)), seed(_mm_set_epi64x(77, 99)), out(out_),
          encSums(groups * depth * 2 * W) {
        for (u64 l = 0; l < depth; ++l) levels[l] = BlockBuf((u64(1) << l) * W);
        reset();
    }
    void reset() { group = 0; parents = 0; beginGroup(); }
    u64 parentsPerGroup() const { return domain - 1; }
    void beginGroup() {
        d = 0; p = 0;
        std::memset(static_cast<void*>(sums), 0, sizeof(sums));
        __m256i r[V];
        for (int v = 0; v < V; ++v) {
            const u64 t = group * W + 2 * v;
            r[v] = _mm256_xor_si256(_mm256_broadcastsi128_si256(seed), _mm256_set_epi64x(0, t + 1, 0, t));
        }
        enc<V>(key, r);
        for (int v = 0; v < V; ++v) _mm256_store_si256(levels[0].y() + v, r[v]);
    }
    OTGEN_NOINLINE void finishLevel() {
        auto* e = encSums.y() + ((group * depth + d) * 2) * V;
        for (int v = 0; v < V; ++v) { _mm256_store_si256(e + v, sums[0][v]); _mm256_store_si256(e + V + v, sums[1][v]); }
        std::memset(static_cast<void*>(sums), 0, sizeof(sums));
        p = 0;
        if (++d == depth) { group = group + 1 == groups ? 0 : group + 1; beginGroup(); }
    }
    OTGEN_INLINE void step(u64 n) {
        parents += n;
        while (n) {
            const u64 width = u64(1) << d;
            const u64 lim = n < width - p ? n : width - p;
            if (d + 1 < depth)
                internalRange<V>(key, levels[d].y() + p * V, levels[d + 1].y() + 2 * p * V, lim, sums[0], sums[1]);
            else
                leafRange<V, L>(leaf0, leaf1, levels[d].y() + p * V, lim, p, out + group * W * domain, domain,
                                sums[0], sums[1]);
            p += lim; n -= lim;
            if (p == width) finishLevel();
        }
    }
    void expandAll() { reset(); step(groups * parentsPerGroup()); }

    // Blocked traversal (output-identical): breadth-first down to level c, then
    // each level-c node's subtree breadth-first in L1-sized scratch.
    BlockBuf local[24], lsum;
    u64 localSub = 0;
    void expandBlocked(u64 sub) {
        if (sub >= depth) { expandAll(); return; }
        const u64 c = depth - sub;
        if (localSub != sub) { // allocate once, outside timed repetitions
            for (u64 j = 1; j < sub; ++j) local[j] = BlockBuf((u64(1) << j) * W);
            lsum = BlockBuf(depth * 2 * W);
            localSub = sub;
        }
        auto* ls = lsum.y();
        for (group = 0; group < groups; ++group) {
            beginGroup();
            std::memset(static_cast<void*>(ls), 0, depth * 2 * W * 16);
            for (u64 l = 0; l < c; ++l)
                internalRange<V>(key, levels[l].y(), levels[l + 1].y(), u64(1) << l, ls + 2 * l * V, ls + (2 * l + 1) * V);
            for (u64 q = 0; q < (u64(1) << c); ++q) {
                const __m256i* root = levels[c].y() + q * V;
                for (u64 j = 0; j + 1 < sub; ++j) {
                    const __m256i* in = j ? local[j].y() : root;
                    internalRange<V>(key, in, local[j + 1].y(), u64(1) << j, ls + 2 * (c + j) * V, ls + (2 * (c + j) + 1) * V);
                }
                const u64 lp = u64(1) << (sub - 1);
                const __m256i* in = sub > 1 ? local[sub - 1].y() : root;
                leafRange<V, L>(leaf0, leaf1, in, lp, q * lp, out + group * W * domain, domain,
                                ls + 2 * (depth - 1) * V, ls + (2 * (depth - 1) + 1) * V);
            }
            std::memcpy(static_cast<void*>(encSums.p + group * depth * 2 * W), ls, depth * 2 * W * 16);
        }
        group = 0;
    }
};

// Scalar reference of the same map (libOTe formulas), for validation.
inline void ggmReference(const Aes& key, const Aes& k0, const Aes& k1, Leaf leaf, __m128i seed,
                         u64 depth, u64 trees, u64 W, __m128i* out, __m128i* encSums) {
    const u64 domain = u64(1) << depth;
    for (u64 t = 0; t < trees; ++t) {
        const u64 g = t / W, j = t % W;
        __m128i* lvl = new __m128i[domain];
        __m128i* nxt = new __m128i[domain];
        lvl[0] = key.enc1(_mm_xor_si128(seed, _mm_set_epi64x(0, t)));
        for (u64 d = 0; d < depth; ++d) {
            __m128i s0 = _mm_setzero_si128(), s1 = _mm_setzero_si128();
            for (u64 p = 0; p < (u64(1) << d); ++p) {
                const __m128i x = lvl[p];
                __m128i c0, c1;
                if (d + 1 < depth || leaf == Leaf::AesRound) {
                    const __m128i a = (d + 1 < depth ? key : k0).enc1(x);
                    c0 = _mm_aesenc_si128(a, x); c1 = _mm_add_epi64(a, x);
                } else {
                    c0 = _mm_xor_si128(k0.enc1(x), x); c1 = _mm_xor_si128(k1.enc1(x), x);
                }
                nxt[2 * p] = c0; nxt[2 * p + 1] = c1;
                s0 = _mm_xor_si128(s0, c0); s1 = _mm_xor_si128(s1, c1);
            }
            // encSums layout of Ggm<W>: [group][level][side][W blocks]
            encSums[((g * depth + d) * 2 + 0) * W + j] = s0;
            encSums[((g * depth + d) * 2 + 1) * W + j] = s1;
            std::swap(lvl, nxt);
        }
        const u64 sg = t / 8, tt = t % 8;
        for (u64 c = 0; c < domain; ++c) out[sg * 8 * domain + 8 * c + tt] = lvl[c];
        delete[] lvl; delete[] nxt;
    }
}

// ---------------------------------------------------------------- hash
// libOTe SilentOtExtSender::hash (ChoiceBitPacking::True), VAES path:
//   r_i &= ~1;  m_i = (H(r_i), H(r_i ^ d)),  H(x) = AES_K(x) ^ x,  d = Delta & ~1.
// Output pairs are adjacent, so one ymm per OT. 20 AES-unit ops per OT.
struct HashCtx {
    Aes key{_mm_set_epi64x(0x4242, 0x4343)};
    __m128i delta = _mm_set_epi64x(0x0123456789abcdefll, 0x0fedcba987654320ll);
};

// Faithful to the libOTe loop: masks r in place, packs pairs through an
// aligned stack batch, hashes 16 blocks, copies the batch out.
OTGEN_NOINLINE inline void hashFaithful(const HashCtx& h, __m128i* r, __m256i* m, u64 n) {
    const __m128i mask = _mm_set_epi64x(-1, -2);
    const __m128i d = _mm_and_si128(h.delta, mask);
    alignas(32) __m128i batch[16];
    for (u64 i = 0; i + 8 <= n; i += 8) {
        for (int j = 0; j < 8; ++j) r[i + j] = _mm_and_si128(r[i + j], mask);
        for (int j = 0; j < 8; ++j) {
            const __m256i pair = _mm256_set_m128i(_mm_xor_si128(r[i + j], d), r[i + j]);
            std::memcpy(batch + 2 * j, &pair, 32);
        }
        __m256i x[8], a[8];
        for (int j = 0; j < 8; ++j) { x[j] = _mm256_load_si256(reinterpret_cast<__m256i*>(batch) + j); a[j] = x[j]; }
        enc<8>(h.key, a);
        for (int j = 0; j < 8; ++j) _mm256_store_si256(reinterpret_cast<__m256i*>(batch) + j, _mm256_xor_si256(a[j], x[j]));
        std::memcpy(static_cast<void*>(m + i), batch, sizeof(batch));
    }
}

// Direct: vbroadcasti128 builds (r, r) on the load port; xor/and form the pair;
// B OTs per iteration; result stored straight to m. Does not write r back.
template<int B, bool Stream>
OTGEN_INLINE void hashDirect(const HashCtx& h, const __m128i* __restrict r, __m256i* __restrict m, u64 n) {
    const __m256i mask = _mm256_set_epi64x(-1, -2, -1, -2);
    const __m256i dpair = _mm256_set_m128i(h.delta, _mm_setzero_si128());
    for (u64 i = 0; i < n; i += B) {
        __m256i x[B], a[B];
#pragma GCC unroll 32
        for (int j = 0; j < B; ++j) {
            x[j] = _mm256_and_si256(_mm256_xor_si256(_mm256_broadcastsi128_si256(_mm_load_si128(r + i + j)), dpair), mask);
            a[j] = x[j];
        }
        enc<B>(h.key, a);
#pragma GCC unroll 32
        for (int j = 0; j < B; ++j) {
            const __m256i o = _mm256_xor_si256(a[j], x[j]);
            if constexpr (Stream) _mm256_stream_si256(m + i + j, o);
            else _mm256_store_si256(m + i + j, o);
        }
    }
}

struct HashCursor {
    HashCtx h;
    const __m128i* r; __m256i* m; u64 n, i = 0, done = 0;
    HashCursor(const __m128i* r_, __m256i* m_, u64 n_) : r(r_), m(m_), n(n_) {}
    OTGEN_INLINE void step(u64 k) { // k multiple of 8
        done += k;
        while (k) {
            const u64 lim = k < n - i ? k : n - i;
            hashDirect<8, false>(h, r + i, m + i, lim);
            i += lim; k -= lim;
            if (i == n) i = 0;
        }
    }
};

inline void hashReference(const HashCtx& h, const __m128i* r, __m128i* m, u64 n) {
    const __m128i mask = _mm_set_epi64x(-1, -2);
    const __m128i d = _mm_and_si128(h.delta, mask);
    for (u64 i = 0; i < n; ++i) {
        const __m128i x0 = _mm_and_si128(r[i], mask), x1 = _mm_xor_si128(x0, d);
        m[2 * i] = _mm_xor_si128(h.key.enc1(x0), x0);
        m[2 * i + 1] = _mm_xor_si128(h.key.enc1(x1), x1);
    }
}

// ---------------------------------------------------------------- clock
inline double measureGhz() { // 4 dependent 1-cycle adds per iteration
    long n = 200000000; std::uint64_t v = 0;
    struct timespec a, b;
    clock_gettime(CLOCK_MONOTONIC, &a);
    asm volatile("1: add $1,%0\n add $1,%0\n add $1,%0\n add $1,%0\n dec %1\n jnz 1b" : "+r"(v), "+r"(n));
    clock_gettime(CLOCK_MONOTONIC, &b);
    const double ns = double(b.tv_sec - a.tv_sec) * 1e9 + double(b.tv_nsec - a.tv_nsec);
    return 4.0 * 200000000 / ns;
}
} // namespace otgen
