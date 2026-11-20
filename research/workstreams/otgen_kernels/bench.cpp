// Standalone OT-generation kernels at libOTe's K = 2^18 SPIN configuration:
// 128 trees x 4096 leaves (N = 2^19 noise blocks), K OTs hashed, SPIN
// transpose of 2K blocks. Validates every variant against scalar references,
// then times each one (median of reps) against its AES-unit floor.
// One process, one thread; run pinned on an idle core, never concurrently.
#include "kernels.h"
#include <spin/Code.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <span>
#include <string>
#include <vector>

using namespace otgen;

template<class F>
static double medianMs(int reps, F&& f) {
    for (int i = 0; i < 3; ++i) f();
    std::vector<double> t(reps);
    for (auto& x : t) {
        const auto a = std::chrono::steady_clock::now();
        f();
        x = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    }
    std::sort(t.begin(), t.end());
    return t[reps / 2];
}

static bool eq(const __m128i* a, const __m128i* b, u64 n) { return !std::memcmp(a, b, n * 16); }
static void fill(__m128i* p, u64 n, u64 s) {
    for (u64 i = 0; i < n; ++i) {
        s += 0x9e3779b97f4a7c15ull; u64 v = s;
        v = (v ^ (v >> 30)) * 0xbf58476d1ce4e5b9ull; v = (v ^ (v >> 27)) * 0x94d049bb133111ebull; v ^= v >> 31;
        p[i] = _mm_set_epi64x(i * 0x5851f42d4c957f2dull ^ v, v);
    }
}
static int failures = 0;
static void check(bool ok, const char* what) { std::printf("  %-44s %s\n", what, ok ? "ok" : "FAIL"); failures += !ok; }

template<int W, Leaf L>
static void validateGgm(u64 depth, u64 trees) {
    const u64 n = trees << depth;
    BlockBuf ref(n), refSums((trees / W) * depth * 2 * W), got(n);
    Ggm<W, L> g(depth, trees, got.p);
    ggmReference(g.key, g.leaf0, g.leaf1, L, g.seed, depth, trees, W, ref.p, refSums.p);
    char name[96];
    g.expandAll();
    std::snprintf(name, sizeof name, "ggm W=%d %s breadth-first", W, L == Leaf::TwoAes ? "2aes" : "aes+round");
    check(eq(ref.p, got.p, n) && eq(refSums.p, g.encSums.p, refSums.n), name);
    std::memset(static_cast<void*>(got.p), 0, n * 16); std::memset(static_cast<void*>(g.encSums.p), 0, refSums.n * 16);
    g.expandBlocked(3);
    std::snprintf(name, sizeof name, "ggm W=%d %s blocked(3)", W, L == Leaf::TwoAes ? "2aes" : "aes+round");
    check(eq(ref.p, got.p, n) && eq(refSums.p, g.encSums.p, refSums.n), name);
    std::memset(static_cast<void*>(got.p), 0, n * 16); std::memset(static_cast<void*>(g.encSums.p), 0, refSums.n * 16);
    g.reset();
    for (u64 i = 0; i < (trees / W) * g.parentsPerGroup(); i += 3) g.step(std::min<u64>(3, (trees / W) * g.parentsPerGroup() - i));
    std::snprintf(name, sizeof name, "ggm W=%d %s resumable step(3)", W, L == Leaf::TwoAes ? "2aes" : "aes+round");
    check(eq(ref.p, got.p, n) && eq(refSums.p, g.encSums.p, refSums.n), name);
}

static void validateHash() {
    const u64 n = 1024;
    HashCtx h;
    BlockBuf r(n), r2(n), ref(2 * n), got(2 * n);
    fill(r.p, n, 5);
    hashReference(h, r.p, ref.p, n);
    std::memcpy(static_cast<void*>(r2.p), r.p, n * 16);
    hashFaithful(h, r2.p, got.y(), n); check(eq(ref.p, got.p, 2 * n), "hash faithful");
    std::memset(static_cast<void*>(got.p), 0, 2 * n * 16);
    hashDirect<8, false>(h, r.p, got.y(), n); check(eq(ref.p, got.p, 2 * n), "hash direct B=8");
    std::memset(static_cast<void*>(got.p), 0, 2 * n * 16);
    hashDirect<16, true>(h, r.p, got.y(), n); _mm_sfence(); check(eq(ref.p, got.p, 2 * n), "hash direct B=16 streaming");
    std::memset(static_cast<void*>(got.p), 0, 2 * n * 16);
    HashCursor c(r.p, got.y(), n);
    for (u64 i = 0; i < n; i += 24) c.step(std::min<u64>(24, n - i));
    check(eq(ref.p, got.p, 2 * n), "hash resumable step(24)");
    const Aes fips(_mm_set_epi64x(0x0f0e0d0c0b0a0908ll, 0x0706050403020100ll));
    const __m128i pt = _mm_set_epi64x(0xffeeddccbbaa9988ll, 0x7766554433221100ll);
    const __m128i ct = _mm_set_epi64x(0x5ac5b47080b7cdd8ll, 0x30047b6ad8e0c469ll);
    __m128i c1 = fips.enc1(pt);
    __m256i y = _mm256_broadcastsi128_si256(pt); enc<1>(fips, &y);
    check(eq(&c1, &ct, 1) && eq(reinterpret_cast<__m128i*>(&y), &ct, 1), "AES-128 FIPS-197 vector (xmm and ymm)");
}

static double aesPeakNsPerBlock() {
    const Aes a(_mm_set_epi64x(1, 2));
    __m256i x[12];
    for (int i = 0; i < 12; ++i) x[i] = _mm256_set1_epi64x(i);
    const long it = 4000000;
    const auto t0 = std::chrono::steady_clock::now();
    for (long i = 0; i < it; ++i) enc<12>(a, x);
    const double ns = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count();
    volatile long long sink = _mm256_extract_epi64(x[0], 0) ^ _mm256_extract_epi64(x[11], 3); (void)sink;
    return ns / (double(it) * 24);
}

int main(int argc, char** argv) {
    const int reps = argc > 1 ? std::atoi(argv[1]) : 21;
    const u64 depth = 12, trees = 128, N = trees << depth, K = N / 2;
    const double ghz = measureGhz();
    const double floorNsPerOp = 1.0 / (4.0 * ghz); // 4 AES block-rounds per cycle
    std::printf("clock %.2f GHz (add chain); AES-unit floor %.3f ns per AES-128 block\n", ghz, 10 * floorNsPerOp);
    std::printf("AES peak, 12 ymm in registers: %.3f ns per block\n", aesPeakNsPerBlock());

    std::printf("validation\n");
    validateGgm<8, Leaf::TwoAes>(6, 32); validateGgm<16, Leaf::TwoAes>(6, 32); validateGgm<32, Leaf::TwoAes>(6, 32);
    validateGgm<8, Leaf::AesRound>(6, 32); validateGgm<16, Leaf::AesRound>(6, 32);
    validateHash();
    if (failures) { std::printf("VALIDATION FAILED (%d)\n", failures); return 1; }

    BlockBuf noise(N);
    std::printf("\nGGM expansion, %llu trees x %llu leaves (N=%llu), median of %d\n",
                (unsigned long long)trees, (unsigned long long)(1ull << depth), (unsigned long long)N, reps);
    std::printf("  %-34s %8s %9s %9s %7s\n", "variant", "ms", "ns/leaf", "floor ms", "%floor");
    auto ggmRow = [&](const char* name, Leaf L, auto&& run) {
        const double ms = medianMs(reps, run);
        const double floorMs = double(N) * aesOpsPerLeaf(L, depth) * floorNsPerOp * 1e-6;
        std::printf("  %-34s %8.3f %9.3f %9.3f %6.1f%%\n", name, ms, ms * 1e6 / double(N), floorMs, 100 * floorMs / ms);
    };
    {
        Ggm<8, Leaf::TwoAes> g(depth, trees, noise.p);
        ggmRow("W=8  2aes breadth-first (libOTe)", Leaf::TwoAes, [&] { g.expandAll(); });
        ggmRow("W=8  2aes blocked sub=6", Leaf::TwoAes, [&] { g.expandBlocked(6); });
        ggmRow("W=8  2aes blocked sub=8", Leaf::TwoAes, [&] { g.expandBlocked(8); });
    }
    {
        Ggm<16, Leaf::TwoAes> g(depth, trees, noise.p);
        ggmRow("W=16 2aes breadth-first", Leaf::TwoAes, [&] { g.expandAll(); });
        ggmRow("W=16 2aes blocked sub=6", Leaf::TwoAes, [&] { g.expandBlocked(6); });
    }
    {
        Ggm<32, Leaf::TwoAes> g(depth, trees, noise.p);
        ggmRow("W=32 2aes breadth-first", Leaf::TwoAes, [&] { g.expandAll(); });
        ggmRow("W=32 2aes blocked sub=5", Leaf::TwoAes, [&] { g.expandBlocked(5); });
    }
    {
        Ggm<8, Leaf::AesRound> g(depth, trees, noise.p);
        ggmRow("W=8  aes+round (UNPROVEN) BF", Leaf::AesRound, [&] { g.expandAll(); });
    }
    {
        Ggm<16, Leaf::AesRound> g(depth, trees, noise.p);
        ggmRow("W=16 aes+round (UNPROVEN) BF", Leaf::AesRound, [&] { g.expandAll(); });
        ggmRow("W=16 aes+round (UNPROVEN) blk6", Leaf::AesRound, [&] { g.expandBlocked(6); });
    }

    std::printf("\nHash, K=%llu OTs (2 AES per OT), median of %d\n", (unsigned long long)K, reps);
    std::printf("  %-34s %8s %9s %9s %7s\n", "variant", "ms", "ns/OT", "floor ms", "%floor");
    {
        HashCtx h;
        BlockBuf r(K), m(2 * K);
        fill(r.p, K, 11);
        const double floorMs = double(K) * 20 * floorNsPerOp * 1e-6;
        auto row = [&](const char* name, auto&& run) {
            const double ms = medianMs(reps, run);
            std::printf("  %-34s %8.3f %9.3f %9.3f %6.1f%%\n", name, ms, ms * 1e6 / double(K), floorMs, 100 * floorMs / ms);
        };
        row("faithful (libOTe VAES loop)", [&] { hashFaithful(h, r.p, m.y(), K); });
        row("direct B=8 cached", [&] { hashDirect<8, false>(h, r.p, m.y(), K); });
        row("direct B=16 cached", [&] { hashDirect<16, false>(h, r.p, m.y(), K); });
        row("direct B=8 streaming", [&] { hashDirect<8, true>(h, r.p, m.y(), K); _mm_sfence(); });
        row("direct B=16 streaming", [&] { hashDirect<16, true>(h, r.p, m.y(), K); _mm_sfence(); });
    }

    std::printf("\nSPIN transpose (library, Automatic backend), 2K -> K, median of %d\n", reps);
    {
        spin::Code code({.message_size = K, .route_seed = 17, .inner_seed = 41});
        auto work = code.make_workspace();
        struct alignas(16) B128 { std::array<std::uint64_t, 2> w; };
        std::vector<B128> buf(2 * K), src(2 * K);
        fill(reinterpret_cast<__m128i*>(src.data()), 2 * K, 3);
        const double ms = medianMs(reps, [&] {
            std::memcpy(buf.data(), src.data(), 2 * K * 16);
            code.transpose_inplace<B128>(std::span<B128>(buf), work);
        });
        const double copyMs = medianMs(reps, [&] { std::memcpy(buf.data(), src.data(), 2 * K * 16); });
        std::printf("  transpose_inplace %.3f ms (refill memcpy %.3f ms subtracted)\n", ms - copyMs, copyMs);
    }
    return 0;
}
