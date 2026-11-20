// Fused feasibility: SPIN's AVX-512 transpose (patched copies of RouteFast.cpp
// and Fast.cpp with an inlined hook per route epoch and per outer group) runs
// chunks of real GGM-expansion or hash work in the gaps. Compares
//   T_fused  : hooked transpose doing `work`
//   T_spin   : same binary, hook mode 0
//   T_work   : the identical chunk sequence run without SPIN
// hidden = T_spin + T_work - T_fused is the time the two overlapped.
// One process, one thread; run pinned on an idle core, never concurrently.
#include "otgen_hook.h"
#include "rs8/Rs8.h"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace otgen;
namespace rs8 = spin::detail::rs8;

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
static void fill(__m128i* p, u64 n, u64 s) {
    for (u64 i = 0; i < n; ++i) {
        s += 0x9e3779b97f4a7c15ull; u64 v = s;
        v = (v ^ (v >> 30)) * 0xbf58476d1ce4e5b9ull; v = (v ^ (v >> 27)) * 0x94d049bb133111ebull; v ^= v >> 31;
        p[i] = _mm_set_epi64x(i ^ v, v);
    }
}

int main(int argc, char** argv) {
    const int reps = argc > 1 ? std::atoi(argv[1]) : 21;
    const u64 depth = 12, trees = 128, N = trees << depth, K = N / 2;
    std::printf("clock %.2f GHz\n", measureGhz());
    rs8::Plan plan(K, 17, 41);
    BlockBuf in(N), src(N), out(K), ref(K), scratch(plan.scratchBlocks()), ggmOut(N), hashIn(K), hashOut(2 * K);
    fill(src.p, N, 3); fill(hashIn.p, K, 9);
    auto* inB = reinterpret_cast<rs8::Block*>(in.p);
    auto* outB = reinterpret_cast<rs8::Block*>(out.p);
    auto* scrB = reinterpret_cast<rs8::Block*>(scratch.p);
    auto run = [&] { rs8::transposeFast(inB, outB, scrB, plan); };

    FusedGgm ggm(depth, trees, ggmOut.p);
    HashCursor hash(hashIn.p, hashOut.y(), K);
    g_hook.ggm = &ggm; g_hook.hash = &hash;

    // Validation: hooked transpose equals the scalar oracle, with and without work.
    std::memcpy(static_cast<void*>(in.p), src.p, N * 16);
    rs8::transposeScalar(reinterpret_cast<rs8::Block*>(in.p), reinterpret_cast<rs8::Block*>(ref.p), scrB, plan);
    for (int mode = 0; mode < 3; ++mode) {
        g_hook = HookState{mode, 4, 16, 0, 0, &ggm, &hash};
        std::memcpy(static_cast<void*>(in.p), src.p, N * 16);
        run();
        const bool ok = !std::memcmp(out.p, ref.p, K * 16);
        std::printf("validation mode %d: transpose %s\n", mode, ok ? "ok" : "FAIL");
        if (!ok) return 1;
    }
    std::memcpy(static_cast<void*>(in.p), src.p, N * 16);

    g_hook = HookState{0, 0, 0, 0, 0, &ggm, &hash};
    run();
    const u64 routeCalls = g_hook.routeCalls, outerCalls = g_hook.outerCalls;
    const double tSpin = medianMs(reps, run);
    std::printf("hooks per transpose: %llu route epochs, %llu outer groups\n",
                (unsigned long long)routeCalls, (unsigned long long)outerCalls);
    std::printf("T_spin (hook mode 0) %.3f ms\n\n", tSpin);

    std::printf("%-6s %6s %6s %12s %9s %9s %9s %9s %8s %8s\n", "work", "route", "outer", "units/call",
                "T_work", "T_fused", "serial", "hidden", "%spin", "%work");
    auto sweep = [&](int mode, const char* name, u64 rc, u64 oc) {
        g_hook = HookState{mode, rc, oc, 0, 0, &ggm, &hash};
        ggm.reset(); hash.i = 0;
        const u64 units0 = mode == 1 ? ggm.parents : hash.done;
        const double tFused = medianMs(reps, run);
        const u64 calls = g_hook.routeCalls / routeCalls;
        const u64 units = ((mode == 1 ? ggm.parents : hash.done) - units0) / calls;
        g_hook.mode = 0;
        ggm.reset(); hash.i = 0;
        const double tWork = medianMs(reps, [&] {
            for (u64 i = 0; i < routeCalls; ++i) { if (mode == 1) ggm.step(rc); else hash.step(rc); }
            for (u64 i = 0; i < outerCalls; ++i) { if (mode == 1) ggm.step(oc); else hash.step(oc); }
        });
        const double serial = tSpin + tWork, hidden = serial - tFused;
        std::printf("%-6s %6llu %6llu %12llu %9.3f %9.3f %9.3f %9.3f %7.1f%% %7.1f%%\n", name,
                    (unsigned long long)rc, (unsigned long long)oc, (unsigned long long)units,
                    tWork, tFused, serial, hidden, 100 * hidden / tSpin, 100 * hidden / tWork);
    };
    // GGM: 16 groups x 4095 parents = 65,520 parents is one full K=2^18 expansion.
    for (u64 rc : {1, 2, 4, 8}) sweep(1, "ggm", rc, 4 * rc);
    // Hash: K = 262,144 OTs per batch.
    for (u64 rc : {8, 16, 32, 64}) sweep(2, "hash", rc, 4 * rc);
    // Diagnostic: identical AES work on a 512-OT (8 KiB in, 16 KiB out) L1-resident
    // buffer. Same execution-unit pressure, almost no memory traffic.
    HashCursor small(hashIn.p, hashOut.y(), 512);
    g_hook.hash = &small;
    auto sweepSmall = [&](u64 rc) {
        g_hook = HookState{2, rc, 4 * rc, 0, 0, &ggm, &small};
        small.i = 0;
        const double tFused = medianMs(reps, run);
        g_hook.mode = 0;
        const double tWork = medianMs(reps, [&] {
            for (u64 i = 0; i < routeCalls; ++i) small.step(rc);
            for (u64 i = 0; i < outerCalls; ++i) small.step(4 * rc);
        });
        const double serial = tSpin + tWork, hidden = serial - tFused;
        std::printf("%-6s %6llu %6llu %12llu %9.3f %9.3f %9.3f %9.3f %7.1f%% %7.1f%%\n", "hashL1",
                    (unsigned long long)rc, (unsigned long long)(4 * rc),
                    (unsigned long long)(rc * routeCalls + 4 * rc * outerCalls),
                    tWork, tFused, serial, hidden, 100 * hidden / tSpin, 100 * hidden / tWork);
    };
    for (u64 rc : {8, 16, 32}) sweepSmall(rc);
    return 0;
}
