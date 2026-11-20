// Hook injected into patched copies of SPIN's AVX-512 transpose loops (route
// epoch loop and outer group loop; see patch_spin.cmake). Mode 0 does nothing,
// so the same binary measures SPIN alone. The kernels inline into SPIN's loops
// so out-of-order execution can overlap them with SPIN's own instructions.
#pragma once
#include "kernels.h"

namespace otgen {
using FusedGgm = Ggm<8, Leaf::TwoAes>;
struct HookState {
    int mode = 0;               // 0 none, 1 GGM parents, 2 hash OTs
    u64 routeChunk = 0, outerChunk = 0;
    u64 routeCalls = 0, outerCalls = 0;
    FusedGgm* ggm = nullptr;
    HashCursor* hash = nullptr;
};
inline HookState g_hook;

// OTGEN_HOOK_OUTLINE=1 keeps the work in a separate function (a call per hook)
// so SPIN's loop keeps its register allocation; 0 inlines it into SPIN's loop.
#ifndef OTGEN_HOOK_OUTLINE
#define OTGEN_HOOK_OUTLINE 0
#endif
#if OTGEN_HOOK_OUTLINE
OTGEN_NOINLINE inline void hookWork(u64 n) {
    if (g_hook.mode == 1) g_hook.ggm->step(n);
    else g_hook.hash->step(n);
}
OTGEN_INLINE void hook(u64 n) { if (__builtin_expect(g_hook.mode != 0, 0)) hookWork(n); }
#else
OTGEN_INLINE void hook(u64 n) {
    if (g_hook.mode == 1) g_hook.ggm->step(n);
    else if (g_hook.mode == 2) g_hook.hash->step(n);
}
#endif
} // namespace otgen

#define OTGEN_HOOK_ROUTE() do { ++::otgen::g_hook.routeCalls; ::otgen::hook(::otgen::g_hook.routeChunk); } while (0)
#define OTGEN_HOOK_OUTER() do { ++::otgen::g_hook.outerCalls; ::otgen::hook(::otgen::g_hook.outerChunk); } while (0)
