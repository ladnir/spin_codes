#include "Cpu.h"
#include <spin/Code.h>
#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86))
#include <intrin.h>
#endif
namespace spin::detail {
namespace {
struct Features { bool avx2=false,gfni=false,packet=false; };
// x86-64 feature detection; other architectures report no x86 features.
Features detect() noexcept {
    Features x;
#if defined(__x86_64__) || defined(_M_X64)
#if defined(_MSC_VER)
    int r[4]; __cpuid(r,0); if(r[0]<7) return x;
    __cpuidex(r,1,0);
    if((r[2]&(1<<27))==0 || (r[2]&(1<<28))==0) return x;
    const auto state=_xgetbv(0);
    if((state&6)!=6) return x;
    __cpuidex(r,7,0); x.avx2=(r[1]&(1<<5))!=0; x.gfni=x.avx2 && (r[2]&(1<<8))!=0;
    if((state&0xe6)!=0xe6) return x;
    // /arch:AVX512 permits DQ/CD as well as the intrinsics' explicit subset.
    x.packet=x.avx2 && (r[1]&(1<<16)) && (r[1]&(1u<<31)) && (r[1]&(1<<17)) && (r[1]&(1<<28)) &&
        (r[1]&(1<<30)) && (r[2]&(1<<1)) && (r[2]&(1<<8));
#else
    __builtin_cpu_init();
    x.avx2=__builtin_cpu_supports("avx2");
    x.gfni=x.avx2 && __builtin_cpu_supports("gfni");
    x.packet=x.avx2 && __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512vl") &&
        __builtin_cpu_supports("avx512dq") && __builtin_cpu_supports("avx512bw") &&
        __builtin_cpu_supports("avx512vbmi") && __builtin_cpu_supports("gfni");
#endif
#endif
    return x;
}
const Features& features() noexcept {static const auto f=detect();return f;}
}
bool cpu_avx2() noexcept {return features().avx2;}
// VEX-encoded GFNI (256-bit) needs AVX with GFNI; AVX2 implies AVX.
bool cpu_avx2gfni() noexcept {return features().gfni;}
bool cpu_packet512() noexcept {return SPIN_AVX512 && features().packet;}
}
namespace spin {
Capabilities capabilities() noexcept {
    return {detail::cpu_avx2(),detail::cpu_avx2gfni(),detail::cpu_packet512()};
}
}
