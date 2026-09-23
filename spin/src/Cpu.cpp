#include "Cpu.h"
#include <spin/Code.h>
#ifdef _MSC_VER
#include <intrin.h>
#endif
namespace spin::detail {
namespace {
struct Features { bool avx2=false,f=false,vl=false,wide=false,masks=false; };
Features detect() noexcept {
    Features x;
#ifdef _MSC_VER
    int r[4]; __cpuid(r,0); if(r[0]<7) return x;
    __cpuidex(r,1,0);
    if((r[2]&(1<<27))==0 || (r[2]&(1<<28))==0) return x;
    const auto state=_xgetbv(0);
    if((state&6)!=6) return x;
    __cpuidex(r,7,0); x.avx2=(r[1]&(1<<5))!=0;
    if((state&0xe6)!=0xe6) return x;
    x.f=(r[1]&(1<<16))!=0; x.vl=(r[1]&(1u<<31))!=0;
    x.wide=x.f && x.vl && (r[1]&(1<<17)) && (r[1]&(1<<30));
    x.masks=x.wide && (r[2]&(1<<14));
#else
    __builtin_cpu_init();
    x.avx2=__builtin_cpu_supports("avx2");
    x.f=__builtin_cpu_supports("avx512f");
    x.vl=__builtin_cpu_supports("avx512vl");
    x.wide=x.f && x.vl && __builtin_cpu_supports("avx512dq") && __builtin_cpu_supports("avx512bw");
    x.masks=x.wide && __builtin_cpu_supports("avx512vpopcntdq");
#endif
    return x;
}
const Features& features() noexcept {static const auto f=detect();return f;}
}
bool cpu_avx2() noexcept {return features().avx2;}
bool cpu_avx512f() noexcept {return SPIN_BCH_AVX512 && features().f;}
bool cpu_avx512vl() noexcept {return SPIN_BCH_AVX512 && features().vl;}
bool cpu_wide512() noexcept {return SPIN_BCH_AVX512 && features().wide;}
bool cpu_mask512() noexcept {return SPIN_BCH_AVX512 && features().masks;}
}
namespace spin {
Capabilities capabilities() noexcept {
    const bool base=detail::cpu_avx2();
    return {base,base && detail::cpu_avx512f() && detail::cpu_avx512vl(),base,
            base && detail::cpu_wide512()};
}
}
