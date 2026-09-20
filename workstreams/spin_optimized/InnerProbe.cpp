// Diagnostic-only entry points, compiled with the production forward ISA flags.
#include "Spin.h"
#include "Inner.h"
using namespace bare_spin;
extern "C" void spin_probe_routed(std::size_t n,const u32* masks,const u32* route,const block* in,block* out) {
    innerForward<Map128S19>(n,masks,[&](std::size_t i){return in[route[i]];},out);
}
extern "C" void spin_probe_contiguous(std::size_t n,const u32* masks,const block* in,block* out) {
    innerForward<Map128S19>(n,masks,[&](std::size_t i){return in[i];},out);
}
extern "C" void spin_probe_gather(std::size_t n,const u32* route,const block* in,block* out) {
    for(std::size_t i=0;i<n;++i) out[i]=in[route[i]];
}
