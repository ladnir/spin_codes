#include "Wide.h"
#include <stdexcept>
#define DECLARE(B) \
extern "C" void* hc_spin_wide_workspace##B(const void*) noexcept; \
extern "C" void hc_spin_wide_destroy##B(void*) noexcept; \
extern "C" std::size_t hc_spin_wide_bytes##B(const void*) noexcept; \
extern "C" int hc_spin_wide_encode##B(const void*,void*,const void*,std::size_t,void*,std::size_t) noexcept;
DECLARE(256)
DECLARE(512)
#undef DECLARE
namespace bare_spin {
bool wideAvailable(unsigned lanes) noexcept {
    if(lanes==2) return __builtin_cpu_supports("avx2");
#if !SPIN_TEST_NO_AVX512
    if(lanes==4) return __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512vl") &&
        __builtin_cpu_supports("avx512bw") && __builtin_cpu_supports("avx512dq");
#endif
    return false;
}
WideWorkspace::WideWorkspace(const Spin& code,unsigned lanes):mCode(code),mLanes(lanes),mWork(nullptr) {
    // Check ISA before even entering an ISA-specific workspace constructor.
    if(!wideAvailable(lanes)) throw std::invalid_argument("unsupported wide ISA or lane count");
    (void)code.wideForwardView();
    mWork=lanes==2?hc_spin_wide_workspace256(&code):hc_spin_wide_workspace512(&code);
    if(!mWork) throw std::runtime_error("wide workspace allocation failed");
}
WideWorkspace::~WideWorkspace() {
    if(mLanes==2) hc_spin_wide_destroy256(mWork); else hc_spin_wide_destroy512(mWork);
}
std::size_t WideWorkspace::bytes() const noexcept {
    return mLanes==2?hc_spin_wide_bytes256(mWork):hc_spin_wide_bytes512(mWork);
}
void WideWorkspace::forward(std::span<const block> in,std::span<block> out) {
    const int error=mLanes==2?
        hc_spin_wide_encode256(&mCode,mWork,in.data(),in.size(),out.data(),out.size()):
        hc_spin_wide_encode512(&mCode,mWork,in.data(),in.size(),out.data(),out.size());
    if(error) throw std::invalid_argument("wide forward: invalid shape, overlap, or workspace");
}
}
