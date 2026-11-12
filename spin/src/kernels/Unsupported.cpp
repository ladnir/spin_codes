// Non-x86 targets: the BCH/IMT families need x86-64 with AVX2. Their public
// entry points already reject them (capabilities().avx2 is false), so these
// definitions only satisfy the linker; reaching one is a logic error.
#include "Spin.h"
#include "BankKernel.h"
#include <cstddef>
#include <stdexcept>

namespace spin::detail::kernel {
[[noreturn]] static void unsupported() {
    throw std::runtime_error("SPIN BCH/IMT families require x86-64 with AVX2");
}
Spin::Workspace::Workspace(const Spin&) { unsupported(); }
std::size_t Spin::Workspace::bytes() const noexcept { return 0; }
Spin::Spin(Configuration,MessageLength,u64,u64,unsigned,BchBackend,bool,std::vector<u32>) { unsupported(); }
std::size_t Spin::setupBytes() const noexcept { return 0; }
void Spin::compact(Layout) { unsupported(); }
void Spin::encodeUnchecked(const block*,block*,Workspace&,Layout) const { unsupported(); }
void Spin::forwardUnchecked(const block*,block*,Workspace&,Layout,bool) const { unsupported(); }
void Spin::forwardBits(const std::uint64_t*,std::size_t,std::uint64_t*,std::size_t,std::uint64_t*,std::size_t) const { unsupported(); }
void bankTranspose(const BankState&,bool,const void*,void*,storage::block*,std::uint32_t*) { unsupported(); }
}
extern "C" void* spin_internal_wide_workspace256(const void*) noexcept { return nullptr; }
extern "C" void spin_internal_wide_destroy256(void*) noexcept {}
extern "C" std::size_t spin_internal_wide_bytes256(const void*) noexcept { return 0; }
extern "C" int spin_internal_wide_encode256(const void*,void*,const void*,std::size_t,void*,std::size_t,bool) noexcept { return 0; }
