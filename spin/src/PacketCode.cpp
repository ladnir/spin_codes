#include <spin/PacketCode.h>
#include "Cpu.h"
#include "MemoryPolicy.h"
#include "packet/PacketForward.h"
#include <optional>
#include "kernels/WorkspaceRouting.h"
#include <stdexcept>
#include <utility>

namespace spin::detail {
struct PacketState {
    PacketSpec spec;
    PacketBackend backend;
    packet::Plan plan;
    std::optional<packet::ForwardPlan> forward;
    PacketState(PacketSpec s,PacketBackend b):spec(s),backend(b),plan(s.message_size,s.seed) {
        if(b==PacketBackend::Avx512Gfni)forward.emplace(plan);
    }
};
struct PacketScratch {
    std::shared_ptr<const PacketState> state;
    Buffer allocation;
    PacketScratch(std::shared_ptr<const PacketState> s,PacketMemory policy)
        :state(std::move(s)),allocation(state->plan.scratchBlocks()*16,policy) {
        // Match the measured routed workspace: preparation advises its owned
        // pages even with normal alignment. The huge-page buffer does this itself.
        if(resolveMemoryPolicy(allocation.bytes().size(),policy)==PacketMemory::Normal)
            kernel::workspace_routing::adviseOwned(allocation.bytes().data(),allocation.allocation_bytes());
    }
};
namespace {
void aligned(const void* p) {
    if(!p || reinterpret_cast<std::uintptr_t>(p)%16)
        throw std::invalid_argument("SPIN packet buffers require 16-byte alignment");
}
bool overlap(const void* a,std::size_t na,const void* b,std::size_t nb) noexcept {
    const auto x=reinterpret_cast<std::uintptr_t>(a),y=reinterpret_cast<std::uintptr_t>(b);
    return x<=y?y-x<na:x-y<nb;
}
void encode(const PacketState& s,const void* in,void* out,void* scratch) {
    auto* input=static_cast<const storage::block*>(in);
    auto* output=static_cast<storage::block*>(out);
    auto* work=static_cast<storage::block*>(scratch);
#if SPIN_BCH_AVX512
    if(s.backend==PacketBackend::Avx512Gfni) {
        packet::transposeFast(input,output,work,s.plan);
        return;
    }
#endif
    packet::transposeScalar(input,output,work,s.plan);
}
}
}
namespace spin {
bool packet_fast_available() noexcept {return detail::cpu_packet512();}
bool valid_packet_message_size(std::size_t k) noexcept {return detail::packet::validMessageSize(k);}
PacketCode::PacketCode(PacketSpec spec,PacketBackend backend) {
    if(!valid_packet_message_size(spec.message_size))
        throw std::invalid_argument("SPIN packet K must be a positive supported multiple of 256");
    switch(backend) {
    case PacketBackend::Automatic:
        backend=packet_fast_available()?PacketBackend::Avx512Gfni:PacketBackend::Portable;
        break;
    case PacketBackend::Portable: break;
    case PacketBackend::Avx512Gfni:
        if(!packet_fast_available()) throw std::runtime_error("SPIN packet AVX512/VBMI/GFNI backend unavailable");
        break;
    default: throw std::invalid_argument("SPIN unknown packet backend");
    }
    state_=std::make_shared<detail::PacketState>(spec,backend);
}
const detail::PacketState& PacketCode::checked_state() const {
    if(!state_) throw std::invalid_argument("SPIN moved-from packet code");
    return *state_;
}
std::size_t PacketCode::message_size() const noexcept {return state_?state_->spec.message_size:0;}
std::size_t PacketCode::code_size() const noexcept {return 2*message_size();}
PacketSpec PacketCode::specification() const noexcept {return state_?state_->spec:PacketSpec{0,0};}
PacketBackend PacketCode::backend() const noexcept {return state_?state_->backend:PacketBackend::Automatic;}
std::size_t PacketCode::setup_bytes() const noexcept {return state_?state_->plan.setupBytes()+(state_->forward?state_->forward->setupBytes():0):0;}
std::array<std::byte,32> PacketCode::descriptor() const noexcept {
    std::array<std::byte,32> out{};
    if(!state_)return out;
    out[0]=std::byte{'S'};out[1]=std::byte{'P'};out[2]=std::byte{'K'};out[3]=std::byte{'P'};
    const auto put=[&](unsigned offset,std::uint64_t v,unsigned bytes) {
        for(unsigned i=0;i<bytes;++i)out[offset+i]=std::byte((v>>(8*i))&255);
    };
    put(4,2,4);put(8,2,4);put(16,message_size(),8);put(24,state_->spec.seed,8);
    return out;
}
PacketCode::Workspace PacketCode::make_workspace(PacketMemory p) const {
    checked_state();return Workspace(state_,p);
}
PacketBuffer PacketCode::make_buffer(PacketMemory p) const {
    checked_state();return PacketBuffer(code_size()*16,p);
}
PacketCode::Workspace::Workspace(std::shared_ptr<const detail::PacketState> s,PacketMemory p)
    :scratch_(std::make_unique<detail::PacketScratch>(std::move(s),p)) {}
PacketCode::Workspace::~Workspace()=default;
PacketCode::Workspace::Workspace(Workspace&&) noexcept=default;
PacketCode::Workspace& PacketCode::Workspace::operator=(Workspace&&) noexcept=default;
std::size_t PacketCode::Workspace::bytes() const noexcept {return scratch_?scratch_->allocation.allocation_bytes():0;}
detail::PacketScratch& PacketCode::checked_workspace(Workspace& w) const {
    checked_state();
    if(!w.scratch_ || w.scratch_->state!=state_)
        throw std::invalid_argument("SPIN packet workspace belongs to another plan or was moved");
    return *w.scratch_;
}
void PacketCode::forward_bytes(std::span<const std::byte> in,std::span<std::byte> out,Workspace& w) const {
    auto& work=checked_workspace(w);
    if(in.size()!=message_size()*16 || out.size()!=code_size()*16)
        throw std::invalid_argument("SPIN packet forward buffer size mismatch");
    detail::aligned(in.data());detail::aligned(out.data());
    if(detail::overlap(in.data(),in.size(),out.data(),out.size()))
        throw std::invalid_argument("SPIN packet forward buffers overlap");
    const auto* input=reinterpret_cast<const detail::storage::block*>(in.data());
    auto* output=reinterpret_cast<detail::storage::block*>(out.data());
    auto* scratch=reinterpret_cast<detail::storage::block*>(work.allocation.bytes().data());
#if SPIN_BCH_AVX512
    if(state_->backend==PacketBackend::Avx512Gfni) {
        detail::packet::forwardSelected(input,output,scratch,state_->plan,*state_->forward);
        return;
    }
#endif
    detail::packet::forwardScalar(input,output,scratch,state_->plan);
}
void PacketCode::transpose_bytes(std::span<const std::byte> in,std::span<std::byte> out,Workspace& w) const {
    auto& work=checked_workspace(w);
    if(in.size()!=code_size()*16 || out.size()!=message_size()*16)
        throw std::invalid_argument("SPIN packet transpose buffer size mismatch");
    detail::aligned(in.data());detail::aligned(out.data());
    if(detail::overlap(in.data(),in.size(),out.data(),out.size()))
        throw std::invalid_argument("SPIN packet out-of-place buffers overlap");
    detail::encode(*state_,in.data(),out.data(),work.allocation.bytes().data());
}
void PacketCode::transpose_inplace_bytes(std::span<std::byte> in,Workspace& w) const {
    auto& work=checked_workspace(w);
    if(in.size()!=code_size()*16) throw std::invalid_argument("SPIN packet inplace buffer size mismatch");
    detail::aligned(in.data());
    detail::encode(*state_,in.data(),in.data(),work.allocation.bytes().data());
}
}
