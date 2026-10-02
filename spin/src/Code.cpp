#include <spin/Code.h>
#include <spin/Generic.h>
#include "Cpu.h"
#include "MemoryPolicy.h"
#include "kernels/Spin.h"
#include "kernels/WorkspaceRouting.h"
#include "packet/PacketForward.h"
#include <limits>
#include <optional>
#include <vector>

#define WIDE_DECL(B) \
extern "C" void* spin_internal_wide_workspace##B(const void*) noexcept; \
extern "C" void spin_internal_wide_destroy##B(void*) noexcept; \
extern "C" std::size_t spin_internal_wide_bytes##B(const void*) noexcept; \
extern "C" int spin_internal_wide_encode##B(const void*,void*,const void*,std::size_t,void*,std::size_t) noexcept;
WIDE_DECL(256)
#if SPIN_BCH_AVX512
WIDE_DECL(512)
#endif
#undef WIDE_DECL

namespace spin::detail::kernel {
struct Access {
    static const std::vector<u32>& route(const Spin& c) {return c.mRoute;}
    static const std::vector<u32>& masks(const Spin& c) {return c.mFieldRows;}
};
}
namespace spin::detail {
static kernel::Configuration config(Parameters p) {
    switch(p) {
    case Parameters::T128S19:return kernel::Configuration::T128S19;
    case Parameters::T64S12:return kernel::Configuration::T64S12;
    case Parameters::T64S12R2:return kernel::Configuration::T64S12R2;
    case Parameters::PacketRsT64S20:break;
    }
    throw std::invalid_argument("SPIN unknown parameter set");
}
static kernel::BchBackend backend(Backend b) {
    switch(b) {
    case Backend::Automatic:return kernel::BchBackend::Auto;
    case Backend::Avx2:return kernel::BchBackend::Avx2;
    case Backend::Avx512:return kernel::BchBackend::Avx512;
    case Backend::Portable:break;
    }
    throw std::invalid_argument("SPIN unknown backend");
}
struct Plan {
    const CodeSpec spec;
    Backend selected;
    std::optional<kernel::Spin> code;
    std::optional<packet::Plan> packets;
    std::optional<packet::ForwardPlan> forward;
    Plan(CodeSpec s,ExecutionOptions o):spec(s),selected(o.backend) {
        if(s.parameters==Parameters::PacketRsT64S20) {
            if(o.tile_rows && o.tile_rows!=8)
                throw std::invalid_argument("SPIN RS packet requires eight-row tiles");
            switch(o.backend) {
            case Backend::Automatic:
                selected=cpu_packet512()?Backend::Avx512:Backend::Portable;
                break;
            case Backend::Portable:break;
            case Backend::Avx512:
                if(!cpu_packet512())
                    throw std::runtime_error("SPIN packet AVX512/VBMI/GFNI backend unavailable");
                break;
            case Backend::Avx2:
                throw std::invalid_argument("SPIN packet has no AVX2 backend; select Portable or Automatic");
            default:throw std::invalid_argument("SPIN unknown backend");
            }
            packets.emplace(s.message_size,s.route_seed,s.inner_seed);
            if(selected==Backend::Avx512)forward.emplace(*packets);
        } else {
            if(o.backend==Backend::Portable)
                throw std::invalid_argument("SPIN Portable backend is available only for PacketRsT64S20");
            if(!capabilities().avx2) throw std::runtime_error("SPIN requires AVX2 with OS support");
            code.emplace(config(s.parameters),kernel::MessageLength{s.message_size},
                s.route_seed,s.inner_seed,o.tile_rows?o.tile_rows:256,backend(o.backend),true);
            code->compact();
            selected=code->bchBackend()==kernel::BchBackend::Avx512?Backend::Avx512:Backend::Avx2;
        }
    }
};
struct Scratch {
    std::shared_ptr<const Plan> plan;
    Width width;
    MemoryPolicy policy;
    std::unique_ptr<kernel::Spin::Workspace> single;
    std::optional<Buffer> packetScratch;
    void* wide=nullptr;
    Scratch(std::shared_ptr<const Plan> p,Width w,MemoryPolicy memory)
        :plan(std::move(p)),width(w),policy(memory) {
        (void)resolveMemoryPolicy(0,memory);
        if(plan->packets) {
            packetScratch.emplace(plan->packets->scratchBlocks()*16,memory);
            // The allocation is owned here, so all advised whole pages belong
            // to this workspace. Foreign input/output buffers are never advised.
            // PreferHugePages already advises its owned allocation in Buffer.
            if(resolveMemoryPolicy(packetScratch->bytes().size(),memory)==MemoryPolicy::Normal)
                kernel::workspace_routing::adviseOwned(packetScratch->bytes().data(),
                                                       packetScratch->allocation_bytes());
        } else if(w==Width::Bits128)
            single=std::make_unique<kernel::Spin::Workspace>(*plan->code);
        else {
            if(w==Width::Bits256) wide=spin_internal_wide_workspace256(&*plan->code);
#if SPIN_BCH_AVX512
            else wide=spin_internal_wide_workspace512(&*plan->code);
#endif
            if(!wide) throw std::bad_alloc();
        }
    }
    ~Scratch() {
        if(!wide) return;
        if(width==Width::Bits256) spin_internal_wide_destroy256(wide);
#if SPIN_BCH_AVX512
        else spin_internal_wide_destroy512(wide);
#endif
    }
    std::size_t bytes() const noexcept {
        if(packetScratch) return packetScratch->allocation_bytes();
        if(single) return single->bytes();
        if(width==Width::Bits256) return spin_internal_wide_bytes256(wide);
#if SPIN_BCH_AVX512
        return spin_internal_wide_bytes512(wide);
#else
        return 0;
#endif
    }
};
static void aligned(const void* p) {
    if(!p || reinterpret_cast<std::uintptr_t>(p)%16)
        throw std::invalid_argument("SPIN SIMD buffers require 16-byte alignment");
}
static bool overlap(const void* a,std::size_t na,const void* b,std::size_t nb) noexcept {
    const auto x=reinterpret_cast<std::uintptr_t>(a),y=reinterpret_cast<std::uintptr_t>(b);
    return x<=y?y-x<na:x-y<nb;
}
static const Plan& checkedPlan(const std::shared_ptr<const Plan>& plan) {
    if(!plan) throw std::invalid_argument("SPIN moved-from code");
    return *plan;
}
static void transposePacket(const Plan& plan,const storage::block* input,
                            storage::block* output,Scratch& scratch) {
    const auto& prepared=*plan.packets;
    auto* work=reinterpret_cast<storage::block*>(scratch.packetScratch->bytes().data());
#if SPIN_BCH_AVX512
    if(plan.selected==Backend::Avx512) {
        packet::transposeFast(input,output,work,prepared);
        return;
    }
#endif
    packet::transposeScalar(input,output,work,prepared);
}
}
namespace spin {
std::size_t message_alignment(Parameters p) {
    switch(p) {
    case Parameters::T128S19:return 16384;
    case Parameters::T64S12:case Parameters::T64S12R2:return 8192;
    case Parameters::PacketRsT64S20:return 256;
    }
    throw std::invalid_argument("SPIN unknown parameter set");
}
bool valid_message_size(Parameters p,std::size_t k) noexcept {
    if(p==Parameters::PacketRsT64S20) return detail::packet::validMessageSize(k);
    const auto unit=p==Parameters::T128S19?16384u:
        (p==Parameters::T64S12 || p==Parameters::T64S12R2)?8192u:0u;
    return unit && k && k<=std::numeric_limits<std::uint32_t>::max()/2 && k%unit==0;
}
Code::Code(CodeSpec s,ExecutionOptions o) {
    if(!valid_message_size(s.parameters,s.message_size))
        throw std::invalid_argument("SPIN invalid parameter set or message length");
    plan_=std::make_shared<detail::Plan>(s,o);
}
std::size_t Code::message_size() const noexcept {return plan_?plan_->spec.message_size:0;}
std::size_t Code::code_size() const noexcept {return 2*message_size();}
CodeSpec Code::specification() const noexcept {return plan_?plan_->spec:CodeSpec{0};}
Backend Code::backend() const noexcept {
    return plan_?plan_->selected:Backend::Automatic;
}
std::size_t Code::setup_bytes() const noexcept {
    if(!plan_) return 0;
    return plan_->packets?plan_->packets->setupBytes()+(plan_->forward?plan_->forward->setupBytes():0):plan_->code->setupBytes();
}
std::array<std::byte,40> Code::descriptor() const noexcept {
    std::array<std::byte,40> out{};
    if(!plan_) return out;
    out[0]=std::byte{'S'};out[1]=std::byte{'P'};out[2]=std::byte{'I'};out[3]=std::byte{'N'};
    const auto put=[&](unsigned offset,std::uint64_t x,unsigned bytes) {
        for(unsigned i=0;i<bytes;++i) out[offset+i]=std::byte((x>>(8*i))&255);
    };
    put(4,1,4);put(8,static_cast<std::uint32_t>(plan_->spec.parameters),4);
    // Bytes 12..15 are reserved, zero. The parameter identifier fixes the construction.
    put(16,message_size(),8);put(24,plan_->spec.route_seed,8);put(32,plan_->spec.inner_seed,8);
    return out;
}
bool Code::supports_forward(Width width) const noexcept {
    if(!plan_) return false;
    if(plan_->packets)return width==Width::Bits128;
    if(width==Width::Bits128) return true;
    if(plan_->spec.parameters==Parameters::T64S12) return false;
    return width==Width::Bits256?capabilities().forward256:
        width==Width::Bits512?capabilities().forward512:false;
}
bool Code::supports_transpose(Width width) const noexcept {
    return plan_ && width==Width::Bits128;
}
bool Code::supports_generic_transpose() const noexcept {
    return plan_ && !plan_->packets;
}
Workspace Code::make_workspace(Width width,MemoryPolicy policy) const {
    detail::checkedPlan(plan_);
    if(!supports_forward(width) && !supports_transpose(width))
        throw std::invalid_argument("SPIN unsupported record width, map, or CPU");
    return Workspace(plan_,width,policy);
}
Buffer Code::make_buffer(Width width,MemoryPolicy policy) const {
    detail::checkedPlan(plan_);
    if(!supports_forward(width) && !supports_transpose(width))
        throw std::invalid_argument("SPIN unsupported record width, map, or CPU");
    const auto element=static_cast<unsigned>(width);
    if(code_size()>std::numeric_limits<std::size_t>::max()/element)
        throw std::length_error("SPIN buffer size overflow");
    return Buffer(code_size()*element,policy);
}
void Code::prepare_workspace(Workspace& w) const {
    detail::checkedPlan(plan_);
    if(w.scratch_ && w.scratch_->plan==plan_)return;
    if(w.scratch_ && w.scratch_->width==Width::Bits128 && plan_->packets &&
       w.scratch_->packetScratch &&
       w.scratch_->packetScratch->bytes().size()==plan_->packets->scratchBlocks()*16) {
        w.scratch_->plan=plan_;
        return;
    }
    if(w.scratch_ && w.scratch_->width==Width::Bits128 && plan_->code && w.scratch_->single &&
       w.scratch_->single->buckets.size()==code_size() &&
       w.scratch_->single->tile.size()==plan_->code->tileBlocks()) {
        w.scratch_->plan=plan_;
        return;
    }
    const auto width=w.scratch_?w.scratch_->width:Width::Bits128;
    const auto policy=w.scratch_?w.scratch_->policy:MemoryPolicy::Automatic;
    w=make_workspace(width,policy);
}
Workspace::Workspace(std::shared_ptr<const detail::Plan> p,Width w,MemoryPolicy policy)
    :scratch_(std::make_unique<detail::Scratch>(std::move(p),w,policy)) {}
Workspace::~Workspace()=default;
Workspace::Workspace(Workspace&&) noexcept=default;
Workspace& Workspace::operator=(Workspace&&) noexcept=default;
std::size_t Workspace::bytes() const noexcept {return scratch_?scratch_->bytes():0;}
Width Workspace::width() const noexcept {return scratch_?scratch_->width:Width::Bits128;}
detail::Scratch& Code::check_workspace(Workspace& w) const {
    detail::checkedPlan(plan_);
    if(!w.scratch_ || w.scratch_->plan!=plan_)
        throw std::invalid_argument("SPIN workspace belongs to another plan or was moved");
    return *w.scratch_;
}
void Code::forward_bytes(std::span<const std::byte> in,std::span<std::byte> out,Workspace& w) const {
    auto& scratch=check_workspace(w);
    if(!supports_forward(scratch.width))
        throw std::invalid_argument("SPIN selected family does not support forward encoding at this width");
    const auto width=static_cast<unsigned>(scratch.width);
    if(in.size()!=message_size()*width || out.size()!=code_size()*width)
        throw std::invalid_argument("SPIN forward buffer size mismatch");
    detail::aligned(in.data());detail::aligned(out.data());
    if(detail::overlap(in.data(),in.size(),out.data(),out.size()))
        throw std::invalid_argument("SPIN forward buffers overlap");
    if(plan_->packets) {
        const auto* input=reinterpret_cast<const detail::storage::block*>(in.data());
        auto* output=reinterpret_cast<detail::storage::block*>(out.data());
        auto* work=reinterpret_cast<detail::storage::block*>(scratch.packetScratch->bytes().data());
#if SPIN_BCH_AVX512
        if(plan_->selected==Backend::Avx512) {
            detail::packet::forwardSelected(input,output,work,*plan_->packets,*plan_->forward);
            return;
        }
#endif
        detail::packet::forwardScalar(input,output,work,*plan_->packets);
        return;
    }
    if(width==16) {
        plan_->code->forwardUnchecked(reinterpret_cast<const detail::storage::block*>(in.data()),
            reinterpret_cast<detail::storage::block*>(out.data()),*scratch.single);
        return;
    }
    int error=1;
    if(width==32) error=spin_internal_wide_encode256(&*plan_->code,scratch.wide,in.data(),in.size()/16,out.data(),out.size()/16);
#if SPIN_BCH_AVX512
    else error=spin_internal_wide_encode512(&*plan_->code,scratch.wide,in.data(),in.size()/16,out.data(),out.size()/16);
#endif
    if(error) throw std::runtime_error("SPIN wide kernel rejected validated buffers");
}
void Code::transpose_bytes(std::span<const std::byte> in,std::span<std::byte> out,Workspace& w) const {
    auto& scratch=check_workspace(w);
    if(scratch.width!=Width::Bits128 || in.size()!=code_size()*16 || out.size()!=message_size()*16)
        throw std::invalid_argument("SPIN transpose requires 128-bit records and exact sizes");
    detail::aligned(in.data());detail::aligned(out.data());
    if(detail::overlap(in.data(),in.size(),out.data(),out.size()))
        throw std::invalid_argument("SPIN transpose buffers overlap; use transpose_inplace");
    const auto* input=reinterpret_cast<const detail::storage::block*>(in.data());
    auto* output=reinterpret_cast<detail::storage::block*>(out.data());
    if(plan_->packets) detail::transposePacket(*plan_,input,output,scratch);
    else plan_->code->encodeUnchecked(input,output,*scratch.single);
}
void Code::transpose_inplace_bytes(std::span<std::byte> buffer,Workspace& w) const {
    auto& scratch=check_workspace(w);
    if(scratch.width!=Width::Bits128 || buffer.size()!=code_size()*16)
        throw std::invalid_argument("SPIN in-place transpose requires 2K 128-bit records");
    detail::aligned(buffer.data());
    auto* p=reinterpret_cast<detail::storage::block*>(buffer.data());
    if(plan_->packets) detail::transposePacket(*plan_,p,p,scratch);
    else plan_->code->encodeUnchecked(p,p,*scratch.single);
}
void Code::forward_bits(std::span<const std::uint64_t> in,std::span<std::uint64_t> out,
                        std::span<std::uint64_t> scratch) const {
    detail::checkedPlan(plan_);
    if(plan_->packets)
        throw std::invalid_argument("SPIN PacketRsT64S20 does not support bit-packed forward encoding");
    plan_->code->forwardBits(in.data(),in.size(),out.data(),out.size(),scratch.data(),scratch.size());
}
GenericTranspose Code::generic_transpose() const {
    detail::checkedPlan(plan_);
    if(!supports_generic_transpose())
        throw std::invalid_argument("SPIN PacketRsT64S20 does not support generic transpose");
    return GenericTranspose(plan_->spec.parameters,message_size(),
        detail::kernel::Access::route(*plan_->code),detail::kernel::Access::masks(*plan_->code));
}
}
