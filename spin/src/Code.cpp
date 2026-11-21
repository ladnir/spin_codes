#include <spin/Code.h>
#include "Cpu.h"
#include "MemoryPolicy.h"
#include "kernels/WorkspaceRouting.h"
#include "Turbo.h"
#include <limits>
#include <mutex>

namespace spin::detail {
Backend selectTurboBackend(Backend requested) {
    switch(requested) {
    case Backend::Automatic:
        return cpu_packet512()?Backend::Avx512:cpu_avx2gfni()?Backend::Avx2Gfni:
               cpu_avx2()?Backend::Avx2:SPIN_AARCH64?Backend::Neon:Backend::Portable;
    case Backend::Portable:case Backend::Reference:return requested;
    case Backend::Avx2Gfni:
        if(!cpu_avx2gfni()) throw std::runtime_error("SPIN Turbo AVX2+GFNI backend unavailable");
        return requested;
    case Backend::Neon:
        if(!SPIN_AARCH64) throw std::runtime_error("SPIN Turbo NEON backend unavailable");
        return requested;
    case Backend::Avx512:
        if(!cpu_packet512()) throw std::runtime_error("SPIN Turbo AVX512/VBMI/GFNI backend unavailable");
        return requested;
    case Backend::Avx2:
        if(!cpu_avx2()) throw std::runtime_error("SPIN Turbo AVX2 backend unavailable");
        return requested;
    }
    throw std::invalid_argument("SPIN unknown backend");
}
certified22::SetupPath turboSetupPath(Backend b) noexcept {
    // An explicit AVX2 or Portable backend also uses that setup (bit-identical).
    return (b==Backend::Avx2 || b==Backend::Avx2Gfni)?certified22::SetupPath::Avx2:
        b==Backend::Portable?certified22::SetupPath::Portable:certified22::SetupPath::Automatic;
}
void turboTranspose(Backend b,const storage::block* input,storage::block* output,
                    storage::block* work,const rs8::Plan& plan) {
#if SPIN_AVX512
    if(b==Backend::Avx512) { rs8::transposeFast(input,output,work,plan); return; }
#endif
#if SPIN_X86
    if(b==Backend::Avx2Gfni) { rs8::transposeAvx2Gfni(input,output,work,plan); return; }
    if(b==Backend::Avx2) { rs8::transposeAvx2(input,output,work,plan); return; }
#endif
#if SPIN_AARCH64
    if(b==Backend::Neon) { rs8::transposeNeon(input,output,work,plan); return; }
#endif
    if(b==Backend::Reference) rs8::transposeScalar(input,output,work,plan);
    else rs8::transposePortable(input,output,work,plan);
}
bool turboBitsFast(Backend b) noexcept {
#if SPIN_AVX512
    return b==Backend::Avx512;
#else
    (void)b;return false;
#endif
}
void turboTransposeBits(Backend b,const std::uint64_t* input,std::uint64_t* output,const rs8::Plan& plan,
                        const rs8::bits::Tables& tables,rs8::bits::Scratch& scratch) {
#if SPIN_AVX512
    if(turboBitsFast(b) && tables.fast) { rs8::bits::transposeFast(input,output,plan,tables,scratch); return; }
#endif
    (void)b;
    rs8::bits::transposePortable(input,output,plan,tables,scratch);
}
struct Plan {
    const CodeSpec spec;
    const Backend selected;
    const rs8::Plan turbo;
    Plan(CodeSpec s,ExecutionOptions o)
        :spec(s),selected(selectTurboBackend(o.backend)),
         turbo(s.message_size,s.route_seed,s.inner_seed,turboSetupPath(selected)) {}
    std::size_t scratchBlocks() const noexcept { return turbo.scratchBlocks(); }
    // Bit tables are built on the first bit workspace and then immutable.
    const rs8::bits::Tables& bitTables() const {
        std::call_once(bitOnce,[&] { bits.build(turbo,turboBitsFast(selected)); });
        return bits;
    }
private:
    mutable std::once_flag bitOnce;
    mutable rs8::bits::Tables bits;
};
struct BitScratch {
    std::shared_ptr<const Plan> plan;
    const rs8::bits::Tables* tables;
    rs8::bits::Scratch scratch;
    explicit BitScratch(std::shared_ptr<const Plan> p)
        :plan(std::move(p)),tables(&plan->bitTables()),scratch(plan->turbo) {}
};
struct Scratch {
    std::shared_ptr<const Plan> plan;
    Width width;
    MemoryPolicy policy;
    Buffer buffer;
    Scratch(std::shared_ptr<const Plan> p,Width w,MemoryPolicy memory)
        :plan(std::move(p)),width(w),policy(memory),
         buffer(plan->scratchBlocks()*static_cast<unsigned>(w),memory) {
        // The allocation is owned here, so all advised whole pages belong
        // to this workspace. Foreign input/output buffers are never advised.
        // PreferHugePages already advises its owned allocation in Buffer.
        if(resolveMemoryPolicy(buffer.bytes().size(),memory)==MemoryPolicy::Normal)
            kernel::workspace_routing::adviseOwned(buffer.bytes().data(),buffer.allocation_bytes());
    }
    storage::block* work() noexcept { return reinterpret_cast<storage::block*>(buffer.bytes().data()); }
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
}
namespace spin {
std::size_t message_alignment(Parameters p) {
    if(p==Parameters::PacketRs8T64) return 128;
    throw std::invalid_argument("SPIN unknown parameter set");
}
bool valid_message_size(Parameters p,std::size_t k) noexcept {
    return p==Parameters::PacketRs8T64 && detail::rs8::validMessageSize(k);
}
Code::Code(CodeSpec s,ExecutionOptions o) {
    if(!valid_message_size(s.parameters,s.message_size))
        throw std::invalid_argument("SPIN invalid parameter set or message length");
    if(o.tile_rows && o.tile_rows!=8)
        throw std::invalid_argument("SPIN Turbo tile geometry is fixed (8 rows)");
    plan_=std::make_shared<detail::Plan>(s,o);
}
std::size_t Code::message_size() const noexcept {return plan_?plan_->spec.message_size:0;}
std::size_t Code::code_size() const noexcept {return 2*message_size();}
CodeSpec Code::specification() const noexcept {return plan_?plan_->spec:CodeSpec{0};}
Backend Code::backend() const noexcept {return plan_?plan_->selected:Backend::Automatic;}
std::size_t Code::setup_bytes() const noexcept {return plan_?plan_->turbo.setupBytes():0;}
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
    return plan_ && (width==Width::Bits128 || width==Width::Bits256 || width==Width::Bits512);
}
bool Code::supports_transpose(Width width) const noexcept {
    return plan_ && width==Width::Bits128;
}
Workspace Code::make_workspace(Width width,MemoryPolicy policy) const {
    detail::checkedPlan(plan_);
    if(!supports_forward(width))
        throw std::invalid_argument("SPIN unsupported record width");
    return Workspace(plan_,width,policy);
}
Buffer Code::make_buffer(Width width,MemoryPolicy policy) const {
    detail::checkedPlan(plan_);
    if(!supports_forward(width))
        throw std::invalid_argument("SPIN unsupported record width");
    const auto element=static_cast<unsigned>(width);
    if(code_size()>std::numeric_limits<std::size_t>::max()/element)
        throw std::length_error("SPIN buffer size overflow");
    return Buffer(code_size()*element,policy);
}
void Code::prepare_workspace(Workspace& w) const {
    detail::checkedPlan(plan_);
    if(w.scratch_ && w.scratch_->plan==plan_)return;
    if(w.scratch_ &&
       w.scratch_->buffer.bytes().size()==plan_->scratchBlocks()*static_cast<unsigned>(w.scratch_->width)) {
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
std::size_t Workspace::bytes() const noexcept {return scratch_?scratch_->buffer.allocation_bytes():0;}
Width Workspace::width() const noexcept {return scratch_?scratch_->width:Width::Bits128;}
detail::Scratch& Code::check_workspace(Workspace& w) const {
    detail::checkedPlan(plan_);
    if(!w.scratch_ || w.scratch_->plan!=plan_)
        throw std::invalid_argument("SPIN workspace belongs to another plan or was moved");
    return *w.scratch_;
}
void Code::forward_bytes(std::span<const std::byte> in,std::span<std::byte> out,Workspace& w,
                         OutputStores stores) const {
    auto& scratch=check_workspace(w);
    if(stores!=OutputStores::Streaming && stores!=OutputStores::Cached)
        throw std::invalid_argument("SPIN unknown forward output store policy");
    const bool stream=stores==OutputStores::Streaming;
    const auto width=static_cast<unsigned>(scratch.width);
    if(in.size()!=message_size()*width || out.size()!=code_size()*width)
        throw std::invalid_argument("SPIN forward buffer size mismatch");
    detail::aligned(in.data());detail::aligned(out.data());
    if(detail::overlap(in.data(),in.size(),out.data(),out.size()))
        throw std::invalid_argument("SPIN forward buffers overlap");
    const auto* input=reinterpret_cast<const detail::storage::block*>(in.data());
    auto* output=reinterpret_cast<detail::storage::block*>(out.data());
    auto* work=scratch.work();
    const auto& turbo=plan_->turbo;
    const unsigned lanes=width/16;
    switch(plan_->selected) {
#if SPIN_AVX512
    case Backend::Avx512:detail::rs8::forwardFast(input,output,work,turbo,lanes,stream);return;
#endif
#if SPIN_X86
    case Backend::Avx2Gfni:detail::rs8::forwardAvx2Gfni(input,output,work,turbo,lanes,stream);return;
    case Backend::Avx2:detail::rs8::forwardAvx2(input,output,work,turbo,lanes,stream);return;
#endif
#if SPIN_AARCH64
    case Backend::Neon:detail::rs8::forwardNeon(input,output,work,turbo,lanes);return;
#endif
    case Backend::Reference:detail::rs8::forwardScalar(input,output,work,turbo,lanes);return;
    default:detail::rs8::forwardPortable(input,output,work,turbo,lanes);return;
    }
}
void Code::transpose_bytes(std::span<const std::byte> in,std::span<std::byte> out,Workspace& w) const {
    auto& scratch=check_workspace(w);
    if(scratch.width!=Width::Bits128 || in.size()!=code_size()*16 || out.size()!=message_size()*16)
        throw std::invalid_argument("SPIN transpose requires 128-bit records and exact sizes");
    detail::aligned(in.data());detail::aligned(out.data());
    if(detail::overlap(in.data(),in.size(),out.data(),out.size()))
        throw std::invalid_argument("SPIN transpose buffers overlap; use transpose_inplace");
    detail::turboTranspose(plan_->selected,reinterpret_cast<const detail::storage::block*>(in.data()),
        reinterpret_cast<detail::storage::block*>(out.data()),scratch.work(),plan_->turbo);
}
BitWorkspace Code::make_bit_workspace() const {
    detail::checkedPlan(plan_);
    return BitWorkspace(plan_);
}
void Code::prepare_bit_workspace(BitWorkspace& w) const {
    detail::checkedPlan(plan_);
    if(w.scratch_ && w.scratch_->plan==plan_)return;
    if(w.scratch_ && w.scratch_->plan->spec.message_size==plan_->spec.message_size) {
        w.scratch_->tables=&plan_->bitTables();
        w.scratch_->plan=plan_;
        return;
    }
    w=make_bit_workspace();
}
void Code::transpose_bits(std::span<const std::uint64_t> in,std::span<std::uint64_t> out,BitWorkspace& w) const {
    detail::checkedPlan(plan_);
    if(!w.scratch_ || w.scratch_->plan!=plan_)
        throw std::invalid_argument("SPIN bit workspace belongs to another plan or was moved");
    if(in.size()!=code_size()/64 || out.size()!=message_size()/64 || !in.data() || !out.data())
        throw std::invalid_argument("SPIN bit transpose requires 2K input bits and K output bits");
    if(detail::overlap(in.data(),in.size()*8,out.data(),out.size()*8))
        throw std::invalid_argument("SPIN bit transpose buffers overlap");
    detail::turboTransposeBits(plan_->selected,in.data(),out.data(),plan_->turbo,*w.scratch_->tables,w.scratch_->scratch);
}
BitWorkspace::BitWorkspace(std::shared_ptr<const detail::Plan> p)
    :scratch_(std::make_unique<detail::BitScratch>(std::move(p))) {}
BitWorkspace::~BitWorkspace()=default;
BitWorkspace::BitWorkspace(BitWorkspace&&) noexcept=default;
BitWorkspace& BitWorkspace::operator=(BitWorkspace&&) noexcept=default;
std::size_t BitWorkspace::bytes() const noexcept {return scratch_?scratch_->scratch.bytes():0;}
void Code::transpose_inplace_bytes(std::span<std::byte> buffer,Workspace& w) const {
    auto& scratch=check_workspace(w);
    if(scratch.width!=Width::Bits128 || buffer.size()!=code_size()*16)
        throw std::invalid_argument("SPIN in-place transpose requires 2K 128-bit records");
    detail::aligned(buffer.data());
    auto* p=reinterpret_cast<detail::storage::block*>(buffer.data());
    detail::turboTranspose(plan_->selected,p,p,scratch.work(),plan_->turbo);
}
}
