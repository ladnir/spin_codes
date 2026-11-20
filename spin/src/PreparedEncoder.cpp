#include <spin/PreparedEncoder.h>
#include "kernels/SetupRandom.h"
#include "kernels/WorkspaceRouting.h"
#include "Cpu.h"
#include "MemoryPolicy.h"
#include "Turbo.h"
#include <cstring>
#include <algorithm>
#include <limits>
#include <numeric>
#include <optional>

namespace spin::detail {
// Banked Turbo (PacketRs8T64) is heuristic. Full setup is certified under an
// ideal-randomness assumption on the seed expander; banked setup instead
// samples a structurally different, correlated distribution (one bank reused
// across codes) that the analysis does not cover. The bank is an exact Turbo
// plan drawn from the bank seed. A code seed selects
//  - route: a uniform order of the 64 bank regions and a rotation of each
//    region's group destinations;
//  - inner: a uniform order of the bank's GL22 epoch updates;
//  - outer: fresh nonzero GF(2^32) scalars, exactly the Full-mode outer.
// Refresh writes the live plan in place: no allocation, O(K) copying. The
// tables are kept in step only for the kernels that read them (see oracle,
// scalars); every kernel reads routeBytes and updates.
struct TurboBank {
    rs8::Plan bank,live;
    Backend kernel;
    // oracle: the Reference kernel reads route and reverseMatrices; scalars:
    // the Reference and Portable kernels read the outer scalars.
    bool oracle,scalars;
    std::array<std::uint32_t,64> regions{};
    std::array<std::uint32_t,64> shifts{};
    std::vector<std::uint32_t> epochs;
    TurboBank(CodeSpec s,std::uint64_t bankSeed,Backend k)
        :bank(s.message_size,bankSeed,bankSeed^0x7475726b62616e6bULL,turboSetupPath(k)),
         live(s.message_size,bankSeed,bankSeed^0x7475726b62616e6bULL,turboSetupPath(k)),kernel(k),oracle(k==Backend::Reference),scalars(k==Backend::Reference || k==Backend::Portable),
         epochs(s.message_size/32) {
        refresh(s.route_seed,s.inner_seed);
    }
    static std::uint32_t below(kernel::setup::Words& words,std::uint64_t n) noexcept {
#ifdef _MSC_VER
        return std::uint32_t(__umulh(words(),n));
#else
        return std::uint32_t((__uint128_t(words())*n)>>64);
#endif
    }
    void refresh(std::uint64_t routeSeed,std::uint64_t innerSeed) noexcept {
        using kernel::setup::Words;
        const auto& b=*bank.certified22Plan;
        auto& l=*live.certified22Plan;
        const auto groups=b.groups;
        Words routeWords(routeSeed^0x62616e6b726f7574ULL);
        std::iota(regions.begin(),regions.end(),0U);
        for(unsigned n=64;n>1;--n) std::swap(regions[n-1],regions[below(routeWords,n)]);
        for(auto& shift:shifts) shift=below(routeWords,groups);
        const auto copyRoute=[&](const std::vector<std::uint32_t>& from,std::vector<std::uint32_t>& to) {
            for(unsigned region=0;region<64;++region) {
                const auto* source=from.data()+regions[region]*groups;
                auto* target=to.data()+region*groups;
                const auto shift=shifts[region];
                std::memcpy(target+shift,source,(groups-shift)*sizeof(std::uint32_t));
                std::memcpy(target,source+groups-shift,shift*sizeof(std::uint32_t));
            }
        };
        copyRoute(b.routeBytes,l.routeBytes);
        if(oracle) copyRoute(b.route,l.route);
        Words innerWords(innerSeed^0x62616e6b696e6e72ULL);
        std::iota(epochs.begin(),epochs.end(),0U);
        for(auto n=epochs.size();n>1;--n) std::swap(epochs[n-1],epochs[below(innerWords,n)]);
        for(std::size_t e=0;e<epochs.size();++e) {
            std::memcpy(l.updates.data()+9*e,b.updates.data()+9*epochs[e],9*sizeof(std::uint64_t));
            if(oracle) l.reverseMatrices[e]=b.reverseMatrices[epochs[e]];
        }
        rs8::refreshOuterTranspose(live,routeSeed,scalars);
    }
    std::size_t bytes() const noexcept {
        return bank.setupBytes()+live.setupBytes()+epochs.capacity()*sizeof(epochs[0]);
    }
};
struct PreparedState {
    CodeSpec spec;
    SetupOptions setup;
    ExecutionOptions execution;
    std::optional<Code> full;
    std::optional<TurboBank> turbo;
    PreparedState(CodeSpec s,SetupOptions o,ExecutionOptions e):spec(s),setup(o),execution(e) {
        if(!valid_message_size(s.parameters,s.message_size))throw std::invalid_argument("SPIN invalid prepared geometry");
        if(o.mode==SetupMode::Full)full.emplace(s,e);
        else if(o.mode==SetupMode::BankedHeuristic) {
            if(e.tile_rows && e.tile_rows!=8)
                throw std::invalid_argument("SPIN Turbo tile geometry is fixed (8 rows)");
            turbo.emplace(s,o.bank_seed,selectTurboBackend(e.backend));
        } else throw std::invalid_argument("SPIN unknown setup mode");
    }
};
struct PreparedScratch {
    std::shared_ptr<PreparedState> owner;
    Width width;
    std::optional<spin::Workspace> full;
    std::optional<Buffer> turbo;
    PreparedScratch(std::shared_ptr<PreparedState> s,Width w,MemoryPolicy policy)
        :owner(std::move(s)),width(w) {
        (void)resolveMemoryPolicy(0,policy);
        if(owner->full)full.emplace(owner->full->make_workspace(w,policy));
        else {
            turbo.emplace(owner->turbo->live.scratchBlocks()*sizeof(storage::block),policy);
            if(resolveMemoryPolicy(turbo->bytes().size(),policy)==MemoryPolicy::Normal)
                kernel::workspace_routing::adviseOwned(turbo->bytes().data(),turbo->allocation_bytes());
        }
    }
};
static const PreparedState& checkedPreparedState(const std::shared_ptr<PreparedState>& state) {
    if(!state) throw std::invalid_argument("SPIN moved-from prepared encoder");
    return *state;
}
}
namespace spin {
PreparedEncoder::PreparedEncoder(CodeSpec s,SetupOptions o,ExecutionOptions e)
    :state_(std::make_shared<detail::PreparedState>(s,o,e)) {}
PreparedEncoder::~PreparedEncoder()=default;
PreparedEncoder::PreparedEncoder(PreparedEncoder&&) noexcept=default;
PreparedEncoder& PreparedEncoder::operator=(PreparedEncoder&&) noexcept=default;
void PreparedEncoder::setCodeSeed(CodeSeed seed) {
    detail::checkedPreparedState(state_);
    if(seed==CodeSeed{state_->spec.route_seed,state_->spec.inner_seed})return;
    if(state_->turbo)state_->turbo->refresh(seed.route,seed.inner);
    else {
        auto spec=state_->spec;spec.route_seed=seed.route;spec.inner_seed=seed.inner;
        // Build before committing, preserving the old encoder on allocation failure.
        state_->full=Code(spec,state_->execution);
    }
    state_->spec.route_seed=seed.route;state_->spec.inner_seed=seed.inner;
}
CodeSpec PreparedEncoder::specification() const noexcept{return state_?state_->spec:CodeSpec{0};}
SetupOptions PreparedEncoder::setup_options() const noexcept{return state_?state_->setup:SetupOptions{};}
std::size_t PreparedEncoder::message_size() const noexcept{return state_?state_->spec.message_size:0;}
std::size_t PreparedEncoder::code_size() const noexcept{return 2*message_size();}
std::size_t PreparedEncoder::setup_bytes() const noexcept {
    if(!state_)return 0;
    return state_->turbo?state_->turbo->bytes():state_->full->setup_bytes();
}
Backend PreparedEncoder::backend() const noexcept {
    if(!state_)return Backend::Automatic;
    return state_->full?state_->full->backend():state_->turbo->kernel;
}
bool PreparedEncoder::supports_transpose(Width width) const noexcept {
    if(!state_)return false;
    return state_->full?state_->full->supports_transpose(width):width==Width::Bits128;
}
std::array<std::byte,56> PreparedEncoder::descriptor() const noexcept {
    std::array<std::byte,56> out{};
    if(!state_)return out;
    out[0]=std::byte{'S'};out[1]=std::byte{'P'};out[2]=std::byte{'I'};out[3]=std::byte{'N'};
    auto put=[&](unsigned pos,std::uint64_t x,unsigned n){for(unsigned i=0;i<n;++i)out[pos+i]=std::byte((x>>(8*i))&255);};
    put(4,2,4);put(8,std::uint32_t(state_->spec.parameters),4);put(12,std::uint32_t(state_->setup.mode),4);
    put(16,message_size(),8);put(24,state_->spec.route_seed,8);put(32,state_->spec.inner_seed,8);
    if(state_->turbo){put(40,state_->setup.bank_seed,8);put(48,1,8);}
    return out;
}
PreparedEncoder::Workspace::Workspace(std::shared_ptr<detail::PreparedState> s,Width width,MemoryPolicy policy)
    :scratch_(std::make_unique<detail::PreparedScratch>(std::move(s),width,policy)) {}
PreparedEncoder::Workspace::~Workspace()=default;
PreparedEncoder::Workspace::Workspace(Workspace&&) noexcept=default;
PreparedEncoder::Workspace& PreparedEncoder::Workspace::operator=(Workspace&&) noexcept=default;
std::size_t PreparedEncoder::Workspace::bytes() const noexcept {
    if(!scratch_)return 0;
    return scratch_->turbo?scratch_->turbo->allocation_bytes():scratch_->full->bytes();
}
Width PreparedEncoder::Workspace::width() const noexcept {
    return scratch_?scratch_->width:Width::Bits128;
}
PreparedEncoder::Workspace PreparedEncoder::make_workspace(Width width,MemoryPolicy policy) const {
    detail::checkedPreparedState(state_);
    if(!supports_transpose(width))
        throw std::invalid_argument("SPIN prepared transpose does not support this record width");
    return Workspace(state_,width,policy);
}
Buffer PreparedEncoder::make_buffer(Width width,MemoryPolicy policy) const {
    detail::checkedPreparedState(state_);
    if(!supports_transpose(width))
        throw std::invalid_argument("SPIN prepared transpose does not support this record width");
    if(state_->full)return state_->full->make_buffer(width,policy);
    const auto element=static_cast<unsigned>(width);
    if(code_size()>std::numeric_limits<std::size_t>::max()/element)
        throw std::length_error("SPIN prepared buffer size overflow");
    return Buffer(code_size()*element,policy);
}
void PreparedEncoder::prepare_workspace(Workspace& w) const {
    detail::checkedPreparedState(state_);
    if(!w.scratch_ || w.scratch_->owner!=state_ || w.scratch_->width!=Width::Bits128)
        throw std::invalid_argument("SPIN prepared workspace mismatch");
    if(state_->full)state_->full->prepare_workspace(*w.scratch_->full);
}
static void turboTranspose(const detail::TurboBank& t,const std::byte* in,std::byte* out,detail::PreparedScratch& s) {
    detail::turboTranspose(t.kernel,reinterpret_cast<const detail::storage::block*>(in),
        reinterpret_cast<detail::storage::block*>(out),
        reinterpret_cast<detail::storage::block*>(s.turbo->bytes().data()),t.live);
}
static bool overlaps(const void* a,std::size_t na,const void* b,std::size_t nb) {
    const auto x=reinterpret_cast<std::uintptr_t>(a),y=reinterpret_cast<std::uintptr_t>(b);
    return x<=y?y-x<na:x-y<nb;
}
void PreparedEncoder::transpose_bytes(std::span<const std::byte> in,std::span<std::byte> out,Workspace& w) const {
    detail::checkedPreparedState(state_);
    if(overlaps(in.data(),in.size(),out.data(),out.size()))
        throw std::invalid_argument("SPIN buffers overlap; use transpose_inplace");
    if(!w.scratch_ || w.scratch_->owner!=state_ || w.scratch_->width!=Width::Bits128)
        throw std::invalid_argument("SPIN prepared workspace mismatch");
    if(in.size()!=code_size()*16 || out.size()!=message_size()*16 || !in.data() || !out.data() ||
       reinterpret_cast<std::uintptr_t>(in.data())%16 || reinterpret_cast<std::uintptr_t>(out.data())%16)
        throw std::invalid_argument("SPIN prepared buffer geometry or alignment mismatch");
    auto& s=*w.scratch_;
    if(state_->full){state_->full->prepare_workspace(*s.full);state_->full->transpose_bytes(in,out,*s.full);}
    else turboTranspose(*state_->turbo,in.data(),out.data(),s);
}
void PreparedEncoder::transpose_inplace_bytes(std::span<std::byte> buf,Workspace& w) const {
    detail::checkedPreparedState(state_);
    if(!w.scratch_ || w.scratch_->owner!=state_ || w.scratch_->width!=Width::Bits128)
        throw std::invalid_argument("SPIN prepared workspace mismatch");
    if(buf.size()!=code_size()*16 || !buf.data() || reinterpret_cast<std::uintptr_t>(buf.data())%16)
        throw std::invalid_argument("SPIN prepared buffer geometry or alignment mismatch");
    auto& s=*w.scratch_;
    if(state_->full){state_->full->prepare_workspace(*s.full);state_->full->transpose_inplace_bytes(buf,*s.full);}
    else turboTranspose(*state_->turbo,buf.data(),buf.data(),s);
}
}
