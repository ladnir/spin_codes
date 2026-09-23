#include <spin/PreparedEncoder.h>
#include "kernels/BankKernel.h"
#include "kernels/SetupRandom.h"
#include "kernels/WorkspaceRouting.h"
#include "Cpu.h"
#include <algorithm>
#include <numeric>
#include <optional>

namespace spin::detail {
BankState::BankState(CodeSpec s,std::uint64_t bankSeed)
    :spec(s),bank_seed(bankSeed),rows(s.message_size/128),bank(2*s.message_size),row_keys(rows),
     masks((s.parameters==Parameters::T128S19?2: s.parameters==Parameters::T64S12?4:8)*s.message_size/64) {
    using kernel::setup::Words;
    using kernel::setup::Divisor;
    std::vector<Divisor> divisors(std::max<std::size_t>(256,rows)+1);
    for(std::size_t n=2;n<divisors.size();++n)divisors[n]=Divisor(n);
    Words rowWords(bankSeed^0x726f7773ULL),regionSeeds(bankSeed^0x72656773ULL);
    std::vector<std::uint8_t> inverse(256*rows);
    std::array<std::uint32_t,256> perm;
    for(std::size_t row=0;row<rows;++row) {
        std::iota(perm.begin(),perm.end(),0);
        for(unsigned n=256;n>1;--n)std::swap(perm[n-1],perm[divisors[n].sample(rowWords,n)]);
        for(unsigned c=0;c<256;++c)inverse[row*256+perm[c]]=std::uint8_t(c);
    }
    for(unsigned g=0;g<256;++g) {
        Words words(regionSeeds());
        auto* table=bank.data()+g*rows;
        std::iota(table,table+rows,0);
        for(std::size_t n=rows;n>1;--n)std::swap(table[n-1],table[divisors[n].sample(words,n)]);
        for(std::size_t p=0;p<rows;++p) {
            const auto row=table[p];table[p]=row*256+inverse[row*256+g];
        }
    }
    if(s.message_size==(1U<<18) && s.parameters==Parameters::T128S19) {
        indexed.resize(bank.size());packed_keys.resize(rows);
        for(std::size_t i=0;i<bank.size();++i) {
            const auto x=bank[i],row=x>>8;
            indexed[i]=(x&~1023U)|((x&255)<<2)|(row&3)|(row<<19);
        }
    }
    if(kernel::workspace_routing::eligible(true,s.message_size)) {
        kernel::workspace_routing::adviseOwned(bank.data(),bank.size()*sizeof(std::uint32_t));
        if(!indexed.empty())kernel::workspace_routing::adviseOwned(indexed.data(),indexed.size()*sizeof(std::uint32_t));
    }
    refresh(s.route_seed,s.inner_seed);
}
void BankState::refresh(std::uint64_t routeSeed,std::uint64_t innerSeed) noexcept {
    using kernel::setup::Words;
    using kernel::setup::Divisor;
    Words words(routeSeed);
    const Divisor rowBound(rows),unusedRotation(11);
    for(unsigned g=0;g<256;++g) {
        // Preserve the K18 prototype's word stream (including unused fields).
        (void)words();(void)words();shift[g]=std::uint32_t(rowBound.sample(words,rows));
        (void)unusedRotation.sample(words,11);
    }
    std::iota(region.begin(),region.end(),0);
    for(unsigned n=256;n>1;--n) {
        const Divisor d(n);std::swap(region[n-1],region[d.sample(words,n)]);
    }
    for(auto& key:row_keys)key=(std::uint32_t(words())&0x7ffffff)|1;
    for(std::size_t i=0;i<packed_keys.size();++i) {
        const auto d=row_keys[i];packed_keys[i]=((d&0xffffff)<<2)|(((d>>24)+2)<<26);
    }
    bool fastMasks=false;
#if SPIN_BCH_AVX512
    if(cpu_mask512())fastMasks=kernel::bankMasks512(masks.data(),masks.size(),spec.parameters==Parameters::T128S19?19:12,innerSeed);
#endif
    if(!fastMasks) {
    Words inner(innerSeed);
    const auto mask=(1U<<(spec.parameters==Parameters::T128S19?19:12))-1;
    for(std::size_t i=0;i<masks.size();i+=2) {
        std::uint32_t u;do{u=std::uint32_t(inner())&mask;}while(!u);
        auto v=std::uint32_t(inner())&mask;
        if(std::popcount(u&v)&1)v^=u&-u;
        masks[i]=u;masks[i+1]=v;
    }
    }
    spec.route_seed=routeSeed;spec.inner_seed=innerSeed;
}
struct PreparedState {
    CodeSpec spec;
    SetupOptions setup;
    ExecutionOptions execution;
    std::optional<Code> full;
    std::shared_ptr<BankState> bank;
    std::optional<GenericTranspose> generic;
    bool four=false;
    PreparedState(CodeSpec s,SetupOptions o,ExecutionOptions e):spec(s),setup(o),execution(e) {
        if(!valid_message_size(s.parameters,s.message_size))throw std::invalid_argument("SPIN invalid prepared geometry");
        if(!capabilities().avx2)throw std::runtime_error("SPIN requires AVX2 with OS support");
        if(o.mode==SetupMode::Full)full.emplace(s,e);
        else if(o.mode==SetupMode::BankedHeuristic) {
            if(e.backend!=Backend::Automatic && e.backend!=Backend::Avx2 && e.backend!=Backend::Avx512)
                throw std::invalid_argument("SPIN unknown backend");
            if(e.backend==Backend::Avx512 && !capabilities().avx512_bch)
                throw std::invalid_argument("SPIN AVX512 backend unavailable");
            four=e.backend!=Backend::Avx2 && capabilities().avx512_bch;
            bank=std::make_shared<BankState>(s,o.bank_seed);
        } else throw std::invalid_argument("SPIN unknown setup mode");
    }
};
struct PreparedScratch {
    std::shared_ptr<PreparedState> owner;
    std::optional<spin::Workspace> full;
    std::vector<storage::block> values;
    std::vector<std::uint32_t> addresses;
    explicit PreparedScratch(std::shared_ptr<PreparedState> s):owner(std::move(s)) {
        if(owner->full)full.emplace(owner->full->make_workspace());
        else {
            values.resize(2*owner->spec.message_size);addresses.resize(owner->bank->rows);
            if(kernel::workspace_routing::eligible(true,owner->spec.message_size))
                kernel::workspace_routing::adviseOwned(values.data(),values.size()*sizeof(storage::block));
        }
    }
};
}
namespace spin {
PreparedEncoder::PreparedEncoder(CodeSpec s,SetupOptions o,ExecutionOptions e)
    :state_(std::make_shared<detail::PreparedState>(s,o,e)) {}
PreparedEncoder::~PreparedEncoder()=default;
PreparedEncoder::PreparedEncoder(PreparedEncoder&&) noexcept=default;
PreparedEncoder& PreparedEncoder::operator=(PreparedEncoder&&) noexcept=default;
void PreparedEncoder::setCodeSeed(CodeSeed seed) {
    if(seed==CodeSeed{state_->spec.route_seed,state_->spec.inner_seed})return;
    if(state_->bank)state_->bank->refresh(seed.route,seed.inner);
    else {
        auto spec=state_->spec;spec.route_seed=seed.route;spec.inner_seed=seed.inner;
        Code next(spec,state_->execution);
        // Build before committing, preserving the old encoder on allocation failure.
        if(state_->generic) {
            auto replacement=next.generic_transpose();
            auto live=std::const_pointer_cast<GenericTranspose::Data>(state_->generic->data_);
            auto fresh=std::const_pointer_cast<GenericTranspose::Data>(replacement.data_);
            live->route=std::move(fresh->route);live->masks=std::move(fresh->masks);
        }
        state_->full=std::move(next);
    }
    state_->spec.route_seed=seed.route;state_->spec.inner_seed=seed.inner;
}
CodeSpec PreparedEncoder::specification() const noexcept{return state_->spec;}
SetupOptions PreparedEncoder::setup_options() const noexcept{return state_->setup;}
std::size_t PreparedEncoder::message_size() const noexcept{return state_->spec.message_size;}
std::size_t PreparedEncoder::code_size() const noexcept{return 2*message_size();}
std::size_t PreparedEncoder::setup_bytes() const noexcept {
    if(state_->bank)return state_->bank->bytes();
    return state_->full->setup_bytes()+
        (state_->generic?state_->generic->setup_bytes():0);
}
std::array<std::byte,56> PreparedEncoder::descriptor() const noexcept {
    std::array<std::byte,56> out{};
    out[0]=std::byte{'S'};out[1]=std::byte{'P'};out[2]=std::byte{'I'};out[3]=std::byte{'N'};
    auto put=[&](unsigned pos,std::uint64_t x,unsigned n){for(unsigned i=0;i<n;++i)out[pos+i]=std::byte((x>>(8*i))&255);};
    put(4,2,4);put(8,std::uint32_t(state_->spec.parameters),4);put(12,std::uint32_t(state_->setup.mode),4);
    put(16,message_size(),8);put(24,state_->spec.route_seed,8);put(32,state_->spec.inner_seed,8);
    if(state_->bank){put(40,state_->setup.bank_seed,8);put(48,1,8);}
    return out;
}
PreparedEncoder::Workspace::Workspace(std::shared_ptr<detail::PreparedState> s)
    :scratch_(std::make_unique<detail::PreparedScratch>(std::move(s))) {}
PreparedEncoder::Workspace::~Workspace()=default;
PreparedEncoder::Workspace::Workspace(Workspace&&) noexcept=default;
PreparedEncoder::Workspace& PreparedEncoder::Workspace::operator=(Workspace&&) noexcept=default;
std::size_t PreparedEncoder::Workspace::bytes() const noexcept {
    if(!scratch_)return 0;
    return scratch_->full?scratch_->full->bytes():
        scratch_->values.capacity()*16+scratch_->addresses.capacity()*4;
}
PreparedEncoder::Workspace PreparedEncoder::make_workspace() const{return Workspace(state_);}
GenericTranspose PreparedEncoder::generic_transpose() const {
    if(!state_->generic) {
        if(state_->full)state_->generic.emplace(state_->full->generic_transpose());
        else {
            GenericTranspose view(state_->spec.parameters,message_size(),{},{});
            std::const_pointer_cast<GenericTranspose::Data>(view.data_)->bank=state_->bank;
            state_->generic.emplace(std::move(view));
        }
    }
    return *state_->generic;
}
bool PreparedEncoder::owns(const GenericTranspose& view) const noexcept {
    return state_->generic && state_->generic->data_==view.data_;
}
static bool overlaps(const void* a,std::size_t na,const void* b,std::size_t nb) {
    const auto x=reinterpret_cast<std::uintptr_t>(a),y=reinterpret_cast<std::uintptr_t>(b);
    return x<=y?y-x<na:x-y<nb;
}
void PreparedEncoder::transpose_bytes(std::span<const std::byte> in,std::span<std::byte> out,Workspace& w) const {
    if(overlaps(in.data(),in.size(),out.data(),out.size()))
        throw std::invalid_argument("SPIN buffers overlap; use transpose_inplace");
    if(!w.scratch_ || w.scratch_->owner!=state_)throw std::invalid_argument("SPIN prepared workspace mismatch");
    if(in.size()!=code_size()*16 || out.size()!=message_size()*16 || !in.data() || !out.data() ||
       reinterpret_cast<std::uintptr_t>(in.data())%16 || reinterpret_cast<std::uintptr_t>(out.data())%16)
        throw std::invalid_argument("SPIN prepared buffer geometry or alignment mismatch");
    auto& s=*w.scratch_;
    if(state_->full){state_->full->prepare_workspace(*s.full);state_->full->transpose_bytes(in,out,*s.full);}
    else detail::kernel::bankTranspose(*state_->bank,state_->four,in.data(),out.data(),s.values.data(),s.addresses.data());
}
void PreparedEncoder::transpose_inplace_bytes(std::span<std::byte> buf,Workspace& w) const {
    if(!w.scratch_ || w.scratch_->owner!=state_)throw std::invalid_argument("SPIN prepared workspace mismatch");
    if(buf.size()!=code_size()*16 || !buf.data() || reinterpret_cast<std::uintptr_t>(buf.data())%16)
        throw std::invalid_argument("SPIN prepared buffer geometry or alignment mismatch");
    auto& s=*w.scratch_;
    if(state_->full){state_->full->prepare_workspace(*s.full);state_->full->transpose_inplace_bytes(buf,*s.full);}
    else detail::kernel::bankTranspose(*state_->bank,state_->four,buf.data(),buf.data(),s.values.data(),s.addresses.data());
}
}
