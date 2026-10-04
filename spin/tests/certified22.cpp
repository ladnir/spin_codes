#include <spin/Code.h>
#include <spin/PacketCode.h>
#include <spin/PreparedEncoder.h>
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <new>
#include <span>
#include <utility>
#ifdef _MSC_VER
#include <malloc.h>
#endif

// Count both ordinary and aligned allocations across the public hot path.
static std::atomic<std::size_t> allocations=0;
void* operator new(std::size_t n) {++allocations;if(auto* p=std::malloc(n?n:1))return p;throw std::bad_alloc();}
void* operator new[](std::size_t n) {return ::operator new(n);}
void operator delete(void* p) noexcept {std::free(p);}
void operator delete[](void* p) noexcept {std::free(p);}
void operator delete(void* p,std::size_t) noexcept {std::free(p);}
void operator delete[](void* p,std::size_t) noexcept {std::free(p);}
void* operator new(std::size_t n,std::align_val_t a) {
    ++allocations;
#ifdef _MSC_VER
    auto* p=_aligned_malloc(n?n:1,std::size_t(a));
#else
    void* p=nullptr;if(posix_memalign(&p,std::size_t(a),n?n:1))p=nullptr;
#endif
    if(p)return p;throw std::bad_alloc();
}
void* operator new[](std::size_t n,std::align_val_t a) {return ::operator new(n,a);}
void operator delete(void* p,std::align_val_t) noexcept {
#ifdef _MSC_VER
    _aligned_free(p);
#else
    std::free(p);
#endif
}
void operator delete[](void* p,std::align_val_t a) noexcept {::operator delete(p,a);}
void operator delete(void* p,std::size_t,std::align_val_t a) noexcept {::operator delete(p,a);}
void operator delete[](void* p,std::size_t,std::align_val_t a) noexcept {::operator delete(p,a);}

namespace {
constexpr auto family=spin::Parameters::PacketRsT64S22;
static_assert(static_cast<std::uint32_t>(family)==7);
struct Block {std::uint64_t words[2]{};bool operator==(const Block&) const=default;};
struct Context {const char* test="metadata";std::size_t k=0;std::uint64_t route=0,inner=0;unsigned width=16,offset=0;} context;
void require(bool value,const char* message) {
    if(value)return;
    std::cerr<<"certified22 context: "<<context.test<<" K="<<context.k<<" seeds="<<context.route<<','<<context.inner
        <<" width="<<8*context.width<<" offset="<<context.offset<<'\n';
    throw std::runtime_error(message);
}
template<class F> void rejects(F f) {
    bool rejected=false;try{f();}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"invalid operation accepted");
}
template<class F> void noAlloc(F f) {
    const auto before=allocations.load();f();require(allocations.load()==before,"encoding/rebind allocated");
}
bool equal(std::span<const std::byte> a,std::span<const std::byte> b) {
    return a.size()==b.size() && !std::memcmp(a.data(),b.data(),a.size());
}
void compare(std::span<const std::byte> a,std::span<const std::byte> b,const char* message) {
    if(equal(a,b))return;
    std::size_t first=0;while(first<std::min(a.size(),b.size()) && a[first]==b[first])++first;
    std::cerr<<"first mismatching byte="<<first<<" record="<<first/context.width<<'\n';require(false,message);
}
void fill(std::span<std::byte> bytes,std::uint64_t seed) {
    require(bytes.size()%8==0,"partial fill word");
    for(std::size_t i=0;i<bytes.size();i+=8) {
        seed+=0x9e3779b97f4a7c15ULL;auto z=seed;
        z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;z=(z^(z>>27))*0x94d049bb133111ebULL;z^=z>>31;
        std::memcpy(bytes.data()+i,&z,8);
    }
}
void guards(std::span<const std::byte> storage,std::size_t offset,std::size_t bytes) {
    const auto marker=[](auto b){return b==std::byte{0xa5};};
    require(std::all_of(storage.begin(),storage.begin()+offset,marker),"prefix guard changed");
    require(std::all_of(storage.begin()+offset+bytes,storage.end(),marker),"suffix guard changed");
}
Block dot(std::span<const std::byte> a,std::span<const std::byte> b) {
    require(a.size()==b.size() && a.size()%16==0,"adjoint shape");Block result;
    for(std::size_t i=0;i<a.size();i+=16) {
        Block x,y;std::memcpy(&x,a.data()+i,16);std::memcpy(&y,b.data()+i,16);
        result.words[0]^=x.words[0]&y.words[0];result.words[1]^=x.words[1]&y.words[1];
    }
    return result;
}
void metadata() {
    require(spin::message_alignment(family)==128,"natural size alignment");
    constexpr std::size_t maximum=(std::numeric_limits<std::uint32_t>::max()/(260ULL*16))*128;
    for(auto k:{std::size_t{128},std::size_t{256},std::size_t{384},std::size_t{2048},std::size_t{65536},
                std::size_t{262144},std::size_t{1048576},maximum})
        noAlloc([&]{require(spin::valid_message_size(family,k),"natural size rejected");});
    for(auto k:{std::size_t{0},std::size_t{1},std::size_t{127},std::size_t{129},maximum+128,
                std::numeric_limits<std::size_t>::max()}) {
        noAlloc([&]{require(!spin::valid_message_size(family,k),"invalid/overflow size accepted");});
        rejects([&]{spin::Code code({k,family,17,43});});
    }
    spin::Code code({384,family,17,43}),scalar(code.specification(),{spin::Backend::Portable});
    require(code.message_size()==384 && code.code_size()==768 && code.setup_bytes()>0,"shape/setup size");
    require(code.descriptor()==scalar.descriptor(),"backend changes map descriptor");
    const auto descriptor=code.descriptor();require(descriptor[8]==std::byte{7},"descriptor family ID");
    for(unsigned i=0;i<8;++i) {
        require(descriptor[16+i]==std::byte((384ULL>>(8*i))&255),"descriptor K");
        require(descriptor[24+i]==std::byte((17ULL>>(8*i))&255),"descriptor route seed");
        require(descriptor[32+i]==std::byte((43ULL>>(8*i))&255),"descriptor inner seed");
    }
    require(code.backend()==(spin::packet_fast_available()?spin::Backend::Avx512:spin::Backend::Portable),"automatic backend");
    for(auto width:{spin::Width::Bits128,spin::Width::Bits256,spin::Width::Bits512})
        require(code.supports_forward(width) && scalar.supports_forward(width),"forward width unsupported");
    require(code.supports_transpose() && !code.supports_transpose(spin::Width::Bits256) &&
        !code.supports_transpose(spin::Width::Bits512) && !code.supports_generic_transpose(),"transpose capabilities");
    rejects([&]{code.generic_transpose();});
    rejects([&]{code.make_workspace(static_cast<spin::Width>(17));});
    std::array<std::uint64_t,6> packedMessage{};
    std::array<std::uint64_t,12> packedOutput{},packedScratch{};
    rejects([&]{code.forward_bits(packedMessage,packedOutput,packedScratch);});
    rejects([&]{spin::Code wrong({384,family,17,43},{spin::Backend::Avx2});});
    if(!spin::packet_fast_available()) {
        bool failed=false;try{spin::Code forced({384,family,17,43},{spin::Backend::Avx512});}
        catch(const std::runtime_error&){failed=true;}require(failed,"unavailable ISA accepted");
    }
}
void transposeCase(std::size_t k,std::uint64_t route,std::uint64_t inner) {
    const auto n=2*k;spin::Code code({k,family,route,inner}),scalar(code.specification(),{spin::Backend::Portable});
    context={"transpose/adjoint",k,route,inner};auto work=code.make_workspace(),sw=scalar.make_workspace();
    spin::Buffer source(n*16),input(n*16+128),output(k*16+128),expected(k*16),message(k*16),encoded(n*16);
    fill(source.bytes(),913);fill(message.bytes(),1949);
    noAlloc([&]{scalar.transpose_bytes(source.bytes(),expected.bytes(),sw);});
    for(unsigned offset:{0U,16U,32U,48U}) {
        context.offset=offset;std::fill(input.bytes().begin(),input.bytes().end(),std::byte{0xa5});
        std::fill(output.bytes().begin(),output.bytes().end(),std::byte{0xa5});
        auto in=input.bytes().subspan(64+offset,n*16),out=output.bytes().subspan(64+offset,k*16);
        std::memcpy(in.data(),source.bytes().data(),in.size());
        noAlloc([&]{code.transpose_bytes(in,out,work);});compare(out,expected.bytes(),"scalar/automatic transpose mismatch");
        compare(in,source.bytes(),"transpose modified input");guards(input.bytes(),64+offset,n*16);guards(output.bytes(),64+offset,k*16);
        noAlloc([&]{code.transpose_inplace_bytes(in,work);});compare(in.first(k*16),expected.bytes(),"in-place prefix mismatch");
        compare(in.subspan(k*16),source.bytes().subspan(k*16),"in-place suffix changed");guards(input.bytes(),64+offset,n*16);
        // The first natural size above the streaming threshold tests every
        // output alignment, including cached fallback for unaligned NT output.
        if(k>65536 && k!=262272)break;
    }
    noAlloc([&]{code.forward_bytes(message.bytes(),encoded.bytes(),work);});
    require(dot(encoded.bytes(),source.bytes())==dot(message.bytes(),expected.bytes()),"128 independent bitwise adjoints fail");
    noAlloc([&]{code.transpose_bytes(source.bytes(),output.bytes().first(k*16),work);});
    compare(output.bytes().first(k*16),expected.bytes(),"forward/transpose dirty-workspace reuse");
    if(k>2048)return;
    rejects([&]{code.transpose_bytes(source.bytes(),source.bytes().first(k*16),work);});
    rejects([&]{code.transpose_bytes(source.bytes(),source.bytes().subspan(16,k*16),work);});
    rejects([&]{code.transpose_bytes(input.bytes().subspan(1,n*16),expected.bytes(),work);});
    rejects([&]{code.transpose_bytes(source.bytes(),output.bytes().subspan(1,k*16),work);});
    rejects([&]{code.transpose_inplace_bytes(input.bytes().subspan(1,n*16),work);});
    rejects([&]{code.transpose_inplace_bytes(source.bytes().first(n*16-1),work);});
    auto wide=code.make_workspace(spin::Width::Bits256);
    rejects([&]{code.transpose_bytes(source.bytes(),expected.bytes(),wide);});
    rejects([&]{code.transpose_inplace_bytes(source.bytes(),wide);});
    for(auto coordinate:{std::size_t{0},std::size_t{63},std::size_t{64},n/2,n-65,n-1}) {
        std::fill(source.bytes().begin(),source.bytes().end(),std::byte{});
        source.bytes()[16*coordinate]=std::byte{1};source.bytes()[16*coordinate+15]=std::byte{128};
        noAlloc([&]{scalar.transpose_bytes(source.bytes(),expected.bytes(),sw);code.transpose_bytes(source.bytes(),output.bytes().first(k*16),work);});
        compare(output.bytes().first(k*16),expected.bytes(),"sparse physical-boundary mismatch");
    }
}
void forwardWidths(spin::Width width) {
    constexpr std::size_t k=384,n=2*k;const auto bytes=static_cast<unsigned>(width);
    spin::Code code({k,family,17,43}),scalar(code.specification(),{spin::Backend::Portable});
    context={"forward widths",k,17,43,bytes};auto work=code.make_workspace(width),sw=scalar.make_workspace(width),narrow=code.make_workspace();
    spin::Buffer input(k*bytes+128),original(k*bytes),output(n*bytes+128),expected(n*bytes),plane(k*16),encoded(n*16);
    auto in=input.bytes().subspan(16,k*bytes);
    for(unsigned pattern=0;pattern<3;++pattern) {
        fill(in,1949);if(pattern)std::fill(in.begin(),in.end(),std::byte{});
        if(pattern==2){in.front()=std::byte{1};in.back()=std::byte{128};}
        std::memcpy(original.bytes().data(),in.data(),in.size());
        noAlloc([&]{scalar.forward_bytes(in,expected.bytes(),sw);});
        if(pattern==1)require(std::all_of(expected.bytes().begin(),expected.bytes().end(),[](auto b){return b==std::byte{};}),"nonzero image of zero");
        for(unsigned offset:{0U,16U,32U,48U})for(auto stores:{spin::OutputStores::Streaming,spin::OutputStores::Cached}) {
            context.offset=offset;std::fill(output.bytes().begin(),output.bytes().end(),std::byte{0xa5});
            auto out=output.bytes().subspan(64+offset,n*bytes);
            noAlloc([&]{code.forward_bytes(in,out,work,stores);});compare(out,expected.bytes(),"forward fallback mismatch");
            compare(in,original.bytes(),"forward modified input");guards(output.bytes(),64+offset,n*bytes);
        }
        if(pattern==0)for(unsigned lane=0;lane<bytes/16;++lane) {
            for(std::size_t i=0;i<k;++i)std::memcpy(plane.bytes().data()+16*i,in.data()+bytes*i+16*lane,16);
            noAlloc([&]{code.forward_bytes(plane.bytes(),encoded.bytes(),narrow);});
            for(std::size_t i=0;i<n;++i)require(!std::memcmp(expected.bytes().data()+bytes*i+16*lane,encoded.bytes().data()+16*i,16),"wide lane differs from narrow map");
        }
    }
    rejects([&]{code.forward_bytes(in.first(in.size()-1),expected.bytes(),work);});
    rejects([&]{code.forward_bytes(in,expected.bytes().first(n*bytes-1),work);});
    rejects([&]{code.forward_bytes(input.bytes().subspan(1,k*bytes),expected.bytes(),work);});
    rejects([&]{code.forward_bytes(in,output.bytes().subspan(1,n*bytes),work);});
    rejects([&]{code.forward_bytes(output.bytes().first(k*bytes),output.bytes().first(n*bytes),work);});
}
void ownershipAndSeeds() {
    constexpr std::size_t k=384,n=2*k;context={"ownership/seeds",k,17,43};
    spin::Code base({k,family,17,43});auto work=base.make_workspace();
    spin::Buffer input(n*16),reference(k*16),actual(k*16),expected(k*16);fill(input.bytes(),913);
    base.transpose_bytes(input.bytes(),reference.bytes(),work);
    spin::Code same(base.specification());rejects([&]{same.transpose_bytes(input.bytes(),actual.bytes(),work);});
    for(auto seeds:std::array<std::pair<std::uint64_t,std::uint64_t>,2>{{{18,43},{17,44}}}) {
        spin::Code changed({k,family,seeds.first,seeds.second});auto fresh=changed.make_workspace();
        changed.transpose_bytes(input.bytes(),expected.bytes(),fresh);
        require(base.descriptor()!=changed.descriptor() && !equal(reference.bytes(),expected.bytes()),"seed has no effect");
        noAlloc([&]{changed.prepare_workspace(work);changed.transpose_bytes(input.bytes(),actual.bytes(),work);});
        compare(actual.bytes(),expected.bytes(),"rebind retained stale setup");
        rejects([&]{base.transpose_bytes(input.bytes(),actual.bytes(),work);});
        auto moved=std::move(changed);auto retained=moved;auto movedWork=std::move(work);
        rejects([&]{changed.transpose_bytes(input.bytes(),actual.bytes(),movedWork);});
        rejects([&]{moved.transpose_bytes(input.bytes(),actual.bytes(),work);});
        noAlloc([&]{retained.transpose_bytes(input.bytes(),actual.bytes(),movedWork);});
        compare(actual.bytes(),expected.bytes(),"move/share invalidates workspace");work=std::move(movedWork);
    }
    spin::Code old({2048,spin::Parameters::PacketRsT64S20,17,43}),small({2048,family,17,43});
    require(old.descriptor()!=small.descriptor(),"s22 aliases larger-outer descriptor");
    auto otherFamily=old.make_workspace();spin::Buffer largerInput(4096*16),largerOutput(2048*16);
    rejects([&]{small.transpose_bytes(largerInput.bytes(),largerOutput.bytes(),otherFamily);});
    // Cross-family rebinding may allocate, but must discard the old geometry.
    small.prepare_workspace(otherFamily);
    noAlloc([&]{small.transpose_bytes(largerInput.bytes(),largerOutput.bytes(),otherFamily);});
    // These profiles share the 260-record group layout. Cross-family rebind
    // should reuse its allocation but must replace the sampled inner entirely.
    spin::Code paired({65536,spin::Parameters::PacketRsT64S15K16,17,43}),certified({65536,family,17,43});
    auto shared=paired.make_workspace(),own=certified.make_workspace();
    spin::Buffer sharedInput(131072*16),sharedOutput(65536*16),sharedExpected(65536*16);fill(sharedInput.bytes(),913);
    certified.transpose_bytes(sharedInput.bytes(),sharedExpected.bytes(),own);
    noAlloc([&]{certified.prepare_workspace(shared);certified.transpose_bytes(sharedInput.bytes(),sharedOutput.bytes(),shared);});
    compare(sharedOutput.bytes(),sharedExpected.bytes(),"paired15 to certified22 rebind map mismatch");
    auto pairedOwn=paired.make_workspace();paired.transpose_bytes(sharedInput.bytes(),sharedExpected.bytes(),pairedOwn);
    noAlloc([&]{paired.prepare_workspace(shared);paired.transpose_bytes(sharedInput.bytes(),sharedOutput.bytes(),shared);});
    compare(sharedOutput.bytes(),sharedExpected.bytes(),"certified22 to paired15 rebind map mismatch");
}
void preparedFull() {
    constexpr std::size_t k=384,n=2*k;
    rejects([&]{spin::PreparedEncoder bad({k,family,17,43},{spin::SetupMode::BankedHeuristic});});
    for(auto backend:{spin::Backend::Portable,spin::Backend::Automatic}) {
        context={"prepared Full",k,17,43};spin::PreparedEncoder prepared({k,family,17,43},{spin::SetupMode::Full},{backend});
        auto work=prepared.make_workspace();const auto bytes=work.bytes();
        require(prepared.supports_transpose() && !prepared.supports_generic_transpose(),"prepared capabilities");
        rejects([&]{prepared.generic_transpose();});rejects([&]{prepared.make_workspace(spin::Width::Bits256);});
        spin::Buffer input(n*16),actual(k*16),expected(k*16),previous(k*16),inplace(n*16);fill(input.bytes(),913);
        bool havePrevious=false;
        for(auto seed:{spin::CodeSeed{17,43},spin::CodeSeed{18,43},spin::CodeSeed{18,44}}) {
            prepared.setCodeSeed(seed);noAlloc([&]{prepared.setCodeSeed(seed);});
            spin::Code fresh({k,family,seed.route,seed.inner},{spin::Backend::Portable});auto fw=fresh.make_workspace();
            fresh.transpose_bytes(input.bytes(),expected.bytes(),fw);
            noAlloc([&]{prepared.transpose_bytes(input.bytes(),actual.bytes(),work);});
            compare(actual.bytes(),expected.bytes(),"prepared refreshed map mismatch");require(work.bytes()==bytes,"refresh resized scratch");
            if(havePrevious)require(!equal(actual.bytes(),previous.bytes()),"refresh did not change map");
            std::memcpy(previous.bytes().data(),actual.bytes().data(),actual.bytes().size());havePrevious=true;
            std::memcpy(inplace.bytes().data(),input.bytes().data(),input.bytes().size());
            noAlloc([&]{prepared.transpose_inplace_bytes(inplace.bytes(),work);});
            compare(inplace.bytes().first(k*16),expected.bytes(),"prepared in-place prefix");
            compare(inplace.bytes().subspan(k*16),input.bytes().subspan(k*16),"prepared in-place suffix");
            const auto descriptor=prepared.descriptor();
            require(descriptor[4]==std::byte{2} && descriptor[8]==std::byte{7} && descriptor[12]==std::byte{},"prepared descriptor");
            for(unsigned i=0;i<8;++i) {
                require(descriptor[24+i]==std::byte((seed.route>>(8*i))&255),"prepared route descriptor");
                require(descriptor[32+i]==std::byte((seed.inner>>(8*i))&255),"prepared inner descriptor");
            }
        }
        prepared.setCodeSeed({91,17});spin::Code fresh({k,family,91,17},{spin::Backend::Portable});auto fw=fresh.make_workspace();
        fresh.transpose_bytes(input.bytes(),expected.bytes(),fw);
        noAlloc([&]{prepared.prepare_workspace(work);prepared.transpose_bytes(input.bytes(),actual.bytes(),work);});
        compare(actual.bytes(),expected.bytes(),"explicit prepared rebind");
        spin::PreparedEncoder other({k,family,91,17},{spin::SetupMode::Full},{backend});auto wrong=other.make_workspace();
        rejects([&]{prepared.prepare_workspace(wrong);});rejects([&]{prepared.transpose_bytes(input.bytes(),actual.bytes(),wrong);});
        auto moved=std::move(prepared);auto movedWork=std::move(work);
        require(!prepared.supports_transpose() && prepared.descriptor()==std::array<std::byte,56>{},"moved-from metadata");
        rejects([&]{prepared.setCodeSeed({17,43});});rejects([&]{prepared.transpose_bytes(input.bytes(),actual.bytes(),movedWork);});
        rejects([&]{moved.prepare_workspace(work);});
        noAlloc([&]{moved.prepare_workspace(movedWork);moved.transpose_bytes(input.bytes(),actual.bytes(),movedWork);});
        compare(actual.bytes(),expected.bytes(),"prepared move lost setup");
    }
}
}
int main() {try {
    metadata();
    for(auto k:{std::size_t{128},std::size_t{256},std::size_t{384},std::size_t{1920},std::size_t{2048},
                std::size_t{2176},std::size_t{65536},std::size_t{262016},std::size_t{262144},
                std::size_t{262272},std::size_t{1048576}})transposeCase(k,17,43);
    transposeCase(2048,1,1);transposeCase(2048,913,1123);
    for(auto width:{spin::Width::Bits128,spin::Width::Bits256,spin::Width::Bits512})forwardWidths(width);
    ownershipAndSeeds();preparedFull();
    std::cout<<"certified22 API PASS: natural sizing, scalar/automatic transpose, wide forward lanes, bitwise adjoints, guards, seeds, ownership, refresh, allocation-free encoding\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
