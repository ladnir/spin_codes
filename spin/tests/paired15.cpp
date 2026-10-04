#include <spin/Code.h>
#include <spin/PacketCode.h>
#include <spin/PreparedEncoder.h>
#include "forward_output_stores.h"
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
#include <vector>
#ifdef _MSC_VER
#include <malloc.h>
#endif

// Encoding must remain allocation-free at every supported record width.
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
constexpr std::size_t k=std::size_t{1}<<16,n=2*k;
constexpr auto family=spin::Parameters::PacketRsT64S15K16;
struct alignas(16) Block {std::uint64_t words[2]{};bool operator==(const Block&) const=default;};
static_assert(static_cast<std::uint32_t>(family)==6);
struct TestContext {
    const char* test="metadata";
    std::uint64_t routeSeed=0,innerSeed=0;
    unsigned width=0,pattern=0,alignment=0;
    spin::Backend backend=spin::Backend::Automatic;
} context;
const char* backendName(spin::Backend backend) {
    switch(backend) {
    case spin::Backend::Automatic:return "Automatic";
    case spin::Backend::Avx2:return "Avx2";
    case spin::Backend::Avx512:return "Avx512";
    case spin::Backend::Portable:return "Portable";
    }
    return "Unknown";
}
void require(bool value,const char* message) {
    if(value)return;
    std::cerr<<"context: test="<<context.test<<" route="<<context.routeSeed<<" inner="<<context.innerSeed
        <<" width_bits="<<8*context.width<<" backend="<<backendName(context.backend)
        <<" pattern="<<context.pattern<<" alignment="<<context.alignment<<'\n';
    throw std::runtime_error(message);
}
template<class F> void rejects(F f) {
    bool rejected=false;try{f();}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"paired15 invalid call accepted");
}
template<class F> void noAlloc(F f) {
    const auto before=allocations.load();f();require(allocations.load()==before,"paired15 encoding allocated");
}
void fill(std::span<std::byte> bytes,std::uint64_t seed=913) {
    require(bytes.size()%8==0,"fill requires complete words");
    for(std::size_t i=0;i<bytes.size();i+=8) {
        seed+=0x9e3779b97f4a7c15ULL;auto z=seed;
        z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;z=(z^(z>>27))*0x94d049bb133111ebULL;z^=z>>31;
        std::memcpy(bytes.data()+i,&z,8);
    }
}
bool equal(std::span<const std::byte> a,std::span<const std::byte> b) {
    return a.size()==b.size() && !std::memcmp(a.data(),b.data(),a.size());
}
void requireEqual(std::span<const std::byte> actual,std::span<const std::byte> expected,const char* message) {
    if(equal(actual,expected))return;
    std::size_t first=0;
    while(first<std::min(actual.size(),expected.size()) && actual[first]==expected[first])++first;
    std::cerr<<"mismatch: first_byte="<<first<<" record="<<(context.width?first/context.width:0)
        <<" actual_bytes="<<actual.size()<<" expected_bytes="<<expected.size();
    if(first<actual.size())std::cerr<<" actual="<<std::to_integer<unsigned>(actual[first]);
    if(first<expected.size())std::cerr<<" expected="<<std::to_integer<unsigned>(expected[first]);
    std::cerr<<'\n';require(false,message);
}
void guards(std::span<const std::byte> storage,std::size_t offset,std::size_t bytes) {
    const auto marker=[](auto x){return x==std::byte{0xa5};};
    require(std::all_of(storage.begin(),storage.begin()+offset,marker),"paired15 prefix guard");
    require(std::all_of(storage.begin()+offset+bytes,storage.end(),marker),"paired15 suffix guard");
}
Block dot(std::span<const std::byte> a,std::span<const std::byte> b) {
    require(a.size()==b.size() && a.size()%16==0,"paired15 adjoint shape");
    Block result{};
    for(std::size_t i=0;i<a.size();i+=16) {
        Block x,y;std::memcpy(&x,a.data()+i,16);std::memcpy(&y,b.data()+i,16);
        for(unsigned j=0;j<2;++j)result.words[j]^=x.words[j]&y.words[j];
    }
    return result;
}
void metadata() {
    require(spin::valid_message_size(family,k),"paired15 K16 rejected");
    for(auto invalid:{std::size_t{0},k-1,k+1,k/2,2*k,4*k,std::numeric_limits<std::size_t>::max()}) {
        require(!spin::valid_message_size(family,invalid),"paired15 unsupported geometry accepted");
        rejects([&]{spin::Code code({invalid,family,17,43});});
    }
    spin::Code code({k,family,17,43}),scalar(code.specification(),{spin::Backend::Portable});
    spin::Code old({k,spin::Parameters::PacketRsT64S20,17,43});
    require(code.message_size()==k && code.code_size()==n && code.setup_bytes()>0,"paired15 shape");
    require(code.descriptor()==scalar.descriptor(),"paired15 backend affects map identifier");
    require(code.descriptor()!=old.descriptor(),"paired15 aliases old profile descriptor");
    const auto descriptor=code.descriptor();
    require(descriptor[8]==std::byte{6},"paired15 public descriptor family");
    for(unsigned i=0;i<8;++i) {
        require(descriptor[16+i]==std::byte((k>>(8*i))&255),"paired15 descriptor K");
        require(descriptor[24+i]==std::byte((17ULL>>(8*i))&255),"paired15 descriptor route seed");
        require(descriptor[32+i]==std::byte((43ULL>>(8*i))&255),"paired15 descriptor inner seed");
    }
    require(code.backend()==(spin::packet_fast_available()?spin::Backend::Avx512:spin::Backend::Portable),"paired15 dispatch");
    for(auto width:{spin::Width::Bits128,spin::Width::Bits256,spin::Width::Bits512})
        require(code.supports_forward(width) && scalar.supports_forward(width),"paired15 forward width");
    require(code.supports_transpose() && !code.supports_transpose(spin::Width::Bits256) &&
        !code.supports_transpose(spin::Width::Bits512) && !code.supports_generic_transpose(),"paired15 capabilities");
    rejects([&]{code.make_workspace(static_cast<spin::Width>(17));});
    rejects([&]{spin::Code wrong({k,family,17,43},{spin::Backend::Avx2});});
    if(!spin::packet_fast_available()) {
        bool failed=false;try{spin::Code forced({k,family,17,43},{spin::Backend::Avx512});}
        catch(const std::runtime_error&){failed=true;}
        require(failed,"paired15 unavailable ISA accepted");
    }
}
void forwardCase(std::uint64_t routeSeed,std::uint64_t innerSeed,spin::Width width) {
    const std::size_t bytes=static_cast<unsigned>(width),lanes=bytes/16;
    spin::Code code({k,family,routeSeed,innerSeed}),scalar(code.specification(),{spin::Backend::Portable});
    context={"forward",routeSeed,innerSeed,unsigned(bytes),0,0,code.backend()};
    auto work=code.make_workspace(width),sw=scalar.make_workspace(width),single=code.make_workspace();
    spin::Buffer input(k*bytes+128),output(n*bytes+128),expected(n*bytes),plane(k*16),encodedPlane(n*16);
    // 16-byte alignment is sufficient even for 256/512-bit records.
    auto in=input.bytes().subspan(16,k*bytes);
    std::vector<std::byte> original(in.size());
    for(unsigned pattern=0;pattern<3;++pattern) {
        context.pattern=pattern;
        fill(in,routeSeed+91);
        if(pattern)std::fill(in.begin(),in.end(),std::byte{});
        if(pattern==2) {in.front()=std::byte{1};in.back()=std::byte{128};}
        std::copy(in.begin(),in.end(),original.begin());
        noAlloc([&]{scalar.forward_bytes(in,expected.bytes(),sw);});
        if(pattern==1)require(std::all_of(expected.bytes().begin(),expected.bytes().end(),[](auto x){return x==std::byte{};}),"paired15 nonzero encoding of zero");
        for(unsigned offset:{0U,16U,32U,48U}) {
            context.alignment=offset;
            std::fill(output.bytes().begin(),output.bytes().end(),std::byte{0xa5});
            auto out=output.bytes().subspan(64+offset,n*bytes);
            noAlloc([&]{code.forward_bytes(in,out,work);});
            requireEqual(out,expected.bytes(),"paired15 scalar vs SIMD forward");
            require(equal(in,original),"paired15 forward input modified");
            guards(output.bytes(),64+offset,n*bytes);
        }
        if(pattern==0 && lanes>1)for(std::size_t lane=0;lane<lanes;++lane) {
            for(std::size_t i=0;i<k;++i)std::memcpy(plane.bytes().data()+16*i,in.data()+bytes*i+16*lane,16);
            noAlloc([&]{code.forward_bytes(plane.bytes(),encodedPlane.bytes(),single);});
            for(std::size_t i=0;i<n;++i)require(!std::memcmp(expected.bytes().data()+bytes*i+16*lane,encodedPlane.bytes().data()+16*i,16),"paired15 wide lane differs from narrow map");
        }
    }
    rejects([&]{code.forward_bytes(in.first(in.size()-1),expected.bytes(),work);});
    rejects([&]{code.forward_bytes(in,expected.bytes().first(n*bytes-1),work);});
    rejects([&]{code.forward_bytes(input.bytes().subspan(1,k*bytes),expected.bytes(),work);});
    rejects([&]{code.forward_bytes(in,output.bytes().subspan(1,n*bytes),work);});
    rejects([&]{code.forward_bytes(output.bytes().first(k*bytes),output.bytes().first(n*bytes),work);});
    rejects([&]{code.forward_bytes(output.bytes().subspan(16,k*bytes),output.bytes().first(n*bytes),work);});
    rejects([&]{code.forward_bytes(output.bytes().first(k*bytes),output.bytes().subspan(16,n*bytes),work);});
    spin::Code sameSpec(code.specification());
    rejects([&]{sameSpec.forward_bytes(in,expected.bytes(),work);});
    spin::Code changed({k,family,routeSeed+1,innerSeed+1});auto cw=changed.make_workspace(width);
    changed.forward_bytes(in,expected.bytes(),cw);
    noAlloc([&]{changed.prepare_workspace(work);changed.forward_bytes(in,output.bytes().first(n*bytes),work);});
    require(equal(output.bytes().first(n*bytes),expected.bytes()),"paired15 rebind left stale workspace");
    rejects([&]{code.forward_bytes(in,expected.bytes(),work);});
    auto moved=std::move(changed);auto retained=moved;auto movedWork=std::move(work);
    rejects([&]{changed.forward_bytes(in,expected.bytes(),movedWork);});
    rejects([&]{moved.forward_bytes(in,expected.bytes(),work);});
    noAlloc([&]{retained.forward_bytes(in,output.bytes().first(n*bytes),movedWork);});
    require(equal(output.bytes().first(n*bytes),expected.bytes()),"paired15 move invalidates workspace");
}
void transposeAndAdjoint(std::uint64_t routeSeed,std::uint64_t innerSeed) {
    spin::Code code({k,family,routeSeed,innerSeed}),scalar(code.specification(),{spin::Backend::Portable});
    context={"transpose/adjoint",routeSeed,innerSeed,16,0,0,code.backend()};
    auto work=code.make_workspace(),sw=scalar.make_workspace();
    spin::Buffer input(n*16+128),output(k*16+128),expected(k*16),message(k*16),encoded(n*16);
    spin::Buffer original(n*16);fill(original.bytes(),192);fill(message.bytes(),997);
    noAlloc([&]{scalar.transpose_bytes(original.bytes(),expected.bytes(),sw);});
    for(unsigned offset:{0U,16U,32U,48U}) {
        context.alignment=offset;
        std::fill(input.bytes().begin(),input.bytes().end(),std::byte{0xa5});
        std::fill(output.bytes().begin(),output.bytes().end(),std::byte{0xa5});
        auto in=input.bytes().subspan(64+offset,n*16),out=output.bytes().subspan(64+offset,k*16);
        std::memcpy(in.data(),original.bytes().data(),in.size());
        noAlloc([&]{code.transpose_bytes(in,out,work);});
        requireEqual(out,expected.bytes(),"paired15 scalar vs SIMD transpose");
        require(equal(in,original.bytes()),"paired15 transpose input modified");
        guards(input.bytes(),64+offset,n*16);guards(output.bytes(),64+offset,k*16);
        noAlloc([&]{code.transpose_inplace_bytes(in,work);});
        requireEqual(in.first(k*16),expected.bytes(),"paired15 in-place transpose prefix");
        require(equal(in.subspan(k*16),original.bytes().subspan(k*16)),"paired15 in-place transpose suffix");
        guards(input.bytes(),64+offset,n*16);
    }
    // Preserve all 128 independent bitwise adjoint identities, not just parity.
    noAlloc([&]{code.forward_bytes(message.bytes(),encoded.bytes(),work);});
    require(dot(encoded.bytes(),original.bytes())==dot(message.bytes(),expected.bytes()),"paired15 binary adjoint identity");
    // Direction changes reuse dirty scratch without initialization.
    noAlloc([&]{code.transpose_bytes(original.bytes(),output.bytes().first(k*16),work);});
    require(equal(output.bytes().first(k*16),expected.bytes()),"paired15 direction reuse");
    rejects([&]{code.transpose_bytes(original.bytes(),original.bytes().first(k*16),work);});
    rejects([&]{code.transpose_bytes(original.bytes(),original.bytes().subspan(16,k*16),work);});
    rejects([&]{code.transpose_inplace_bytes(input.bytes().subspan(1,n*16),work);});
    rejects([&]{code.transpose_inplace_bytes(original.bytes().first(n*16-1),work);});
    auto wide=code.make_workspace(spin::Width::Bits256);
    rejects([&]{code.transpose_bytes(original.bytes(),expected.bytes(),wide);});
    rejects([&]{code.transpose_inplace_bytes(original.bytes(),wide);});
    // Sparse probes independently exercise both ends and physical boundaries.
    for(auto coordinate:{std::size_t{0},std::size_t{63},std::size_t{64},std::size_t{1023},n-65,n-1}) {
        std::fill(original.bytes().begin(),original.bytes().end(),std::byte{});
        original.bytes()[16*coordinate]=std::byte{1};original.bytes()[16*coordinate+15]=std::byte{128};
        noAlloc([&]{scalar.transpose_bytes(original.bytes(),expected.bytes(),sw);code.transpose_bytes(original.bytes(),output.bytes().first(k*16),work);});
        requireEqual(output.bytes().first(k*16),expected.bytes(),"paired15 sparse transpose boundary");
    }
}
void independentSeeds() {
    context={"seed independence",17,43,16,0,0,spin::Backend::Automatic};
    spin::Buffer message(k*16),first(n*16),second(n*16);fill(message.bytes(),19);
    spin::Code base({k,family,17,43});auto work=base.make_workspace();
    base.forward_bytes(message.bytes(),first.bytes(),work);
    for(auto seed:std::array<std::pair<std::uint64_t,std::uint64_t>,2>{{{18,43},{17,44}}}) {
        spin::Code changed({k,family,seed.first,seed.second});auto other=changed.make_workspace();
        require(base.descriptor()!=changed.descriptor(),"paired15 seed missing from descriptor");
        changed.forward_bytes(message.bytes(),second.bytes(),other);
        require(!equal(first.bytes(),second.bytes()),"paired15 setup seed has no effect on code");
    }
    spin::Code old({k,spin::Parameters::PacketRsT64S20,17,43});auto oldWork=old.make_workspace();
    rejects([&]{base.forward_bytes(message.bytes(),second.bytes(),oldWork);});
    rejects([&]{old.forward_bytes(message.bytes(),second.bytes(),work);});
}
void preparedFull() {
    rejects([&]{spin::PreparedEncoder invalid({k,family,17,43},{spin::SetupMode::BankedHeuristic});});
    for(auto backend:{spin::Backend::Portable,spin::Backend::Automatic}) {
        context={"prepared Full",17,43,16,0,0,backend};
        spin::PreparedEncoder prepared({k,family,17,43},{spin::SetupMode::Full},{backend});
        auto work=prepared.make_workspace();const auto scratchBytes=work.bytes();
        require(prepared.supports_transpose() && !prepared.supports_generic_transpose(),"paired15 prepared capabilities");
        rejects([&]{prepared.generic_transpose();});
        rejects([&]{prepared.make_workspace(spin::Width::Bits256);});
        spin::Buffer input(n*16),actual(k*16),expected(k*16),previous(k*16),inplace(n*16);
        spin::Buffer message(k*16),encoded(n*16);fill(input.bytes(),131);fill(message.bytes(),997);
        bool havePrevious=false;
        for(auto seed:{spin::CodeSeed{17,43},spin::CodeSeed{18,43},spin::CodeSeed{18,44}}) {
            context.routeSeed=seed.route;context.innerSeed=seed.inner;
            prepared.setCodeSeed(seed);
            noAlloc([&]{prepared.setCodeSeed(seed);});
            spin::Code fresh({k,family,seed.route,seed.inner},{spin::Backend::Portable});
            auto fw=fresh.make_workspace();fresh.transpose_bytes(input.bytes(),expected.bytes(),fw);
            // Full mode automatically rebinds an existing workspace after refresh.
            noAlloc([&]{prepared.transpose_bytes(input.bytes(),actual.bytes(),work);});
            requireEqual(actual.bytes(),expected.bytes(),"paired15 prepared refresh mismatch");
            require(work.bytes()==scratchBytes,"paired15 prepared refresh resized scratch");
            if(havePrevious)require(!equal(actual.bytes(),previous.bytes()),"paired15 prepared seed refresh has no effect");
            std::memcpy(previous.bytes().data(),actual.bytes().data(),actual.bytes().size());havePrevious=true;
            std::memcpy(inplace.bytes().data(),input.bytes().data(),input.bytes().size());
            noAlloc([&]{prepared.transpose_inplace_bytes(inplace.bytes(),work);});
            require(equal(inplace.bytes().first(k*16),expected.bytes()),"paired15 prepared in-place refresh");
            require(equal(inplace.bytes().subspan(k*16),input.bytes().subspan(k*16)),"paired15 prepared suffix");
            // PreparedEncoder is transpose-only: compare its refreshed map with
            // the immutable Code forward API under exactly the same specification.
            noAlloc([&]{fresh.forward_bytes(message.bytes(),encoded.bytes(),fw);});
            require(dot(encoded.bytes(),input.bytes())==dot(message.bytes(),actual.bytes()),"paired15 prepared/fresh forward adjoint");
            const auto descriptor=prepared.descriptor();
            require(descriptor[4]==std::byte{2} && descriptor[8]==std::byte{6} && descriptor[12]==std::byte{},"paired15 prepared descriptor");
            for(unsigned i=0;i<8;++i) {
                require(descriptor[24+i]==std::byte((seed.route>>(8*i))&255),"paired15 prepared route descriptor");
                require(descriptor[32+i]==std::byte((seed.inner>>(8*i))&255),"paired15 prepared inner descriptor");
            }
            const auto spec=prepared.specification();
            require(spec.parameters==family && spec.route_seed==seed.route && spec.inner_seed==seed.inner,"paired15 prepared specification");
        }
        // Explicit rebind has the same allocation-free contract as automatic rebind.
        prepared.setCodeSeed({91,17});
        spin::Code fresh({k,family,91,17},{spin::Backend::Portable});auto fw=fresh.make_workspace();
        fresh.transpose_bytes(input.bytes(),expected.bytes(),fw);
        noAlloc([&]{prepared.prepare_workspace(work);prepared.transpose_bytes(input.bytes(),actual.bytes(),work);});
        require(equal(actual.bytes(),expected.bytes()) && work.bytes()==scratchBytes,"paired15 explicit prepared rebind");
        spin::PreparedEncoder other({k,family,91,17},{spin::SetupMode::Full},{backend});auto wrong=other.make_workspace();
        rejects([&]{prepared.prepare_workspace(wrong);});
        rejects([&]{prepared.transpose_bytes(input.bytes(),actual.bytes(),wrong);});
        auto moved=std::move(prepared);auto movedWork=std::move(work);
        require(!prepared.supports_transpose() && prepared.descriptor()==std::array<std::byte,56>{},"paired15 moved prepared metadata");
        rejects([&]{prepared.setCodeSeed({17,43});});
        rejects([&]{prepared.transpose_bytes(input.bytes(),actual.bytes(),movedWork);});
        rejects([&]{moved.prepare_workspace(work);});
        noAlloc([&]{moved.prepare_workspace(movedWork);moved.transpose_bytes(input.bytes(),actual.bytes(),movedWork);});
        require(equal(actual.bytes(),expected.bytes()),"paired15 prepared move loses map");
    }
}
void knownAnswers() {
    // Independently frozen by paired15_research.cpp against the original mode52
    // setup and literal research oracles, not either production backend.
    // Each payload word is SplitMix64 in memory order: forward message seed1949,
    // transpose input seed913. Fold little-endian 64-bit output words below.
    struct Answer {std::uint64_t route,inner,forward,transpose;};
    constexpr Answer answers[] = {
        {1,1,0x13e0e84053fe613aULL,0x96f538d4aaa15e09ULL},
        {17,43,0x2a2b1fd8cbe63e27ULL,0x50e2f0fb80a44631ULL},
        {913,1123,0x54a88ac92075ff50ULL,0xbb2dddd2feab1d05ULL}};
    const auto hash=[](std::span<const std::byte> bytes) {
        std::uint64_t result=0;
        for(std::size_t i=0;i<bytes.size();i+=8) {
            std::uint64_t word;std::memcpy(&word,bytes.data()+i,8);
            result=(result^word)*0x100000001b3ULL;
        }
        return result;
    };
    spin::Buffer message(k*16),input(n*16),forward(n*16),transpose(k*16);
    fill(message.bytes(),1949);fill(input.bytes(),913);
    for(const auto& answer:answers)for(auto backend:{spin::Backend::Portable,spin::Backend::Automatic}) {
        spin::Code code({k,family,answer.route,answer.inner},{backend});auto work=code.make_workspace();
        context={"frozen research known answers",answer.route,answer.inner,16,0,0,code.backend()};
        noAlloc([&]{code.forward_bytes(message.bytes(),forward.bytes(),work);code.transpose_bytes(input.bytes(),transpose.bytes(),work);});
        const auto forwardHash=hash(forward.bytes()),transposeHash=hash(transpose.bytes());
        if(forwardHash!=answer.forward || transposeHash!=answer.transpose)
            std::cerr<<std::hex<<"known answers: forward actual=0x"<<forwardHash<<" expected=0x"<<answer.forward
                <<" transpose actual=0x"<<transposeHash<<" expected=0x"<<answer.transpose<<std::dec<<'\n';
        require(forwardHash==answer.forward,"paired15 forward differs from frozen research known answer");
        require(transposeHash==answer.transpose,"paired15 transpose differs from frozen research known answer");
    }
}
}
int main() {try {
    metadata();
    for(auto seeds:std::array<std::pair<std::uint64_t,std::uint64_t>,2>{{{1,1},{17,43}}}) {
        for(auto width:{spin::Width::Bits128,spin::Width::Bits256,spin::Width::Bits512})forwardCase(seeds.first,seeds.second,width);
        transposeAndAdjoint(seeds.first,seeds.second);
    }
    independentSeeds();preparedFull();knownAnswers();
    for(auto backend:{spin::Backend::Portable,spin::Backend::Automatic}) {
        spin::Code stores({k,family,17,43},{backend});
        spin::test::forwardOutputStores(stores,[](auto operation){noAlloc(operation);});
    }
    std::cout<<"paired15 API PASS: K16 profile, scalar/SIMD, 128/256/512 forward, bitwise adjoints, sparse boundaries, guards, descriptors, seeds, workspace ownership/reuse, prepared Full refresh, frozen research known answers, no allocation\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
