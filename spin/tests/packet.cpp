#include <spin/PacketCode.h>
#include <spin/PreparedEncoder.h>
#include "../src/packet/PacketPlan.h"
#include <algorithm>
#include <array>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <future>
#include <iostream>
#include <limits>
#include <new>
#include <utility>
#include <vector>
#ifdef _MSC_VER
#include <malloc.h>
#endif

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

struct alignas(16) Block {std::uint64_t words[2]{};bool operator==(const Block&) const=default;};
static void require(bool value,const char* message) {if(!value)throw std::runtime_error(message);}
template<class F> static void rejects(F f) {
    bool rejected=false;try{f();}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"invalid packet call accepted");
}
template<class F> static void noAlloc(F f) {
    const auto before=allocations.load();f();require(allocations.load()==before,"packet encoding allocated");
}
static void fill(std::span<Block> values,std::uint64_t seed=913) {
    const auto next=[&]() {seed+=0x9e3779b97f4a7c15ULL;auto z=seed;
        z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;z=(z^(z>>27))*0x94d049bb133111ebULL;return z^(z>>31);};
    for(auto& b:values){b.words[0]=next();b.words[1]=next();}
}
static void test(std::size_t k,std::uint64_t seed) {
    spin::PacketCode reference({k,seed},spin::PacketBackend::Portable), code({k,seed});
    auto work=code.make_workspace();auto rw=reference.make_workspace();
    require(reference.descriptor()==code.descriptor(),"backend affects descriptor");
    require(code.message_size()==k && code.code_size()==2*k && code.setup_bytes()>0,"packet shape");
    require(code.specification().message_size==k,"packet construction rounded K");
    require(code.backend()==(spin::packet_fast_available()?spin::PacketBackend::Avx512Gfni:
        spin::PacketBackend::Portable),"natural length disabled fast packet backend");
    std::vector<Block> input(2*k),expected(k),actual(k);fill(input);
    const auto original=input;
    noAlloc([&]{reference.transpose<Block>(input,expected,rw);code.transpose<Block>(input,actual,work);});
    require(input==original && actual==expected,"packet portable vs fast transpose");
    noAlloc([&]{code.transpose_inplace<Block>(input,work);});
    require(std::equal(expected.begin(),expected.end(),input.begin()),"packet inplace output");
    require(std::equal(input.begin()+k,input.end(),original.begin()+k),"packet suffix changed");
    // All allowed cache-line alignments, guards, separate and in-place outputs.
    std::vector<Block> guarded(2*k+12),outguard(k+12);
    auto* base=reinterpret_cast<Block*>((reinterpret_cast<std::uintptr_t>(guarded.data()+4)+63)&~std::uintptr_t(63));
    auto* outbase=reinterpret_cast<Block*>((reinterpret_cast<std::uintptr_t>(outguard.data()+4)+63)&~std::uintptr_t(63));
    const Block marker{{0xab84fe052987134fULL,0xff59c47186180239ULL}};
    for(unsigned offset=0;offset<4;++offset) {
        auto* x=base+offset;auto* y=outbase+offset;
        x[-1]=x[2*k]=y[-1]=y[k]=marker;
        std::copy(original.begin(),original.end(),x);
        noAlloc([&]{code.transpose<Block>({x,2*k},{y,k},work);});
        require(std::equal(expected.begin(),expected.end(),y),"unaligned cache-line output");
        require(x[-1]==marker && x[2*k]==marker && y[-1]==marker && y[k]==marker,"out-of-place guard");
        noAlloc([&]{code.transpose_inplace<Block>({x,2*k},work);});
        require(std::equal(expected.begin(),expected.end(),x),"unaligned cache-line inplace");
        require(std::equal(original.begin()+k,original.end(),x+k),"unaligned suffix");
        require(x[-1]==marker && x[2*k]==marker,"inplace guard");
    }
    // Exercise first/last physical steps independently of pseudorandom payloads.
    if(k<=2560)for(auto coordinate:{std::size_t{0},std::size_t{63},std::size_t{64},2*k-65,2*k-1}) {
        std::fill(input.begin(),input.end(),Block{});input[coordinate]={{1,0x843234543123abcdULL}};
        reference.transpose<Block>(input,expected,rw);code.transpose<Block>(input,actual,work);
        require(expected==actual,"sparse packet boundary");
    }
    rejects([&]{code.transpose<Block>(std::span<const Block>(original).first(2*k-1),actual,work);});
    rejects([&]{code.transpose<Block>(original,std::span<Block>(actual).first(k-1),work);});
    rejects([&]{code.transpose<Block>(input,std::span<Block>(input).first(k),work);});
    rejects([&]{code.transpose<Block>(input,std::span<Block>(input).subspan(1,k),work);});
    rejects([&]{code.transpose_bytes(std::as_bytes(std::span<const Block>(guarded)).subspan(16,2*k*16),
        std::as_writable_bytes(std::span(guarded)).first(k*16),work);});
    auto bytes=std::as_writable_bytes(std::span(guarded));
    rejects([&]{code.transpose_inplace_bytes(bytes.subspan(1,2*k*16),work);});
    spin::PacketCode other({k,seed});auto wrong=other.make_workspace();
    rejects([&]{code.transpose_inplace<Block>(input,wrong);});
    auto movedWork=std::move(work);require(work.bytes()==0,"moved workspace size");
    rejects([&]{code.transpose_inplace<Block>(input,work);});
    auto copy=code;auto moved=std::move(code);
    require(code.code_size()==0 && code.setup_bytes()==0,"moved code queries");
    rejects([&]{code.make_workspace();});rejects([&]{code.make_buffer();});
    noAlloc([&]{copy.transpose<Block>(input,actual,movedWork);});
    auto buffer=moved.make_buffer(spin::PacketMemory::PreferHugePages);
    require(buffer.bytes().size()==2*k*16 && !(reinterpret_cast<std::uintptr_t>(buffer.bytes().data())%(2*1024*1024)),"owned huge alignment");
    require(std::all_of(buffer.bytes().begin(),buffer.bytes().end(),[](auto b){return b==std::byte{};}),"owned initialization");
    std::memcpy(buffer.bytes().data(),input.data(),input.size()*16);
    auto hugeWork=moved.make_workspace(spin::PacketMemory::PreferHugePages);
    noAlloc([&]{moved.transpose_inplace_bytes(buffer.bytes(),hugeWork);});
    require(!std::memcmp(buffer.bytes().data(),actual.data(),actual.size()*16),"owned huge policy changes output");
    auto movedBuffer=std::move(buffer);require(buffer.bytes().empty() && buffer.allocation_bytes()==0,"moved buffer");
    rejects([&]{moved.make_workspace(static_cast<spin::PacketMemory>(99));});
}
static void concurrent() {
    spin::PacketCode code({1536,17});auto work=code.make_workspace();
    std::vector<Block> in(3072),expected(1536);fill(in);code.transpose<Block>(in,expected,work);
    auto run=[code,&in,&expected] {auto w=code.make_workspace();std::vector<Block> out(expected.size());
        code.transpose<Block>(in,out,w);require(out==expected,"concurrent packet encoding");};
    auto first=std::async(std::launch::async,run),second=std::async(std::launch::async,run);first.get();second.get();
}
static void knownAnswers() {
    // Frozen against research/workstreams/k16_design/implementation at
    // ad4e7711, not either library backend. RsBorderSetup/RsBorderScalar and
    // RsWideSetup/RsWideScalar, Plan(K, routeSeed, innerSeed, 20), scalar
    // transpose; the RS oracle uses literal Lagrange interpolation. Hash the
    // K output records, low then high payload word, using the fold below.
    struct Answer {std::size_t k;std::uint64_t routeSeed,innerSeed,hash;};
    constexpr Answer answers[]={{4096,1,18,0x326c6ee2a5f2493fULL},
        {4096,17,34,0xc38644b2f45f8158ULL},{12288,1,18,0x6c01adf670e2fd3bULL},
        {12288,17,34,0xe91168212580ef87ULL}};
    for(const auto& a:answers) {
        spin::Code code({a.k,spin::Parameters::PacketRsT64S20,a.routeSeed,a.innerSeed});
        auto work=code.make_workspace();
        std::vector<Block> input(code.code_size()),expected(a.k);fill(input);
        // Research block(hi,lo) draws the high half first.
        for(auto& b:input)std::swap(b.words[0],b.words[1]);
        spin::Code portable(code.specification(),{spin::Backend::Portable});
        auto pw=portable.make_workspace();portable.transpose<Block>(input,expected,pw);
        code.transpose_inplace<Block>(input,work);
        std::uint64_t hash=0;
        for(const auto& b:std::span<const Block>(input).first(a.k))
            for(auto w:b.words)hash=(hash^w)*0x100000001b3ULL;
        require(hash==a.hash,"research-anchored packet known answer");
        require(std::equal(expected.begin(),expected.end(),input.begin()),"research known answer portable backend");
        spin::PacketCode packet({a.k,a.routeSeed});
        auto descriptor=packet.descriptor();
        require(descriptor[0]==std::byte{'S'} && descriptor[1]==std::byte{'P'} &&
            descriptor[2]==std::byte{'K'} && descriptor[3]==std::byte{'P'} &&
            descriptor[4]==std::byte{2} && descriptor[8]==std::byte{2},"packet descriptor family");
        for(unsigned i=0;i<8;++i)require(descriptor[16+i]==std::byte((a.k>>(8*i))&255) &&
            descriptor[24+i]==std::byte((a.routeSeed>>(8*i))&255),"packet descriptor K/seed");
    }
}
static void coordinateBasis() {
    if(!spin::packet_fast_available())return;
    spin::PacketCode code({256,43}),reference({256,43},spin::PacketBackend::Portable);
    auto work=code.make_workspace(),rw=reference.make_workspace();
    std::vector<Block> input(512),expected(256),actual(256);
    // All input coordinates of the smallest natural geometry, through both
    // payload halves. The fallback uses literal physical maps throughout.
    for(unsigned i=0;i<512;++i) {
        if(i)input[i-1]={};input[i]={{0x948763ab1352fcedULL,0x0360d45b87abc921ULL}};
        reference.transpose<Block>(input,expected,rw);code.transpose<Block>(input,actual,work);
        require(expected==actual,"packet complete coordinate basis");
    }
}
static void binaryAdjoint() {
    namespace packet=spin::detail::packet;
    const auto dot=[](std::span<const Block> a,std::span<const Block> b) {
        require(a.size()==b.size(),"adjoint vector size");Block result{};
        for(std::size_t i=0;i<a.size();++i)for(unsigned word=0;word<2;++word)
            result.words[word]^=a[i].words[word]&b[i].words[word];
        return result;
    };
    for(auto k:{std::size_t{256},std::size_t{768}}) {
        const packet::Plan plan(k,17,43);
        spin::Code code({k,spin::Parameters::PacketRsT64S20,17,43});
        auto work=code.make_workspace();
        std::vector<Block> message(k),input(2*k),encoded(2*k),transposed(k);
        std::vector<packet::Block> oracleMessage(k),oracleEncoded(2*k);
        for(auto seed:{1ULL,2947ULL,79897ULL}) {
            fill(message,seed);fill(input,seed+17);
            std::memcpy(oracleMessage.data(),message.data(),k*sizeof(Block));
            packet::forwardScalar(oracleMessage.data(),oracleEncoded.data(),plan);
            std::memcpy(encoded.data(),oracleEncoded.data(),2*k*sizeof(Block));
            noAlloc([&]{code.transpose<Block>(input,transposed,work);});
            // Check all 128 independent payload-bit identities, rather than
            // reducing them to a single parity that could hide cancellation.
            require(dot(encoded,input)==dot(message,transposed),"literal forward/binary transpose adjoint");
        }
    }
}
static void sizing() {
    constexpr auto unit=spin::packet_message_alignment();
    static_assert(unit==256);
    // Validate the complete requested engineering range without allocating
    // giant code plans. Execution coverage is separate below.
    noAlloc([&]{
        for(std::size_t k=unit;k<=(std::size_t{1}<<26);k+=unit) {
            require(spin::valid_packet_message_size(k),"natural size query rejected");
            require(!spin::valid_packet_message_size(k-1),"unaligned size query accepted");
            require(!spin::valid_packet_message_size(k+1),"unaligned size query accepted");
            for(auto value:{k-1,k,k+1})require(
                spin::valid_message_size(spin::Parameters::PacketRsT64S20,value)==spin::valid_packet_message_size(value),
                "unified natural size query disagrees with PacketCode");
        }
    });
    // Padded route offsets use 32 bits. Exercise the bound without allocating.
    constexpr auto maxK=(std::size_t{std::numeric_limits<std::uint32_t>::max()}/516)*unit;
    require(spin::valid_packet_message_size(maxK),"last representable size rejected");
    require(!spin::valid_packet_message_size(maxK+unit),"padded route overflow accepted");
    require(spin::valid_message_size(spin::Parameters::PacketRsT64S20,maxK) &&
        !spin::valid_message_size(spin::Parameters::PacketRsT64S20,maxK+unit),"unified padded route limit");
    rejects([&]{spin::PacketCode code({maxK+unit,1});});
    // Around both the inner/region alignment boundary (16 groups) and the
    // measured powers of two. These execute the same selected backend.
    for(auto k:{15*unit,16*unit,17*unit,(std::size_t{1}<<16)-unit,
        (std::size_t{1}<<16)+unit,(std::size_t{1}<<18)-unit,(std::size_t{1}<<18)+unit})
        test(k,43);
}

static void largePacketCase(std::size_t k) {
    constexpr auto family=spin::Parameters::PacketRsT64S20;
    constexpr std::uint64_t seed=43;
    const spin::CodeSpec spec{k,family,seed,seed};
    spin::Code reference(spec,{spin::Backend::Portable}),code(spec,{spin::Backend::Avx512});
    auto rw=reference.make_workspace(),work=code.make_workspace();
    require(code.backend()==spin::Backend::Avx512 && reference.descriptor()==code.descriptor(),
        "large packet backend or descriptor");
    const auto n=2*k;
    std::vector<Block> original(n),expected(k);fill(original,seed);
    // Force the legal16-byte but non-cache-line-aligned API boundary. Routing
    // scratch remains owned and64-byte aligned, as required by streaming stores.
    spin::Buffer input((n+2)*sizeof(Block)),output((k+2)*sizeof(Block));
    auto* x=reinterpret_cast<Block*>(input.bytes().data())+1;
    auto* y=reinterpret_cast<Block*>(output.bytes().data())+1;
    require(reinterpret_cast<std::uintptr_t>(x)%64==16 &&
        reinterpret_cast<std::uintptr_t>(y)%64==16,"large packet alignment fixture");
    const Block marker{{0x764ac51efe0328b9ULL,0x84c0d973123afed1ULL}};
    x[-1]=x[n]=y[-1]=y[k]=marker;
    std::copy(original.begin(),original.end(),x);
    noAlloc([&]{reference.transpose<Block>(original,expected,rw);});
    noAlloc([&]{code.transpose<Block>({x,n},{y,k},work);});
    require(std::equal(expected.begin(),expected.end(),y),"routing boundary scalar vs fast");
    require(std::equal(original.begin(),original.end(),x),"routing boundary changed separate input");
    require(x[-1]==marker && x[n]==marker && y[-1]==marker && y[k]==marker,
        "routing boundary separate guard");
    noAlloc([&]{code.transpose_inplace_bytes(std::as_writable_bytes(std::span<Block>(x,n)),work);});
    require(std::equal(expected.begin(),expected.end(),x),"routing boundary inplace result");
    require(std::equal(original.begin()+k,original.end(),x+k),"routing boundary inplace suffix");
    require(x[-1]==marker && x[n]==marker,"routing boundary inplace guard");
    // The compatibility facade must reach the same large-size map.
    if(k==(std::size_t{1}<<20)) {
        spin::PacketCode legacy({k,seed},spin::PacketBackend::Avx512Gfni);
        auto legacyWork=legacy.make_workspace();
        std::copy(original.begin(),original.end(),x);
        noAlloc([&]{legacy.transpose_inplace<Block>({x,n},legacyWork);});
        require(std::equal(expected.begin(),expected.end(),x),"K20 PacketCode differs from Code");
        require(std::equal(original.begin()+k,original.end(),x+k),"K20 PacketCode suffix changed");
        require(x[-1]==marker && x[n]==marker,"K20 PacketCode guard");
    }
}

static void largePacket() {
    // Exercise the certified geometry without timing assertions. Natural
    // non-power-of-two lengths and physical-step tails are covered by sizing().
    if(spin::packet_fast_available())largePacketCase(std::size_t{1}<<20);
}

static void retiredFamily() {
    constexpr auto retired=static_cast<spin::Parameters>(4);
    require(!spin::valid_message_size(retired,512),"retired packet family accepted by size query");
    rejects([]{spin::message_alignment(retired);});
    rejects([]{spin::Code code({512,retired,17,43});});
    rejects([]{spin::PreparedEncoder code({512,retired,17,43});});
}

static void commonBuffers() {
    noAlloc([]{
        spin::Buffer empty(0), hugeEmpty(0,spin::MemoryPolicy::PreferHugePages);
        require(empty.bytes().empty() && !empty.bytes().data() && empty.allocation_bytes()==0,"empty common buffer");
        require(hugeEmpty.bytes().empty() && !hugeEmpty.bytes().data() && hugeEmpty.allocation_bytes()==0,"empty huge common buffer");
    });
    spin::Buffer normal(65);
    require(normal.bytes().size()==65 && normal.allocation_bytes()==128,"common buffer rounded size");
    require(!(reinterpret_cast<std::uintptr_t>(normal.bytes().data())%64),"common buffer alignment");
    require(std::all_of(normal.bytes().begin(),normal.bytes().end(),[](auto b){return b==std::byte{};}),"common buffer initialized");
    normal.bytes()[64]=std::byte{91};const auto* address=normal.bytes().data();
    auto moved=std::move(normal);
    require(normal.bytes().empty() && normal.allocation_bytes()==0,"common moved-from buffer");
    require(moved.bytes().data()==address && std::as_const(moved).bytes()[64]==std::byte{91},"common buffer move contents");
    spin::Buffer assigned(1);assigned=std::move(moved);
    require(moved.bytes().empty() && assigned.bytes().data()==address,"common buffer move assignment");
    rejects([]{spin::Buffer invalid(0,static_cast<spin::MemoryPolicy>(99));});
    for(auto policy:{spin::MemoryPolicy::Normal,spin::MemoryPolicy::PreferHugePages,spin::MemoryPolicy::Automatic}) {
        noAlloc([&]{
            bool failed=false;
            try{spin::Buffer overflow(std::numeric_limits<std::size_t>::max(),policy);}
            catch(const std::bad_alloc&){failed=true;}
            require(failed,"common buffer overflow accepted");
        });
    }
}

static constexpr std::size_t memoryAlignment(std::size_t bytes,spin::MemoryPolicy policy) {
    if(policy==spin::MemoryPolicy::PreferHugePages)return 2*1024*1024;
#if defined(__linux__)
    if(policy==spin::MemoryPolicy::Automatic && bytes>=16*1024*1024)return 2*1024*1024;
#endif
    return 64;
}
static constexpr std::size_t memoryBytes(std::size_t bytes,spin::MemoryPolicy policy) {
    const auto a=memoryAlignment(bytes,policy);
    return (bytes+a-1)&~(a-1);
}
static void automaticMemory() {
    using spin::MemoryPolicy;
    constexpr std::size_t threshold=16*1024*1024;
    static_assert(int(MemoryPolicy::Normal)==0 && int(MemoryPolicy::PreferHugePages)==1);
    for(auto bytes:{std::size_t{65},threshold-1,threshold,threshold+1}) {
        // The omitted argument must select Automatic, including at the exact
        // logical-byte boundary; explicit policies override it at every size.
        spin::Buffer automatic(bytes);
        require(automatic.allocation_bytes()==memoryBytes(bytes,MemoryPolicy::Automatic),"automatic buffer rounding");
        require(reinterpret_cast<std::uintptr_t>(automatic.bytes().data())%
            memoryAlignment(bytes,MemoryPolicy::Automatic)==0,"automatic buffer alignment");
        require(automatic.bytes().size()==bytes && automatic.bytes().front()==std::byte{} &&
            automatic.bytes().back()==std::byte{},"automatic buffer logical size/initialization");
        for(auto policy:{MemoryPolicy::Normal,MemoryPolicy::PreferHugePages}) {
            spin::Buffer explicitBuffer(bytes,policy);
            require(explicitBuffer.allocation_bytes()==memoryBytes(bytes,policy),"explicit memory policy rounding");
            require(reinterpret_cast<std::uintptr_t>(explicitBuffer.bytes().data())%memoryAlignment(bytes,policy)==0,
                "explicit memory policy alignment");
        }
    }
    constexpr auto family=spin::Parameters::PacketRsT64S20;
    // Just above K19 distinguishes the factory defaults by allocation size:
    // Normal rounds by 64 bytes, Automatic/Linux rounds by 2 MiB.
    constexpr std::size_t smallK=256,bigK=(std::size_t{1}<<19)+256;
    const auto scratchBytes=[](std::size_t k){return (k/256)*516*16;};
    // Padded scratch crosses first: the input still fits below 16 MiB at
    // K=2033*256. Resolve both from their own logical bytes, not from K alone.
    for(auto k:{std::size_t{2032*256},std::size_t{2033*256}}) {
        spin::Code boundary({k,family,17,43});
        auto bw=boundary.make_workspace();auto bb=boundary.make_buffer();
        require(bw.bytes()==memoryBytes(scratchBytes(k),MemoryPolicy::Automatic),"scratch-specific automatic boundary");
        require(bb.allocation_bytes()==memoryBytes(32*k,MemoryPolicy::Automatic),"input-specific automatic boundary");
    }
    spin::Code small({smallK,family,17,43}),big({bigK,family,17,43});
    auto work=small.make_workspace();
    require(work.bytes()==memoryBytes(scratchBytes(smallK),MemoryPolicy::Automatic),"small automatic workspace");
    big.prepare_workspace(work);
    require(work.bytes()==memoryBytes(scratchBytes(bigK),MemoryPolicy::Automatic),"automatic workspace resize up");
    auto normal=big.make_workspace(spin::Width::Bits128,MemoryPolicy::Normal);
    std::vector<Block> input(2*bigK),expected(bigK),actual(bigK);fill(input);
    noAlloc([&]{big.transpose<Block>(input,expected,normal);big.transpose<Block>(input,actual,work);});
    require(expected==actual,"automatic policy changed code map");
    noAlloc([&]{big.prepare_workspace(work);});
    small.prepare_workspace(work);
    require(work.bytes()==memoryBytes(scratchBytes(smallK),MemoryPolicy::Automatic),"automatic workspace resize down");
    small.prepare_workspace(normal);big.prepare_workspace(normal);
    require(normal.bytes()==memoryBytes(scratchBytes(bigK),MemoryPolicy::Normal),"resize lost explicit Normal");
    auto buffer=big.make_buffer();
    require(buffer.allocation_bytes()==memoryBytes(32*bigK,MemoryPolicy::Automatic),"Code buffer default policy");
    spin::PacketCode packet({bigK,17});auto pw=packet.make_workspace();auto pb=packet.make_buffer();
    require(pw.bytes()==memoryBytes(scratchBytes(bigK),MemoryPolicy::Automatic) &&
        pb.allocation_bytes()==memoryBytes(32*bigK,MemoryPolicy::Automatic),"PacketCode default policy");
    spin::PreparedEncoder prepared({bigK,family,17,43});auto fw=prepared.make_workspace();auto fb=prepared.make_buffer();
    require(fw.bytes()==memoryBytes(scratchBytes(bigK),MemoryPolicy::Automatic) &&
        fb.allocation_bytes()==memoryBytes(32*bigK,MemoryPolicy::Automatic),"PreparedEncoder default policy");
    rejects([&]{big.make_workspace(spin::Width::Bits128,static_cast<MemoryPolicy>(99));});
    rejects([&]{prepared.make_workspace(spin::Width::Bits128,static_cast<MemoryPolicy>(99));});
}

static void unifiedCode() {
    constexpr std::size_t k=1536;
    constexpr auto family=spin::Parameters::PacketRsT64S20;
    require(spin::message_alignment(family)==256 && spin::valid_message_size(family,k),"unified packet geometry");
    require(!spin::valid_message_size(family,k+1),"unified packet geometry accepted misalignment");
    for(auto seeds:{spin::CodeSeed{17,17},spin::CodeSeed{43,17},spin::CodeSeed{43,91}}) {
        const spin::CodeSpec spec{k,family,seeds.route,seeds.inner};
        spin::Code code(spec),reference(spec,{spin::Backend::Portable});
        auto work=code.make_workspace(),rw=reference.make_workspace();
        require(code.message_size()==k && code.code_size()==2*k && code.setup_bytes()>0,"unified packet shape");
        require(work.width()==spin::Width::Bits128,"unified packet workspace width");
        require(code.backend()==(spin::packet_fast_available()?spin::Backend::Avx512:spin::Backend::Portable),"unified packet backend");
        require(reference.backend()==spin::Backend::Portable && reference.descriptor()==code.descriptor(),"unified portable descriptor");
        require(code.supports_transpose() && code.supports_forward() && !code.supports_generic_transpose(),"unified packet capabilities");
        require(!code.supports_transpose(spin::Width::Bits256) && !code.supports_forward(spin::Width::Bits512),"unified packet wide capability");
        std::vector<Block> input(2*k),expected(k),actual(k);fill(input);
        const auto original=input;
        noAlloc([&]{reference.transpose<Block>(input,expected,rw);code.transpose<Block>(input,actual,work);});
        require(input==original && actual==expected,"unified independent seeds portable vs fast");
        noAlloc([&]{code.transpose_inplace<Block>(input,work);});
        require(std::equal(expected.begin(),expected.end(),input.begin()),"unified inplace result");
        require(std::equal(original.begin()+k,original.end(),input.begin()+k),"unified suffix preserved");
        if(seeds.route==seeds.inner) {
            spin::PacketCode legacy({k,seeds.route});auto lw=legacy.make_workspace();
            noAlloc([&]{legacy.transpose<Block>(original,actual,lw);});
            require(actual==expected,"unified equal seeds differ from PacketCode");
            const auto d=legacy.descriptor();
            std::array<std::byte,32> frozen{};
            frozen[0]=std::byte{'S'};frozen[1]=std::byte{'P'};frozen[2]=std::byte{'K'};frozen[3]=std::byte{'P'};
            frozen[4]=std::byte{2};frozen[8]=std::byte{2};
            for(unsigned i=0;i<8;++i) {frozen[16+i]=std::byte((k>>(8*i))&255);frozen[24+i]=std::byte((seeds.route>>(8*i))&255);}
            require(d==frozen,"PacketCode descriptor family/version");
        }
        std::array<std::byte,40> descriptor{};
        descriptor[0]=std::byte{'S'};descriptor[1]=std::byte{'P'};descriptor[2]=std::byte{'I'};descriptor[3]=std::byte{'N'};
        descriptor[4]=std::byte{1};descriptor[8]=std::byte{5};
        for(unsigned i=0;i<8;++i) {
            descriptor[16+i]=std::byte((k>>(8*i))&255);
            descriptor[24+i]=std::byte((seeds.route>>(8*i))&255);
            descriptor[32+i]=std::byte((seeds.inner>>(8*i))&255);
        }
        require(code.descriptor()==descriptor,"unified packet descriptor or seed fields");
        std::vector<Block> forwardExpected(2*k);
        noAlloc([&]{reference.forward<Block>(expected,forwardExpected,rw);code.forward<Block>(expected,input,work);});
        require(input==forwardExpected,"unified forward result");
        rejects([&]{code.generic_transpose();});
        std::vector<std::uint64_t> bitInput(k/64),bitOutput(2*k/64),bitScratch(2*k/64);
        rejects([&]{code.forward_bits(bitInput,bitOutput,bitScratch);});
        for(auto width:{spin::Width::Bits256,spin::Width::Bits512,static_cast<spin::Width>(17)}) {
            rejects([&]{code.make_workspace(width);});rejects([&]{code.make_buffer(width);});
        }
        spin::Code other(spec);auto wrong=other.make_workspace();
        rejects([&]{code.transpose<Block>(original,actual,wrong);});
        const auto oldBytes=wrong.bytes();
        noAlloc([&]{code.prepare_workspace(wrong);code.transpose<Block>(original,actual,wrong);});
        require(wrong.bytes()==oldBytes && actual==expected,"same K workspace rebind");
        rejects([&]{other.transpose<Block>(original,actual,wrong);});
        spin::Code differentSeed({k,family,seeds.route+1,seeds.inner+1});
        auto dw=differentSeed.make_workspace();differentSeed.transpose<Block>(original,expected,dw);
        require(expected!=actual,"independent seed changes did not change map");
        noAlloc([&]{differentSeed.prepare_workspace(wrong);differentSeed.transpose<Block>(original,actual,wrong);});
        require(actual==expected,"rebound workspace retained old packet setup");
        auto owned=code.make_buffer();auto ownedWork=code.make_workspace();
        require(owned.bytes().size()==code.code_size()*16,"unified owned buffer size");
        std::memcpy(owned.bytes().data(),original.data(),original.size()*16);
        noAlloc([&]{code.transpose_inplace_bytes(owned.bytes(),ownedWork);});
        reference.transpose<Block>(original,expected,rw);
        require(!std::memcmp(owned.bytes().data(),expected.data(),expected.size()*16),"unified owned buffer result");
    }
    rejects([&]{spin::Code c({k,family},{spin::Backend::Avx2});});
    rejects([&]{spin::Code c({k,family},{static_cast<spin::Backend>(99)});});
    rejects([&]{spin::Code c({k,family},{spin::Backend::Portable,4});});
    spin::Code fixedTile({k,family},{spin::Backend::Portable,8});
    require(fixedTile.backend()==spin::Backend::Portable,"explicit fixed packet tile");
    if(!spin::packet_fast_available()) {
        bool failed=false;
        try{spin::Code unavailable({k,family},{spin::Backend::Avx512});}catch(const std::runtime_error&){failed=true;}
        require(failed,"unified forced unavailable packet ISA");
    }
}

static void preparedPacket() {
    constexpr std::size_t k=1536;
    constexpr auto family=spin::Parameters::PacketRsT64S20;
    spin::PreparedEncoder prepared({k,family,17,17});
    auto work=prepared.make_workspace(spin::Width::Bits128,spin::MemoryPolicy::Normal);
    require(prepared.supports_transpose() && !prepared.supports_transpose(spin::Width::Bits256) &&
        !prepared.supports_generic_transpose(),"prepared packet capabilities");
    require(prepared.backend()==(spin::packet_fast_available()?spin::Backend::Avx512:spin::Backend::Portable),"prepared packet backend");
    rejects([&]{prepared.generic_transpose();});
    rejects([&]{prepared.make_workspace(spin::Width::Bits256);});
    rejects([&]{spin::PreparedEncoder invalid({k,family,17,17},{spin::SetupMode::BankedHeuristic});});
    std::vector<Block> input(2*k),actual(k),expected(k);fill(input);
    const auto original=input;const auto scratchBytes=work.bytes();
    std::vector<Block> previous;
    for(auto seeds:{spin::CodeSeed{17,17},spin::CodeSeed{43,17},spin::CodeSeed{43,91}}) {
        prepared.setCodeSeed(seeds);
        noAlloc([&]{prepared.setCodeSeed(seeds);});
        spin::Code fresh({k,family,seeds.route,seeds.inner},{spin::Backend::Portable});auto fw=fresh.make_workspace();
        fresh.transpose<Block>(original,expected,fw);
        noAlloc([&]{prepared.transpose<Block>(original,actual,work);});
        require(actual==expected && work.bytes()==scratchBytes,"prepared packet refresh does not match fresh Code");
        if(!previous.empty())require(previous!=actual,"route-only or inner-only refresh did not change output");
        previous=actual;input=original;
        noAlloc([&]{prepared.transpose_inplace<Block>(input,work);});
        require(std::equal(expected.begin(),expected.end(),input.begin()),"prepared packet inplace refresh");
        require(std::equal(original.begin()+k,original.end(),input.begin()+k),"prepared packet suffix changed");
        const auto spec=prepared.specification();
        require(spec.route_seed==seeds.route && spec.inner_seed==seeds.inner,"prepared packet seed fields");
        const auto d=prepared.descriptor();
        require(d[4]==std::byte{2} && d[8]==std::byte{5} && d[12]==std::byte{},"prepared packet descriptor family/mode");
        for(unsigned i=0;i<8;++i)require(d[24+i]==std::byte((seeds.route>>(8*i))&255) &&
            d[32+i]==std::byte((seeds.inner>>(8*i))&255),"prepared packet descriptor seed refresh");
    }
    spin::PreparedEncoder other({k,family,43,91});auto wrong=other.make_workspace();
    rejects([&]{prepared.transpose<Block>(original,actual,wrong);});
    auto huge=prepared.make_workspace(spin::Width::Bits128,spin::MemoryPolicy::PreferHugePages);
    noAlloc([&]{prepared.transpose<Block>(original,actual,huge);});
    require(actual==expected,"prepared packet huge workspace changed output");
    spin::PreparedEncoder portable({k,family,43,91},{spin::SetupMode::Full},{spin::Backend::Portable});
    auto portableWork=portable.make_workspace();
    require(portable.backend()==spin::Backend::Portable && portable.supports_transpose(),"prepared portable packet controller");
    noAlloc([&]{portable.transpose<Block>(original,actual,portableWork);});
    require(actual==expected,"prepared portable packet result");
    portable.setCodeSeed({91,43});
    spin::Code refreshed({k,family,91,43},{spin::Backend::Portable});auto freshWork=refreshed.make_workspace();
    refreshed.transpose<Block>(original,expected,freshWork);
    noAlloc([&]{portable.transpose<Block>(original,actual,portableWork);});
    require(actual==expected,"prepared portable packet refresh");

    // Callers can rebind after refresh before entering their timed operation.
    // The earlier loop separately tests the automatic first-encode rebind.
    prepared.setCodeSeed({91,17});
    spin::Code explicitlyPrepared({k,family,91,17},{spin::Backend::Portable});
    auto explicitReference=explicitlyPrepared.make_workspace();
    explicitlyPrepared.transpose<Block>(original,expected,explicitReference);
    noAlloc([&]{prepared.prepare_workspace(work);});
    noAlloc([&]{prepared.transpose<Block>(original,actual,work);});
    require(actual==expected && work.bytes()==scratchBytes,"explicit prepared workspace refresh");
    rejects([&]{prepared.prepare_workspace(wrong);});
    auto retainedWork=std::move(work);
    require(work.bytes()==0,"moved prepared workspace query");
    rejects([&]{prepared.prepare_workspace(work);});
    noAlloc([&]{prepared.prepare_workspace(retainedWork);prepared.transpose<Block>(original,actual,retainedWork);});
    require(actual==expected,"explicit prepared workspace move");
}

static void workspaceConversions() {
    constexpr auto family=spin::Parameters::PacketRsT64S20;
    constexpr std::size_t smallK=1536;
    spin::Code first({smallK,family,1,17}),second({smallK,family,43,91});
    auto huge=first.make_workspace(spin::Width::Bits128,spin::MemoryPolicy::PreferHugePages);
    const auto bytes=huge.bytes();
    require(bytes && bytes%(2*1024*1024)==0,"huge common workspace allocation");
    std::vector<Block> in(2*smallK),expected(smallK),actual(smallK);fill(in);
    auto fresh=second.make_workspace();second.transpose<Block>(in,expected,fresh);
    rejects([&]{second.transpose<Block>(in,actual,huge);});
    noAlloc([&]{second.prepare_workspace(huge);second.transpose<Block>(in,actual,huge);});
    require(huge.bytes()==bytes && actual==expected,"huge workspace rebind changed allocation or output");
    rejects([&]{first.transpose<Block>(in,actual,huge);});
    noAlloc([&]{second.prepare_workspace(huge);});

    if(!spin::capabilities().avx2)return;
    constexpr std::size_t k=16384;
    spin::Code old({k,spin::Parameters::T128S19,17,43}),packet({k,family,17,43});
    auto shared=old.make_workspace(),oldFresh=old.make_workspace(),packetFresh=packet.make_workspace();
    std::vector<Block> input(2*k),oldExpected(k),packetExpected(k),output(k);fill(input);
    old.transpose<Block>(input,oldExpected,oldFresh);packet.transpose<Block>(input,packetExpected,packetFresh);
    rejects([&]{packet.transpose<Block>(input,output,shared);});
    packet.prepare_workspace(shared);
    noAlloc([&]{packet.transpose<Block>(input,output,shared);});
    require(output==packetExpected && shared.width()==spin::Width::Bits128,"old-to-packet workspace conversion");
    rejects([&]{old.transpose<Block>(input,output,shared);});
    old.prepare_workspace(shared);
    noAlloc([&]{old.transpose<Block>(input,output,shared);});
    require(output==oldExpected && shared.width()==spin::Width::Bits128,"packet-to-old workspace conversion");
    rejects([&]{packet.transpose<Block>(input,output,shared);});
}

static void movedUnifiedHandles() {
    constexpr std::size_t k=512;
    constexpr auto family=spin::Parameters::PacketRsT64S20;
    std::vector<Block> input(2*k),expected(k),actual(k);fill(input);
    spin::Code source({k,family,17,43},{spin::Backend::Portable});auto work=source.make_workspace();
    source.transpose<Block>(input,expected,work);
    auto moved=std::move(source);
    require(source.message_size()==0 && source.code_size()==0 && source.setup_bytes()==0 &&
        source.specification().message_size==0 && source.backend()==spin::Backend::Automatic,"moved Code getters");
    require(!source.supports_forward() && !source.supports_transpose() && !source.supports_generic_transpose() &&
        source.descriptor()==std::array<std::byte,40>{},"moved Code capabilities/descriptor");
    rejects([&]{source.make_workspace();});rejects([&]{source.make_buffer();});
    rejects([&]{source.prepare_workspace(work);});rejects([&]{source.generic_transpose();});
    rejects([&]{source.transpose_bytes({},{},work);});rejects([&]{source.transpose_inplace_bytes({},work);});
    rejects([&]{source.forward_bytes({},{},work);});rejects([&]{source.forward_bits({},{},{});});
    noAlloc([&]{moved.transpose<Block>(input,actual,work);});
    require(actual==expected,"Code move invalidated retained workspace");

    spin::PreparedEncoder prepared({k,family,17,43},{spin::SetupMode::Full},{spin::Backend::Portable});
    auto pw=prepared.make_workspace();auto movedPrepared=std::move(prepared);
    require(prepared.message_size()==0 && prepared.code_size()==0 && prepared.setup_bytes()==0 &&
        prepared.specification().message_size==0 && prepared.backend()==spin::Backend::Automatic,"moved Prepared getters");
    require(prepared.setup_options().mode==spin::SetupMode::Full && !prepared.supports_transpose() &&
        !prepared.supports_generic_transpose() && prepared.descriptor()==std::array<std::byte,56>{},"moved Prepared capabilities/descriptor");
    rejects([&]{prepared.setCodeSeed({43,91});});rejects([&]{prepared.make_workspace();});
    rejects([&]{prepared.make_buffer();});rejects([&]{prepared.generic_transpose();});
    rejects([&]{prepared.prepare_workspace(pw);});
    rejects([&]{prepared.transpose_bytes({},{},pw);});rejects([&]{prepared.transpose_inplace_bytes({},pw);});
    noAlloc([&]{movedPrepared.transpose<Block>(input,actual,pw);});
    require(actual==expected,"Prepared move invalidated retained workspace");
}


static void forwardApi() {
    constexpr auto family=spin::Parameters::PacketRsT64S20;
    for(auto k:{256U,768U,65536U,262144U,262400U,1048576U}) {
        const std::size_t n=2*std::size_t(k);
        spin::Code code({k,family,17,43}),reference({k,family,17,43},{spin::Backend::Portable});
        auto work=code.make_workspace(),rw=reference.make_workspace();
        spin::Buffer input(k*16+64),output(n*16+128),expected(n*16);
        auto in=input.bytes().subspan(16,k*16); // Foreign minimum alignment.
        auto typed=std::span<Block>(reinterpret_cast<Block*>(in.data()),k);
        auto inputBefore=std::vector<std::byte>(in.size());
        for(unsigned mode=0;mode<3;++mode) {
            fill(typed,mode+77);
            if(mode)std::fill(in.begin(),in.end(),std::byte{});
            if(mode==2) {in.front()=std::byte{1};in.back()=std::byte{0x80};}
            std::copy(in.begin(),in.end(),inputBefore.begin());
            noAlloc([&]{reference.forward_bytes(in,expected.bytes(),rw);});
            for(unsigned offset=0;offset<4;++offset) {
                std::fill(output.bytes().begin(),output.bytes().end(),std::byte{0xa5});
                auto out=output.bytes().subspan(64+16*offset,n*16);
                noAlloc([&]{code.forward_bytes(in,out,work);});
                require(!std::memcmp(out.data(),expected.bytes().data(),out.size()),"forward API alignment/reference");
                require(std::equal(in.begin(),in.end(),inputBefore.begin()),"forward modified input");
                require(std::all_of(output.bytes().begin(),output.bytes().begin()+64+16*offset,[](auto x){return x==std::byte{0xa5};}),"forward prefix guard");
                require(std::all_of(out.end(),output.bytes().end(),[](auto x){return x==std::byte{0xa5};}),"forward suffix guard");
            }
        }
        // Switch directions and rebind dirty workspace to different seeds.
        noAlloc([&]{code.transpose_bytes(expected.bytes(),in,work);code.forward_bytes(in,output.bytes().first(n*16),work);});
        spin::Code different({k,family,43,17});auto dw=different.make_workspace();
        different.forward_bytes(in,expected.bytes(),dw);
        noAlloc([&]{different.prepare_workspace(work);different.forward_bytes(in,output.bytes().first(n*16),work);});
        require(!std::memcmp(output.bytes().data(),expected.bytes().data(),n*16),"forward workspace rebind");
        rejects([&]{code.forward_bytes(in,output.bytes().first(n*16),work);});
        code.prepare_workspace(work);
        rejects([&]{code.forward_bytes(in.first(in.size()-1),output.bytes().first(n*16),work);});
        rejects([&]{code.forward_bytes(in,output.bytes().first(n*16-1),work);});
        rejects([&]{code.forward_bytes(in,output.bytes().subspan(1,n*16),work);});
        rejects([&]{code.forward_bytes(input.bytes().subspan(1,k*16),output.bytes().first(n*16),work);});
        rejects([&]{code.forward_bytes(output.bytes().first(k*16),output.bytes().first(n*16),work);});
        rejects([&]{code.forward_bytes(output.bytes().subspan(32,k*16),output.bytes().first(n*16),work);});
        rejects([&]{code.forward_bytes(output.bytes().first(k*16),output.bytes().subspan(16,n*16),work);});
        auto saved=code, moved=std::move(code);auto movedWork=std::move(work);
        require(!code.supports_forward(),"moved forward capability");
        rejects([&]{code.forward_bytes(in,output.bytes().first(n*16),movedWork);});
        rejects([&]{saved.forward_bytes(in,output.bytes().first(n*16),work);});
        noAlloc([&]{saved.forward_bytes(in,output.bytes().first(n*16),movedWork);});
        auto huge=moved.make_buffer(spin::Width::Bits128,spin::MemoryPolicy::PreferHugePages);
        auto hw=moved.make_workspace(spin::Width::Bits128,spin::MemoryPolicy::PreferHugePages);
        noAlloc([&]{moved.forward_bytes(in,huge.bytes(),hw);});
        require(!std::memcmp(huge.bytes().data(),output.bytes().data(),n*16),"forward memory policy");
    }
    for(auto backend:{spin::PacketBackend::Portable,spin::PacketBackend::Automatic}) {
        constexpr std::size_t k=768;
        spin::PacketCode code({k,17},backend);
        spin::Code common({k,family,17,17});auto cw=common.make_workspace();
        auto w=code.make_workspace();std::vector<Block> in(k),out(2*k),expected(2*k);fill(in);
        require(code.supports_forward() && !code.supports_forward(spin::Width::Bits256),"PacketCode forward capability");
        common.forward<Block>(in,expected,cw);
        noAlloc([&]{code.forward<Block>(in,out,w);});require(out==expected,"PacketCode forward map");
        rejects([&]{code.forward<Block>(std::span<const Block>(in).first(k-1),out,w);});
        rejects([&]{code.forward_bytes(std::as_bytes(std::span(in)),std::as_writable_bytes(std::span(out)).first(16),w);});
        rejects([&]{code.forward<Block>(std::span<const Block>(out).first(k),out,w);});
        spin::PacketCode other({k,17});auto wrong=other.make_workspace();
        rejects([&]{code.forward<Block>(in,out,wrong);});
        auto movedWork=std::move(w);auto saved=code,moved=std::move(code);
        require(!code.supports_forward(),"moved PacketCode forward");
        rejects([&]{code.forward<Block>(in,out,movedWork);});
        rejects([&]{saved.forward<Block>(in,out,w);});
        noAlloc([&]{saved.forward<Block>(in,out,movedWork);});require(out==expected,"copied PacketCode forward");
        // Shared immutable setup, one workspace per concurrent call.
        auto run=[saved,&in,&expected] {auto w=saved.make_workspace();std::vector<Block> out(2*k);
            saved.forward<Block>(in,out,w);require(out==expected,"concurrent forward");};
        auto a=std::async(std::launch::async,run),b=std::async(std::launch::async,run);a.get();b.get();
    }
}

int main() {try {
    for(auto k:{std::size_t{0},std::size_t{1},std::size_t{255},std::size_t{257},std::size_t{511},std::size_t{513},std::numeric_limits<std::size_t>::max()}) {
        require(!spin::valid_packet_message_size(k),"invalid K query accepted");rejects([&]{spin::PacketCode c({k,1});});
    }
    rejects([]{spin::PacketCode c({512,1},static_cast<spin::PacketBackend>(99));});
    if(!spin::packet_fast_available()) {
        bool failed=false;try{spin::PacketCode c({512,1},spin::PacketBackend::Avx512Gfni);}catch(const std::runtime_error&){failed=true;}
        require(failed,"forced unavailable packet ISA");
    }
    for(auto k:{256U,512U,768U,1024U,1536U,2560U,4096U,12288U,65536U,262144U})for(auto seed:{1ULL,17ULL})test(k,seed);
    sizing();largePacket();retiredFamily();knownAnswers();coordinateBasis();binaryAdjoint();concurrent();commonBuffers();unifiedCode();preparedPacket();
    workspaceConversions();movedUnifiedHandles();automaticMemory();forwardApi();
    std::cout<<"packet API PASS: scalar/fast, research known answers, coordinate basis, natural sizes and K20, retired family rejection, unified/prepared APIs, buffers, boundaries, alignment, ownership, no allocation\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
