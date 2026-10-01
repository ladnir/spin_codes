#include <spin/PacketCode.h>
#include <spin/PreparedEncoder.h>
#include <algorithm>
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
    // Frozen against the retained research encoder, not generated from either
    // library backend: K, seed, hash of all 2K records after one inplace call.
    struct Answer {std::size_t k;std::uint64_t seed,hash;};
    constexpr Answer answers[]={{512,1,0xa011c1a36862a018ULL},{512,17,0x6d4c655324880a31ULL},
        {1536,1,0x1af105fb37d235efULL},{1536,17,0x2032ac91043db7b6ULL},
        {65536,1,0xa611428b3f82c409ULL},{262144,17,0x7867d2fb90feb72cULL}};
    for(const auto& a:answers) {
        spin::PacketCode code({a.k,a.seed});auto work=code.make_workspace();
        std::vector<Block> input(code.code_size());fill(input);
        // Research block(hi,lo) draws the high half first.
        for(auto& b:input)std::swap(b.words[0],b.words[1]);
        code.transpose_inplace<Block>(input,work);
        std::uint64_t hash=0;
        for(const auto& b:input)for(auto w:b.words)hash=(hash^w)*0x100000001b3ULL;
        require(hash==a.hash,"research-anchored packet known answer");
        auto descriptor=code.descriptor();
        require(descriptor[0]==std::byte{'S'} && descriptor[1]==std::byte{'P'} &&
            descriptor[2]==std::byte{'K'} && descriptor[3]==std::byte{'P'} &&
            descriptor[4]==std::byte{1} && descriptor[8]==std::byte{1},"packet descriptor family");
        for(unsigned i=0;i<8;++i)require(descriptor[16+i]==std::byte((a.k>>(8*i))&255) &&
            descriptor[24+i]==std::byte((a.seed>>(8*i))&255),"packet descriptor K/seed");
    }
}
static void coordinateBasis() {
    if(!spin::packet_fast_available())return;
    spin::PacketCode code({512,43}),reference({512,43},spin::PacketBackend::Portable);
    auto work=code.make_workspace(),rw=reference.make_workspace();
    std::vector<Block> input(1024),expected(512),actual(512);
    // All input coordinates of the smallest natural geometry, through both
    // payload halves. The fallback uses literal physical maps throughout.
    for(unsigned i=0;i<1024;++i) {
        if(i)input[i-1]={};input[i]={{0x948763ab1352fcedULL,0x0360d45b87abc921ULL}};
        reference.transpose<Block>(input,expected,rw);code.transpose<Block>(input,actual,work);
        require(expected==actual,"packet complete coordinate basis");
    }
}
static void sizing() {
    constexpr auto unit=spin::packet_message_alignment();
    static_assert(unit==4*128);
    // Validate the complete requested engineering range without allocating
    // giant code plans. Execution coverage is separate below.
    noAlloc([&]{
        for(std::size_t k=unit;k<=(std::size_t{1}<<26);k+=unit) {
            require(spin::valid_packet_message_size(k),"natural size query rejected");
            require(!spin::valid_packet_message_size(k-1),"unaligned size query accepted");
            require(!spin::valid_packet_message_size(k+1),"unaligned size query accepted");
            for(auto value:{k-1,k,k+1})require(
                spin::valid_message_size(spin::Parameters::PacketT64S16,value)==spin::valid_packet_message_size(value),
                "unified natural size query disagrees with PacketCode");
        }
    });
    // Padded route offsets use 32 bits. Exercise the bound without allocating.
    constexpr auto maxK=(std::size_t{std::numeric_limits<std::uint32_t>::max()}/1028)*unit;
    require(spin::valid_packet_message_size(maxK),"last representable size rejected");
    require(!spin::valid_packet_message_size(maxK+unit),"padded route overflow accepted");
    require(spin::valid_message_size(spin::Parameters::PacketT64S16,maxK) &&
        !spin::valid_message_size(spin::Parameters::PacketT64S16,maxK+unit),"unified padded route limit");
    rejects([&]{spin::PacketCode code({maxK+unit,1});});
    // Around both the inner/region alignment boundary (16 groups) and the
    // measured powers of two. These execute the same selected backend.
    for(auto k:{15*unit,16*unit,17*unit,(std::size_t{1}<<16)-unit,
        (std::size_t{1}<<16)+unit,(std::size_t{1}<<18)-unit,(std::size_t{1}<<18)+unit})
        test(k,43);
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
    for(auto policy:{spin::MemoryPolicy::Normal,spin::MemoryPolicy::PreferHugePages}) {
        noAlloc([&]{
            bool failed=false;
            try{spin::Buffer overflow(std::numeric_limits<std::size_t>::max(),policy);}
            catch(const std::bad_alloc&){failed=true;}
            require(failed,"common buffer overflow accepted");
        });
    }
}

static void unifiedCode() {
    constexpr std::size_t k=1536;
    constexpr auto family=spin::Parameters::PacketT64S16;
    require(spin::message_alignment(family)==512 && spin::valid_message_size(family,k),"unified packet geometry");
    require(!spin::valid_message_size(family,k+1),"unified packet geometry accepted misalignment");
    for(auto seeds:{spin::CodeSeed{17,17},spin::CodeSeed{43,17},spin::CodeSeed{43,91}}) {
        const spin::CodeSpec spec{k,family,seeds.route,seeds.inner};
        spin::Code code(spec),reference(spec,{spin::Backend::Portable});
        auto work=code.make_workspace(),rw=reference.make_workspace();
        require(code.message_size()==k && code.code_size()==2*k && code.setup_bytes()>0,"unified packet shape");
        require(work.width()==spin::Width::Bits128,"unified packet workspace width");
        require(code.backend()==(spin::packet_fast_available()?spin::Backend::Avx512:spin::Backend::Portable),"unified packet backend");
        require(reference.backend()==spin::Backend::Portable && reference.descriptor()==code.descriptor(),"unified portable descriptor");
        require(code.supports_transpose() && !code.supports_forward() && !code.supports_generic_transpose(),"unified packet capabilities");
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
            frozen[4]=std::byte{1};frozen[8]=std::byte{1};
            for(unsigned i=0;i<8;++i) {frozen[16+i]=std::byte((k>>(8*i))&255);frozen[24+i]=std::byte((seeds.route>>(8*i))&255);}
            require(d==frozen,"PacketCode descriptor changed during API unification");
        }
        std::array<std::byte,40> descriptor{};
        descriptor[0]=std::byte{'S'};descriptor[1]=std::byte{'P'};descriptor[2]=std::byte{'I'};descriptor[3]=std::byte{'N'};
        descriptor[4]=std::byte{1};descriptor[8]=std::byte{4};
        for(unsigned i=0;i<8;++i) {
            descriptor[16+i]=std::byte((k>>(8*i))&255);
            descriptor[24+i]=std::byte((seeds.route>>(8*i))&255);
            descriptor[32+i]=std::byte((seeds.inner>>(8*i))&255);
        }
        require(code.descriptor()==descriptor,"unified packet descriptor or seed fields");
        rejects([&]{code.forward<Block>(expected,input,work);});
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
    rejects([&]{spin::Code c({k,family},{spin::Backend::Portable,8});});
    spin::Code fixedTile({k,family},{spin::Backend::Portable,4});
    require(fixedTile.backend()==spin::Backend::Portable,"explicit fixed packet tile");
    if(!spin::packet_fast_available()) {
        bool failed=false;
        try{spin::Code unavailable({k,family},{spin::Backend::Avx512});}catch(const std::runtime_error&){failed=true;}
        require(failed,"unified forced unavailable packet ISA");
    }
}

static void preparedPacket() {
    constexpr std::size_t k=1536;
    constexpr auto family=spin::Parameters::PacketT64S16;
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
        require(d[4]==std::byte{2} && d[8]==std::byte{4} && d[12]==std::byte{},"prepared packet descriptor family/mode");
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
    constexpr auto family=spin::Parameters::PacketT64S16;
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
    constexpr auto family=spin::Parameters::PacketT64S16;
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

int main() {try {
    for(auto k:{std::size_t{0},std::size_t{1},std::size_t{511},std::size_t{513},std::numeric_limits<std::size_t>::max()}) {
        require(!spin::valid_packet_message_size(k),"invalid K query accepted");rejects([&]{spin::PacketCode c({k,1});});
    }
    rejects([]{spin::PacketCode c({512,1},static_cast<spin::PacketBackend>(99));});
    if(!spin::packet_fast_available()) {
        bool failed=false;try{spin::PacketCode c({512,1},spin::PacketBackend::Avx512Gfni);}catch(const std::runtime_error&){failed=true;}
        require(failed,"forced unavailable packet ISA");
    }
    for(auto k:{512U,1024U,1536U,2560U,65536U,262144U})for(auto seed:{1ULL,17ULL})test(k,seed);
    sizing();knownAnswers();coordinateBasis();concurrent();commonBuffers();unifiedCode();preparedPacket();
    workspaceConversions();movedUnifiedHandles();
    std::cout<<"packet API PASS: scalar/fast, research known answers, coordinate basis, natural sizes, unified/prepared APIs, buffers, boundaries, alignment, ownership, no allocation\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
