#include <spin/Code.h>
#include <spin/Generic.h>
#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <future>
#include <iostream>
#include <new>
#include <vector>
#include "../src/kernels/Spin.h"
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

template<unsigned L> struct alignas(16) Record {
    std::array<std::uint64_t,2*L> words{};
    bool operator==(const Record&) const=default;
    friend Record operator^(const Record& a,const Record& b) {
        Record c;for(unsigned i=0;i<2*L;++i)c.words[i]=a.words[i]^b.words[i];return c;
    }
};
using Block=Record<1>;
static void require(bool good,const char* why) {if(!good)throw std::runtime_error(why);}
template<class F> static void rejects(F f) {
    bool rejected=false;try{f();}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"invalid call accepted");
}
template<class F> static void no_alloc(F f) {
    const auto before=allocations.load();f();require(allocations.load()==before,"encoding allocated");
}
static std::uint64_t random_word(std::uint64_t& s) {s^=s<<13;s^=s>>7;s^=s<<17;return s;}
template<unsigned L> static void fill(std::vector<Record<L>>& x) {
    std::uint64_t s=139;for(auto& v:x)for(auto& u:v.words)u=random_word(s);
}
template<unsigned L> static void wide_test(const spin::Code& c) {
    const auto width=L==2?spin::Width::Bits256:spin::Width::Bits512;
    if(!c.supports_forward(width))return;
    auto w=c.make_workspace(width);auto one=c.make_workspace();
    std::vector<Record<L>> x(c.message_size()),y(c.code_size());fill(x);
    no_alloc([&]{c.forward<Record<L>>(x,y,w);});
    std::vector<Block> a(c.message_size()),b(c.code_size());
    for(unsigned lane=0;lane<L;++lane) {
        for(std::size_t i=0;i<a.size();++i)for(unsigned j=0;j<2;++j)a[i].words[j]=x[i].words[2*lane+j];
        c.forward<Block>(a,b,one);
        for(std::size_t i=0;i<b.size();++i)for(unsigned j=0;j<2;++j)
            require(b[i].words[j]==y[i].words[2*lane+j],"wide lane differs");
    }
}
static void test(spin::Parameters p,std::size_t k,spin::Backend backend,unsigned tile) {
    spin::Code c({k,p,17,29},{backend,tile});auto w=c.make_workspace();
    require(c.message_size()==k && c.code_size()==2*k,"shape");
    std::vector<Block> x(k),y(2*k),in(2*k),tr(k);fill(x);fill(in);
    no_alloc([&]{c.forward<Block>(x,y,w);c.transpose<Block>(in,tr,w);});
    Block left{},right{};
    for(std::size_t i=0;i<k;++i)for(unsigned j=0;j<2;++j)left.words[j]^=x[i].words[j]&tr[i].words[j];
    for(std::size_t i=0;i<2*k;++i)for(unsigned j=0;j<2;++j)right.words[j]^=y[i].words[j]&in[i].words[j];
    require(left==right,"adjoint identity");
    auto inplace=in;
    no_alloc([&]{c.transpose_inplace<Block>(inplace,w);});
    require(std::equal(tr.begin(),tr.end(),inplace.begin()),"inplace result");
    require(std::equal(in.begin()+k,in.end(),inplace.begin()+k),"inplace changed suffix");

    // Independent dense oracle, not a call to the production recurrence.
    namespace kernel=spin::detail::kernel;
    const auto cfg=p==spin::Parameters::T128S19?kernel::Configuration::T128S19:
        p==spin::Parameters::T64S12?kernel::Configuration::T64S12:kernel::Configuration::T64S12R2;
    kernel::Spin reference(cfg,kernel::MessageLength{k},17,29,256,kernel::BchBackend::Avx2);
    std::vector<kernel::block> ri(2*k),ro(2*k);
    std::memcpy(ri.data(),x.data(),k*16);reference.forwardReference(ri.data(),ro.data());
    require(!std::memcmp(y.data(),ro.data(),2*k*16),"forward dense oracle");
    std::memcpy(ri.data(),in.data(),2*k*16);reference.reference(ri.data(),ro.data());
    require(!std::memcmp(tr.data(),ro.data(),k*16),"transpose dense oracle");

    auto generic=c.generic_transpose();auto gw=generic.make_workspace<std::uint8_t>();
    std::vector<std::uint8_t> gi(2*k),go(k);
    for(std::size_t i=0;i<2*k;++i)gi[i]=std::uint8_t(in[i].words[0]);
    no_alloc([&]{generic.transpose<std::uint8_t>(gi,go,gw);});
    for(std::size_t i=0;i<k;++i)require(go[i]==std::uint8_t(tr[i].words[0]),"generic byte mismatch");
    no_alloc([&]{generic.transpose_inplace<std::uint8_t>(gi,gw);});
    require(std::equal(go.begin(),go.end(),gi.begin()),"generic inplace");
    auto g64=generic.make_workspace<std::uint64_t>();
    std::vector<std::uint64_t> i64(2*k),o64(k);
    for(std::size_t i=0;i<2*k;++i)i64[i]=in[i].words[1];
    no_alloc([&]{generic.transpose<std::uint64_t>(i64,o64,g64);});
    for(std::size_t i=0;i<k;++i)require(o64[i]==tr[i].words[1],"generic uint64 mismatch");

    std::vector<std::uint64_t> bits(k/64),encoded(2*k/64),scratch(2*k/64);
    std::uint64_t seed=18;for(auto& b:bits)b=random_word(seed);
    for(std::size_t i=0;i<k;++i)x[i]={{(bits[i/64]>>(i%64))&1,0}};
    c.forward<Block>(x,y,w);
    no_alloc([&]{c.forward_bits(bits,encoded,scratch);});
    for(std::size_t i=0;i<2*k;++i)require(((encoded[i/64]>>(i%64))&1)==y[i].words[0],"packed bits");

    rejects([&]{c.forward<Block>(std::span<const Block>(x).first(k-1),y,w);});
    rejects([&]{c.forward_bytes(std::as_bytes(std::span<const Block>(in)).first(k*16),std::as_writable_bytes(std::span(in)),w);});
    rejects([&]{c.transpose<Block>(in,std::span<Block>(in).subspan(1,k),w);});
    spin::Code other({k,p,17,29});auto wrong=other.make_workspace();
    rejects([&]{c.forward<Block>(x,y,wrong);});
    auto otherGeneric=other.generic_transpose();auto wrongGeneric=otherGeneric.make_workspace<std::uint8_t>();
    rejects([&]{generic.transpose<std::uint8_t>(gi,go,wrongGeneric);});
    auto moved=std::move(w);rejects([&]{c.forward<Block>(x,y,w);});
    auto shared=c;auto relocated=std::move(c);
    shared.forward<Block>(x,y,moved);relocated.forward<Block>(x,y,moved);
    require(shared.descriptor()==other.descriptor(),"backend/tile affected descriptor");
    wide_test<2>(shared);wide_test<4>(shared);
    std::cout<<"parameters="<<unsigned(p)<<" K="<<k<<" backend="<<unsigned(backend)<<" PASS\n";
}
int main() {try {
    if(!spin::capabilities().avx2)return 77;
    std::cout<<"AVX512 BCH="<<spin::capabilities().avx512_bch
             <<" wide512="<<spin::capabilities().forward512<<'\n';
    for(auto p:{spin::Parameters::T128S19,spin::Parameters::T64S12,spin::Parameters::T64S12R2}) {
        const auto unit=spin::message_alignment(p);
        require(spin::valid_message_size(p,std::size_t{1}<<30),"representation boundary");
        require(!spin::valid_message_size(p,std::size_t{1}<<31),"overflow accepted");
        rejects([&]{spin::Code c({unit-1,p});});
        for(auto b:{spin::Backend::Avx2,spin::Backend::Automatic}) {
            test(p,3*unit,b,1024);test(p,65536,b,256);
        }
    }
    // Exercise range-direct without a forward-direct table, then tiled routing.
    test(spin::Parameters::T128S19,393216,spin::Backend::Automatic,256);
    test(spin::Parameters::T128S19,540672,spin::Backend::Automatic,1024);
    rejects([]{spin::Code c({16384,static_cast<spin::Parameters>(99)});});
    spin::Code c({16384});auto a=c.make_workspace(),b=c.make_workspace();
    rejects([&]{auto bad=c.make_workspace(static_cast<spin::Width>(99));});
    if(!spin::capabilities().avx512_bch)
        rejects([]{spin::Code forced({16384},{spin::Backend::Avx512,256});});
    const auto descriptor=c.descriptor();
    require(descriptor[0]==std::byte{'S'} && descriptor[4]==std::byte{1} &&
            descriptor[8]==std::byte{1} && descriptor[16]==std::byte{0} &&
            descriptor[17]==std::byte{0x40} && descriptor[24]==std::byte{1} &&
            descriptor[32]==std::byte{2},"descriptor encoding");
    require(c.descriptor()!=spin::Code({16384,spin::Parameters::T128S19,2,2}).descriptor(),"route seed identity");
    require(c.descriptor()!=spin::Code({16384,spin::Parameters::T128S19,1,3}).descriptor(),"inner seed identity");
    std::vector<Block> in(16384),x(32768),y(32768);fill(in);
    std::vector<std::byte> unaligned(in.size()*16+1);
    rejects([&]{c.forward_bytes(std::span<const std::byte>(unaligned).subspan(1),std::as_writable_bytes(std::span(x)),a);});
    auto f=std::async(std::launch::async,[&]{c.forward<Block>(in,x,a);});
    c.forward<Block>(in,y,b);f.get();require(x==y,"concurrent workspaces");
    auto survivor=[] {spin::Code local({16384});return local.make_workspace();}();
    require(survivor.bytes()>0,"workspace failed to retain owner");
    auto generic=[] {spin::Code local({16384});return local.generic_transpose();}();
    auto gw=generic.make_workspace<std::uint8_t>();std::vector<std::uint8_t> gi(32768),go(16384);
    generic.transpose<std::uint8_t>(gi,go,gw);
    std::cout<<"API, dense references, ownership, concurrency, no-allocation PASS\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
