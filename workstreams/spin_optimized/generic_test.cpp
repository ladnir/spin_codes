#include "GenericSpin.h"
#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <new>
using namespace bare_spin;
static std::atomic<std::size_t> allocations=0;
void* operator new(std::size_t n) {++allocations;if(auto p=std::malloc(n?n:1)) return p;throw std::bad_alloc();}
void* operator new[](std::size_t n) {return ::operator new(n);}
void operator delete(void* p) noexcept {std::free(p);}
void operator delete[](void* p) noexcept {std::free(p);}
void operator delete(void* p,std::size_t) noexcept {std::free(p);}
void operator delete[](void* p,std::size_t) noexcept {std::free(p);}
void* operator new(std::size_t n,std::align_val_t a) {
    ++allocations;void* p=nullptr;
    if(posix_memalign(&p,std::size_t(a),n?n:1)) throw std::bad_alloc();return p;
}
void* operator new[](std::size_t n,std::align_val_t a) {return ::operator new(n,a);}
void operator delete(void* p,std::align_val_t) noexcept {std::free(p);}
void operator delete[](void* p,std::align_val_t) noexcept {std::free(p);}
void operator delete(void* p,std::size_t,std::align_val_t) noexcept {std::free(p);}
void operator delete[](void* p,std::size_t,std::align_val_t) noexcept {std::free(p);}
template<unsigned W> struct alignas(W*8) Wide {
    std::array<u64,W> words;
    friend Wide operator^(const Wide& a,const Wide& b) {
        Wide c;for(unsigned i=0;i<W;++i)c.words[i]=a.words[i]^b.words[i];return c;
    }
};
static_assert(XorElement<u8> && XorElement<u64> && XorElement<Wide<8>>);
static_assert(!XorElement<bool>);
static void check(bool b,const char* why) {if(!b)throw std::runtime_error(why);}
template<class F> static void rejects(F f) {
    bool caught=false;try {f();}catch(const std::invalid_argument&){caught=true;}
    check(caught,"invalid generic arguments accepted");
}
template<XorElement E> void compare(const Spin& reference,const GenericTranspose& generic,bool checks=false) {
    const auto n=reference.codeBlocks(),k=reference.messageBlocks();
    std::vector<E> in(n),out(k),expected(k);u64 seed=77;
    for(auto& x:in) {
        auto* bytes=reinterpret_cast<unsigned char*>(&x);
        for(std::size_t j=0;j<sizeof(E);++j)bytes[j]=u8(splitmix(seed));
    }
    Spin::Workspace rw(reference);std::vector<block> packed(n),answer(k);
    for(std::size_t offset=0;offset<sizeof(E);offset+=16) {
        const auto bytes=std::min<std::size_t>(16,sizeof(E)-offset);
        for(std::size_t i=0;i<n;++i) {packed[i]=block(0,0);std::memcpy(&packed[i],reinterpret_cast<const unsigned char*>(&in[i])+offset,bytes);}
        reference.encode(packed.data(),n,answer.data(),k,rw);
        for(std::size_t i=0;i<k;++i)std::memcpy(reinterpret_cast<unsigned char*>(&expected[i])+offset,&answer[i],bytes);
    }
    GenericTranspose::Workspace<E> w(generic);
    const auto before=allocations.load();
    generic.encode<E>(in,out,w);
    check(allocations.load()==before,"generic encode allocated");
    check(!std::memcmp(expected.data(),out.data(),k*sizeof(E)),"generic lane mismatch");
    auto inplace=in;
    const auto beforeInplace=allocations.load();
    generic.encodeInplace<E>(inplace,w);
    check(allocations.load()==beforeInplace,"generic inplace allocated");
    check(!std::memcmp(expected.data(),inplace.data(),k*sizeof(E)),"generic inplace mismatch");
    check(!std::memcmp(in.data()+k,inplace.data()+k,k*sizeof(E)),"generic suffix changed");
    if(checks) {
        rejects([&]{generic.encode<E>(std::span<const E>(in.data(),n-1),out,w);});
        rejects([&]{generic.encode<E>(in,std::span<E>(out.data(),k-1),w);});
        rejects([&]{generic.encode<E>(in,std::span<E>(in.data()+1,k),w);});
        rejects([&]{generic.encode<E>(std::span<const E>(w.values),out,w);});
        w.values.pop_back();rejects([&]{generic.encode<E>(in,out,w);});
    }
}
int main() {try {
    for(auto cfg:{Configuration::T128S19,Configuration::T64S12,Configuration::T64S12R2}) {
        const std::size_t unit=cfg==Configuration::T128S19?16384:8192;
        Spin reference(cfg,MessageLength{3*unit},17,29,0,Outer::Bch256x128,BchBackend::Avx2);
        // Snapshot owns its data and survives destruction of the source object.
        auto generic=[&]{Spin c(cfg,MessageLength{3*unit},17,29);return c.generic();}();
        compare<u8>(reference,generic,true);compare<std::byte>(reference,generic);
        compare<std::uint32_t>(reference,generic);compare<u64>(reference,generic);
        compare<block>(reference,generic);compare<Wide<4>>(reference,generic);compare<Wide<8>>(reference,generic);
        // Direct K16/K18, runtime-length direct, full tiles, and partial tiles.
        for(const auto k:{unit,3*unit,std::size_t{65536},std::size_t{262144},67*unit}) {
            Spin ref(cfg,MessageLength{k},31,47,0,Outer::Bch256x128,BchBackend::Avx2);
            for(auto backend:{BchBackend::Auto,BchBackend::Avx2})
             for(auto layout:{Layout::Packed24,Layout::Indices32}) {
                Spin c(cfg,MessageLength{k},31,47,0,Outer::Bch256x128,backend);
                c.compact(layout);auto g=c.generic();compare<u64>(ref,g);
            }
        }
        std::cout<<reference.name()<<" generic widths and compacted routing PASS\n";
    }
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
