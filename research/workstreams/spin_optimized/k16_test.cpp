#include "Spin.h"
#include <array>
#include <cstring>
#include <iostream>
#include <stdexcept>
using namespace bare_spin;
static void check(bool v,const char* s) {if(!v) throw std::runtime_error(s);}
static bool equal(const std::vector<block>& a,const std::vector<block>& b) {return !std::memcmp(a.data(),b.data(),a.size()*16);}
static Configuration selected=Configuration::T64S12;
static u64 coefficient=2;
static Spin make(unsigned m,u64 seed,unsigned tile,BchBackend backend) {
    return Spin(selected,m,seed,coefficient,tile,
#ifndef SPIN_K16_BIDIRECTIONAL
        Outer::Bch256x128,
#endif
        backend);
}
template<class F> static void rejects(F f) {
    bool rejected=false;try {f();} catch(const std::invalid_argument&) {rejected=true;}
    check(rejected,"invalid configuration was accepted");
}
static __m128i dot(const std::vector<block>& a,const std::vector<block>& b) {
    auto x=_mm_setzero_si128();for(std::size_t i=0;i<a.size();++i) x=_mm_xor_si128(x,_mm_and_si128(a[i].mData,b[i].mData));return x;
}
int main() {
 try {
   for(auto cfg:{Configuration::T64S12,Configuration::T64S12R2}) {
    selected=cfg;
   for(u64 cs:{2ULL,23ULL}) {
    coefficient=cs;
#if !SPIN_GENERAL_LENGTHS
    for(unsigned m:{14U,15U,17U,18U,20U}) rejects([&]{make(m,1,0,BchBackend::Auto);});
#endif
    rejects([]{make(16,1,1,BchBackend::Auto);});
    rejects([]{make(16,1,2,BchBackend::Avx512);});
    if(!bchAvx512Available()) rejects([]{make(16,1,0,BchBackend::Avx512);});
    for(u64 routeSeed:{1ULL,17ULL}) {
        auto oracle=make(16,routeSeed,0,BchBackend::Avx2);
        std::vector<block> q(oracle.codeBlocks()),z(oracle.messageBlocks()),ez(z.size());
        u64 seed=123;for(auto& v:q) {const auto lo=splitmix(seed),hi=splitmix(seed);v=block(hi,lo);}
        oracle.reference(q.data(),ez.data());
#ifdef SPIN_K16_BIDIRECTIONAL
        std::vector<block> x(z.size()),y(q.size()),ey(q.size());
        for(auto& v:x) {const auto lo=splitmix(seed),hi=splitmix(seed);v=block(hi,lo);}
        oracle.forwardReference(x.data(),ey.data());
#endif
        std::vector<BchBackend> backends{BchBackend::Auto,BchBackend::Avx2};
        if(bchAvx512Available()) backends.push_back(BchBackend::Avx512);
        for(auto layout:{Layout::Packed24,Layout::Indices32}) for(unsigned tile:{0U,2U,4U,16U}) for(auto backend:backends) {
            if(tile==2 && backend==BchBackend::Avx512) continue;
            auto code=make(16,routeSeed,tile,backend);code.validateSetup();
            check(code.step()==64 && code.state()==12,"wrong configuration dimensions");
            check(code.routeHash()==oracle.routeHash(),"route changed with backend/tile");
#ifdef SPIN_K16_BIDIRECTIONAL
            rejects([&]{code.wideView();});
#endif
            Spin::Workspace w(code);
            for(unsigned compact=0;compact<2;++compact) {
                if(compact) code.compact(layout);
                code.encode(q.data(),q.size(),z.data(),z.size(),w,layout);
                check(equal(z,ez),"transpose dense oracle mismatch");
                auto inplace=q;code.encodeInplace(inplace.data(),inplace.size(),w,layout);
                check(!std::memcmp(inplace.data(),ez.data(),ez.size()*16),"inplace output mismatch");
                check(!std::memcmp(inplace.data()+ez.size(),q.data()+ez.size(),ez.size()*16),"inplace suffix changed");
#ifdef SPIN_K16_BIDIRECTIONAL
                code.forward(x.data(),x.size(),y.data(),y.size(),w,layout);
                check(equal(y,ey),"forward dense oracle mismatch");
                check(_mm_movemask_epi8(_mm_cmpeq_epi8(dot(y,q),dot(x,z)))==0xffff,"adjoint mismatch");
#endif
            }
#ifdef SPIN_K16_BIDIRECTIONAL
            std::vector<u64> bits(x.size()/64),out(q.size()/64),scratch(out.size());
            std::vector<block> bx(x.size()),by(q.size());
            for(auto& v:bits) v=splitmix(seed);
            for(std::size_t i=0;i<bx.size();++i) bx[i]=block(0,(bits[i/64]>>(i%64))&1);
            code.forwardBits(bits.data(),bits.size(),out.data(),out.size(),scratch.data(),scratch.size());
            code.forward(bx.data(),bx.size(),by.data(),by.size(),w,layout);
            for(std::size_t i=0;i<by.size();++i) check((u64(_mm_cvtsi128_si64(by[i].mData))&1)==((out[i/64]>>(i%64))&1),"packed bit forward mismatch");
            rejects([&]{code.forwardBits(bits.data(),bits.size(),scratch.data(),out.size(),scratch.data(),scratch.size());});
#endif
            rejects([&]{code.compact(layout==Layout::Packed24?Layout::Indices32:Layout::Packed24);});
        }
        u64 hash=0xcbf29ce484222325ULL;
        for(const auto& v:ez) {std::array<u64,2> words;std::memcpy(words.data(),&v,16);for(auto word:words) {hash^=word;hash*=0x100000001b3ULL;}}
        if(coefficient==2) {
            const u64 expected=selected==Configuration::T64S12 ?
                (routeSeed==1?0xe3b2dd281c937281ULL:0xc6686f2d018544ddULL) :
                (routeSeed==1?0xbfb9cacfde1d6b2fULL:0x988f0b98f31f77dfULL);
            check(hash==expected,"frozen r1/prototype-r2 output changed");
        }
        std::cout<<oracle.name()<<" seed="<<routeSeed<<" coefficient="<<coefficient<<" transpose_hash="<<std::hex<<hash<<std::dec<<" PASS\n";
    }
   }
   }
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
