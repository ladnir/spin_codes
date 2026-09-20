#include "Spin.h"
#include <algorithm>
#include <bit>
#include <cstring>
#include <iostream>
#include <stdexcept>
using namespace bare_spin;
static void check(bool b,const char* m) {if(!b) throw std::runtime_error(m);}
template<class F> static void rejects(F f) {
    bool rejected=false;try {f();} catch(const std::invalid_argument&) {rejected=true;}
    check(rejected,"invalid geometry accepted");
}
int main(int argc,char** argv) {
 try {
    check(length_geometry::valid(std::size_t{1}<<27,128),"old cap remains");
    check(length_geometry::valid(std::size_t{1}<<30,128),"32-bit geometry rejected");
    check(!length_geometry::valid(std::size_t{1}<<31,128),"32-bit overflow accepted");
    // Optional large test crosses Packed24's exact boundary, without a huge
    // dense oracle. Exercise both tail dispatch and automatic 32-bit routing.
    if(argc>1) {
        const std::size_t k=std::stoull(argv[1]);
        Spin c(Configuration::T128S19,MessageLength{k});
        check(!c.packed24Available() && c.preferredLayout()==Layout::Indices32,"large layout");
        rejects([&]{c.compact(Layout::Packed24);});
        c.validateSetup();c.compact();Spin::Workspace w(c);
        std::vector<block> in(2*k,block(0,0)),out(k,block(1,1));
        c.encode(in.data(),in.size(),out.data(),out.size(),w);
        for(const auto& v:out) check(_mm_testz_si128(v.mData,v.mData),"large zero map");
        // Nonzero, then AVX2 comparison with identical routing and coefficients.
        in[0]=block(3,5);in[k+7]=block(11,13);in.back()=block(17,19);
        c.encode(in.data(),in.size(),out.data(),out.size(),w);
        Spin fallback(Configuration::T128S19,MessageLength{k},1,2,0,Outer::Bch256x128,BchBackend::Avx2);
        fallback.compact();Spin::Workspace fw(fallback);std::vector<block> expected(k);
        fallback.encode(in.data(),in.size(),expected.data(),expected.size(),fw);
        check(!std::memcmp(out.data(),expected.data(),k*16),"large backend mismatch");
        std::cout<<"large 32-bit indices, partial tile, auto/AVX2 PASS\n";return 0;
    }
    for(auto cfg:{Configuration::T128S19,Configuration::T64S12,Configuration::T64S12R2}) {
        const std::size_t unit=cfg==Configuration::T128S19?16384:8192;
        for(auto k:{std::size_t(0),unit-128,unit+128,Spin::maxMessageBlocks+unit})
            rejects([&]{Spin c(cfg,MessageLength{k});});
        for(unsigned mult:{1U,3U,5U,8U,17U,32U,67U}) {
            const auto k=unit*mult;
            Spin oracle(cfg,MessageLength{k},17,29,0,Outer::Bch256x128,BchBackend::Avx2);
            oracle.validateSetup();
            std::vector<block> input(2*k),expected(k),actual(k);u64 seed=123;
            for(auto& v:input) {auto lo=splitmix(seed),hi=splitmix(seed);v=block(hi,lo);}
            oracle.reference(input.data(),expected.data());
            for(auto backend:{BchBackend::Auto,BchBackend::Avx2}) for(unsigned tile:{0U,2U,64U,256U})
             for(auto layout:{Layout::Packed24,Layout::Indices32}) {
                Spin c(cfg,MessageLength{k},17,29,tile,Outer::Bch256x128,backend);
                check(c.routeHash()==oracle.routeHash(),"tile changed permutation");
                c.validateSetup();Spin::Workspace w(c);
                c.encode(input.data(),input.size(),actual.data(),actual.size(),w,layout);
                check(!std::memcmp(actual.data(),expected.data(),k*16),"aligned oracle mismatch");
                c.compact(layout);
                auto inplace=input;c.encodeInplace(inplace.data(),inplace.size(),w);
                check(!std::memcmp(inplace.data(),expected.data(),k*16),"compacted auto/inplace mismatch");
                check(!std::memcmp(inplace.data()+k,input.data()+k,k*16),"inplace suffix changed");
            }
            if(std::has_single_bit(k)) {
                Spin legacy(cfg,unsigned(std::countr_zero(k)),17,29);
                check(legacy.routeHash()==oracle.routeHash(),"legacy length changed route");
            }
            std::cout<<oracle.name()<<" K="<<k<<" PASS\n";
        }
    }
    rejects([] {Spin c(Configuration::T128S19,31);});
    rejects([] {Spin c(Configuration::T128S19,64);});
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
