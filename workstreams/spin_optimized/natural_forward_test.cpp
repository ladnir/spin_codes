#include "Spin.h"
#ifdef SPIN_TEST_WIDE
#include "Wide.h"
#endif
#include <bit>
#include <cstring>
#include <iostream>
#include <stdexcept>
using namespace bare_spin;
static void check(bool value,const char* message) {if(!value) throw std::runtime_error(message);}
template<class F> static void rejects(F f) {
    bool caught=false;try {f();} catch(const std::invalid_argument&) {caught=true;}
    check(caught,"invalid geometry accepted");
}
static void fill(std::vector<block>& a) {
    u64 seed=73;for(auto& v:a) {auto lo=splitmix(seed),hi=splitmix(seed);v=block(hi,lo);}
}
static void single(Configuration cfg,std::size_t k,bool large=false) {
    Spin oracle(cfg,MessageLength{k},17,29,0,BchBackend::Avx2);
    std::vector<block> in(k),expected(2*k),actual(2*k),dual(2*k),et(k),at(k);
    fill(in);fill(dual);
    if(large) {
        oracle.compact();Spin::Workspace w(oracle);
        oracle.forward(in.data(),k,expected.data(),2*k,w);
        oracle.encode(dual.data(),2*k,et.data(),k,w);
    } else {
        oracle.forwardReference(in.data(),expected.data());oracle.reference(dual.data(),et.data());
    }
    for(auto backend:{BchBackend::Auto,BchBackend::Avx2})
     for(unsigned tile:{2U,256U,1024U}) for(auto layout:{Layout::Packed24,Layout::Indices32}) {
        if(large && (backend==BchBackend::Avx2 || tile!=1024 || layout==Layout::Packed24)) continue;
        Spin c(cfg,MessageLength{k},17,29,tile,backend);c.validateSetup();
        check(c.routeHash()==oracle.routeHash(),"length/tile changed permutation");
        c.compact(layout);Spin::Workspace w(c);
        c.forward(in.data(),k,actual.data(),2*k,w);
        check(!std::memcmp(actual.data(),expected.data(),32*k),"natural forward mismatch");
        c.encode(dual.data(),2*k,at.data(),k,w);
        check(!std::memcmp(at.data(),et.data(),16*k),"natural bidirectional transpose mismatch");
    }
    if(std::has_single_bit(k)) {
        Spin legacy(cfg,unsigned(std::countr_zero(k)),17,29);
        check(legacy.routeHash()==oracle.routeHash(),"legacy constructor changed setup");
    }
}
#ifdef SPIN_TEST_WIDE
static void wide(Configuration cfg,std::size_t k,bool large=false) {
    for(unsigned lanes:{2U,4U}) {
        if(!wideAvailable(lanes)) continue;
        for(unsigned tile:{256U,1024U}) for(auto layout:{Layout::Packed24,Layout::Indices32}) {
            if(large && (tile!=1024 || layout==Layout::Packed24)) continue;
            Spin c(cfg,MessageLength{k},17,29,tile);c.validateSetup();
            if(large) {
                check(!c.packed24Available() && c.preferredLayout()==Layout::Indices32,"32-bit auto selection");
                rejects([&]{c.compact(Layout::Packed24);});
                c.compact();
            } else c.compact(layout);
            WideWorkspace w(c,lanes);
            std::vector<block> x(k*lanes),y(2*k*lanes),plane(k),expected(2*k);fill(x);
            w.forward(x,y);
            Spin control(cfg,MessageLength{k},17,29,0,BchBackend::Avx2);
            control.compact();Spin::Workspace cw(control);
            for(unsigned lane=0;lane<lanes;++lane) {
                for(std::size_t i=0;i<k;++i) plane[i]=x[i*lanes+lane];
                control.forward(plane.data(),k,expected.data(),2*k,cw);
                for(std::size_t i=0;i<2*k;++i)
                    check(!std::memcmp(&expected[i],&y[i*lanes+lane],16),"natural wide mismatch");
            }
        }
    }
}
#endif
int main(int argc,char** argv) {
 try {
    if(argc>1) {
        const auto k=std::stoull(argv[1]);
#ifdef SPIN_TEST_WIDE
        wide(Configuration::T128S19,k,true);
#else
        single(Configuration::T128S19,k,true);
#endif
        std::cout<<"large natural length PASS\n";return 0;
    }
    // Allocation-free boundary checks, including beyond the old 2^26 cap.
    check(length_geometry::valid(std::size_t{1}<<27,128),"old cap remains");
    check(length_geometry::valid(std::size_t{1}<<30,128),"32-bit range rejected");
    check(!length_geometry::valid(std::size_t{1}<<31,128),"32-bit range overflow");
    check(length_geometry::valid((std::size_t{1}<<31)-16384,128),"last aligned length rejected");
    for(auto cfg:{Configuration::T128S19,Configuration::T64S12,Configuration::T64S12R2}) {
        const std::size_t unit=cfg==Configuration::T128S19?16384:8192;
        for(auto k:{std::size_t{0},unit-128,unit+128,Spin::maxMessageBlocks+unit})
            rejects([&]{Spin c(cfg,MessageLength{k});});
        for(unsigned multiple:{1U,3U,5U,8U,17U,67U}) {
            const auto k=unit*multiple;
#ifdef SPIN_TEST_WIDE
            if(cfg==Configuration::T64S12) continue;
            wide(cfg,k);
#else
            single(cfg,k);
#endif
            std::cout<<"configuration="<<unsigned(cfg)<<" K="<<k<<" PASS\n";
        }
    }
    rejects([] {Spin c(Configuration::T128S19,31);});
    rejects([] {Spin c(Configuration::T128S19,64);});
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
