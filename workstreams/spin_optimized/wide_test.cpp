#include "Wide.h"
#include <cstring>
#include <iostream>
#include <stdexcept>
using namespace bare_spin;
static void check(bool b,const char* s) {if(!b) throw std::runtime_error(s);}
template<class F> static void rejects(F f) {
    bool rejected=false;try {f();} catch(const std::invalid_argument&) {rejected=true;}
    check(rejected,"invalid request accepted");
}
int main() {
 try {
    rejects([] {Spin c(Configuration::T128S19,14);WideWorkspace w(c,3);});
    for(unsigned lanes:{2U,4U}) {
        if(!wideAvailable(lanes)) {
            rejects([&] {Spin c(Configuration::T128S19,14);WideWorkspace w(c,lanes);});
            continue;
        }
        rejects([&] {Spin c(Configuration::T64S12,16);WideWorkspace w(c,lanes);});
        rejects([&] {Spin c(Configuration::T128S19,14);c.compact(Layout::Indices32);WideWorkspace w(c,lanes);});
        for(auto cfg:{Configuration::T128S19,Configuration::T64S12R2})
         for(unsigned m:{14U,16U,20U}) for(unsigned tile:{0U,2U,64U,128U,512U}) for(bool compact:{false,true}) {
            if(cfg==Configuration::T64S12R2 && m!=16) continue;
            if(m==20 && ((tile!=64 && tile!=512) || !compact)) continue;
            Spin c(cfg,m,17,29,tile);
            if(compact) c.compact();
            const auto k=c.messageBlocks(),n=c.codeBlocks();
            WideWorkspace wide(c,lanes);
            // Independent default tile and AVX2 backend: tile tuning is a
            // schedule change, never a change to the encoded linear map.
            Spin control(cfg,m,17,29,0,BchBackend::Avx2);control.compact();
            Spin::Workspace controlWork(control);
            std::vector<block> x(k), y(n), input(k*lanes+4), output(n*lanes+4);
            // Deliberately 16-byte aligned but not 32/64-byte aligned.
            auto* in=input.data();while((reinterpret_cast<std::uintptr_t>(in)&63)!=16) ++in;
            auto* out=output.data();while((reinterpret_cast<std::uintptr_t>(out)&63)!=16) ++out;
            u64 seed=123;
            for(std::size_t i=0;i<k*lanes;++i) {auto lo=splitmix(seed),hi=splitmix(seed);in[i]=block(hi,lo);}
            wide.forward({in,k*lanes},{out,n*lanes});
            for(unsigned l=0;l<lanes;++l) {
                for(std::size_t i=0;i<k;++i) x[i]=in[i*lanes+l];
                if(cfg==Configuration::T64S12R2 && !compact) c.forwardReference(x.data(),y.data());
                else control.forward(x.data(),k,y.data(),n,controlWork);
                for(std::size_t i=0;i<n;++i) check(!std::memcmp(&out[i*lanes+l],&y[i],16),"wide plane differs");
            }
            rejects([&] {wide.forward({in,k*lanes-1},{out,n*lanes});});
            rejects([&] {wide.forward({out,k*lanes},{out,n*lanes});});
            std::cout<<"lanes="<<lanes<<" config="<<c.name()<<" m="<<m<<" tile="<<tile<<" compact="<<compact<<" PASS\n";
        }
    }
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
