#include "Spin.h"
#include "BenchmarkLock.h"
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
using namespace bare_spin;
extern "C" void spin_probe_routed(std::size_t,const u32*,const u32*,const block*,block*);
extern "C" void spin_probe_contiguous(std::size_t,const u32*,const block*,block*);
extern "C" void spin_probe_gather(std::size_t,const u32*,const block*,block*);
static volatile u64 sink;
int main(int argc,char** argv) {
 try {
    const unsigned m=argc>1?std::stoul(argv[1]):18,trials=argc>2?std::stoul(argv[2]):101;
    const unsigned mode=argc>3?std::stoul(argv[3]):0;
    if((m!=18 && m!=20) || trials<3 || !(trials&1) || mode>2 || argc>4) throw std::invalid_argument("m=18|20, odd trials>=3, mode=0|1|2");
    lockBenchmarks();
    if(!bchAvx512Available()) throw std::runtime_error("probe requires AVX-512");
    Spin code(Configuration::T128S19,m,1,2);code.compact();
    const auto v=code.wideView();const auto n=v.n;
    const auto unpack=[](const u8* p){u32 x;std::memcpy(&x,p,4);return x&0xffffff;};
    std::vector<u32> route(n);
    for(std::size_t i=0;i<n;++i) {
        const auto slot=unpack(v.slots+3*i);
        const auto pos=(slot/v.tile)*v.tile+unpack(v.offsets+3*slot);
        route[i]=(pos&~1023U)|((pos&255U)<<2)|((pos>>8)&3U);
    }
    std::vector<block> input(n),linear(n),out(n),expected(n);u64 seed=123;
    for(auto& x:input) {const auto lo=splitmix(seed),hi=splitmix(seed);x=block(hi,lo);}
    for(std::size_t i=0;i<n;++i) linear[i]=input[route[i]];
    spin_probe_contiguous(n,v.fieldRows,linear.data(),expected.data());
    std::cout<<"m,trial,mode,ms\n"<<std::setprecision(10);
    for(int trial=-4;trial<int(trials);++trial) {
        const auto start=std::chrono::steady_clock::now();
        if(mode==0) spin_probe_routed(n,v.fieldRows,route.data(),input.data(),out.data());
        else if(mode==1) spin_probe_contiguous(n,v.fieldRows,linear.data(),out.data());
        else spin_probe_gather(n,route.data(),input.data(),out.data());
        const auto ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
        if(trial==-4 || trial==int(trials)-1) {
            const auto* check=mode==2?linear.data():expected.data();
            if(std::memcmp(out.data(),check,n*sizeof(block))) throw std::runtime_error("probe mismatch");
        }
        u64 word;std::memcpy(&word,out.data(),8);sink=word;
        if(trial>=0) std::cout<<m<<','<<trial<<','<<(mode==0?"routed":mode==1?"contiguous":"gather")<<','<<ms<<'\n';
    }
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
