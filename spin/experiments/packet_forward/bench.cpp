#include "Forward.h"
#include <spin/Code.h>
#include <spin/PacketCode.h>
#include <algorithm>
#include <chrono>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <stdexcept>
using namespace spin::detail::packet;
void fill(std::span<std::byte> bytes,std::uint64_t state=913) {
    for(std::size_t i=0;i<bytes.size();i+=8) {
        state+=0x9e3779b97f4a7c15ULL;auto v=state;v=(v^(v>>30))*0xbf58476d1ce4e5b9ULL;
        v=(v^(v>>27))*0x94d049bb133111ebULL;v^=v>>31;
        std::memcpy(bytes.data()+(i^8),&v,8);
    }
}
Block* ptr(spin::Buffer& b) {return reinterpret_cast<Block*>(b.bytes().data());}
std::uint64_t hash(const Block* p,std::size_t n) {
    std::uint64_t h=0;for(std::size_t i=0;i<2*n;++i) {
        std::uint64_t x;std::memcpy(&x,reinterpret_cast<const char*>(p)+8*i,8);h=(h^x)*0x100000001b3ULL;
    }return h;
}
void verify() {
    unsigned cases=0;
    for(auto k:{256U,768U,16384U,65536U,262144U,1048576U})for(auto seed:{1U,17U}) {
        Plan p(k,seed,seed+3);ForwardPlan f(p);
        spin::Buffer input(16*p.k),output(16*p.n+64),ref(16*p.n),scratch(16*p.scratchBlocks()),trans(16*p.k);
        for(unsigned mode=0;mode<3;++mode) {
            fill(input.bytes(),913+mode);
            if(mode!=0)std::memset(input.bytes().data(),0,input.bytes().size());
            if(mode==2)input.bytes()[19]=std::byte{1};
            std::memset(scratch.bytes().data(),0xa5,scratch.bytes().size());
            std::memset(output.bytes().data(),0x93,output.bytes().size());
            forwardScalar(ptr(input),ptr(ref),p);
            // Exercise minimum 16-byte output alignment as well as 64-byte.
            auto* dest=ptr(output)+(mode==1?1:0);
            forwardPacket(ptr(input),dest,ptr(scratch),p,f);
            if(std::memcmp(dest,ptr(ref),16*p.n))throw std::runtime_error("forward mismatch K="+std::to_string(k)+" mode="+std::to_string(mode));
            forwardOuterGather(ptr(input),ptr(scratch),p,f);forwardInnerGather(ptr(scratch),dest,p,f);
            if(std::memcmp(dest,ptr(ref),16*p.n))throw std::runtime_error("gather forward mismatch");
            forwardOuter(ptr(input),ptr(output),p,f);forwardInner(ptr(output),ptr(output),p,f);
            if(std::memcmp(ptr(output),ptr(ref),16*p.n))throw std::runtime_error("inplace forward mismatch");
            forwardSelected(ptr(input),ptr(output),ptr(scratch),p,f);
            if(std::memcmp(ptr(output),ptr(ref),16*p.n))throw std::runtime_error("selected forward mismatch");
            const auto* guard=reinterpret_cast<const unsigned char*>(dest+p.n);
            for(unsigned j=0;j<16;++j)if(guard[j]!=0x93)throw std::runtime_error("output guard");
            ++cases;
        }
        fill(ref.bytes(),47);
        forwardPacket(ptr(input),ptr(output),ptr(scratch),p,f);
        transposeFast(ptr(ref),ptr(trans),ptr(scratch),p);
        __m128i left{},right{};
        for(std::size_t i=0;i<p.n;++i)left=_mm_xor_si128(left,_mm_and_si128(ptr(output)[i].mData,ptr(ref)[i].mData));
        for(std::size_t i=0;i<p.k;++i)right=_mm_xor_si128(right,_mm_and_si128(ptr(input)[i].mData,ptr(trans)[i].mData));
        if(_mm_movemask_epi8(_mm_cmpeq_epi8(left,right))!=65535)throw std::runtime_error("adjoint mismatch");
    }
    std::cout<<"Verified "<<cases<<" independent literal forward comparisons plus adjoint identities.\n";
}
int main(int argc,char** argv) {try {
    if(!spin::packet_fast_available()){std::cerr<<"AVX512/VBMI/GFNI required\n";return 77;}
    if(argc==2 && std::string(argv[1])=="verify"){verify();return 0;}
    if(argc<3)throw std::invalid_argument("MODE K [seed] [calls] [normal|huge]");
    const std::string mode=argv[1];const auto k=std::stoull(argv[2]);
    const auto seed=argc>3?std::stoull(argv[3]):1,calls=argc>4?std::stoull(argv[4]):501;
    const std::string policy=argc>5?argv[5]:"normal";
    const auto memory=policy=="huge"?spin::MemoryPolicy::PreferHugePages:spin::MemoryPolicy::Normal;
    if(!calls || (mode!="selected" && mode!="forward" && mode!="gather" && mode!="transpose" && mode!="outer" && mode!="inner"))throw std::invalid_argument("mode/calls");
    Plan p(k,seed);ForwardPlan f(p);
    spin::Buffer input(16*p.n,memory),output(16*p.n,memory),scratch(16*p.scratchBlocks(),memory);
    fill(input.bytes());fill(output.bytes());fill(scratch.bytes());
    // Dispatch outside measured calls; kernels allocate no storage.
    std::vector<double> times(calls);
    auto measure=[&](auto&& operation) {
        for(unsigned i=0;i<5;++i)operation();
        for(auto& t:times){const auto start=std::chrono::steady_clock::now();operation();t=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();}
    };
    if(mode=="forward")measure([&]{forwardOuter(ptr(input),ptr(output),p,f);forwardInner(ptr(output),ptr(output),p,f);});
    else if(mode=="selected")measure([&]{forwardSelected(ptr(input),ptr(output),ptr(scratch),p,f);});
    else if(mode=="gather")measure([&]{forwardOuterGather(ptr(input),ptr(scratch),p,f);forwardInnerGather(ptr(scratch),ptr(output),p,f);});
    else if(mode=="transpose")measure([&]{transposeFast(ptr(input),ptr(input),ptr(scratch),p);});
    else if(mode=="outer")measure([&]{forwardOuter(ptr(input),ptr(scratch),p,f);});
    else measure([&]{forwardInner(ptr(scratch),ptr(output),p,f);});
    const auto checksum=mode=="transpose"?hash(ptr(input),k):mode=="outer"?hash(ptr(scratch),p.scratchBlocks()):hash(ptr(output),p.n);
    if(const auto* path=std::getenv("SPIN_FORWARD_RAW")) {
        std::ofstream raw(path);if(!raw)throw std::runtime_error("raw output");
        raw<<"trial,ms\n"<<std::setprecision(12);
        for(std::size_t i=0;i<times.size();++i)raw<<i<<','<<times[i]<<'\n';
    }
    std::sort(times.begin(),times.end());
    std::cout<<"mode,K,seed,calls,memory,median_ms,p10_ms,p90_ms,checksum\n"<<mode<<','<<k<<','<<seed<<','<<calls<<','<<policy<<','<<std::fixed<<std::setprecision(6)<<times[calls/2]<<','<<times[calls/10]<<','<<times[9*calls/10]<<','<<std::hex<<checksum<<'\n';
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
