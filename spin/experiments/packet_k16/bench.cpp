#include <spin/Code.h>
#include <spin/PacketCode.h>
#include "../../src/packet/PacketPlan.h"
#include <algorithm>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace p=spin::detail::packet;
namespace spin::detail::packet {
void reverseRouteCached(const Block*,Block*,const Plan&);
void outerFastCached(const Block*,Block*,const Plan&);
}
using Clock=std::chrono::steady_clock;
static p::Block* blocks(spin::Buffer& b) {return reinterpret_cast<p::Block*>(b.bytes().data());}
static void fill(spin::Buffer& b) {
    std::uint64_t state=913;
    for(std::size_t i=0;i<b.bytes().size();i+=8) {
        state+=0x9e3779b97f4a7c15ULL;auto v=state;
        v=(v^(v>>30))*0xbf58476d1ce4e5b9ULL;v=(v^(v>>27))*0x94d049bb133111ebULL;v^=v>>31;
        std::memcpy(b.bytes().data()+(i^8),&v,8);
    }
}
static std::uint64_t hash(const p::Block* x,std::size_t k) {
    std::uint64_t value=0;
    for(std::size_t i=0;i<k;++i) {
        std::uint64_t w[2];std::memcpy(w,x+i,16);
        for(auto v:w)value=(value^v)*0x100000001b3ULL;
    }
    return value;
}
template<bool CachedRoute,bool CachedOuter>
static void encode(p::Block* x,p::Block* work,const p::Plan& plan) {
    if constexpr(CachedRoute)p::reverseRouteCached(x,work,plan);else p::reverseRoute(x,work,plan);
    if constexpr(CachedOuter)p::outerFastCached(work,x,plan);else p::outerFast(work,x,plan);
}
static double median(std::vector<double>& values) {
    std::sort(values.begin(),values.end());return values[values.size()/2];
}
template<bool R,bool O>
static void run(std::size_t k,std::uint64_t seed,unsigned calls,bool phases) {
    const p::Plan plan(k,seed);
    spin::Buffer input(2*k*16),scratch(plan.scratchBlocks()*16);
    fill(input);auto* x=blocks(input);auto* work=blocks(scratch);
    for(unsigned i=0;i<5;++i)encode<R,O>(x,work,plan);
    std::vector<double> total(calls),route(calls),outer(calls);
    for(unsigned i=0;i<calls;++i) {
        const auto begin=Clock::now();
        if(phases) {
            if constexpr(R)p::reverseRouteCached(x,work,plan);else p::reverseRoute(x,work,plan);
            const auto middle=Clock::now();
            if constexpr(O)p::outerFastCached(work,x,plan);else p::outerFast(work,x,plan);
            const auto end=Clock::now();
            route[i]=std::chrono::duration<double,std::micro>(middle-begin).count();
            outer[i]=std::chrono::duration<double,std::micro>(end-middle).count();
            total[i]=std::chrono::duration<double,std::micro>(end-begin).count();
        } else {
            encode<R,O>(x,work,plan);
            total[i]=std::chrono::duration<double,std::micro>(Clock::now()-begin).count();
        }
    }
    std::cout<<"K,seed,route_cached,outer_cached,calls,total_us,route_us,outer_us,checksum\n"
        <<k<<','<<seed<<','<<R<<','<<O<<','<<calls<<','<<std::fixed<<std::setprecision(4)
        <<median(total)<<','<<median(route)<<','<<median(outer)<<','<<std::hex<<hash(x,k)<<'\n';
}
template<bool R,bool O>
static void check(std::size_t k,std::uint64_t seed) {
    const p::Plan plan(k,seed);
    spin::Buffer input(2*k*16),expected(k*16),scratch(plan.scratchBlocks()*16);
    fill(input);
    p::transposeScalar(blocks(input),blocks(expected),blocks(scratch),plan);
    encode<R,O>(blocks(input),blocks(scratch),plan);
    if(std::memcmp(input.bytes().data(),expected.bytes().data(),k*16))throw std::runtime_error("store-policy mismatch");
}
int main(int argc,char** argv) {try {
    if(!spin::packet_fast_available())return 77;
    const std::string mode=argc>1?argv[1]:"check";
    if(mode=="check") {
        for(auto k:{256U,768U,4096U,65280U,65536U,65792U})for(auto seed:{1ULL,17ULL}) {
            check<false,false>(k,seed);check<true,false>(k,seed);
            check<false,true>(k,seed);check<true,true>(k,seed);
        }
        std::cout<<"PASS four store policies versus scalar, natural lengths, two seeds\n";return 0;
    }
    if(argc<3 || argc>6)throw std::invalid_argument("usage: MODE K [seed=1] [calls=1001] [phases=0]");
    const auto k=std::stoull(argv[2]),seed=argc>3?std::stoull(argv[3]):1;
    const auto calls=argc>4?std::stoul(argv[4]):1001;
    const bool phases=argc>5?std::stoul(argv[5])!=0:false;
    if(!calls || calls>1000000)throw std::invalid_argument("invalid calls");
    if(mode=="nn")run<false,false>(k,seed,calls,phases);
    else if(mode=="cn")run<true,false>(k,seed,calls,phases);
    else if(mode=="nc")run<false,true>(k,seed,calls,phases);
    else if(mode=="cc")run<true,true>(k,seed,calls,phases);
    else throw std::invalid_argument("unknown mode");
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
