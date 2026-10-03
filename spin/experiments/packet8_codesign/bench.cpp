#include "Packet8.h"
#include <spin/Code.h>
#include <spin/PacketCode.h>
#include <algorithm>
#include <bit>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace p8=spin::research::packet8;
using Block=p8::Block;
using Clock=std::chrono::steady_clock;
static Block* blocks(spin::Buffer& b){return reinterpret_cast<Block*>(b.bytes().data());}
static void fill(spin::Buffer& b,std::uint64_t s=913) {
    for(std::size_t i=0;i<b.bytes().size();i+=8) {
        s+=0x9e3779b97f4a7c15ULL;auto v=s;
        v=(v^(v>>30))*0xbf58476d1ce4e5b9ULL;v=(v^(v>>27))*0x94d049bb133111ebULL;v^=v>>31;
        std::memcpy(b.bytes().data()+i,&v,8);
    }
}
static unsigned dot(const Block* a,const Block* b,std::size_t n) {
    unsigned result=0;
    for(std::size_t i=0;i<n;++i) {
        std::uint64_t x[2],y[2];std::memcpy(x,a+i,16);std::memcpy(y,b+i,16);
        result^=std::popcount(x[0]&y[0])^std::popcount(x[1]&y[1]);
    }
    return result&1;
}
static void equalScratch(const Block* a,const Block* b,const p8::Plan& p,const char* what) {
    for(std::size_t g=0;g<p.outer.groups;++g)
        if(std::memcmp(a+g*spin::research::rs::groupStride,b+g*spin::research::rs::groupStride,256*16))
            throw std::runtime_error(what);
}
static void checks(const p8::Plan& p) {
    spin::Buffer input((p.n()+12)*16),packed(p.scratchBlocks()*16),raw(p.scratchBlocks()*16),
        expected(p.k()*16),actual((p.k()+12)*16),reference(p.scratchBlocks()*16),
        message(p.k()*16),encoded(p.n()*16);
    fill(input);fill(message,1949);
    auto* x=blocks(input)+4;
    p8::transposeScalar(x,blocks(expected),blocks(raw),p);
    p8::packScratchScalar(blocks(raw),blocks(reference),p);
    const std::vector<std::byte> saved(input.bytes().begin(),input.bytes().end());
    for(unsigned offset=0;offset<4;++offset) {
        x=blocks(input)+4+offset;std::memcpy(x,saved.data()+64,p.n()*16);
        std::fill(packed.bytes().begin(),packed.bytes().end(),std::byte{0xa5});
        p8::reverseRouteFast(x,blocks(packed),p);
        equalScratch(blocks(packed),blocks(reference),p,"packed inner differs from literal scalar");
        for(auto run:{p8::reverseRoutePrefetch8,p8::reverseRoutePrefetch16,p8::reverseRoutePrefetch32}) {
            run(x,blocks(packed),p);
            equalScratch(blocks(packed),blocks(reference),p,"prefetched inner differs from literal scalar");
        }
        for(std::size_t g=0;g<p.outer.groups;++g)
            for(std::size_t b=256*16;b<spin::research::rs::groupStride*16;++b)
                if(packed.bytes()[g*spin::research::rs::groupStride*16+b]!=std::byte{0xa5})
                    throw std::runtime_error("scratch padding changed");
        std::fill(actual.bytes().begin(),actual.bytes().end(),std::byte{0xa5});
        auto* y=blocks(actual)+4+offset;
        p8::outerFast(blocks(packed),y,p);
        if(std::memcmp(y,blocks(expected),p.k()*16))throw std::runtime_error("outer scalar mismatch");
        p8::outerFastShared(blocks(packed),y,p);
        if(std::memcmp(y,blocks(expected),p.k()*16))throw std::runtime_error("shared outer scalar mismatch");
        p8::outerFastHalf(blocks(packed),y,p);
        if(std::memcmp(y,blocks(expected),p.k()*16))throw std::runtime_error("half-payload outer scalar mismatch");
        if(std::memcmp(x,saved.data()+64,p.n()*16))throw std::runtime_error("input changed");
        for(std::size_t b=0;b<actual.bytes().size();++b)
            if((b<(4+offset)*16||b>=(4+offset+p.k())*16)&&actual.bytes()[b]!=std::byte{0xa5})
                throw std::runtime_error("output guard changed");
        p8::transposeFast(x,x,blocks(packed),p);
        if(std::memcmp(x,blocks(expected),p.k()*16)||
            std::memcmp(x+p.k(),saved.data()+64+p.k()*16,(p.n()-p.k())*16))
            throw std::runtime_error("in-place transpose mismatch");
        std::memcpy(x,saved.data()+64,p.n()*16);
        p8::reverseRouteFast(x,blocks(packed),p);p8::outerFastShared(blocks(packed),x,p);
        if(std::memcmp(x,blocks(expected),p.k()*16)||
            std::memcmp(x+p.k(),saved.data()+64+p.k()*16,(p.n()-p.k())*16))
            throw std::runtime_error("in-place shared transpose mismatch");
        std::memcpy(x,saved.data()+64,p.n()*16);
        p8::reverseRouteFast(x,blocks(packed),p);p8::outerFastHalf(blocks(packed),x,p);
        if(std::memcmp(x,blocks(expected),p.k()*16)||
            std::memcmp(x+p.k(),saved.data()+64+p.k()*16,(p.n()-p.k())*16))
            throw std::runtime_error("in-place half-payload transpose mismatch");
    }
    std::memcpy(input.bytes().data(),saved.data(),saved.size());x=blocks(input)+4;
    p8::forwardScalar(blocks(message),blocks(encoded),p);
    if(dot(blocks(message),blocks(expected),p.k())!=dot(blocks(encoded),x,p.n()))
        throw std::runtime_error("full forward/transpose adjoint mismatch");
    // Independently compare route-only paired stores against a literal route
    // followed by scalar bit packing. No old four-bit route is involved.
    for(std::size_t packet=0;packet<p.n()/8;++packet)
        std::memcpy(blocks(raw)+p.route[packet],x+8*packet,128);
    p8::packScratchScalar(blocks(raw),blocks(reference),p);
    p8::routeOnlyFast(x,blocks(packed),p);
    equalScratch(blocks(packed),blocks(reference),p,"route-only paired-store mismatch");
    for(auto run:{p8::routeOnlyPrefetch8,p8::routeOnlyPrefetch16,p8::routeOnlyPrefetch32}) {
        run(x,blocks(packed),p);
        equalScratch(blocks(packed),blocks(reference),p,"prefetched route-only mismatch");
    }
    p8::unpackScratchScalar(blocks(packed),blocks(reference),p);
    equalScratch(blocks(raw),blocks(reference),p,"scalar packing roundtrip mismatch");
}
int main(int argc,char** argv) {try {
    if(!spin::packet_fast_available()||!spin::capabilities().forward512)return 77;
    if(argc>6)throw std::invalid_argument("[K=65536] [seed=1] [calls=0] [mode:0..4 baseline,5..7 shared+prefetch8/16/32,8..10 route+prefetch8/16/32,11..12 half-payload full/outer] [phases=0]");
    const std::size_t k=argc>1?std::stoull(argv[1]):65536;
    const std::uint64_t seed=argc>2?std::stoull(argv[2]):1;
    const auto requestedCalls=argc>3?std::stoull(argv[3]):0;
    if(requestedCalls>1000000)throw std::invalid_argument("calls must be at most1000000");
    const unsigned calls=unsigned(requestedCalls),mode=argc>4?std::stoul(argv[4]):0;
    const bool phases=argc>5&&std::stoul(argv[5]);
    if(mode>12)throw std::invalid_argument("mode must be0..12");
    const bool outerOnly=mode==2||mode==4||mode==12;
    const bool routeOnly=mode==1||(mode>=8&&mode<=10);
    p8::Plan plan(k,seed);checks(plan);
    if(!calls){std::cout<<"PASS packet8 K="<<k<<" seed="<<seed<<'\n';return 0;}
    spin::Buffer input(plan.n()*16),scratch(plan.scratchBlocks()*16);
    fill(input);auto* x=blocks(input);auto* s=blocks(scratch);
    p8::reverseRouteFast(x,s,plan);
    auto route=[&]{switch(mode) {
        case 1:p8::routeOnlyFast(x,s,plan);break;
        case 5:p8::reverseRoutePrefetch8(x,s,plan);break;
        case 6:p8::reverseRoutePrefetch16(x,s,plan);break;
        case 7:p8::reverseRoutePrefetch32(x,s,plan);break;
        case 8:p8::routeOnlyPrefetch8(x,s,plan);break;
        case 9:p8::routeOnlyPrefetch16(x,s,plan);break;
        case 10:p8::routeOnlyPrefetch32(x,s,plan);break;
        default:p8::reverseRouteFast(x,s,plan);
    }};
    auto outerStage=[&]{if(mode>=11)p8::outerFastHalf(s,x,plan);
        else if(mode>=3)p8::outerFastShared(s,x,plan);else p8::outerFast(s,x,plan);};
    for(unsigned i=0;i<5;++i){if(!outerOnly)route();if(!routeOnly)outerStage();}
    std::vector<double> total(calls),inner(calls),outer(calls);
    for(unsigned i=0;i<calls;++i) {
        auto a=Clock::now();if(!outerOnly)route();auto b=a;
        if(phases)b=Clock::now();
        if(!routeOnly)outerStage();auto c=Clock::now();
        total[i]=std::chrono::duration<double,std::micro>(c-a).count();
        inner[i]=phases?std::chrono::duration<double,std::micro>(b-a).count():0.;
        outer[i]=phases?std::chrono::duration<double,std::micro>(c-b).count():0.;
    }
    for(auto* v:{&total,&inner,&outer})std::sort(v->begin(),v->end());
    const auto* bytes=reinterpret_cast<const std::byte*>(routeOnly?s:x);
    const auto length=routeOnly?plan.scratchBlocks()*16:plan.k()*16;
    std::uint64_t checksum=0;
    for(std::size_t i=0;i<length;i+=8){std::uint64_t v;std::memcpy(&v,bytes+i,8);checksum=(checksum^v)*0x100000001b3ULL;}
    std::cout<<"mode,K,seed,calls,phases,median_us,p10_us,p90_us,inner_us,outer_us,checksum\n"
        <<mode<<','<<k<<','<<seed<<','<<calls<<','<<phases<<','<<std::fixed<<std::setprecision(4)
        <<total[calls/2]<<','<<total[calls/10]<<','<<total[9*calls/10]<<','
        <<inner[calls/2]<<','<<outer[calls/2]<<','<<std::hex<<checksum<<'\n';
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
