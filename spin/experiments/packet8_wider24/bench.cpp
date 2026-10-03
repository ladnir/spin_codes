#include "Packet8Wide24.h"
#include <spin/Code.h>
#include <spin/PacketCode.h>
#include <algorithm>
#include <bit>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace candidate=spin::research::packet8wide24;
using Block=candidate::Block;
using Clock=std::chrono::steady_clock;
static Block* blocks(spin::Buffer& b){return reinterpret_cast<Block*>(b.bytes().data());}
static void fill(spin::Buffer& b,std::uint64_t s=913) {
    for(std::size_t i=0;i<b.bytes().size();i+=8) {
        s+=0x9e3779b97f4a7c15ULL;auto v=s;
        v=(v^(v>>30))*0xbf58476d1ce4e5b9ULL;v=(v^(v>>27))*0x94d049bb133111ebULL;v^=v>>31;
        std::memcpy(b.bytes().data()+i,&v,8);
    }
}
static std::array<std::uint64_t,2> dotLanes(const Block* a,const Block* b,std::size_t n) {
    std::array<std::uint64_t,2> result{};
    for(std::size_t i=0;i<n;++i) {
        std::uint64_t x[2],y[2];std::memcpy(x,a+i,16);std::memcpy(y,b+i,16);
        result[0]^=x[0]&y[0];result[1]^=x[1]&y[1];
    }
    return result;
}
static void equalScratch(const Block* a,const Block* b,const candidate::Plan& p,const char* what) {
    for(std::size_t g=0;g<p.groups;++g)
        if(std::memcmp(a+g*candidate::groupStride,b+g*candidate::groupStride,512*16))
            throw std::runtime_error(what);
}
static void checkSetup(const candidate::Plan& p) {
    if(p.groups*256!=p.k||p.n!=2*p.k||p.route.size()!=p.n/8||p.updates.size()!=p.n/64)
        throw std::runtime_error("setup geometry mismatch");
    std::vector<unsigned> seen(p.groups*64),regionSeen(p.groups);
    for(unsigned region=0;region<64;++region) {
        std::fill(regionSeen.begin(),regionSeen.end(),0);
        for(std::size_t slot=0;slot<p.groups;++slot) {
            const auto address=p.route[region*p.groups+slot];
            const auto group=address/candidate::groupStride,offset=address%candidate::groupStride;
            if(group>=p.groups||offset>=512||offset%8||++seen[64*group+offset/8]!=1||++regionSeen[group]!=1)
                throw std::runtime_error("route not a bijection with one packet per group per region");
        }
    }
    if(std::any_of(seen.begin(),seen.end(),[](auto n){return n!=1;}))
        throw std::runtime_error("route omitted a packet");
    for(const auto& u:p.updates) {
        if(!u.scalar||u.scalar>0xffffffU)throw std::runtime_error("invalid scalar update");
        const auto a=std::uint8_t(u.scalar),b=std::uint8_t(u.scalar>>8),c=std::uint8_t(u.scalar>>16);
        if(u.coefficients!=std::array<std::uint8_t,6>{a,b,c,std::uint8_t(a^b),std::uint8_t(a^c),std::uint8_t(b^c)})
            throw std::runtime_error("update coefficients differ from literal scalar");
    }
}
static void checks(const candidate::Plan& p) {
    checkSetup(p);
    spin::Buffer input((p.n+12)*16),packed(p.scratchBlocks()*16),raw(p.scratchBlocks()*16),
        expected(p.k*16),actual((p.k+12)*16),reference(p.scratchBlocks()*16),message(p.k*16),encoded(p.n*16);
    fill(input);fill(message,1949);
    auto* x=blocks(input)+4;
    candidate::transposeScalar(x,blocks(expected),blocks(raw),p);
    candidate::packScratchScalar(blocks(raw),blocks(reference),p);
    const std::vector<std::byte> saved(input.bytes().begin(),input.bytes().end());
    for(unsigned offset=0;offset<4;++offset) {
        x=blocks(input)+4+offset;std::memcpy(x,saved.data()+64,p.n*16);
        std::fill(packed.bytes().begin(),packed.bytes().end(),std::byte{0xa5});
        candidate::reverseRouteFast(x,blocks(packed),p);
        equalScratch(blocks(packed),blocks(reference),p,"packed inner differs from literal scalar");
        for(std::size_t g=0;g<p.groups;++g)
            for(std::size_t b=512*16;b<candidate::groupStride*16;++b)
                if(packed.bytes()[g*candidate::groupStride*16+b]!=std::byte{0xa5})
                    throw std::runtime_error("scratch padding changed");
        std::fill(actual.bytes().begin(),actual.bytes().end(),std::byte{0xa5});
        auto* y=blocks(actual)+4+offset;
        candidate::outerFast(blocks(packed),y,p);
        if(std::memcmp(y,blocks(expected),p.k*16))throw std::runtime_error("outer scalar mismatch");
        if(std::memcmp(x,saved.data()+64,p.n*16))throw std::runtime_error("input changed");
        for(std::size_t b=0;b<actual.bytes().size();++b)
            if((b<(4+offset)*16||b>=(4+offset+p.k)*16)&&actual.bytes()[b]!=std::byte{0xa5})
                throw std::runtime_error("output guard changed");
        candidate::transposeFast(x,x,blocks(packed),p);
        if(std::memcmp(x,blocks(expected),p.k*16)||std::memcmp(x+p.k,saved.data()+64+p.k*16,(p.n-p.k)*16))
            throw std::runtime_error("in-place transpose mismatch");
    }
    std::memcpy(input.bytes().data(),saved.data(),saved.size());x=blocks(input)+4;
    candidate::forwardScalar(blocks(message),blocks(encoded),p);
    if(dotLanes(blocks(message),blocks(expected),p.k)!=dotLanes(blocks(encoded),x,p.n))
        throw std::runtime_error("full forward/transpose adjoint mismatch");
    for(std::size_t packet=0;packet<p.n/8;++packet)
        std::memcpy(blocks(raw)+p.route[packet],x+8*packet,128);
    candidate::packScratchScalar(blocks(raw),blocks(reference),p);
    candidate::routeOnlyFast(x,blocks(packed),p);
    equalScratch(blocks(packed),blocks(reference),p,"route-only scalar mismatch");
    candidate::unpackScratchScalar(blocks(packed),blocks(reference),p);
    equalScratch(blocks(raw),blocks(reference),p,"scalar packing roundtrip mismatch");
}
int main(int argc,char** argv) {try {
    if(!spin::packet_fast_available()||!spin::capabilities().forward512)return 77;
    if(argc>6)throw std::invalid_argument("[K=65536] [seed=1] [calls=0] [mode:0 full,1 route,2 outer,3 inner+route] [phases=0]");
    const std::size_t k=argc>1?std::stoull(argv[1]):65536;
    const std::uint64_t seed=argc>2?std::stoull(argv[2]):1;
    const auto requested=argc>3?std::stoull(argv[3]):0;
    if(requested>1000000)throw std::invalid_argument("calls must be at most1000000");
    const unsigned calls=unsigned(requested),mode=argc>4?std::stoul(argv[4]):0;
    const bool phases=argc>5&&std::stoul(argv[5]);
    if(mode>3)throw std::invalid_argument("mode must be0..3");
    for(const auto invalid:std::array<std::size_t,5>{0,1,2047,2049,std::size_t(-1)}) {
        bool rejected=false;
        try {candidate::Plan bad(invalid,seed);}catch(const std::invalid_argument&){rejected=true;}
        if(!rejected)throw std::runtime_error("invalid dimension accepted");
    }
    candidate::Plan plan(k,seed);checks(plan);
    if(!calls){std::cout<<"PASS packet8wide24 K="<<k<<" seed="<<seed<<'\n';return 0;}
    spin::Buffer input(plan.n*16,spin::MemoryPolicy::Normal),scratch(plan.scratchBlocks()*16,spin::MemoryPolicy::Normal);
    fill(input);auto* x=blocks(input);auto* s=blocks(scratch);
    candidate::reverseRouteFast(x,s,plan);
    auto route=[&]{if(mode==1)candidate::routeOnlyFast(x,s,plan);else candidate::reverseRouteFast(x,s,plan);};
    auto outer=[&]{candidate::outerFast(s,x,plan);};
    for(unsigned i=0;i<5;++i){if(mode!=2)route();if(mode==0||mode==2)outer();}
    std::vector<double> total(calls),inner(calls),outerTimes(calls);
    for(unsigned i=0;i<calls;++i) {
        const auto a=Clock::now();if(mode!=2)route();auto b=a;
        if(phases)b=Clock::now();
        if(mode==0||mode==2)outer();const auto c=Clock::now();
        total[i]=std::chrono::duration<double,std::micro>(c-a).count();
        inner[i]=phases?std::chrono::duration<double,std::micro>(b-a).count():0.;
        outerTimes[i]=phases?std::chrono::duration<double,std::micro>(c-b).count():0.;
    }
    for(auto* v:{&total,&inner,&outerTimes})std::sort(v->begin(),v->end());
    const auto* bytes=reinterpret_cast<const std::byte*>((mode==1||mode==3)?s:x);
    const auto length=(mode==1||mode==3)?plan.scratchBlocks()*16:plan.k*16;
    std::uint64_t checksum=0;
    for(std::size_t i=0;i<length;i+=8){std::uint64_t v;std::memcpy(&v,bytes+i,8);checksum=(checksum^v)*0x100000001b3ULL;}
    std::cout<<"mode,K,seed,calls,phases,median_us,p10_us,p90_us,inner_us,outer_us,checksum\n"
        <<mode<<','<<k<<','<<seed<<','<<calls<<','<<phases<<','<<std::fixed<<std::setprecision(4)
        <<total[calls/2]<<','<<total[calls/10]<<','<<total[9*calls/10]<<','
        <<inner[calls/2]<<','<<outerTimes[calls/2]<<','<<std::hex<<checksum<<'\n';
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
