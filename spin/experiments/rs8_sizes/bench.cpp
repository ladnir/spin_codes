#include <cstdlib>
#include "Kernels.h"
#include "../k20_two_ms/rs8/Rs8.h"
#include "../../src/certified22/productionMap.h"
#include <spin/Code.h>
#include <spin/PacketCode.h>
#include <algorithm>
#include <bit>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace code=spin::detail::certified22;
namespace p15=spin::detail::paired15;
namespace local=spin::research::rs8sizes;
namespace rs8=spin::research::twoms::rs8;
using Block=code::Block;
using Clock=std::chrono::steady_clock;

static Block* blocks(spin::Buffer& b) { return reinterpret_cast<Block*>(b.bytes().data()); }
static void fill(spin::Buffer& b,std::uint64_t state=913) {
    for(std::size_t i=0;i<b.bytes().size();i+=8) {
        state+=0x9e3779b97f4a7c15ULL;auto v=state;
        v=(v^(v>>30))*0xbf58476d1ce4e5b9ULL;
        v=(v^(v>>27))*0x94d049bb133111ebULL;v^=v>>31;
        std::memcpy(b.bytes().data()+i,&v,8);
    }
}
static void paint(spin::Buffer& b,std::byte value) {
    std::fill(b.bytes().begin(),b.bytes().end(),value);
}
static std::uint64_t hash(const void* pointer,std::size_t bytes) {
    const auto* p=static_cast<const std::byte*>(pointer);std::uint64_t h=0;
    for(std::size_t i=0;i<bytes;i+=8) {
        std::uint64_t v;std::memcpy(&v,p+i,8);h=(h^v)*0x100000001b3ULL;
    }
    return h;
}
static void bytesEqual(const void* a,const void* b,std::size_t bytes,const char* failure) {
    if(std::memcmp(a,b,bytes)) throw std::runtime_error(failure);
}
static void guards(spin::Buffer& b,std::size_t offset,std::size_t count,std::byte value) {
    const auto guard=[=](std::byte x){return x==value;};
    if(!std::all_of(b.bytes().begin(),b.bytes().begin()+16*offset,guard) ||
       !std::all_of(b.bytes().begin()+16*(offset+count),b.bytes().end(),guard))
        throw std::runtime_error("buffer guard changed");
}

// Literal coordinates from production Paired15Scalar.cpp. This oracle uses
// literal binary matrices, independent of the fast slot layout/GFNI updates.
constexpr std::uint16_t columns15[64]={
    0x1,0x3,0x5,0x87,0x9,0xb,0x400d,0x408f,0x11,0x113,0x815,0x997,0x419,0x51b,0x4c1d,0x4d9f,
    0x21,0x1023,0x25,0x10a7,0x829,0x182b,0x482d,0x58af,0x2031,0x3133,0x2835,0x39b7,0x2c39,0x3d3b,0x643d,0x75bf,
    0x41,0x443,0x2045,0x24c7,0x1049,0x144b,0x704d,0x74cf,0x4051,0x4553,0x6855,0x6dd7,0x5459,0x515b,0x3c5d,0x39df,
    0x261,0x1663,0x2265,0x36e7,0x1a69,0xe6b,0x7a6d,0x6eef,0x6271,0x7773,0x4a75,0x5ff7,0x7e79,0x6b7b,0x167d,0x3ff};

template<bool S15,class Plan>
static void forwardInnerScalar(const Block* input,Block* output,const Plan& p) {
    constexpr unsigned stateSize=S15?16:22;
    alignas(16) __m128i state[stateSize]{},next[stateSize],feedback[stateSize];
    for(std::size_t epoch=0;epoch<p.n/64;++epoch) {
        for(auto& value:feedback) value=_mm_setzero_si128();
        for(unsigned j=0;j<64;++j) {
            const auto i=64*epoch+j,routed=p.route[i/4]+(i&3);
            const auto raw=input[routed].mData;auto value=raw;
            const auto column=[&] {
                if constexpr(S15) return std::uint32_t(columns15[j]);
                else return code::map::oldColumns[j];
            }();
            for(auto mask=column;mask;mask&=mask-1) {
                const auto c=std::countr_zero(mask);
                value=_mm_xor_si128(value,state[c]);
                feedback[c]=_mm_xor_si128(feedback[c],raw);
            }
            output[i]=Block(value);
        }
        if(epoch+1==p.n/64) break;
        for(unsigned j=0;j<stateSize;++j) {
            auto value=feedback[j];
            for(unsigned c=0;c<stateSize;++c)
                if((p.reverseMatrices[epoch][c]>>j)&1U)
                    value=_mm_xor_si128(value,state[c]);
            next[j]=value;
        }
        std::memcpy(state,next,sizeof(state));
    }
}

template<class Route,class Outer,class TransposeScalar,class ForwardScalar>
static void verify(const char* mode,const code::Plan& p,std::uint64_t seed,
                   Route&& route,Outer&& outer,TransposeScalar&& scalar,ForwardScalar&& forward) {
    const auto scratchCount=p.scratchBlocks();
    spin::Buffer original(p.n*16),input((p.n+12)*16),output((p.k+12)*16),
        work((scratchCount+8)*16),referenceWork((scratchCount+8)*16),expected(p.k*16),
        message(p.k*16),encoded(p.n*16),forwardWork(scratchCount*16);
    fill(original,seed+913);fill(message,seed+0xdef931);
    paint(referenceWork,std::byte{0x3c});
    scalar(blocks(original),blocks(expected),blocks(referenceWork)+4);
    guards(referenceWork,4,scratchCount,std::byte{0x3c});
    forward(blocks(message),blocks(encoded),blocks(forwardWork));
    auto lhs=_mm_setzero_si128(),rhs=_mm_setzero_si128();
    for(std::size_t i=0;i<p.n;++i)
        lhs=_mm_xor_si128(lhs,_mm_and_si128(blocks(original)[i].mData,blocks(encoded)[i].mData));
    for(std::size_t i=0;i<p.k;++i)
        rhs=_mm_xor_si128(rhs,_mm_and_si128(blocks(message)[i].mData,blocks(expected)[i].mData));
    if(!_mm_testz_si128(_mm_xor_si128(lhs,rhs),_mm_set1_epi32(-1)))
        throw std::runtime_error("complete forward/transpose binary adjoint differs");

    for(unsigned offset=0;offset<4;++offset) {
        paint(input,std::byte{0x91});paint(output,std::byte{0xa5});paint(work,std::byte{0x3c});
        auto* x=blocks(input)+4+offset;auto* y=blocks(output)+4+offset;auto* w=blocks(work)+4;
        std::memcpy(x,blocks(original),p.n*16);
        const auto saved=std::vector<std::byte>(input.bytes().begin(),input.bytes().end());
        route(x,w);
        bytesEqual(work.bytes().data(),referenceWork.bytes().data(),work.bytes().size(),
              "route, scratch padding, or scratch guard differs from independent scalar");
        outer(w,y);
        bytesEqual(y,blocks(expected),p.k*16,"fast outer/full transpose differs from scalar");
        bytesEqual(work.bytes().data(),referenceWork.bytes().data(),work.bytes().size(),"outer modified scratch");
        bytesEqual(input.bytes().data(),saved.data(),saved.size(),"out-of-place input changed");
        guards(output,4+offset,p.k,std::byte{0xa5});

        // Poison only the four padding blocks; no outer may consume them.
        for(std::size_t g=0;g<p.groups;++g)
            std::memset(w+g*code::groupStride+256,0xe7,4*sizeof(Block));
        outer(w,y);
        bytesEqual(y,blocks(expected),p.k*16,"outer depends on scratch padding");
        guards(output,4+offset,p.k,std::byte{0xa5});

        paint(work,std::byte{0x3c});route(x,w);outer(w,x);
        bytesEqual(x,blocks(expected),p.k*16,"in-place result differs");
        bytesEqual(input.bytes().data(),saved.data(),16*(4+offset),"in-place prefix guard changed");
        const auto tail=16*(4+offset+p.k);
        bytesEqual(input.bytes().data()+tail,saved.data()+tail,saved.size()-tail,
              "in-place input suffix or end guard changed");
        bytesEqual(work.bytes().data(),referenceWork.bytes().data(),work.bytes().size(),
              "in-place scratch or scratch guard changed");
    }
    std::cout<<"PASS mode="<<mode<<" K="<<p.k<<" seed="<<seed
             <<" scalar-route/full-adjoint/offsets0-3/guards/in-place/padding hash="
             <<std::hex<<hash(blocks(expected),p.k*16)<<std::dec<<'\n';
}

static void verifyRs8Outer(const code::Plan& p,const rs8::Tables32& tables) {
    if(!rs8::selfTest()) throw std::runtime_error("RS8 flattened field arithmetic failed");
    spin::Buffer message(p.k*16),word(p.scratchBlocks()*16),
        forward(p.scratchBlocks()*16),transposed(p.k*16);
    fill(message,0x719def);fill(word,0xb706c5a81924ULL);
    rs8::outerForwardScalar(blocks(message),blocks(forward),p,tables);
    rs8::outerScalar(blocks(word),blocks(transposed),p,tables);
    auto lhs=_mm_setzero_si128(),rhs=_mm_setzero_si128();
    for(std::size_t g=0;g<p.groups;++g) for(unsigned j=0;j<256;++j) {
        const auto i=g*code::groupStride+j;
        lhs=_mm_xor_si128(lhs,_mm_and_si128(blocks(word)[i].mData,blocks(forward)[i].mData));
    }
    for(std::size_t i=0;i<p.k;++i)
        rhs=_mm_xor_si128(rhs,_mm_and_si128(blocks(message)[i].mData,blocks(transposed)[i].mData));
    if(!_mm_testz_si128(_mm_xor_si128(lhs,rhs),_mm_set1_epi32(-1)))
        throw std::runtime_error("RS8 outer binary adjoint differs");
}

template<unsigned Phase,class Route,class Outer>
static void measure(const char* mode,const code::Plan& p,std::uint64_t seed,unsigned calls,
                    Route&& route,Outer&& outer) {
    // SPIN_BENCH_MEMORY=huge|normal overrides the size-based default (research harness only).
    const char* memoryEnv=std::getenv("SPIN_BENCH_MEMORY");
    const bool huge=memoryEnv?std::string(memoryEnv)=="huge":p.k>65536;
    const auto policy=huge?spin::MemoryPolicy::PreferHugePages:spin::MemoryPolicy::Normal;
    spin::Buffer input((p.n+12)*16,policy),scratch(p.scratchBlocks()*16,policy);
    fill(input);paint(scratch,std::byte{0});auto* x=blocks(input)+4;auto* w=blocks(scratch);
    route(x,w);
    const auto encode=[&] {
        if constexpr(Phase!=3) route(x,w);
        if constexpr(Phase!=2) outer(w,x);
    };
    for(unsigned i=0;i<5;++i) encode();
    std::vector<double> all(calls),routes(calls),outers(calls);
    for(unsigned i=0;i<calls;++i) {
        const auto start=Clock::now();
        if constexpr(Phase==1) {
            route(x,w);const auto middle=Clock::now();outer(w,x);const auto end=Clock::now();
            all[i]=std::chrono::duration<double,std::micro>(end-start).count();
            routes[i]=std::chrono::duration<double,std::micro>(middle-start).count();
            outers[i]=std::chrono::duration<double,std::micro>(end-middle).count();
        } else {
            encode();all[i]=std::chrono::duration<double,std::micro>(Clock::now()-start).count();
        }
    }
    std::sort(all.begin(),all.end());std::sort(routes.begin(),routes.end());std::sort(outers.begin(),outers.end());
    std::cout<<"mode,K,seed,calls,phase,memory,total_us,p10_us,p90_us,route_us,outer_us,checksum\n"
        <<mode<<','<<p.k<<','<<seed<<','<<calls<<','<<Phase<<','<<(huge?"huge-preferred":"normal")<<','
        <<std::fixed<<std::setprecision(4)<<all[calls/2]<<','<<all[calls/10]<<','<<all[9*calls/10]<<','
        <<routes[calls/2]<<','<<outers[calls/2]<<','<<std::hex
        <<hash(Phase==2?static_cast<void*>(w):static_cast<void*>(x),Phase==2?scratch.bytes().size():p.k*16)
        <<std::dec<<'\n';
}

template<class Route,class Outer,class TransposeScalar,class ForwardScalar>
static void run(const char* mode,const code::Plan& p,std::uint64_t seed,unsigned calls,unsigned phase,
                Route&& route,Outer&& outer,TransposeScalar&& scalar,ForwardScalar&& forward) {
    if(!calls) { verify(mode,p,seed,route,outer,scalar,forward);return; }
    switch(phase) {
    case 0:measure<0>(mode,p,seed,calls,route,outer);break;
    case 1:measure<1>(mode,p,seed,calls,route,outer);break;
    case 2:measure<2>(mode,p,seed,calls,route,outer);break;
    case 3:measure<3>(mode,p,seed,calls,route,outer);break;
    }
}

int main(int argc,char** argv) { try {
    if(!spin::packet_fast_available()) return 77;
    if(argc<3 || argc>6)
        throw std::invalid_argument("baseline15|rs8-15|baseline22|rs8-22 K [seed41] [calls301;0=verify] [phase0:full,1:timed,2:route,3:outer]");
    const std::string mode=argv[1];const auto k=std::stoull(argv[2]);
    const auto seed=argc>3?std::stoull(argv[3]):41;
    const auto calls=argc>4?std::stoul(argv[4]):301;
    const auto phase=argc>5?std::stoul(argv[5]):0;
    if(k<2048 || phase>3 || calls>1000000) throw std::invalid_argument("bad K/phase/calls");
    if(mode!="baseline15" && mode!="rs8-15" && mode!="baseline22" && mode!="rs8-22")
        throw std::invalid_argument("unknown mode");
    code::Plan geometry(k,seed,seed);
    if(mode=="baseline15" || mode=="rs8-15") {
        p15::Plan plan(65536,seed,seed);local::resizePaired15(plan,geometry,seed);
        const auto route=[&](const Block*x,Block*w){local::reverseRoute15(x,w,plan);};
        if(mode=="baseline15") {
            run(mode.c_str(),geometry,seed,unsigned(calls),unsigned(phase),route,
                [&](const Block*w,Block*y){local::outer15(w,y,plan);},
                [&](const Block*x,Block*y,Block*w){p15::transposeScalar(x,y,w,plan);},
                [&](const Block*x,Block*y,Block*w){p15::forwardScalar(x,y,w,plan,1);});
        } else {
            rs8::Tables32 tables;rs8::setupField32(geometry,seed,tables);
            if(!calls) verifyRs8Outer(geometry,tables);
            run(mode.c_str(),geometry,seed,unsigned(calls),unsigned(phase),route,
                [&](const Block*w,Block*y){rs8::outerFast(w,y,geometry,tables);},
                [&](const Block*x,Block*y,Block*w){
                    p15::transposeScalar(x,y,w,plan);rs8::outerScalar(w,y,geometry,tables);
                },
                [&](const Block*x,Block*y,Block*w){
                    rs8::outerForwardScalar(x,w,geometry,tables);forwardInnerScalar<true>(w,y,plan);
                });
        }
    } else {
        const auto route=[&](const Block*x,Block*w){code::reverseRouteFast(x,w,geometry);};
        if(mode=="baseline22") {
            run(mode.c_str(),geometry,seed,unsigned(calls),unsigned(phase),route,
                [&](const Block*w,Block*y){
                    if(k<=262144) code::outerFastCached(w,y,geometry);
                    else code::outerFastStreaming(w,y,geometry);
                },
                [&](const Block*x,Block*y,Block*w){code::transposeScalar(x,y,w,geometry);},
                [&](const Block*x,Block*y,Block*w){code::forwardScalar(x,y,w,geometry,1);});
        } else {
            rs8::Tables32 tables;rs8::setupField32(geometry,seed,tables);
            if(!calls) verifyRs8Outer(geometry,tables);
            run(mode.c_str(),geometry,seed,unsigned(calls),unsigned(phase),route,
                [&](const Block*w,Block*y){rs8::outerFast(w,y,geometry,tables);},
                [&](const Block*x,Block*y,Block*w){
                    code::transposeScalar(x,y,w,geometry);rs8::outerScalar(w,y,geometry,tables);
                },
                [&](const Block*x,Block*y,Block*w){
                    rs8::outerForwardScalar(x,w,geometry,tables);forwardInnerScalar<false>(w,y,geometry);
                });
        }
    }
    return 0;
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;} }
