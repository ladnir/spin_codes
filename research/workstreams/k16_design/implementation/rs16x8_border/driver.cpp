#include "RsBorder.h"
#include "RsBorderCandidates.h"
#include "../rs16x8/RsWideOuterVariants.h"
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

namespace rs=spin::research::rsborder;
namespace outer=spin::research::rswide;
using Block=rs::Block;
static Block* blocks(spin::Buffer& b){return reinterpret_cast<Block*>(b.bytes().data());}
static void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
static void fill(spin::Buffer& b,std::uint64_t state=913) {
    for(std::size_t i=0;i<b.bytes().size();i+=8) {
        state+=0x9e3779b97f4a7c15ULL;auto v=state;
        v=(v^(v>>30))*0xbf58476d1ce4e5b9ULL;v=(v^(v>>27))*0x94d049bb133111ebULL;v^=v>>31;
        std::memcpy(b.bytes().data()+(i^8),&v,8);
    }
}
static unsigned dot(const Block* a,const Block* b,std::size_t n) {
    unsigned parity=0;
    for(std::size_t i=0;i<n;++i) {
        std::uint64_t x[2],y[2];std::memcpy(x,a+i,16);std::memcpy(y,b+i,16);
        parity^=std::popcount(x[0]&y[0])^std::popcount(x[1]&y[1]);
    }
    return parity&1;
}
static std::uint64_t checksum(const Block* a,std::size_t n) {
    std::uint64_t value=0;
    for(std::size_t i=0;i<n;++i) {
        std::uint64_t x[2];std::memcpy(x,a+i,16);
        value=(value^x[0])*0x100000001b3ULL;value=(value^x[1])*0x100000001b3ULL;
    }
    return value;
}
template<class Word,std::size_t Width> static bool fullRank(std::array<Word,Width> rows) {
    for(unsigned c=0;c<Width;++c) {
        unsigned p=c;while(p<Width && !(rows[p]&(Word{1}<<c)))++p;
        if(p==Width)return false;std::swap(rows[p],rows[c]);
        for(unsigned j=c+1;j<Width;++j)if(rows[j]&(Word{1}<<c))rows[j]^=rows[c];
    }
    return true;
}
static void check(std::size_t k,std::uint64_t seed,unsigned stateBits) {
    const rs::Plan plan(k,seed,seed+17,stateBits);const auto& base=plan.outer;const auto n=base.n;
    for(auto rows:plan.reverseMatrices) {
        for(unsigned c=0;c<stateBits;++c) {
            unsigned p=c;while(p<stateBits && !(rows[p]&(1U<<c)))++p;
            require(p<stateBits,"singular extended state matrix");
            std::swap(rows[p],rows[c]);
            for(unsigned j=c+1;j<stateBits;++j)if(rows[j]&(1U<<c))rows[j]^=rows[c];
        }
    }
    std::vector<bool> seen(128*base.groups);
    require(base.route.size()==n/4,"route size");
    for(unsigned region=0;region<128;++region) {
        std::vector<bool> groups(base.groups);
        for(std::size_t slot=0;slot<base.groups;++slot) {
            const auto routeBase=base.route[region*base.groups+slot];
            const auto group=routeBase/outer::groupStride,offset=routeBase%outer::groupStride;
            require(group<base.groups && offset<512 && offset%4==0,"route address");
            require(!groups[group],"region repeats group");groups[group]=true;
            require(!seen[128*group+offset/4],"repeated packet");seen[128*group+offset/4]=true;
        }
    }
    require(std::all_of(seen.begin(),seen.end(),[](bool v){return v;}),"missing packet");
    require(base.reverseMatrices.size()==n/64 && base.denseUpdates.size()==n/64,"inner count");
    require(base.outerMatrices.size()==16*base.groups,"outer count");
    for(const auto& matrix:base.reverseMatrices)require(fullRank(matrix),"singular inner matrix");
    for(const auto& matrix:base.outerMatrices)require(fullRank(matrix),"singular outer matrix");
    spin::Buffer input((n+8)*16),output((k+8)*16),expected(k*16),scratch(base.scratchBlocks()*16),message(k*16),encoded(n*16);
    fill(input);fill(message,2947);
    const std::vector<std::byte> saved(input.bytes().begin(),input.bytes().end());
    rs::transposeScalar(blocks(input)+4,blocks(expected),blocks(scratch),plan);
    for(unsigned offset=0;offset<4;++offset) {
        std::fill(output.bytes().begin(),output.bytes().end(),std::byte{0xa5});
        std::memmove(blocks(input)+4+offset,saved.data()+64,n*16);
        auto* y=blocks(output)+4+offset;
        rs::transposeFast(blocks(input)+4+offset,y,blocks(scratch),plan);
        require(std::memcmp(y,blocks(expected),k*16)==0,"fast versus scalar mismatch");
        rs::reverseRouteFused(blocks(input)+4+offset,blocks(scratch),plan);
        outer::outerFastFusedPairNt(blocks(scratch),y,base);
        require(std::memcmp(y,blocks(expected),k*16)==0,"NT outer versus scalar mismatch");
        require(std::all_of(output.bytes().begin(),output.bytes().begin()+(4+offset)*16,[](std::byte v){return v==std::byte{0xa5};}),"output prefix guard");
        require(std::all_of(output.bytes().begin()+(4+offset+k)*16,output.bytes().end(),[](std::byte v){return v==std::byte{0xa5};}),"output suffix guard");
        require(std::memcmp(blocks(input)+4+offset,saved.data()+64,n*16)==0,"input was changed");
    }
    std::memcpy(input.bytes().data(),saved.data(),saved.size());
    rs::reverseRoute(blocks(input)+4,blocks(scratch),plan);
    for(auto variant:{outer::outerFastPair,outer::outerFastNoInline,outer::outerFastPairNoInline,outer::outerFastPlaneNoInline,outer::outerFastFusedPair}) {
        variant(blocks(scratch),blocks(output)+4,base);
        require(std::memcmp(blocks(output)+4,blocks(expected),k*16)==0,"outer variant mismatch");
    }
    rs::transposeFastFused(blocks(input)+4,blocks(output)+4,blocks(scratch),plan);
    require(std::memcmp(blocks(output)+4,blocks(expected),k*16)==0,"fused refresh mismatch");
    rs::transposeFastTernary(blocks(input)+4,blocks(output)+4,blocks(scratch),plan);
    require(std::memcmp(blocks(output)+4,blocks(expected),k*16)==0,"ternary refresh mismatch");
    rs::transposeFastEarly(blocks(input)+4,blocks(output)+4,blocks(scratch),plan);
    require(std::memcmp(blocks(output)+4,blocks(expected),k*16)==0,"early refresh mismatch");
    auto* x=blocks(input)+4;rs::transposeFast(x,x,blocks(scratch),plan);
    require(std::memcmp(x,blocks(expected),k*16)==0,"inplace mismatch");
    require(std::memcmp(x+k,saved.data()+64+k*16,k*16)==0,"inplace suffix changed");
    require(std::memcmp(blocks(input),saved.data(),64)==0,"inplace prefix guard");
    require(std::memcmp(x+n,saved.data()+(n+4)*16,64)==0,"inplace suffix guard");
    std::memcpy(input.bytes().data(),saved.data(),saved.size());
    rs::reverseRouteFused(x,blocks(scratch),plan);outer::outerFastFusedPairNt(blocks(scratch),x,base);
    require(std::memcmp(x,blocks(expected),k*16)==0,"NT inplace mismatch");
    require(std::memcmp(x+k,saved.data()+64+k*16,k*16)==0,"NT inplace suffix changed");
    std::memcpy(input.bytes().data(),saved.data(),saved.size());
    rs::forwardScalar(blocks(message),blocks(encoded),plan);
    require(dot(blocks(encoded),x,n)==dot(blocks(message),blocks(expected),k),"binary adjoint identity");
    std::cout<<"PASS rs16x8-bordered s="<<stateBits<<" K="<<k<<" seed="<<seed<<" checksum="<<std::hex<<checksum(blocks(expected),k)<<std::dec<<'\n';
}
template<class Call> static void timeCalls(const std::string& mode,std::size_t k,std::uint64_t seed,std::size_t count,Call&& call,const Block* output) {
    for(unsigned i=0;i<5;++i)call();
    std::vector<double> times(count);
    for(auto& time:times) {
        const auto begin=std::chrono::steady_clock::now();call();
        time=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-begin).count();
    }
    std::sort(times.begin(),times.end());
    std::cout<<mode<<','<<k<<','<<seed<<','<<count<<','<<std::fixed<<std::setprecision(6)<<times[count/2]<<','<<times[count/10]<<','<<times[9*count/10]<<','<<std::hex<<checksum(output,k)<<std::dec<<'\n';
}
int main(int argc,char** argv){try {
    if(!spin::packet_fast_available() || !spin::capabilities().forward512){std::cerr<<"AVX512/DQ/VBMI/GFNI required\n";return 77;}
    if(argc>=2 && std::string(argv[1])=="check") {
        if(argc==2)for(auto k:{4096ULL,12288ULL})for(auto seed:{1ULL,7ULL})
            for(unsigned bits=16;bits<=20;++bits)check(k,seed,bits);
        else {
            require(argc<=5,"usage: check [K [seed [stateBits]]]");
            if(argc==5)check(std::stoull(argv[2]),std::stoull(argv[3]),unsigned(std::stoul(argv[4])));
            else for(unsigned bits=16;bits<=20;++bits)check(std::stoull(argv[2]),argc>3?std::stoull(argv[3]):1,bits);
        }
        return 0;
    }
    require(argc>=4 && argc<=7,"usage: {run|phases} K stateBits [seed=1] [calls=301] [normal|huge]");
    const std::string mode=argv[1];require(mode=="run" || mode=="phases" || mode=="outer-variants" || mode=="inner-variants" || mode=="final-variants" || mode=="joint-variants" || mode=="fused-nt" || mode=="fused-pair" || mode=="fused-nt-phases","unknown mode");
    const auto k=std::stoull(argv[2]),seed=argc>4?std::stoull(argv[4]):1,calls=argc>5?std::stoull(argv[5]):301;
    const auto bits=unsigned(std::stoul(argv[3]));require(calls && calls<=1000000,"invalid call count");
    const std::string memory=argc>6?argv[6]:"normal";require(memory=="normal" || memory=="huge","invalid memory policy");
    const auto policy=memory=="huge"?spin::MemoryPolicy::PreferHugePages:spin::MemoryPolicy::Normal;
    const rs::Plan plan(k,seed,seed,bits);
    spin::Buffer input(2*k*16,policy),scratch(plan.outer.scratchBlocks()*16,policy);fill(input);
    auto* x=blocks(input);auto* work=blocks(scratch);
    const auto label=std::string("rs16x8-byte-bordered-s")+std::to_string(bits)+"-"+memory;
    std::cout<<"mode,K,seed,calls,median_ms,p10_ms,p90_ms,checksum\n";
    if(mode=="fused-nt" || mode=="fused-nt-phases") {
        timeCalls(label+"-fused-pair-nt",k,seed,calls,[&]{rs::reverseRouteFused(x,work,plan);outer::outerFastFusedPairNt(work,x,plan.outer);},x);
        if(mode=="fused-nt-phases") {
            timeCalls(label+"-fused-route",k,seed,calls,[&]{rs::reverseRouteFused(x,work,plan);},work);
            timeCalls(label+"-pair-nt-outer",k,seed,calls,[&]{outer::outerFastFusedPairNt(work,x,plan.outer);},x);
        }
        return 0;
    }
    if(mode=="fused-pair") {
        timeCalls(label+"-fused-pair",k,seed,calls,[&]{rs::reverseRouteFused(x,work,plan);outer::outerFastFusedPair(work,x,plan.outer);},x);
        return 0;
    }
    timeCalls(label,k,seed,calls,[&]{rs::transposeFast(x,x,work,plan);},x);
    if(mode=="outer-variants") {
        timeCalls(label+"-noinline",k,seed,calls,[&]{rs::reverseRoute(x,work,plan);outer::outerFastNoInline(work,x,plan.outer);},x);
        timeCalls(label+"-pair-noinline",k,seed,calls,[&]{rs::reverseRoute(x,work,plan);outer::outerFastPairNoInline(work,x,plan.outer);},x);
        timeCalls(label+"-control",k,seed,calls,[&]{rs::transposeFast(x,x,work,plan);},x);
    }
    if(mode=="inner-variants") {
        timeCalls(label+"-noinline",k,seed,calls,[&]{rs::reverseRoute(x,work,plan);outer::outerFastNoInline(work,x,plan.outer);},x);
        timeCalls(label+"-fused",k,seed,calls,[&]{rs::transposeFastFused(x,x,work,plan);},x);
        timeCalls(label+"-fused-noinline",k,seed,calls,[&]{rs::reverseRouteFused(x,work,plan);outer::outerFastNoInline(work,x,plan.outer);},x);
        timeCalls(label+"-control",k,seed,calls,[&]{rs::transposeFast(x,x,work,plan);},x);
    }
    if(mode=="final-variants") {
        timeCalls(label+"-fused-noinline",k,seed,calls,[&]{rs::reverseRouteFused(x,work,plan);outer::outerFastNoInline(work,x,plan.outer);},x);
        timeCalls(label+"-fused-plane-noinline",k,seed,calls,[&]{rs::reverseRouteFused(x,work,plan);outer::outerFastPlaneNoInline(work,x,plan.outer);},x);
        timeCalls(label+"-fused-outer-pair",k,seed,calls,[&]{rs::reverseRouteFused(x,work,plan);outer::outerFastFusedPair(work,x,plan.outer);},x);
        timeCalls(label+"-ternary-noinline",k,seed,calls,[&]{rs::reverseRouteTernary(x,work,plan);outer::outerFastNoInline(work,x,plan.outer);},x);
        timeCalls(label+"-fused-noinline-control",k,seed,calls,[&]{rs::reverseRouteFused(x,work,plan);outer::outerFastNoInline(work,x,plan.outer);},x);
    }
    if(mode=="joint-variants") {
        timeCalls(label+"-fused-pair",k,seed,calls,[&]{rs::reverseRouteFused(x,work,plan);outer::outerFastFusedPair(work,x,plan.outer);},x);
        timeCalls(label+"-early-pair",k,seed,calls,[&]{rs::reverseRouteEarly(x,work,plan);outer::outerFastFusedPair(work,x,plan.outer);},x);
        timeCalls(label+"-fused-pair-nt",k,seed,calls,[&]{rs::reverseRouteFused(x,work,plan);outer::outerFastFusedPairNt(work,x,plan.outer);},x);
        timeCalls(label+"-early-pair-nt",k,seed,calls,[&]{rs::reverseRouteEarly(x,work,plan);outer::outerFastFusedPairNt(work,x,plan.outer);},x);
        timeCalls(label+"-fused-pair-control",k,seed,calls,[&]{rs::reverseRouteFused(x,work,plan);outer::outerFastFusedPair(work,x,plan.outer);},x);
    }
    if(mode=="phases") {
        timeCalls(label+"-reverse-route",k,seed,calls,[&]{rs::reverseRoute(x,work,plan);},work);
        timeCalls(label+"-outer",k,seed,calls,[&]{outer::outerFast(work,x,plan.outer);},x);
    }
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
