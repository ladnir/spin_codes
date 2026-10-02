#include "RsWide64.h"
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

namespace rs=spin::research::rswide64;
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
static void check(std::size_t k,std::uint64_t seed,bool fast=true) {
    const rs::Plan plan(k,seed,seed+17);const auto n=plan.n;
    std::vector<bool> seen(256*plan.groups);
    require(plan.route.size()==n/4,"route size");
    for(unsigned region=0;region<256;++region) {
        std::vector<bool> groups(plan.groups);
        for(std::size_t slot=0;slot<plan.groups;++slot) {
            const auto base=plan.route[region*plan.groups+slot];
            const auto group=base/rs::groupStride,offset=base%rs::groupStride;
            require(group<plan.groups && offset<1024 && offset%4==0,"route address");
            require(!groups[group],"region repeats group");groups[group]=true;
            require(!seen[256*group+offset/4],"repeated packet");seen[256*group+offset/4]=true;
        }
    }
    require(std::all_of(seen.begin(),seen.end(),[](bool v){return v;}),"missing packet");
    require(plan.reverseMatrices.size()==n/64 && plan.denseUpdates.size()==n/64,"inner count");
    require(plan.outerMatrices.size()==16*plan.groups,"outer count");
    for(const auto& matrix:plan.reverseMatrices)require(fullRank(matrix),"singular inner matrix");
    for(const auto& matrix:plan.outerMatrices)require(fullRank(matrix),"singular outer matrix");
    spin::Buffer input((n+8)*16),output((k+8)*16),expected(k*16),scratch(plan.scratchBlocks()*16),message(k*16),encoded(n*16);
    fill(input);fill(message,2947);
    const std::vector<std::byte> saved(input.bytes().begin(),input.bytes().end());
    rs::transposeScalar(blocks(input)+4,blocks(expected),blocks(scratch),plan);
    auto* x=blocks(input)+4;
    if(fast) {
    for(unsigned offset=0;offset<4;++offset) {
        std::fill(output.bytes().begin(),output.bytes().end(),std::byte{0xa5});
        std::memmove(blocks(input)+4+offset,saved.data()+64,n*16);
        auto* y=blocks(output)+4+offset;
        rs::transposeFast(blocks(input)+4+offset,y,blocks(scratch),plan);
        require(std::memcmp(y,blocks(expected),k*16)==0,"fast versus scalar mismatch");
        require(std::all_of(output.bytes().begin(),output.bytes().begin()+(4+offset)*16,[](std::byte v){return v==std::byte{0xa5};}),"output prefix guard");
        require(std::all_of(output.bytes().begin()+(4+offset+k)*16,output.bytes().end(),[](std::byte v){return v==std::byte{0xa5};}),"output suffix guard");
        require(std::memcmp(blocks(input)+4+offset,saved.data()+64,n*16)==0,"input was changed");
    }
    std::memcpy(input.bytes().data(),saved.data(),saved.size());
    rs::transposeFast(x,x,blocks(scratch),plan);
    require(std::memcmp(x,blocks(expected),k*16)==0,"inplace mismatch");
    require(std::memcmp(x+k,saved.data()+64+k*16,k*16)==0,"inplace suffix changed");
    require(std::memcmp(blocks(input),saved.data(),64)==0,"inplace prefix guard");
    require(std::memcmp(x+n,saved.data()+(n+4)*16,64)==0,"inplace suffix guard");
    }
    std::memcpy(input.bytes().data(),saved.data(),saved.size());
    rs::forwardScalar(blocks(message),blocks(encoded),plan);
    require(dot(blocks(encoded),x,n)==dot(blocks(message),blocks(expected),k),"binary adjoint identity");
    std::cout<<"PASS rs16x16 fast="<<fast<<" K="<<k<<" seed="<<seed<<" checksum="<<std::hex<<checksum(blocks(expected),k)<<std::dec<<'\n';
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
    if(argc>=2 && std::string(argv[1])=="scalar-check") {
        require(argc<=4,"usage: scalar-check [K [seed]]");
        check(argc>2?std::stoull(argv[2]):8192,argc>3?std::stoull(argv[3]):1,false);return 0;
    }
    if(!spin::packet_fast_available() || !spin::capabilities().forward512){std::cerr<<"AVX512/DQ/VBMI/GFNI required\n";return 77;}
    if(argc>=2 && std::string(argv[1])=="check") {
        if(argc==2)for(auto k:{8192ULL,24576ULL,65536ULL})for(auto seed:{1ULL,7ULL})check(k,seed);
        else {require(argc<=4,"usage: check [K [seed]]");check(std::stoull(argv[2]),argc>3?std::stoull(argv[3]):1);}
        return 0;
    }
    require(argc>=3 && argc<=6,"usage: {run|phases} K [seed=1] [calls=301] [normal|huge]");
    const std::string mode=argv[1];require(mode=="run" || mode=="phases","unknown mode");
    const auto k=std::stoull(argv[2]),seed=argc>3?std::stoull(argv[3]):1,calls=argc>4?std::stoull(argv[4]):301;
    require(calls && calls<=1000000,"invalid call count");
    const std::string memory=argc>5?argv[5]:"normal";require(memory=="normal" || memory=="huge","invalid memory policy");
    const auto policy=memory=="huge"?spin::MemoryPolicy::PreferHugePages:spin::MemoryPolicy::Normal;
    const rs::Plan plan(k,seed,seed);
    spin::Buffer input(2*k*16,policy),scratch(plan.scratchBlocks()*16,policy);fill(input);
    auto* x=blocks(input);auto* work=blocks(scratch);
    const auto label=std::string("rs16x16-byte-stream-")+memory;
    std::cout<<"mode,K,seed,calls,median_ms,p10_ms,p90_ms,checksum\n";
    timeCalls(label,k,seed,calls,[&]{rs::transposeFast(x,x,work,plan);},x);
    if(mode=="phases") {
        timeCalls(label+"-reverse-route",k,seed,calls,[&]{rs::reverseRoute(x,work,plan);},work);
        timeCalls(label+"-outer",k,seed,calls,[&]{rs::outerFast(work,x,plan);},x);
    }
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
