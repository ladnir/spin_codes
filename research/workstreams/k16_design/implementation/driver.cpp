#include "RsPrototype.h"
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

namespace rs = spin::research::rs;
using Block = rs::Block;

static Block* blocks(spin::Buffer& b) {
    return reinterpret_cast<Block*>(b.bytes().data());
}
static void require(bool ok, const char* message) {
    if(!ok) throw std::runtime_error(message);
}
static void fill(spin::Buffer& b, std::uint64_t state=913) {
    for(std::size_t i=0; i<b.bytes().size(); i+=8) {
        state += 0x9e3779b97f4a7c15ULL;
        auto v=state;
        v=(v^(v>>30))*0xbf58476d1ce4e5b9ULL;
        v=(v^(v>>27))*0x94d049bb133111ebULL;
        v^=v>>31;
        std::memcpy(b.bytes().data()+(i^8), &v, 8);
    }
}
static unsigned dot(const Block* a,const Block* b,std::size_t n) {
    unsigned parity=0;
    for(std::size_t i=0;i<n;++i) {
        std::uint64_t x[2],y[2];
        std::memcpy(x,a+i,16);std::memcpy(y,b+i,16);
        parity^=std::popcount(x[0]&y[0])^std::popcount(x[1]&y[1]);
    }
    return parity&1;
}
static std::uint64_t checksum(const Block* data,std::size_t n) {
    std::uint64_t result=0;
    for(std::size_t i=0;i<n;++i) {
        std::uint64_t x[2];std::memcpy(x,data+i,16);
        result=(result^x[0])*0x100000001b3ULL;
        result=(result^x[1])*0x100000001b3ULL;
    }
    return result;
}

template<class Word,std::size_t Width>
static bool fullRank(std::array<Word,Width> rows) {
    unsigned rank=0;
    for(unsigned column=0;column<Width;++column) {
        unsigned pivot=rank;
        while(pivot<Width && !(rows[pivot]&(Word{1}<<column)))++pivot;
        if(pivot==Width)return false;
        std::swap(rows[pivot],rows[rank]);
        for(unsigned j=rank+1;j<Width;++j)
            if(rows[j]&(Word{1}<<column))rows[j]^=rows[rank];
        ++rank;
    }
    return true;
}

static void checkSetup(const rs::Plan& plan) {
    require(plan.route.size()==plan.n/4,"RS route length");
    std::vector<bool> seen(plan.groups*64);
    for(unsigned region=0;region<64;++region) {
        std::vector<bool> groups(plan.groups);
        for(std::size_t slot=0;slot<plan.groups;++slot) {
            const auto base=plan.route[region*plan.groups+slot];
            const auto group=base/rs::groupStride,offset=base%rs::groupStride;
            require(group<plan.groups && offset<256 && offset%4==0,"RS invalid route address");
            require(!groups[group],"RS region repeats a group");groups[group]=true;
            const auto column=offset/4;
            require(!seen[64*group+column],"RS column permutation repeats a column");
            seen[64*group+column]=true;
        }
    }
    require(std::all_of(seen.begin(),seen.end(),[](bool v){return v;}),"RS route misses columns");
    require(plan.reverseMatrices.size()==plan.n/64,"RS update count");
    if(plan.innerKernel==rs::InnerKernel::Cached)
        require(plan.composedUpdates.size()==4*(plan.n/64) && plan.denseUpdates.empty(),"RS cached table selection");
    else
        require(plan.denseUpdates.size()==plan.n/64 && plan.composedUpdates.empty(),"RS streaming table selection");
    for(const auto& rows:plan.reverseMatrices)require(fullRank(rows),"RS update matrix is singular");
    if(plan.variant==rs::Variant::Rs8Gf256) {
        require(plan.outerMatrices.size()==8*plan.groups && plan.outerMatrices16.empty(),"RS8 outer count");
        for(const auto& rows:plan.outerMatrices)require(fullRank(rows),"RS8 outer matrix is singular");
    } else {
        require(plan.outerMatrices16.size()==16*plan.groups && plan.outerMatrices.empty(),"RS16 outer count");
        for(const auto& rows:plan.outerMatrices16)require(fullRank(rows),"RS16 outer matrix is singular");
    }
}

static void check(std::size_t k,std::uint64_t seed,rs::Variant variant,rs::InnerKernel kernel) {
    const rs::Plan plan(k,seed,seed+17,variant,kernel);
    checkSetup(plan);
    if(kernel==rs::InnerKernel::RetainedStreaming) {
        const rs::Plan control(k,seed,seed+17,variant);
        require(plan.route==control.route && plan.reverseMatrices==control.reverseMatrices &&
            plan.outerMatrices==control.outerMatrices && plan.outerMatrices16==control.outerMatrices16,
            "RS kernel selection changed sampled map");
        const auto coefficientCount=(variant==rs::Variant::Rs16Gf16?64:128)*plan.groups;
        require(std::equal(plan.compactCoefficients(),plan.compactCoefficients()+coefficientCount,
            control.compactCoefficients()),"RS kernel selection changed outer coefficients");
        for(std::size_t epoch=0;epoch<plan.n/64;++epoch)
            for(unsigned word=0;word<4;++word) {
                std::uint64_t reversed=0;
                for(unsigned byte=0;byte<8;++byte)
                    reversed|=((control.composedUpdates[4*epoch+word]>>(8*byte))&255ULL)<<(8*(7-byte));
                require(plan.denseUpdates[epoch].matrix[word]==reversed,"RS retained inner matrix repacking");
            }
    }
    const auto encode=[&](const Block* input,Block* output,Block* scratch) {
        if(kernel==rs::InnerKernel::RetainedStreaming)rs::transposeLarge(input,output,scratch,plan);
        else rs::transposeFast(input,output,scratch,plan);
    };
    const auto n=2*k;
    spin::Buffer input((n+8)*16),expected(k*16),output((k+8)*16);
    spin::Buffer scratch(plan.scratchBlocks()*16),message(k*16),encoded(n*16);
    fill(input);fill(message,2947);
    const std::vector<std::byte> saved(input.bytes().begin(),input.bytes().end());
    auto* x=blocks(input)+4;
    rs::transposeScalar(x,blocks(expected),blocks(scratch),plan);
    for(unsigned offset=0;offset<4;++offset) {
        std::fill(output.bytes().begin(),output.bytes().end(),std::byte{0xa5});
        auto* y=blocks(output)+4+offset;
        // Change input alignment without changing the mathematical input.
        std::memmove(blocks(input)+4+offset,saved.data()+64,n*16);
        encode(blocks(input)+4+offset,y,blocks(scratch));
        require(std::memcmp(y,blocks(expected),k*16)==0,"RS fast vs scalar mismatch");
        require(std::all_of(output.bytes().begin(),output.bytes().begin()+(4+offset)*16,
            [](std::byte v){return v==std::byte{0xa5};}),"RS output prefix guard");
        require(std::all_of(output.bytes().begin()+(4+offset+k)*16,output.bytes().end(),
            [](std::byte v){return v==std::byte{0xa5};}),"RS output suffix guard");
        require(std::memcmp(blocks(input)+4+offset,saved.data()+64,n*16)==0,
            "RS separate transpose changed input");
    }
    std::memcpy(input.bytes().data(),saved.data(),saved.size());
    encode(x,x,blocks(scratch));
    require(std::memcmp(x,blocks(expected),k*16)==0,"RS inplace mismatch");
    require(std::memcmp(x+k,saved.data()+64+k*16,k*16)==0,"RS inplace suffix changed");
    require(std::memcmp(blocks(input),saved.data(),64)==0,"RS inplace prefix guard");
    require(std::memcmp(x+n,saved.data()+(n+4)*16,64)==0,"RS inplace suffix guard");
    std::memcpy(input.bytes().data(),saved.data(),saved.size());
    rs::forwardScalar(blocks(message),blocks(encoded),plan);
    require(dot(blocks(encoded),x,n)==dot(blocks(message),blocks(expected),k),
        "RS binary adjoint identity failed");
    std::cout<<"PASS variant="<<int(variant)<<" kernel="<<int(kernel)<<" K="<<k<<" route_seed="<<seed<<" inner_seed="<<seed+17
             <<" checksum="<<std::hex<<checksum(blocks(expected),k)<<std::dec<<'\n';
}

template<class Encode>
static void timeCalls(const std::string& label,std::size_t k,std::uint64_t seed,
                      std::size_t calls,Encode&& encode,const Block* output) {
    for(unsigned i=0;i<5;++i)encode();
    std::vector<double> times(calls);
    for(auto& time:times) {
        const auto begin=std::chrono::steady_clock::now();
        encode();
        time=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-begin).count();
    }
    std::sort(times.begin(),times.end());
    std::cout<<label<<','<<k<<','<<seed<<','<<calls<<','<<std::fixed<<std::setprecision(6)
             <<times[calls/2]<<','<<times[calls/10]<<','<<times[9*calls/10]<<','
             <<std::hex<<checksum(output,k)<<std::dec<<'\n';
}

int main(int argc,char** argv) {try {
    if(!spin::packet_fast_available() || !spin::capabilities().forward512) {
        // Research sources also allow AVX512DQ. The library's packet-only
        // GCC gate does not require it; forward512 includes that extra check.
        std::cerr<<"AVX512/DQ/VBMI/GFNI required for this research comparison\n";return 77;
    }
    if(argc>=2 && std::string(argv[1])=="check") {
        if(argc==2) {
            for(auto variant:{rs::Variant::Rs8Gf256,rs::Variant::Rs16Gf16})
                for(auto k:{std::size_t(4096),std::size_t(65536)})
                    for(auto seed:{1ULL,7ULL})
                        for(auto kernel:{rs::InnerKernel::Cached,rs::InnerKernel::RetainedStreaming})
                            check(k,seed,variant,kernel);
        } else {
            if(argc>4)throw std::invalid_argument("usage: spin_k16_rs check [K [seed=1]]");
            const auto k=std::stoull(argv[2]),seed=argc>3?std::stoull(argv[3]):1;
            for(auto variant:{rs::Variant::Rs8Gf256,rs::Variant::Rs16Gf16})
                for(auto kernel:{rs::InnerKernel::Cached,rs::InnerKernel::RetainedStreaming})
                    check(k,seed,variant,kernel);
        }
        return 0;
    }
    if(argc<3 || argc>6)throw std::invalid_argument("usage: spin_k16_rs check [K [seed=1]] | {rs|rs16}[-stream][-phases]|bch K [seed=1] [calls=301] [normal|huge]");
    const std::string mode=argv[1];
    const auto k=std::stoull(argv[2]),seed=argc>3?std::stoull(argv[3]):1;
    const auto calls=argc>4?std::stoull(argv[4]):301;
    if(!calls || calls>1000000)throw std::invalid_argument("calls must be in [1,1000000]");
    const std::string memory=argc>5?argv[5]:"normal";
    if(memory!="normal" && memory!="huge")throw std::invalid_argument("memory must be normal or huge");
    const auto policy=memory=="huge"?spin::MemoryPolicy::PreferHugePages:spin::MemoryPolicy::Normal;
    const std::string suffix=memory=="huge"?"-huge":"";
    std::cout<<"mode,K,seed,calls,median_ms,p10_ms,p90_ms,checksum\n";
    if(mode=="bch") {
        spin::Code code({k,spin::Parameters::PacketT64S16,seed,seed});
        auto scratch=code.make_workspace(spin::Width::Bits128,policy);
        auto input=code.make_buffer(spin::Width::Bits128,policy);fill(input);
        timeCalls(mode+suffix,k,seed,calls,[&]{code.transpose_inplace_bytes(input.bytes(),scratch);},blocks(input));
    } else {
        const bool phases=mode.ends_with("-phases");
        auto family=phases?mode.substr(0,mode.size()-7):mode;
        const bool streaming=family.ends_with("-stream");
        if(streaming)family.resize(family.size()-7);
        if(family!="rs" && family!="rs16")throw std::invalid_argument("unknown mode");
        const bool narrow=family=="rs16";
        const std::string name=family+(streaming?"-stream":"")+suffix;
        const rs::Plan plan(k,seed,seed,narrow?rs::Variant::Rs16Gf16:rs::Variant::Rs8Gf256,
            streaming?rs::InnerKernel::RetainedStreaming:rs::InnerKernel::Cached);
        spin::Buffer input(2*k*16,policy),scratch(plan.scratchBlocks()*16,policy);fill(input);
        auto* x=blocks(input);auto* work=blocks(scratch);
        if(streaming)timeCalls(name,k,seed,calls,[&]{rs::transposeLarge(x,x,work,plan);},x);
        else timeCalls(name,k,seed,calls,[&]{rs::transposeFast(x,x,work,plan);},x);
        if(phases) {
            if(streaming)timeCalls(name+"-reverse-route",k,seed,calls,[&]{rs::reverseRouteLarge(x,work,plan);},work);
            else timeCalls(name+"-reverse-route",k,seed,calls,[&]{rs::reverseRouteFast(x,work,plan);},work);
            timeCalls(name+"-outer",k,seed,calls,[&]{rs::outerFast(work,x,plan);},x);
        }
    }
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
