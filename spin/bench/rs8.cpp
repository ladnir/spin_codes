// Precomputed, single-threaded 128-bit transpose. Run each mode in a separate
// process under the same affinity; never run benchmark processes concurrently.
#include <spin/Code.h>
#include <spin/PacketCode.h>
#ifndef SPIN_RS8_RESEARCH_ORACLE
#define SPIN_RS8_RESEARCH_ORACLE 0
#endif
#if SPIN_RS8_RESEARCH_ORACLE
#include "../experiments/rs8_sizes/Kernels.h"
#include "../experiments/k20_two_ms/rs8/Rs8.h"
#endif
#include <algorithm>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#if SPIN_RS8_RESEARCH_ORACLE
namespace inner22=spin::detail::certified22;
namespace inner15=spin::detail::paired15;
namespace sizes=spin::research::rs8sizes;
namespace frozen=spin::research::twoms::rs8;
using Block=inner22::Block;
#endif
using Clock=std::chrono::steady_clock;

namespace {
#if SPIN_RS8_RESEARCH_ORACLE
Block* blocks(spin::Buffer& buffer) {
    return reinterpret_cast<Block*>(buffer.bytes().data());
}
#endif
void fill(spin::Buffer& buffer) {
    std::uint64_t state=913;
    for(std::size_t i=0;i<buffer.bytes().size();i+=8) {
        state+=0x9e3779b97f4a7c15ULL;auto word=state;
        word=(word^(word>>30))*0xbf58476d1ce4e5b9ULL;
        word=(word^(word>>27))*0x94d049bb133111ebULL;word^=word>>31;
        std::memcpy(buffer.bytes().data()+i,&word,8);
    }
}
std::uint64_t checksum(const std::byte* bytes,std::size_t count) {
    std::uint64_t result=0;
    for(std::size_t i=0;i<count;i+=8) {
        std::uint64_t word;std::memcpy(&word,bytes+i,8);
        result=(result^word)*0x100000001b3ULL;
    }
    return result;
}
template<class Encode>
void measure(const std::string& mode,std::size_t k,std::uint64_t seed,unsigned calls,
             std::span<std::byte> input,Encode&& encode) {
    // Setup, allocation, initial fill, and the five warmups are all excluded.
    // Inputs are not reset between calls, identically to the research campaign.
    for(unsigned i=0;i<5;++i)encode();
    std::vector<double> samples(calls);
    for(unsigned i=0;i<calls;++i) {
        const auto start=Clock::now();encode();
        samples[i]=std::chrono::duration<double,std::micro>(Clock::now()-start).count();
    }
    std::sort(samples.begin(),samples.end());
    std::cout<<"mode,K,seed,calls,memory,median_us,p10_us,p90_us,checksum\n"
        <<mode<<','<<k<<','<<seed<<','<<calls<<','<<(k<=65536?"normal":"huge-preferred")<<','
        <<std::fixed<<std::setprecision(4)<<samples[calls/2]<<','<<samples[calls/10]<<','
        <<samples[9*calls/10]<<','<<std::hex<<checksum(input.data(),k*16)<<std::dec<<'\n';
}
}

int main(int argc,char** argv) {try {
    if(!spin::packet_fast_available())return 77;
    if(argc<3 || argc>5)
        throw std::invalid_argument("usage: spin_rs8_bench public|frozen K [seed=113] [calls=301]");
    const std::string mode=argv[1];const auto k=std::stoull(argv[2]);
    const std::uint64_t seed=argc>3?std::stoull(argv[3]):113;
    const auto count=argc>4?std::stoull(argv[4]):301;
    if((mode!="public" && mode!="frozen") || k<2048 ||
       !spin::valid_message_size(spin::Parameters::PacketRs8T64,k) || count<1 || count>1000000)
        throw std::invalid_argument("invalid mode, natural size, or call count");
#if !SPIN_RS8_RESEARCH_ORACLE
    if(mode=="frozen")throw std::invalid_argument("frozen mode requires SPIN_BUILD_RESEARCH_ORACLE_TESTS");
#endif
    const auto calls=static_cast<unsigned>(count);
    const auto policy=k<=65536?spin::MemoryPolicy::Normal:spin::MemoryPolicy::PreferHugePages;
    // Match the research harness's four-record prefix and twelve-record extra
    // allocation; the 64-byte-aligned payload permits the same store policy.
    spin::Buffer buffer((2*k+12)*16,policy);fill(buffer);
    auto input=buffer.bytes().subspan(64,2*k*16);
    if(mode=="public") {
        spin::Code code({k,spin::Parameters::PacketRs8T64,seed,seed});
        if(code.backend()!=spin::Backend::Avx512)return 77;
        auto work=code.make_workspace(spin::Width::Bits128,policy);
        measure(mode,k,seed,calls,input,[&]{code.transpose_inplace_bytes(input,work);});
    }
#if SPIN_RS8_RESEARCH_ORACLE
    else {
        inner22::Plan geometry(k,seed,seed);
        frozen::Tables32 tables;frozen::setupField32(geometry,seed,tables);
        spin::Buffer scratch(geometry.scratchBlocks()*16,policy);
        auto* x=reinterpret_cast<Block*>(input.data());auto* work=blocks(scratch);
        if(k<=65536) {
            inner15::Plan plan(65536,seed,seed);sizes::resizePaired15(plan,geometry,seed);
            measure(mode,k,seed,calls,input,[&]{
                sizes::reverseRoute15(x,work,plan);frozen::outerFast(work,x,geometry,tables);
            });
        } else {
            measure(mode,k,seed,calls,input,[&]{
                inner22::reverseRouteFast(x,work,geometry);frozen::outerFast(work,x,geometry,tables);
            });
        }
    }
#endif
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
