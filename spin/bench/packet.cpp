#include <spin/Code.h>
#include <algorithm>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// Serial, precomputed transpose benchmark. Setup and memory preparation are
// outside timing; repeated calls reuse the in-place buffer, as the research run.
int main(int argc,char** argv) {try {
    if(argc<2 || argc>5)throw std::invalid_argument("usage: spin_packet_bench K [seed=1] [calls=301] [normal|huge]");
    const auto k=std::stoull(argv[1]);const auto seed=argc>2?std::stoull(argv[2]):1;
    const auto calls=argc>3?std::stoull(argv[3]):301;
    if(!calls || calls>1000000)throw std::invalid_argument("calls must be in [1,1000000]");
    const std::string policy=argc>4?argv[4]:"normal";
    if(policy!="normal" && policy!="huge")throw std::invalid_argument("unknown memory policy");
    const auto memory=policy=="huge"?spin::MemoryPolicy::PreferHugePages:spin::MemoryPolicy::Normal;
    spin::Code code({k,spin::Parameters::PacketT64S16,seed,seed});
    auto work=code.make_workspace(spin::Width::Bits128,memory);
    auto buffer=code.make_buffer(spin::Width::Bits128,memory);
    std::uint64_t state=913;
    for(std::size_t i=0;i<buffer.bytes().size();i+=8) {
        state+=0x9e3779b97f4a7c15ULL;auto v=state;v=(v^(v>>30))*0xbf58476d1ce4e5b9ULL;
        v=(v^(v>>27))*0x94d049bb133111ebULL;v^=v>>31;
        // Retained research block(hi,lo) draws hi first.
        std::memcpy(buffer.bytes().data()+(i^8),&v,8);
    }
    for(unsigned i=0;i<5;++i)code.transpose_inplace_bytes(buffer.bytes(),work);
    std::vector<double> times(calls);
    for(auto& t:times) {const auto begin=std::chrono::steady_clock::now();
        code.transpose_inplace_bytes(buffer.bytes(),work);
        t=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-begin).count();}
    std::sort(times.begin(),times.end());
    std::uint64_t hash=0;
    for(std::size_t i=0;i<buffer.bytes().size();i+=8) {
        std::uint64_t v;std::memcpy(&v,buffer.bytes().data()+i,8);hash=(hash^v)*0x100000001b3ULL;
    }
    std::cout<<"K,seed,calls,backend,memory,median_ms,p10_ms,p90_ms,checksum\n"
        <<k<<','<<seed<<','<<calls<<','<<int(code.backend())<<','<<policy<<','<<std::fixed<<std::setprecision(6)
        <<times[calls/2]<<','<<times[calls/10]<<','<<times[9*calls/10]<<','<<std::hex<<hash<<'\n';
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
