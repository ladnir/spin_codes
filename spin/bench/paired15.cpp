#include <spin/Code.h>
#include <algorithm>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// One process measures one case. Invoke serially in both orders for comparisons.
// Setup, allocations, input generation and checksum are outside the timer.
int main(int argc,char** argv) {try {
    if(argc!=7)throw std::invalid_argument(
        "usage: spin_paired15_bench old|paired15 forward|transpose 128|256|512 seed calls normal|huge");
    const std::string family=argv[1],direction=argv[2],policy=argv[6];
    if(family!="old" && family!="paired15")throw std::invalid_argument("unknown profile");
    if(direction!="forward" && direction!="transpose")throw std::invalid_argument("unknown direction");
    const auto bits=std::stoul(argv[3]);
    if(bits!=128 && bits!=256 && bits!=512)throw std::invalid_argument("unknown width");
    if(direction=="transpose" && bits!=128)throw std::invalid_argument("transpose is 128-bit only");
    const auto seed=std::stoull(argv[4]),calls=std::stoull(argv[5]);
    if(!calls || calls>1000000)throw std::invalid_argument("invalid calls");
    if(policy!="normal" && policy!="huge")throw std::invalid_argument("unknown memory policy");
    const auto width=static_cast<spin::Width>(bits/8);
    const auto memory=policy=="normal"?spin::MemoryPolicy::Normal:spin::MemoryPolicy::PreferHugePages;
    const auto parameter=family=="old"?spin::Parameters::PacketRsT64S20:spin::Parameters::PacketRsT64S15K16;
    constexpr std::size_t k=65536;
    spin::Code code({k,parameter,seed,seed});
    auto work=code.make_workspace(width,memory);
    auto input=code.make_buffer(width,memory),output=code.make_buffer(width,memory);
    std::uint64_t state=913;
    for(std::size_t i=0;i<input.bytes().size();i+=8) {
        state+=0x9e3779b97f4a7c15ULL;auto v=state;v=(v^(v>>30))*0xbf58476d1ce4e5b9ULL;
        v=(v^(v>>27))*0x94d049bb133111ebULL;v^=v>>31;
        std::memcpy(input.bytes().data()+i,&v,8);
    }
    const auto forward=direction=="forward";
    const auto in=input.bytes().first((forward?k:2*k)*(bits/8));
    const auto out=output.bytes().first((forward?2*k:k)*(bits/8));
    const auto run=[&] {
        if(forward)code.forward_bytes(in,out,work);
        else code.transpose_bytes(in,out,work);
    };
    for(unsigned i=0;i<8;++i)run();
    std::vector<double> times(calls);
    for(auto& t:times) {
        const auto begin=std::chrono::steady_clock::now();run();
        t=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-begin).count();
    }
    std::sort(times.begin(),times.end());
    std::uint64_t hash=0;
    for(std::size_t i=0;i<out.size();i+=8) {
        std::uint64_t v;std::memcpy(&v,out.data()+i,8);hash=(hash^v)*0x100000001b3ULL;
    }
    std::cout<<"profile,direction,bits,seed,calls,backend,memory,setup_bytes,scratch_bytes,median_ms,p10_ms,p90_ms,checksum\n"
        <<family<<','<<direction<<','<<bits<<','<<seed<<','<<calls<<','<<int(code.backend())<<','<<policy<<','
        <<code.setup_bytes()<<','<<work.bytes()<<','<<std::fixed<<std::setprecision(6)
        <<times[calls/2]<<','<<times[calls/10]<<','<<times[9*calls/10]<<','<<std::hex<<hash<<'\n';
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
