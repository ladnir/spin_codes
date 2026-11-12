// Warm Turbo (PacketRs8T64) encode time per backend: 128-bit transpose and
// forward at 128/256/512-bit records with cached and streaming (s) output
// stores, median of reps. One backend per process;
// run serially, never concurrently with other benchmarks.
#include <spin/Code.h>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using Clock=std::chrono::steady_clock;
static double median(std::vector<double> x){std::sort(x.begin(),x.end());return x[x.size()/2];}
template<class F> static double time(unsigned reps,F f) {
    std::vector<double> t;f();
    for(unsigned r=0;r<reps;++r){const auto a=Clock::now();f();t.push_back(std::chrono::duration<double,std::milli>(Clock::now()-a).count());}
    return median(t);
}
int main(int argc,char** argv){try{
    if(argc!=4) throw std::invalid_argument("usage: spin_turbo_backends_bench auto|avx2|avx512|portable|reference|neon K reps");
    const std::string b=argv[1];const auto k=std::stoull(argv[2]);const auto reps=unsigned(std::stoul(argv[3]));
    const spin::ExecutionOptions exec{b=="avx2"?spin::Backend::Avx2:b=="avx512"?spin::Backend::Avx512:
        b=="portable"?spin::Backend::Portable:b=="reference"?spin::Backend::Reference:b=="neon"?spin::Backend::Neon:spin::Backend::Automatic};
    spin::Code code({k,spin::Parameters::PacketRs8T64,17,43},exec);
    std::cout<<b<<','<<k<<std::fixed<<std::setprecision(4);
    {
        auto w=code.make_workspace();
        spin::Buffer in(2*k*16),out(k*16);
        for(std::size_t i=0;i<in.bytes().size();++i) in.bytes()[i]=std::byte(i*131+7);
        std::cout<<",transpose_ms="<<time(reps,[&]{code.transpose_bytes(in.bytes(),out.bytes(),w);});
    }
    for(auto width:{spin::Width::Bits128,spin::Width::Bits256,spin::Width::Bits512}) {
        const auto bytes=static_cast<unsigned>(width);
        auto w=code.make_workspace(width);
        spin::Buffer in(k*bytes),out(2*k*bytes);
        for(std::size_t i=0;i<in.bytes().size();++i) in.bytes()[i]=std::byte(i*131+7);
        std::cout<<",forward"<<8*bytes<<"_ms="<<time(reps,[&]{code.forward_bytes(in.bytes(),out.bytes(),w,spin::OutputStores::Cached);})
                 <<",forward"<<8*bytes<<"s_ms="<<time(reps,[&]{code.forward_bytes(in.bytes(),out.bytes(),w,spin::OutputStores::Streaming);});
    }
    std::cout<<'\n';
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
