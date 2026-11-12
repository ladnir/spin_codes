// Stage timing for the tuned Turbo RS8/s22 forward encoder: outer only,
// inner only, and both, per lane count. One process per invocation; run
// serially. Setup, allocation and input generation are outside the timer.
#include "../src/rs8/Rs8.h"
#include <algorithm>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <optional>
#include <vector>

namespace rs8=spin::detail::rs8;
namespace c22=spin::detail::certified22;
using Block=rs8::Block;

int main(int argc,char** argv) {try {
    if(argc!=5) throw std::invalid_argument("usage: spin_rs8_stages_bench K lanes outer|inner|both calls");
    const auto k=std::stoull(argv[1]);
    const auto lanes=unsigned(std::stoul(argv[2]));
    const std::string stage=argv[3];
    const auto calls=std::stoull(argv[4]);
    if(lanes!=1 && lanes!=2 && lanes!=4) throw std::invalid_argument("lanes");
    rs8::Plan plan(k,1,2);
    const auto scratchBlocks=plan.scratchBlocks()*lanes;
    auto alloc=[](std::size_t blocks) {
        auto* p=static_cast<Block*>(::operator new[](blocks*sizeof(Block),std::align_val_t{64}));
        std::memset(static_cast<void*>(p),0,blocks*sizeof(Block));
        return p;
    };
    auto* input=alloc(k*lanes);auto* output=alloc(2*k*lanes);auto* scratch=alloc(scratchBlocks);
    std::uint64_t state=913;
    auto* bytes=reinterpret_cast<unsigned char*>(input);
    for(std::size_t i=0;i<k*lanes*16;i+=8) {
        state+=0x9e3779b97f4a7c15ULL;auto v=state;v=(v^(v>>30))*0xbf58476d1ce4e5b9ULL;
        v=(v^(v>>27))*0x94d049bb133111ebULL;v^=v>>31;std::memcpy(bytes+i,&v,8);
    }
    rs8::outerForwardFast(input,scratch,plan,lanes);
    std::optional<spin::detail::paired15::Plan> p15;
    if(stage.rfind("p15",0)==0) p15.emplace(k,1,2);
    const auto run=[&] {
        if(stage=="p15reverse") { spin::detail::paired15::reverseRouteFast(output,scratch,*p15); return; }
        if(stage=="p15transpose") { spin::detail::paired15::transposeFast(output,input,scratch,*p15); return; }
        if(stage=="outer") rs8::outerForwardFast(input,scratch,plan,lanes);
        else if(stage=="inner") c22::forwardInnerFast(scratch,output,*plan.certified22Plan,lanes,false,plan.scratchBlocks());
        else if(stage=="treverse") c22::reverseRouteFast(output,scratch,*plan.certified22Plan);
        else if(stage=="touter") rs8::outerFast(scratch,input,plan);
        else if(stage=="transpose") rs8::transposeFast(output,input,scratch,plan);
        else rs8::forwardFast(input,output,scratch,plan,lanes,false);
    };
    for(unsigned i=0;i<8;++i) run();
    std::vector<double> times(calls);
    for(auto& t:times) {
        const auto begin=std::chrono::steady_clock::now();run();
        t=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-begin).count();
    }
    std::sort(times.begin(),times.end());
    std::cout<<"rs8-stage,"<<k<<','<<lanes<<','<<stage<<','<<calls<<','<<std::fixed<<std::setprecision(6)
        <<times[calls/2]<<','<<times[calls/10]<<','<<times[9*calls/10]<<'\n';
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
