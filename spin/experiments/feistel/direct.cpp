#include "Direct.h"
#include <algorithm>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
using namespace spin::experimental::feistel;
using Clock=std::chrono::steady_clock;
static double ms(Clock::time_point a,Clock::time_point b){return std::chrono::duration<double,std::milli>(b-a).count();}
static double median(std::vector<double> x){std::sort(x.begin(),x.end());return x[x.size()/2];}
template<unsigned Rounds,bool Materialized> void run(bool reuse) {
    constexpr unsigned N=1U<<19;
    std::vector<kernel::block> input(N),scratch;
    // Validate this fused path against the full-schedule implementation before
    // timing it. The latter is checked separately against the dense oracle.
    {
        std::vector<std::uint32_t> route;
        if constexpr(Rounds)route=Routing<11,Rounds>(17).materialize();
        kernel::Spin control(kernel::Configuration::T128S19,kernel::MessageLength{N/2},17,29,256,kernel::BchBackend::Auto,true,std::move(route));
        control.compact();kernel::Spin::Workspace w(control);
        kernelWords words(913);for(auto& x:input){auto a=words(),b=words();x=kernel::block(a,b);}
        auto expected=input;control.encodeInplace(expected.data(),N,w);
        Direct<Rounds,Materialized> candidate(17,29);scratch.resize(N);candidate.encode(input.data(),scratch);
        if(std::memcmp(input.data(),expected.data(),N*16))throw std::runtime_error("direct path reference mismatch");
        std::vector<kernel::block>().swap(scratch);
    }
    std::vector<double> plan,work,first,total,warm;std::uint64_t checksum=0;std::size_t bytes=0;
    for(unsigned rep=0;rep<9;++rep) {
        for(unsigned i=0;i<N;++i)input[i]=kernel::block(i*0x9e3779b97f4a7c15ULL+rep,i^0xa531b78420ULL);
        const auto a=Clock::now();Direct<Rounds,Materialized> code(17+rep,29+rep);const auto b=Clock::now();
        if(scratch.empty())scratch.resize(N);
        const auto c=Clock::now();code.encode(input.data(),scratch);const auto d=Clock::now();
        for(unsigned j=0;j<3;++j)code.encode(input.data(),scratch);
        const auto e=Clock::now();
        if(rep>=2){plan.push_back(ms(a,b));work.push_back(ms(b,c));first.push_back(ms(c,d));total.push_back(ms(a,d));warm.push_back(ms(d,e)/3);}
        bytes=code.bytes();for(const auto& v:input){std::uint64_t lanes[2];std::memcpy(lanes,&v,16);for(auto x:lanes)checksum=(checksum^x)*0x100000001b3ULL;}
        if(!reuse)std::vector<kernel::block>().swap(scratch);
    }
    std::cout<<std::fixed<<std::setprecision(6)<<Rounds<<','<<Materialized<<','<<reuse<<','<<median(plan)<<','<<median(work)
             <<','<<median(first)<<','<<median(total)<<','<<median(warm)<<','<<bytes<<','<<std::hex<<checksum<<std::dec<<'\n';
}
template<unsigned Rounds> void mode(bool materialized,bool reuse){if(materialized)run<Rounds,true>(reuse);else run<Rounds,false>(reuse);}
int main(int argc,char** argv){try {
    const unsigned rounds=argc>1?unsigned(std::stoul(argv[1])):0;
    const bool materialized=argc>2 && std::string(argv[2])=="materialized";
    const bool reuse=argc>3 && std::string(argv[3])=="--reuse";
    std::cout<<"rounds,materialized,reuse,plan_ms,workspace_ms,first_ms,total_ms,warm_ms,setup_bytes,checksum\n";
    switch(rounds){case 0:run<0,true>(reuse);break;case 4:mode<4>(materialized,reuse);break;case 6:mode<6>(materialized,reuse);break;case 8:mode<8>(materialized,reuse);break;default:throw std::invalid_argument("round count");}
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
