#include "Feistel.h"
#include "../../src/kernels/Spin.h"
#include <algorithm>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
namespace k=spin::detail::kernel;
using namespace spin::experimental::feistel;
using Clock=std::chrono::steady_clock;
static double ms(Clock::time_point a,Clock::time_point b){return std::chrono::duration<double,std::milli>(b-a).count();}
static double median(std::vector<double> x){std::sort(x.begin(),x.end());return x[x.size()/2];}
template<unsigned Bits,unsigned Rounds> void run(bool reuse) {
    constexpr unsigned K=(1U<<Bits)*128,N=2*K;
    std::vector<k::block> input(N);
    std::optional<k::Spin::Workspace> workspace;
    std::vector<double> preparation,plan,scratch,first,total,warm;
    std::uint64_t digest=0;std::size_t setupBytes=0;
    for(unsigned rep=0;rep<9;++rep) {
        for(unsigned i=0;i<N;++i)input[i]=k::block(i*0x9e3779b97f4a7c15ULL+rep,i^0xa531b78420ULL);
        const auto start=Clock::now();
        std::vector<std::uint32_t> route;
        if constexpr(Rounds) {Routing<Bits,Rounds> routing(17+rep);route=routing.materialize();}
        const auto prepared=Clock::now();
        k::Spin code(K==65536?k::Configuration::T64S12R2:k::Configuration::T128S19,
                     k::MessageLength{K},17+rep,29+rep,256,k::BchBackend::Auto,true,std::move(route));
        code.compact();const auto planned=Clock::now();
        if(!workspace)workspace.emplace(code);
        const auto ready=Clock::now();code.encodeInplace(input.data(),N,*workspace);const auto encoded=Clock::now();
        for(unsigned i=0;i<3;++i)code.encodeInplace(input.data(),N,*workspace);
        const auto warmed=Clock::now();
        if(rep>=2){preparation.push_back(ms(start,prepared));plan.push_back(ms(start,planned));
            scratch.push_back(ms(planned,ready));first.push_back(ms(ready,encoded));total.push_back(ms(start,encoded));warm.push_back(ms(encoded,warmed)/3);}
        for(const auto& v:input){std::uint64_t lanes[2];std::memcpy(lanes,&v,16);for(auto x:lanes)digest=(digest^x)*0x100000001b3ULL;}
        setupBytes=code.setupBytes();if(!reuse)workspace.reset();
    }
    std::cout<<std::fixed<<std::setprecision(6)<<Bits+7<<','<<Rounds<<','<<reuse<<','<<median(preparation)<<','<<median(plan)
             <<','<<median(scratch)<<','<<median(first)<<','<<median(total)<<','<<median(warm)<<','<<setupBytes<<','<<std::hex<<digest<<std::dec<<'\n';
}
template<unsigned Rounds> void sizes(bool reuse){run<9,Rounds>(reuse);run<11,Rounds>(reuse);run<13,Rounds>(reuse);}
int main(int argc,char** argv){try {
    const unsigned rounds=argc>1?unsigned(std::stoul(argv[1])):0;
    const bool reuse=argc>2 && std::string(argv[2])=="--reuse";
    std::cout<<"log_k,rounds,reuse,feistel_prep_ms,plan_ms,workspace_ms,first_encode_ms,total_ms,warm_ms,setup_bytes,checksum\n";
    switch(rounds){case 0:sizes<0>(reuse);break;case 4:sizes<4>(reuse);break;case 6:sizes<6>(reuse);break;case 8:sizes<8>(reuse);break;
        default:throw std::invalid_argument("rounds must be 0, 4, 6, or 8");}
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
