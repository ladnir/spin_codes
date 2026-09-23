#include "Routing.h"
#include "../feistel/Feistel.h"
#include <chrono>
#include <atomic>
#include <iomanip>
#include <iostream>
#include <string>
using Clock=std::chrono::steady_clock;
static double ms(Clock::time_point a,Clock::time_point b){return std::chrono::duration<double,std::milli>(b-a).count();}
static double median(std::vector<double> v){std::sort(v.begin(),v.end());return v[v.size()/2];}
constexpr unsigned N=1U<<19;
template<bool Batch,class Route> void fill(const Route& route,std::vector<unsigned>& out) {
    if constexpr(Batch)for(unsigned i=N;i;i-=16)route.outerBatch(i-16,out.data()+i-16);
    else for(unsigned i=N;i--;)out[i]=route.outer(i);
}
template<bool Batch,class Factory> void run(const char* name,Factory factory,double bankMs,std::size_t bankBytes) {
    std::vector<unsigned> out(N);std::vector<double> setup,first,warm;
    std::uint64_t checksum=0;std::size_t instanceBytes=0;
    for(unsigned rep=0;rep<11;++rep) {
        const auto a=Clock::now();auto route=factory(17+rep);const auto b=Clock::now();
        fill<Batch>(route,out);const auto c=Clock::now();
        for(unsigned j=0;j<4;++j){fill<Batch>(route,out);std::atomic_signal_fence(std::memory_order_seq_cst);}const auto d=Clock::now();
        if(rep>=2){setup.push_back(ms(a,b));first.push_back(ms(b,c));warm.push_back(ms(c,d)/4);}
        instanceBytes=route.bytes();
        // Outside timings: every output must occur once and batching must
        // equal the scalar path, not merely have a stable checksum.
        std::vector<bool> seen(N);
        for(unsigned i=0;i<N;++i) {
            if(out[i]>=N || seen[out[i]] || out[i]!=route.outer(i))throw std::runtime_error("routing mismatch");
            seen[out[i]]=true;checksum=(checksum^out[i])*0x100000001b3ULL;
        }
    }
    std::cout<<name<<','<<Batch<<','<<bankMs<<','<<bankBytes<<','<<instanceBytes<<','<<median(setup)<<','<<median(first)<<','<<median(warm)<<','<<std::hex<<checksum<<std::dec<<'\n';
}
template<unsigned Count,bool Batch> void bankRun() {
    const auto a=Clock::now();spin::experimental::bank::Tables<Count> tables(913);const auto b=Clock::now();
    run<Batch>(Count==16?"bank16":"bank64",[&](unsigned seed){return spin::experimental::bank::Routing<Count>(tables,seed);},ms(a,b),tables.bytes());
}
int main(int argc,char** argv){try {
    const std::string mode=argc>1?argv[1]:"bank16";const bool batch=argc<3 || std::string(argv[2])=="batch";
    std::cout<<std::fixed<<std::setprecision(6)<<"family,batch,bank_ms,bank_bytes,instance_bytes,setup_ms,first_route_ms,warm_route_ms,checksum\n";
    if(mode=="bank16"){if(batch)bankRun<16,true>();else bankRun<16,false>();}
    else if(mode=="bank64"){if(batch)bankRun<64,true>();else bankRun<64,false>();}
    else if(mode=="feistel6") {
        const auto factory=[](unsigned seed){return spin::experimental::feistel::Routing<11,6>(seed);};
        if(batch)run<true>("feistel6",factory,0,0);else run<false>("feistel6",factory,0,0);
    }else throw std::invalid_argument("family");
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
