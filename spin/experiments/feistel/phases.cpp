#include "Direct.h"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
using namespace spin::experimental::feistel;
using Clock=std::chrono::steady_clock;
static double ms(Clock::time_point a,Clock::time_point b){return std::chrono::duration<double,std::milli>(b-a).count();}
static double median(std::vector<double> x){std::sort(x.begin(),x.end());return x[x.size()/2];}
template<unsigned Rounds> void run(bool materialized) {
    std::vector<double> tables,expand,masks,total;
    std::uint64_t checksum=0;
    for(unsigned rep=0;rep<35;++rep) {
        const auto a=Clock::now();
        // Match the production order: route, then masks. Retain all outputs
        // through the timed region, and consume them afterwards.
        if constexpr(Rounds) {
            Routing<11,Rounds> routing(17+rep);
            const auto b=Clock::now();
            auto route=materialized?routing.materialize():std::vector<std::uint32_t>{};
            const auto c=Clock::now();
            auto mask=makeMasksK18(29+rep);
            const auto d=Clock::now();
            if(rep>=4){tables.push_back(ms(a,b));expand.push_back(ms(b,c));masks.push_back(ms(c,d));total.push_back(ms(a,d));}
            for(auto x:route)checksum=(checksum^x)*0x100000001b3ULL;
            for(auto x:mask)checksum=(checksum^x)*0x100000001b3ULL;
            for(unsigned i=0;i<(1U<<19);i+=127)checksum=(checksum^routing.outer(i))*0x100000001b3ULL;
        }else {
            auto route=uniformK18(17+rep);
            const auto b=Clock::now();
            auto mask=makeMasksK18(29+rep);
            const auto c=Clock::now();
            if(rep>=4){tables.push_back(ms(a,b));expand.push_back(0);masks.push_back(ms(b,c));total.push_back(ms(a,c));}
            for(auto x:route)checksum=(checksum^x)*0x100000001b3ULL;
            for(auto x:mask)checksum=(checksum^x)*0x100000001b3ULL;
        }
    }
    std::cout<<std::fixed<<std::setprecision(6)<<Rounds<<','<<materialized<<','<<median(tables)<<','<<median(expand)<<','<<median(masks)<<','<<median(total)<<','<<std::hex<<checksum<<std::dec<<'\n';
}
int main(int argc,char** argv){try {
    const unsigned rounds=argc>1?unsigned(std::stoul(argv[1])):0;
    const bool materialized=argc>2 && std::string(argv[2])=="materialized";
    std::cout<<"rounds,materialized,permutation_ms,materialize_ms,masks_ms,setup_ms,checksum\n";
    switch(rounds){case 0:run<0>(true);break;case 4:run<4>(materialized);break;case 6:run<6>(materialized);break;case 8:run<8>(materialized);break;default:throw std::invalid_argument("round count");}
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
