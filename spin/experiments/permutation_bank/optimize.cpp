#include "Rolling.h"
#include "Composed.h"
#include "../../src/kernels/WorkspaceRouting.h"
#include <spin/Code.h>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <optional>
using namespace spin::experimental::bank;
using Clock=std::chrono::steady_clock;
#ifndef SPIN_BANK_TIMING_REPS
#define SPIN_BANK_TIMING_REPS 15
#endif
#ifndef SPIN_BANK_TIMING_WARMUP
#define SPIN_BANK_TIMING_WARMUP 2
#endif
static_assert(SPIN_BANK_TIMING_REPS>SPIN_BANK_TIMING_WARMUP);
static double ms(Clock::time_point a,Clock::time_point b){return std::chrono::duration<double,std::milli>(b-a).count();}
static double median(std::vector<double> x){std::sort(x.begin(),x.end());return x[x.size()/2];}
static void setupCosts() {
    std::vector<double> bankTimes,rowTimes,maskTimes,fixedTimes;
    std::uint64_t checksum=0;std::size_t bankBytes=0,planBytes=0,workspaceBytes=0;
    for(unsigned rep=0;rep<15;++rep) {
        auto a=Clock::now();ComposedBank<1> bank(913+rep);auto b=Clock::now();
        bankBytes=bank.bytes();if(rep>=2)bankTimes.push_back(ms(a,b));
        checksum^=bank.routes[rep];
        a=Clock::now();ComposedRouting<1,true,true,true> route(bank,17+rep);b=Clock::now();
        if(rep>=2)rowTimes.push_back(ms(a,b));checksum^=route.outer(rep);
        a=Clock::now();auto masks=fastMasks(29+rep);b=Clock::now();
        if(rep>=2)maskTimes.push_back(ms(a,b));checksum^=masks[rep];
        a=Clock::now();PrecomputedCode fixedCode(17+rep,29+rep);b=Clock::now();
        planBytes=fixedCode.planBytes();workspaceBytes=fixedCode.workspaceBytes();
        if(rep>=2)fixedTimes.push_back(ms(a,b));
    }
    std::cout<<"component,median_ms,retained_bytes\n"
      <<"reusable_composed_bank,"<<median(bankTimes)<<','<<bankBytes<<'\n'
      <<"fresh_affine_parameters,"<<median(rowTimes)<<','<<sizeof(ComposedRouting<1,true,true,true>)<<'\n'
      <<"fresh_imt_masks,"<<median(maskTimes)<<','<<2*(CodeSize/128)*sizeof(std::uint32_t)<<'\n'
      <<"original_fixed_plan_and_workspace,"<<median(fixedTimes)<<','<<planBytes+workspaceBytes<<'\n';
    std::cerr<<"setup checksum="<<std::hex<<checksum<<'\n';
}
static void fixed(bool bank) {
    BankFlow<16> flow(913);
    PrecomputedCode code(17,29,bank?flow.materialize(17):std::vector<std::uint32_t>{});
    std::vector<kernel::block> input(CodeSize);std::vector<double> times;std::uint64_t checksum=0;
    for(unsigned rep=0;rep<SPIN_BANK_TIMING_REPS;++rep) {
        Words words(913+rep);for(auto& v:input){auto a=words(),b=words();v=kernel::block(a,b);}
        const auto a=Clock::now();code.encode(input.data());const auto b=Clock::now();
        if(rep>=SPIN_BANK_TIMING_WARMUP)times.push_back(ms(a,b));
        for(unsigned i=0;i<MessageSize;++i){std::uint64_t lanes[2];std::memcpy(lanes,&input[i],16);for(auto x:lanes)checksum=(checksum^x)*0x100000001b3ULL;}
    }
    std::cout<<(bank?"fixed-bank":"fixed")<<",0.000000,"<<median(times)<<','<<median(times)<<','<<std::hex<<checksum<<std::dec<<'\n';
}
template<class TablesType=Tables<16>,class RouteType=Routing<16>,class Factory> void run(const char* mode,Factory make,bool testOnly) {
    TablesType tables(913);std::vector<kernel::block> input(CodeSize),scratch(CodeSize);
    if constexpr(requires{tables.routes;}) {
        const auto bankAdvice=kernel::workspace_routing::adviseOwned(tables.routes.data(),tables.routes.size()*sizeof(unsigned));
        std::cerr<<"bank hugepage hint="<<bankAdvice.hinted<<" collapsed="<<bankAdvice.collapsed<<'\n';
    }
    const auto advice=kernel::workspace_routing::adviseOwned(scratch.data(),scratch.size()*sizeof(kernel::block));
    std::cerr<<"scratch hugepage hint="<<advice.hinted<<" collapsed="<<advice.collapsed<<'\n';
    for(unsigned seed:{17U,191U}) {
        RouteType routing(tables,seed);std::vector<unsigned> route(CodeSize);
        for(unsigned i=0;i<CodeSize;++i)route[i]=routing.outer(i);
        std::vector<bool> seen(CodeSize);std::array<unsigned,2048> regionCounts{};
        for(unsigned i=0;i<CodeSize;++i){const auto x=route[i];if(x>=CodeSize || seen[x])throw std::runtime_error("route not bijective");seen[x]=true;if(!(i&2047))regionCounts.fill(0);if(++regionCounts[x>>8]!=1)throw std::runtime_error("multiple row bits in a region");}
        PrecomputedCode reference(seed,seed+12,std::move(route));auto code=make(tables,seed,seed+12);
        Words words(791);for(auto& v:input){auto a=words(),b=words();v=kernel::block(a,b);}
        auto expected=input;reference.encode(expected.data());code.encode(input.data(),scratch);
        // The Split schedule is allowed to overwrite the unused suffix.
        if(std::memcmp(input.data(),expected.data(),MessageSize*sizeof(kernel::block)))throw std::runtime_error("optimized output mismatch");
    }
    if(testOnly){std::cout<<mode<<" reference PASS\n";return;}
    std::vector<double> setup,encode,total;std::uint64_t checksum=0;
    for(unsigned rep=0;rep<SPIN_BANK_TIMING_REPS;++rep) {
        Words words(913+rep);for(auto& v:input){auto a=words(),b=words();v=kernel::block(a,b);}
        const auto a=Clock::now();auto code=std::make_optional(make(tables,17+rep,29+rep));const auto b=Clock::now();
        code->encode(input.data(),scratch);const auto c=Clock::now();code.reset();const auto d=Clock::now();
        if(rep>=SPIN_BANK_TIMING_WARMUP){setup.push_back(ms(a,b));encode.push_back(ms(b,c));total.push_back(ms(a,d));}
        for(unsigned i=0;i<MessageSize;++i){std::uint64_t lanes[2];std::memcpy(lanes,&input[i],16);for(auto x:lanes)checksum=(checksum^x)*0x100000001b3ULL;}
    }
    std::cout<<mode<<','<<median(setup)<<','<<median(encode)<<','<<median(total)<<','<<std::hex<<checksum<<std::dec<<'\n';
}
int main(int argc,char** argv){try {
    if(!spin::capabilities().avx2)return 77;
    const std::string mode=argc>1?argv[1]:"rolling";const bool test=argc>2 && std::string(argv[2])=="verify";
    if(mode=="setup-composed"){std::cout<<std::fixed<<std::setprecision(6);setupCosts();return 0;}
    if(!test)std::cout<<std::fixed<<std::setprecision(6)<<"mode,setup_ms,encode_ms,total_ms,checksum\n";
    if(mode=="fixed" || mode=="fixed-bank")fixed(mode=="fixed-bank");
    else if(mode=="row-rotate1")run<ComposedBank<1>,ComposedRouting<1,true,true,true,true>>("row-rotate1",[](const auto& t,unsigned s,unsigned m){return RollingCode<1,false,false,0,true,ComposedRouting<1,true,true,true,true>>(ComposedRouting<1,true,true,true,true>(t,s),m);},test);
    else if(mode=="row-affine1")run<ComposedBank<1>,ComposedRouting<1,true,true,true>>("row-affine1",[](const auto& t,unsigned s,unsigned m){return RollingCode<1,false,false,0,true,ComposedRouting<1,true,true,true>>(ComposedRouting<1,true,true,true>(t,s),m);},test);
    else if(mode=="row-shift4-epoch")run<ComposedBank<4>,ComposedRouting<4,true,true>>("row-shift4-epoch",[](const auto& t,unsigned s,unsigned m){return RollingCode<4,false,false,0,true,ComposedRouting<4,true,true>>(ComposedRouting<4,true,true>(t,s),m);},test);
    else if(mode=="row-shift1-epoch")run<ComposedBank<1>,ComposedRouting<1,true,true>>("row-shift1-epoch",[](const auto& t,unsigned s,unsigned m){return RollingCode<1,false,false,0,true,ComposedRouting<1,true,true>>(ComposedRouting<1,true,true>(t,s),m);},test);
    else if(mode=="row-shift4")run<ComposedBank<4>,ComposedRouting<4,true,true>>("row-shift4",[](const auto& t,unsigned s,unsigned m){return RollingCode<4,false,false,0,false,ComposedRouting<4,true,true>>(ComposedRouting<4,true,true>(t,s),m);},test);
    else if(mode=="row-shift1")run<ComposedBank<1>,ComposedRouting<1,true,true>>("row-shift1",[](const auto& t,unsigned s,unsigned m){return RollingCode<1,false,false,0,false,ComposedRouting<1,true,true>>(ComposedRouting<1,true,true>(t,s),m);},test);
    else if(mode=="shift16")run<ComposedBank<16>,ComposedRouting<16,true>>("shift16",[](const auto& t,unsigned s,unsigned m){return RollingCode<16,false,false,0,false,ComposedRouting<16,true>>(ComposedRouting<16,true>(t,s),m);},test);
    else if(mode=="shift4")run<ComposedBank<4>,ComposedRouting<4,true>>("shift4",[](const auto& t,unsigned s,unsigned m){return RollingCode<4,false,false,0,false,ComposedRouting<4,true>>(ComposedRouting<4,true>(t,s),m);},test);
    else if(mode=="shift1")run<ComposedBank<1>,ComposedRouting<1,true>>("shift1",[](const auto& t,unsigned s,unsigned m){return RollingCode<1,false,false,0,false,ComposedRouting<1,true>>(ComposedRouting<1,true>(t,s),m);},test);
    else if(mode=="composed-direct")run<ComposedBank<16>,ComposedRouting<16>>("composed-direct",[](const auto& t,unsigned s,unsigned m){return DirectComposed<16>(t,s,m);},test);
    else if(mode=="composed4")run<ComposedBank<4>,ComposedRouting<4>>("composed4",[](const auto& t,unsigned s,unsigned m){return RollingCode<4,false,false,0,false,ComposedRouting<4>>(ComposedRouting<4>(t,s),m);},test);
    else if(mode=="composed2")run<ComposedBank<2>,ComposedRouting<2>>("composed2",[](const auto& t,unsigned s,unsigned m){return RollingCode<2,false,false,0,false,ComposedRouting<2>>(ComposedRouting<2>(t,s),m);},test);
    else if(mode=="composed1")run<ComposedBank<1>,ComposedRouting<1>>("composed1",[](const auto& t,unsigned s,unsigned m){return RollingCode<1,false,false,0,false,ComposedRouting<1>>(ComposedRouting<1>(t,s),m);},test);
    else if(mode=="composed")run<ComposedBank<16>,ComposedRouting<16>>("composed",[](const auto& t,unsigned s,unsigned m){return RollingCode<16,false,false,0,false,ComposedRouting<16>>(ComposedRouting<16>(t,s),m);},test);
    else if(mode=="composed-epoch")run<ComposedBank<16>,ComposedRouting<16>>("composed-epoch",[](const auto& t,unsigned s,unsigned m){return RollingCode<16,false,false,0,true,ComposedRouting<16>>(ComposedRouting<16>(t,s),m);},test);
    else if(mode=="epoch")run("epoch",[](const auto& t,unsigned s,unsigned m){return RollingCode<16,false,true,0,true>(t,s,m);},test);
    else if(mode=="rolling")run("rolling",[](const auto& t,unsigned s,unsigned m){return RollingCode<16,false>(t,s,m);},test);
    else if(mode=="vector")run("vector",[](const auto& t,unsigned s,unsigned m){return RollingCode<16,false,true>(t,s,m);},test);
    else if(mode=="vector-split")run("vector-split",[](const auto& t,unsigned s,unsigned m){return RollingCode<16,true,true>(t,s,m);},test);
    else if(mode=="prefetch")run("prefetch",[](const auto& t,unsigned s,unsigned m){return RollingCode<16,false,true,32>(t,s,m);},test);
    else if(mode=="prefetch-split")run("prefetch-split",[](const auto& t,unsigned s,unsigned m){return RollingCode<16,true,true,32>(t,s,m);},test);
    else if(mode=="split")run("split",[](const auto& t,unsigned s,unsigned m){return RollingCode<16,true>(t,s,m);},test);
    else if(mode=="direct")run("direct",[](const auto& t,unsigned s,unsigned m){return BankFlow<16>::FreshCode(Routing<16>(t,s),m);},test);
    else throw std::invalid_argument("mode");
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
