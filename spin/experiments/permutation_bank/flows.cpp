#include "Flows.h"
#include <spin/Code.h>
#include <atomic>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
using namespace spin::experimental::bank;
using Clock=std::chrono::steady_clock;
static double ms(Clock::time_point a,Clock::time_point b){return std::chrono::duration<double,std::milli>(b-a).count();}
static double median(std::vector<double> v){std::sort(v.begin(),v.end());return v[v.size()/2];}
static void input(std::vector<kernel::block>& x,unsigned seed) {
    Words words(seed);for(auto& v:x){auto a=words(),b=words();v=kernel::block(a,b);}
}
static void hash(const std::vector<kernel::block>& x,std::uint64_t& value) {
    for(const auto& v:x){std::uint64_t lanes[2];std::memcpy(lanes,&v,16);for(auto lane:lanes)value=(value^lane)*0x100000001b3ULL;}
}
template<unsigned Count> void verify() {
    BankFlow<Count> bank(913);std::vector<kernel::block> x(CodeSize),scratch(CodeSize);
    for(unsigned seed:{17U,191U,2029U}) {
        auto dynamic=bank.fresh(seed,seed+12);PrecomputedCode fixed(seed,seed+12,bank.materialize(seed));
        for(unsigned trial=0;trial<3;++trial) {
            input(x,seed*7+trial);if(!trial)std::fill(x.begin(),x.end(),kernel::block(0,0));
            auto expected=x;fixed.encode(expected.data());dynamic.encode(x.data(),scratch);
            if(std::memcmp(x.data(),expected.data(),CodeSize*sizeof(kernel::block)))throw std::runtime_error("precomputed/fresh full encoder mismatch");
        }
    }
}
static void output(const char* mode,unsigned count,double bankMs,double onceMs,double workspaceMs,std::size_t bankBytes,std::size_t planBytes,std::size_t workBytes,
                   const std::vector<double>& setup,const std::vector<double>& encode,const std::vector<double>& total,std::uint64_t digest) {
    std::cout<<mode<<','<<count<<','<<bankMs<<','<<onceMs<<','<<workspaceMs<<','<<bankBytes<<','<<planBytes<<','<<workBytes<<','<<median(setup)<<','<<median(encode)<<','<<median(total)<<','<<std::hex<<digest<<std::dec<<'\n';
}
template<unsigned Count> void fixed() {
    std::vector<kernel::block> x(CodeSize);
    double bankMs=0;std::size_t bankBytes=0;
    std::vector<std::uint32_t> route;
    Clock::time_point a;
    // Bank materialization is part of one-time code preprocessing, never
    // part of an encode. For exact shuffles, Spin constructs its own route.
    if constexpr(Count) {
        const auto start=Clock::now();BankFlow<Count> bank(913);const auto bankReady=Clock::now();
        bankMs=ms(start,bankReady);a=Clock::now();route=bank.materialize(17);
    }else a=Clock::now();
    PrecomputedCode code(17,29,std::move(route));const auto b=Clock::now();
    std::vector<double> setup,encoding,total;std::uint64_t digest=0;
    for(unsigned rep=0;rep<13;++rep) {
        input(x,913+rep);
        const auto c=Clock::now();code.encode(x.data());std::atomic_signal_fence(std::memory_order_seq_cst);const auto d=Clock::now();
        if(rep>=2){setup.push_back(0);encoding.push_back(ms(c,d));total.push_back(ms(c,d));}
        hash(x,digest);
    }
    output("precomputed",Count,bankMs,ms(a,b),0,bankBytes,code.planBytes(),code.workspaceBytes(),setup,encoding,total,digest);
}
template<unsigned Count> void fresh() {
    std::vector<kernel::block> x(CodeSize);
    const auto a=Clock::now();BankFlow<Count> bank(913);const auto b=Clock::now();
    std::vector<kernel::block> scratch(CodeSize);const auto c=Clock::now();
    std::vector<double> setup,encoding,total;std::uint64_t digest=0;std::size_t bytes=0;
    for(unsigned rep=0;rep<13;++rep) {
        input(x,913+rep);
        const auto start=Clock::now();
        auto code=std::make_optional(bank.fresh(17+rep,29+rep));const auto ready=Clock::now();
        code->encode(x.data(),scratch);std::atomic_signal_fence(std::memory_order_seq_cst);const auto done=Clock::now();
        bytes=code->bytes();code.reset();const auto end=Clock::now();
        if(rep>=2){setup.push_back(ms(start,ready));encoding.push_back(ms(ready,done));total.push_back(ms(start,end));}
        hash(x,digest);
    }
    output("fresh",Count,ms(a,b),0,ms(b,c),bank.bytes(),bytes,CodeSize*sizeof(kernel::block),setup,encoding,total,digest);
}
int main(int argc,char** argv){try {
    if(!spin::capabilities().avx2)return 77;
    const std::string mode=argc>1?argv[1]:"verify";
    if(mode=="verify"){verify<16>();verify<64>();std::cout<<"full precomputed/fresh bank encoder references PASS\n";return 0;}
    std::cout<<std::fixed<<std::setprecision(6)<<"flow,bank_size,bank_once_ms,code_and_workspace_once_ms,workspace_once_ms,bank_bytes,plan_bytes,workspace_bytes,fresh_setup_ms,encode_ms,total_ms,checksum\n";
    if(mode=="fixed-exact")fixed<0>();else if(mode=="fixed16"){verify<16>();fixed<16>();}else if(mode=="fixed64"){verify<64>();fixed<64>();}
    else if(mode=="fresh16"){verify<16>();fresh<16>();}else if(mode=="fresh64"){verify<64>();fresh<64>();}
    else throw std::invalid_argument("flow");
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
