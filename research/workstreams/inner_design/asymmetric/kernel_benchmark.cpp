#include "Spin.h"
#include "Inner.h"
#if ASYM_BASELINE
#include "BalancedMap.h"
#include "../MixerInner.h"
using SelectedMap=bare_spin::BalancedMap;
#else
#include "AsymmetricMap.h"
#include "AsymmetricInner.h"
using SelectedMap=bare_spin::AsymmetricMap;
#endif
#include <algorithm>
#include <chrono>
#include <cctype>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <fcntl.h>
#include <sched.h>
#include <sys/file.h>
#include <unistd.h>
using namespace bare_spin;
static volatile u64 sink=0;
static u64 randomWord(u64& state) {
    auto v=(state+=0x9e3779b97f4a7c15ULL);
    v=(v^(v>>30))*0xbf58476d1ce4e5b9ULL;
    v=(v^(v>>27))*0x94d049bb133111ebULL;
    return v^(v>>31);
}
static void guard() {
    for(const auto& e:std::filesystem::directory_iterator("/proc")) {
        auto name=e.path().filename().string();
        if(name.empty() || !std::all_of(name.begin(),name.end(),[](unsigned char c){return std::isdigit(c);}) ||
           std::stoul(name)==static_cast<unsigned long>(getpid())) continue;
        std::error_code error;
        auto file=std::filesystem::read_symlink(e.path()/"exe",error).filename().string();
        if(error) continue;
        std::transform(file.begin(),file.end(),file.begin(),[](unsigned char c){return char(std::tolower(c));});
        if(file.find("bench")!=std::string::npos) throw std::runtime_error("another benchmark is active: "+file);
    }
}
// A real store callback retains the streaming cost but excludes routing/outer work.
__attribute__((noinline)) static void run(const block* input,block* output,std::size_t n,const u32* masks) {
    auto emit=[output](std::size_t i,block v){output[i]=v;};
#if ASYM_BASELINE
    mixerReverse<SelectedMap,true,true>(input,n,masks,emit);
#else
    asymmetricReverse<SelectedMap,true,false>(input,n,masks,emit);
#endif
}
static void check(const block* input,const block* output,std::size_t n,const u32* rawMasks) {
    block state[19]{},next[19]{};
    for(std::size_t e=n/128;e-->0;) {
        for(unsigned p=0;p<128;++p) {
            auto v=input[e*128+p];
            for(unsigned j=0;j<19;++j) if((SelectedMap::columns[p]>>j)&1) v^=state[j];
            if(v!=output[e*128+p]) throw std::runtime_error("inner oracle mismatch");
        }
        block dot{};
        for(unsigned j=0;j<19;++j) if((rawMasks[2*e]>>j)&1) dot^=state[j];
        for(unsigned j=0;j<19;++j) {
            next[j]=state[j];
            if((rawMasks[2*e+1]>>j)&1) next[j]^=dot;
            for(unsigned p=0;p<128;++p) {
#if ASYM_BASELINE
                const auto column=SelectedMap::columns[p];
#else
                const auto column=SelectedMap::feedbackColumns[p];
#endif
                if((column>>j)&1) next[j]^=input[e*128+p];
            }
        }
        std::copy(std::begin(next),std::end(next),state);
    }
}
int main() {
    try {
        const int lock=open("/tmp/bare-spin-benchmark.lock",O_CREAT|O_RDWR,0600);
        if(lock<0 || flock(lock,LOCK_EX|LOCK_NB)) throw std::runtime_error("benchmark lock unavailable");
        cpu_set_t cpus;CPU_ZERO(&cpus);CPU_SET(15,&cpus);
        if(sched_setaffinity(0,sizeof(cpus),&cpus)) throw std::runtime_error("cannot pin CPU 15");
        guard();
        for(std::size_t n:{std::size_t(32768),std::size_t(4194304)}) {
            u64 seed=2;std::vector<u32> masks(n/64),rawMasks(n/64);
            for(std::size_t e=0;e<n/128;++e) {
                u32 u;do u=u32(randomWord(seed))&((1U<<19)-1);while(!u);
                auto v=u32(randomWord(seed))&((1U<<19)-1);
                if(std::popcount(u&v)&1) v^=u&-u;
                rawMasks[2*e]=u;rawMasks[2*e+1]=v;
                u32 grouped=0;
                for(unsigned j=0;j<19;++j) grouped|=((u>>SelectedMap::groupOrder[j])&1)<<j;
                masks[2*e]=grouped;masks[2*e+1]=v;
            }
            std::vector<block> input(n),output(n);seed=123;
            for(auto& v:input) {auto lo=randomWord(seed);auto hi=randomWord(seed);v=block(hi,lo);}
            run(input.data(),output.data(),n,masks.data());
            // Full dense oracle on the smaller diagnostic; full-encoder tests
            // separately validate the large workload and fused routing path.
            if(n==32768) check(input.data(),output.data(),n,rawMasks.data());
            for(unsigned i=0;i<3;++i) run(input.data(),output.data(),n,masks.data());
            std::vector<double> samples;samples.reserve(101);
            for(unsigned i=0;i<101;++i) {
                auto begin=std::chrono::steady_clock::now();
                run(input.data(),output.data(),n,masks.data());
                samples.push_back(std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-begin).count());
                sink=output[(997*i)%n].get<u64>()[0];
            }
            std::sort(samples.begin(),samples.end());
            u64 hash=0xcbf29ce484222325ULL;
            for(auto v:output) for(auto x:v.get<u64>()) {hash^=x;hash*=0x100000001b3ULL;}
            std::cout<<std::setprecision(10)<<"{\"variant\":\""<<ASYM_NAME<<"\",\"blocks\":"<<n
                <<",\"trials\":101,\"median_ms\":"<<samples[50]<<",\"p10_ms\":"<<samples[10]
                <<",\"p90_ms\":"<<samples[90]<<",\"output_hash\":\""<<std::hex<<hash<<std::dec<<"\"}\n"<<std::flush;
        }
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
