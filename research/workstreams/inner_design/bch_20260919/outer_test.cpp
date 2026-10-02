#include "Spin.h"
#include "generated/BchCircuit.h"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <filesystem>
#include <cctype>
#include <fcntl.h>
#include <sys/file.h>
#include <sched.h>
#include <unistd.h>
using namespace bare_spin;
namespace bare_spin {void bchTranspose4(const block*,block*);}
#ifndef BCH_FOUR
#define BCH_FOUR 0
#endif
#ifndef BCH_PACKED
#define BCH_PACKED 0
#endif
constexpr unsigned lanes=BCH_FOUR?4:2;
static unsigned position(unsigned row,unsigned col) {return BCH_PACKED?4*col+row:256*row+col;}
static void invoke(const block* in,block* out) {
#if BCH_FOUR
    bchTranspose4(in,out);
#else
    bchTranspose2(in,in+256,out,out+128);
#endif
}
static void require(bool b) {if(!b) throw std::runtime_error("BCH scalar oracle mismatch");}
static u64 randomWord(u64& seed) {seed^=seed<<13;seed^=seed>>7;seed^=seed<<17;return seed;}
static void tests() {
    std::vector<block> in(256*lanes),out(128*lanes+2),expected(out.size());
    const block guard(123,456);
    out.front()=out.back()=expected.front()=expected.back()=guard;
    // Exhaust the binary matrix basis, separately in every SIMD lane.
    for(unsigned row=0;row<lanes;++row) for(unsigned col=0;col<256;++col) {
        std::fill(in.begin(),in.end(),block(0,0));
        std::fill(expected.begin()+1,expected.end()-1,block(0,0));
        const block value(0xa35c9dULL,0x1234567ULL);
        in[position(row,col)]=value;
        for(unsigned j=0;j<128;++j)
            if((BchRows[j][col/64]>>(col%64))&1) expected[1+row*128+j]=value;
        invoke(in.data(),out.data()+1);require(out==expected);
    }
    u64 seed=713;
    for(unsigned trial=0;trial<12;++trial) {
        for(auto& x:in) {const auto lo=randomWord(seed),hi=randomWord(seed);x=block(hi,lo);}
        for(unsigned row=0;row<lanes;++row) for(unsigned j=0;j<128;++j) {
            block sum(0,0);
            for(unsigned col=0;col<256;++col)
                if((BchRows[j][col/64]>>(col%64))&1) sum^=in[position(row,col)];
            expected[1+row*128+j]=sum;
        }
        invoke(in.data(),out.data()+1);require(out==expected);
    }
    std::cout<<"BCH basis, all SIMD lanes, random dense oracle, canaries PASS\n";
}
static void guardBenchmark() {
    for(const char* path:{"/tmp/prindal-addition-encoder-benchmark.lock","/tmp/bare-spin-benchmark.lock"}) {
        int fd=open(path,O_CREAT|O_RDWR,0600);
        if(fd<0 || flock(fd,LOCK_EX|LOCK_NB)) throw std::runtime_error("benchmark lock unavailable");
    }
    for(const auto& e:std::filesystem::directory_iterator("/proc")) {
        const auto name=e.path().filename().string();
        if(name.empty() || !std::all_of(name.begin(),name.end(),[](unsigned char c){return std::isdigit(c);}) ||
           std::stoul(name)==static_cast<unsigned long>(getpid())) continue;
        std::error_code error;
        auto file=std::filesystem::read_symlink(e.path()/"exe",error).filename().string();
        if(error) continue;
        std::transform(file.begin(),file.end(),file.begin(),[](unsigned char c){return std::tolower(c);});
        if(file.find("bench")!=std::string::npos) throw std::runtime_error("another benchmark is active");
    }
    cpu_set_t cpus;CPU_ZERO(&cpus);CPU_SET(15,&cpus);
    if(sched_setaffinity(0,sizeof(cpus),&cpus)) throw std::runtime_error("affinity failed");
}
int main(int argc,char** argv) {
    try {
        if(argc==1) {tests();return 0;}
        guardBenchmark();
        const bool hot=std::string(argv[1])=="hot";
        constexpr unsigned rows=8192,trials=31;
        std::vector<block> in(hot?256*lanes:rows*256),out(hot?128*lanes:rows*128);
        u64 seed=123;
        for(auto& x:in) {const auto lo=randomWord(seed),hi=randomWord(seed);x=block(hi,lo);}
        std::vector<double> samples;
        volatile u64 sink=0;
        for(unsigned trial=0;trial<trials+3;++trial) {
            const auto start=std::chrono::steady_clock::now();
            for(unsigned row=0;row<rows;row+=lanes)
                invoke(in.data()+(hot?0:row*256),out.data()+(hot?0:row*128));
            const double elapsed=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
            if(trial>=3) samples.push_back(elapsed);
            sink=out[trial%out.size()].get<u64>()[0];
        }
        auto sorted=samples;std::sort(sorted.begin(),sorted.end());
        std::cout<<std::setprecision(10)<<"{\"hot\":"<<(hot?"true":"false")<<",\"rows\":"<<rows
                 <<",\"median_ms\":"<<sorted[trials/2]<<",\"samples_ms\":[";
        for(unsigned i=0;i<samples.size();++i) std::cout<<(i?",":"")<<samples[i];
        std::cout<<"]}\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
