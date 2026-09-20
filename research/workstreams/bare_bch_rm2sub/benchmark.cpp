#include "Spin.h"
#include <algorithm>
#include <chrono>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef __linux__
#include <sched.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/file.h>
#endif
using namespace bare_spin;
static volatile u64 sink=0;
static void guard() {
#ifdef __linux__
    for(const auto& e:std::filesystem::directory_iterator("/proc")) {
        const auto name=e.path().filename().string();
        if(name.empty() || !std::all_of(name.begin(),name.end(),[](unsigned char c){return std::isdigit(c);}) ||
           std::stoul(name)==static_cast<unsigned long>(getpid())) continue;
        std::error_code error;
        auto file=std::filesystem::read_symlink(e.path()/"exe",error).filename().string();
        if(error) continue;
        std::transform(file.begin(),file.end(),file.begin(),[](unsigned char c){return char(std::tolower(c));});
        if(file.find("bench")!=std::string::npos) throw std::runtime_error("another benchmark is active: "+file);
    }
#endif
}
static double elapsed(std::chrono::steady_clock::time_point begin) {
    return std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-begin).count();
}
int main(int argc,char** argv) {
    try {
        const unsigned trials=argc>1?std::stoul(argv[1]):31;
        const unsigned tile=argc>2?std::stoul(argv[2]):0;
        const unsigned layout=argc>3?std::stoul(argv[3]):0;
        const unsigned onlyM=argc>4?std::stoul(argv[4]):0;
        const unsigned quarter=argc>5?std::stoul(argv[5]):0;
        const unsigned onlyC=argc>6?std::stoul(argv[6]):4;
        const unsigned inplace=argc>7?std::stoul(argv[7]):0;
        if(trials<3 || !(trials&1) || layout>1 || quarter>1 || onlyC>4 || inplace>1 || (onlyM && onlyM!=16 && onlyM!=18 && onlyM!=20))
            throw std::invalid_argument("usage: benchmark odd_trials>=3 tile_rows layout[0,1] [m=16,18,20] [quarter=0,1] [config=0..3,4=all] [inplace=0,1]");
#ifdef __linux__
        const int lock=open("/tmp/bare-spin-benchmark.lock",O_CREAT|O_RDWR,0600);
        if(lock<0 || flock(lock,LOCK_EX|LOCK_NB)) throw std::runtime_error("benchmark lock unavailable");
        cpu_set_t cpus; CPU_ZERO(&cpus);CPU_SET(15,&cpus);
        if(sched_setaffinity(0,sizeof(cpus),&cpus)) throw std::runtime_error("cannot pin CPU 15");
#endif
        guard();
        for(unsigned m:{16U,18U,20U}) {
            if(onlyM && m!=onlyM) continue;
            for(unsigned c=0;c<4;++c) {
                if((quarter && c!=2) || (onlyC!=4 && onlyC!=c)) continue;
                guard();
                auto begin=std::chrono::steady_clock::now();
                Spin code(static_cast<Configuration>(c),m,1,2,tile,quarter?Outer::Bch128x32:Outer::Bch256x128);
                const double setup=elapsed(begin);
                const auto diagnosticBytes=code.setupBytes();
                const auto route=layout?Layout::Indices32:Layout::Packed24;
                code.compact(route);
                Spin::Workspace w(code);
                std::vector<block> in(code.codeBlocks()),out(inplace?0:code.messageBlocks());
                block* output=inplace?in.data():out.data();
                u64 seed=123;
                for(auto& v:in) {auto lo=splitmix(seed);auto hi=splitmix(seed);v=block(hi,lo);}
                for(unsigned i=0;i<3;++i) code.encodeUnchecked(in.data(),output,w,route);
                std::vector<double> samples; samples.reserve(trials);
                for(unsigned i=0;i<trials;++i) {
                    begin=std::chrono::steady_clock::now();
                    code.encodeUnchecked(in.data(),output,w,route);
                    samples.push_back(elapsed(begin));
                    sink=output[(i*997)%code.messageBlocks()].get<u64>()[0];
                }
                std::sort(samples.begin(),samples.end());
                u64 hash=0xcbf29ce484222325ULL;
                for(std::size_t j=0;j<code.messageBlocks();++j) for(auto x:output[j].get<u64>()) {hash^=x;hash*=0x100000001b3ULL;}
                std::cout<<std::setprecision(10)<<"{\"configuration\":\""<<code.name()<<"\",\"m\":"<<m
                    <<",\"outer_length\":"<<code.outerLength()<<",\"outer_dimension\":"<<code.outerDimension()
                    <<",\"inplace\":"<<(inplace?"true":"false")
                    <<",\"tile_rows\":"<<code.tileBlocks()/code.outerLength()<<",\"layout\":\""<<(layout?"indices32":"packed24")
                    <<"\",\"trials\":"<<trials<<",\"median_ms\":"<<samples[trials/2]
                    <<",\"p10_ms\":"<<samples[trials/10]<<",\"p90_ms\":"<<samples[(trials*9)/10]
                    <<",\"setup_ms\":"<<setup<<",\"setup_bytes_with_oracles\":"<<diagnosticBytes
                    <<",\"retained_setup_bytes\":"<<code.setupBytes()
                    <<",\"workspace_bytes\":"<<w.bytes()<<",\"output_hash\":\""<<std::hex<<hash<<std::dec<<"\"}\n"<<std::flush;
            }
        }
    } catch(const std::exception& e) {std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}
}
