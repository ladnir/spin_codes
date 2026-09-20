// One 128-bit element stream per call; setup and validation are not timed.
#include "Spin.h"
#include <algorithm>
#include <chrono>
#include <cctype>
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <fcntl.h>
#include <sched.h>
#include <sys/file.h>
#include <unistd.h>
using namespace bare_spin;
static volatile u64 sink;
static void lockBenchmarks() {
    for(const char* path:{"/tmp/prindal-addition-encoder-benchmark.lock","/tmp/bare-spin-benchmark.lock"}) {
        const int fd=open(path,O_CREAT|O_RDWR,0600);
        if(fd<0 || flock(fd,LOCK_EX|LOCK_NB)) throw std::runtime_error("benchmark lock unavailable");
    }
    for(const auto& e:std::filesystem::directory_iterator("/proc")) {
        const auto id=e.path().filename().string();
        if(id.empty() || !std::all_of(id.begin(),id.end(),[](unsigned char c){return std::isdigit(c);}) ||
           std::stoul(id)==static_cast<unsigned long>(getpid())) continue;
        std::error_code error;
        auto name=std::filesystem::read_symlink(e.path()/"exe",error).filename().string();
        if(error) continue;
        std::transform(name.begin(),name.end(),name.begin(),[](unsigned char c){return std::tolower(c);});
        if(name.find("bench")!=std::string::npos) throw std::runtime_error("another benchmark is active: "+name);
    }
    cpu_set_t cpus;CPU_ZERO(&cpus);CPU_SET(15,&cpus);
    if(sched_setaffinity(0,sizeof(cpus),&cpus)) throw std::runtime_error("cannot pin CPU 15");
}
int main(int argc,char** argv) {
 try {
    if(argc<3 || argc>5) throw std::invalid_argument("usage: single_benchmark m forward|transpose|inplace [odd trials=101] [tile=0]");
    const unsigned m=std::stoul(argv[1]),trials=argc>3?std::stoul(argv[3]):101,tile=argc>4?std::stoul(argv[4]):0;
    const std::string direction=argv[2];
    if((m!=16 && m!=18 && m!=20) || trials<3 || !(trials&1) ||
       (direction!="forward" && direction!="transpose" && direction!="inplace")) throw std::invalid_argument("invalid benchmark argument");
#if SPIN_SINGLE_TRANSPOSE
    if(direction=="forward") throw std::invalid_argument("standalone transpose target has no forward encoder");
#endif
    lockBenchmarks();
    Spin code(m==16?Configuration::T64S12R2:Configuration::T128S19,m,1,2,tile);
    const auto ni=direction=="forward"?code.messageBlocks():code.codeBlocks();
    const auto no=direction=="forward"?code.codeBlocks():code.messageBlocks();
    std::vector<block> input(ni),output(no),expected(no);u64 seed=123;
    for(auto& v:input) {const auto lo=splitmix(seed),hi=splitmix(seed);v=block(hi,lo);}
#if !SPIN_SINGLE_TRANSPOSE
    if(direction=="forward") code.forwardReference(input.data(),expected.data()); else
#endif
    code.reference(input.data(),expected.data());
    code.compact();Spin::Workspace work(code);
    auto run=[&] {
#if !SPIN_SINGLE_TRANSPOSE
        if(direction=="forward") {code.forwardUnchecked(input.data(),output.data(),work,Layout::Packed24);return;}
#endif
        code.encodeUnchecked(input.data(),direction=="inplace"?input.data():output.data(),work,Layout::Packed24);
    };
    run();
    if(std::memcmp(direction=="inplace"?input.data():output.data(),expected.data(),no*sizeof(block)))
        throw std::runtime_error("dense reference mismatch");
    // In-place mode intentionally continues evolving the buffer, matching the
    // established transpose harness. Out-of-place modes reuse a fixed input.
    for(unsigned i=0;i<3;++i) run();
    std::vector<double> samples(trials);
    for(unsigned i=0;i<trials;++i) {
        const auto start=std::chrono::steady_clock::now();run();
        samples[i]=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
        u64 word;std::memcpy(&word,(direction=="inplace"?input.data():output.data())+(i*997)%no,8);sink=word;
    }
    if(direction!="inplace" && std::memcmp(output.data(),expected.data(),no*sizeof(block)))
        throw std::runtime_error("timed output mismatch");
    auto sorted=samples;std::sort(sorted.begin(),sorted.end());
    std::cout<<std::setprecision(10)<<"{\"m\":"<<m<<",\"direction\":\""<<direction<<"\",\"configuration\":\""<<code.name()
        <<"\",\"tile_rows\":"<<code.tileBlocks()/256<<",\"median_ms\":"<<sorted[trials/2]
        <<",\"workspace_bytes\":"<<work.bytes()<<",\"setup_bytes\":"<<code.setupBytes()<<",\"samples_ms\":[";
    for(unsigned i=0;i<trials;++i) std::cout<<(i?",":"")<<samples[i];
    std::cout<<"]}\n";
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
