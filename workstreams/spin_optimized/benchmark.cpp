#include "Spin.h"
#include <algorithm>
#include <chrono>
#include <array>
#include <bit>
#include <cstring>
#include <cctype>
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
static std::array<u64,2> words(const block& v) {std::array<u64,2> a;std::memcpy(a.data(),&v,16);return a;}
static void lockBenchmarks() {
    // Held until process exit. Acquire in the same order as companion experiments.
    for(const char* path:{"/tmp/prindal-addition-encoder-benchmark.lock","/tmp/bare-spin-benchmark.lock"}) {
        const int fd=open(path,O_CREAT|O_RDWR,0600);
        if(fd<0 || flock(fd,LOCK_EX|LOCK_NB)) throw std::runtime_error("benchmark lock unavailable");
    }
    for(const auto& e:std::filesystem::directory_iterator("/proc")) {
        auto id=e.path().filename().string();
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
    if(argc<3) throw std::invalid_argument("usage: spin_benchmark m backend[auto,avx2,avx512] [odd trials=101] [seed=1] [tile=0] [layout=0] [configuration=12819|6412|6412r2] [K override, generalized transpose only]");
    const unsigned m=std::stoul(argv[1]),trials=argc>3?std::stoul(argv[3]):101;
    const std::string which=argv[2];const u64 routeSeed=argc>4?std::stoull(argv[4]):1;
    const unsigned tile=argc>5?std::stoul(argv[5]):0,layout=argc>6?std::stoul(argv[6]):0;
    const std::string configuration=argc>7?argv[7]:"12819";
#if !SPIN_GENERAL_LENGTHS
    if(argc>8) throw std::invalid_argument("K override requires the generalized transpose target");
#endif
    if(configuration!="12819" && configuration!="6412" && configuration!="6412r2") throw std::invalid_argument("configuration must be 12819, 6412, or 6412r2");
    if((m!=16 && m!=18 && m!=20) || trials<3 || !(trials&1) ||
       layout>1 || (which!="auto" && which!="avx2" && which!="avx512")) throw std::invalid_argument("invalid benchmark argument");
    lockBenchmarks();
    const auto backend=which=="auto"?BchBackend::Auto:which=="avx2"?BchBackend::Avx2:BchBackend::Avx512;
    Spin code(configuration=="6412r2"?Configuration::T64S12R2:configuration=="6412"?Configuration::T64S12:Configuration::T128S19,
#if SPIN_GENERAL_LENGTHS
        MessageLength{argc>8?std::stoull(argv[8]):std::size_t{1}<<m},
#else
        m,
#endif
        routeSeed,2,tile,
#ifndef SPIN_BIDIRECTIONAL_BENCH
        Outer::Bch256x128,
#endif
        backend);
    const auto route=layout?Layout::Indices32:Layout::Packed24;
    code.compact(route);Spin::Workspace work(code);
    std::vector<block> data(code.codeBlocks());u64 seed=123;
    for(auto& v:data) {const auto lo=splitmix(seed),hi=splitmix(seed);v=block(hi,lo);}
    for(unsigned i=0;i<3;++i) code.encodeUnchecked(data.data(),data.data(),work,route);
    std::vector<double> samples(trials);
    for(unsigned i=0;i<trials;++i) {
        const auto start=std::chrono::steady_clock::now();
        code.encodeUnchecked(data.data(),data.data(),work,route);
        samples[i]=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
        sink=words(data[(i*997)%code.messageBlocks()])[0];
    }
    auto sorted=samples;std::sort(sorted.begin(),sorted.end());u64 hash=0xcbf29ce484222325ULL;
    for(std::size_t j=0;j<code.messageBlocks();++j) for(auto v:words(data[j])) {hash^=v;hash*=0x100000001b3ULL;}
    std::cout<<std::setprecision(10)<<"{\"configuration\":\""<<code.name()<<"\",\"m\":";
    if(std::has_single_bit(code.messageBlocks())) std::cout<<std::countr_zero(code.messageBlocks());
    else std::cout<<"null";
    std::cout<<",\"backend\":\""<<(code.bchBackend()==BchBackend::Avx512?"avx512":"avx2")
        <<"\",\"route_seed\":"<<routeSeed<<",\"tile_rows\":"<<code.tileBlocks()/256<<",\"layout\":"<<layout<<",\"median_ms\":"<<sorted[trials/2]
        <<",\"K\":"<<code.messageBlocks()<<",\"retained_setup_bytes\":"<<code.setupBytes()<<",\"workspace_bytes\":"<<work.bytes()
        <<",\"output_hash\":\""<<std::hex<<hash<<std::dec<<"\",\"samples_ms\":[";
    for(unsigned i=0;i<trials;++i) std::cout<<(i?",":"")<<samples[i];
    std::cout<<"]}\n";
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
