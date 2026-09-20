#include "Wide.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <fcntl.h>
#include <sched.h>
#include <sys/file.h>
#include <unistd.h>
using namespace bare_spin;
extern "C" void hc_pcs_assemble512_stream(const void*,void*,std::size_t,std::size_t,std::size_t,std::size_t);
extern "C" void spin_columns256(const block*,block*,std::size_t);
extern "C" void spin_columns512(const block*,block*,std::size_t);
static void lockBenchmarks() {
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
struct Buffer {
    block* p=nullptr;
    explicit Buffer(std::size_t n) {
        if(posix_memalign(reinterpret_cast<void**>(&p),64,n*16)) throw std::bad_alloc();
        // Fault in pages before timing, including output and scratch.
        std::memset(p,0,n*16);
    }
    ~Buffer() {std::free(p);}
    Buffer(const Buffer&)=delete;
    Buffer& operator=(const Buffer&)=delete;
};
template<unsigned W> __attribute__((noinline)) static void pack(const block* in,block* out,std::size_t k) {
    for(std::size_t j=0;j<k;++j) for(unsigned l=0;l<W;++l) out[j*W+l]=in[l*k+j];
}
using Clock=std::chrono::steady_clock;
static double ms(Clock::time_point a,Clock::time_point b) {return std::chrono::duration<double,std::milli>(b-a).count();}
template<unsigned W> static void run(unsigned m,int trial,bool native,const Spin& code,Spin::Workspace& base,
    WideWorkspace* wide,const block* input,const block* prepared,block* scratch,block* encoded,block* out,const block* reference) {
    const auto k=code.messageBlocks(),n=code.codeBlocks();
    double packing=0,encoding=0;
    const auto begin=Clock::now();
    for(unsigned g=0;g<16;g+=W) {
        const block* src=native?prepared+g*k:input+g*k;
        if constexpr(W>1) {
            if(!native) {auto t=Clock::now();pack<W>(src,scratch,k);packing+=ms(t,Clock::now());src=scratch;}
        }
        auto t=Clock::now();
        if constexpr(W==1) code.forward(src,k,encoded+g*n,n,base);
        else wide->forward({src,W*k},{encoded+g*n,W*n});
        encoding+=ms(t,Clock::now());
    }
    const auto columnStart=Clock::now();
    if constexpr(W==1) hc_pcs_assemble512_stream(encoded,out,n,16,16,n);
    if constexpr(W==2) spin_columns256(encoded,out,n);
    if constexpr(W==4) spin_columns512(encoded,out,n);
    const auto end=Clock::now();
    if(std::memcmp(out,reference,16*n*16)) throw std::runtime_error("canonical column output mismatch");
    if(trial>=0) std::cout<<m<<','<<trial<<','<<W<<','<<(native?"native":"planar")<<','<<packing<<','<<encoding<<','
        <<ms(columnStart,end)<<','<<ms(begin,end)<<','<<(W==1?base.bytes():wide->bytes())<<','
        <<code.name()<<','<<code.tileBlocks()/256<<'\n';
}
int main(int argc,char** argv) {
 try {
    const unsigned m=argc>1?std::stoul(argv[1]):16;
    const int trials=argc>2?std::stoi(argv[2]):7;
    const std::string config=argc>3?argv[3]:"12819";
    const unsigned tile=argc>4?std::stoul(argv[4]):0,only=argc>5?std::stoul(argv[5]):0;
    if(m<14 || m>20 || trials<1 || (config!="12819" && config!="6412r2") ||
       (only!=0 && only!=1 && only!=2 && only!=4))
        throw std::invalid_argument("usage: spin_wide_benchmark m[14..20] trials [12819|6412r2] [tile=0] [lanes=0(all)|1|2|4]");
    if(!wideAvailable(4)) {std::cerr<<"Column benchmark requires AVX-512\n";return 77;}
    lockBenchmarks();
    Spin code(config=="6412r2"?Configuration::T64S12R2:Configuration::T128S19,m,1,2,tile);code.compact();
    Spin::Workspace base(code);WideWorkspace w2(code,2),w4(code,4);
    const auto k=code.messageBlocks(),n=code.codeBlocks();
    Buffer input(16*k),prepared2(16*k),prepared4(16*k),scratch(4*k),encoded(16*n),output(16*n),reference(16*n);
    u64 seed=123;
    for(std::size_t i=0;i<16*k;++i) {auto lo=splitmix(seed),hi=splitmix(seed);input.p[i]=block(hi,lo);}
    for(unsigned g=0;g<16;g+=2) pack<2>(input.p+g*k,prepared2.p+g*k,k);
    for(unsigned g=0;g<16;g+=4) pack<4>(input.p+g*k,prepared4.p+g*k,k);
    // Reference assembly is scalar and independent of all three timed kernels.
    for(unsigned g=0;g<16;++g) {
        code.forward(input.p+g*k,k,encoded.p,n,base);
        for(std::size_t j=0;j<n;++j) reference.p[j*16+g]=encoded.p[j];
    }
    std::cout<<std::setprecision(10)<<"m,trial,lanes,input_layout,pack_ms,encode_ms,column_ms,total_ms,workspace_bytes,configuration,tile_rows\n";
    for(int t=-1;t<trials;++t) {
        // Alternate order; one untimed warmup of every case. All allocations and
        // input generation/prepacking are outside trials; planar mode pays repacking.
        std::array<unsigned,5> order={0,1,2,3,4};if(t&1) std::reverse(order.begin(),order.end());
        for(auto c:order) {
            const unsigned width=c==0?1:c<3?2:4;
            if(only && only!=width) continue;
            if(c==0) run<1>(m,t,true,code,base,nullptr,input.p,input.p,scratch.p,encoded.p,output.p,reference.p);
            if(c==1 || c==2) run<2>(m,t,c==1,code,base,&w2,input.p,prepared2.p,scratch.p,encoded.p,output.p,reference.p);
            if(c==3 || c==4) run<4>(m,t,c==3,code,base,&w4,input.p,prepared4.p,scratch.p,encoded.p,output.p,reference.p);
        }
    }
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
