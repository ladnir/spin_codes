#include "Spin.h"
#include "QuarterCircuit.h"
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cctype>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <numeric>
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
static u64 below(u64& seed,u64 bound) {
    const auto threshold=(0-bound)%bound;
    for(;;) {const auto v=randomWord(seed);if(v>=threshold) return v%bound;}
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
struct Cascade {
    std::size_t n;
    unsigned rounds;
    bool bucketed;
    std::vector<std::vector<u32>> permutations;
    std::vector<std::vector<u32>> slots,offsets;
    std::vector<block> a,b;
    Cascade(std::size_t size,unsigned r,bool bucket=false):n(size),rounds(r),bucketed(bucket),permutations(r),slots(bucket?r:0),offsets(bucket?r:0),a(n),b(n) {
        if(n%256 || (r!=2 && r!=3)) throw std::invalid_argument("cascade geometry");
        u64 seed=12345;
        for(auto& p:permutations) {
            p.resize(n);std::iota(p.begin(),p.end(),u32(0));
            for(std::size_t i=n;i>1;--i) std::swap(p[i-1],p[below(seed,i)]);
        }
        if(bucketed) for(unsigned r=0;r<rounds;++r) {
            const auto tile=std::min<std::size_t>(n,1<<19);
            slots[r].resize(n);offsets[r].resize(n);
            std::vector<u32> counts(n/tile,u32(tile));
            for(std::size_t i=n;i-->0;) {
                const auto outer=permutations[r][i],bucket=outer/u32(tile);
                const auto slot=bucket*u32(tile)+--counts[bucket];
                slots[r][i]=slot;offsets[r][slot]=outer&u32(tile-1);
            }
        }
    }
    // Public precomputed indices; fixed eight-way batches and no hot allocation.
    static void stage(const block* in,block* out,const u32* p,std::size_t n) {
        block sum(0,0);
        for(std::size_t end=n;end;end-=8) {
            const auto i=end-8;
            const block v7=sum^in[i+7],v6=v7^in[i+6],v5=v6^in[i+5],v4=v5^in[i+4];
            const block v3=v4^in[i+3],v2=v3^in[i+2],v1=v2^in[i+1],v0=v1^in[i];
            sum=v0;
            out[p[i+7]]=v7;out[p[i+6]]=v6;out[p[i+5]]=v5;out[p[i+4]]=v4;
            out[p[i+3]]=v3;out[p[i+2]]=v2;out[p[i+1]]=v1;out[p[i]]=v0;
        }
    }
    void encode(const block* input,block* output) {
        const block* current=input;block* destination=a.data();
        for(unsigned r=rounds;r-->0;) {
            if(bucketed) {
                stage(current,a.data(),slots[r].data(),n);
                const auto tile=std::min<std::size_t>(n,1<<19);
                for(std::size_t base=0;base<n;base+=tile) for(std::size_t j=0;j<tile;++j) {
                    if(j+32<tile) _mm_prefetch(reinterpret_cast<const char*>(b.data()+base+offsets[r][base+j+32]),_MM_HINT_T0);
                    b[base+offsets[r][base+j]]=a[base+j];
                }
                current=b.data();
            } else {
                stage(current,destination,permutations[r].data(),n);
                current=destination;destination=destination==a.data()?b.data():a.data();
            }
        }
        for(std::size_t j=0;j<n;j+=256) quarterTranspose2(current+j,current+j+128,output+j/4,output+j/4+32);
    }
    std::vector<block> forward(const std::vector<block>& input) const {
        std::vector<block> current(n),next(n);
        for(std::size_t row=0;row<n/128;++row) for(unsigned c=0;c<128;++c) {
            auto value=block(0,0);
            for(unsigned j=0;j<32;++j) if((QuarterRows[j][c/64]>>(c%64))&1) value^=input[row*32+j];
            current[row*128+c]=value;
        }
        for(unsigned r=0;r<rounds;++r) {
            block sum(0,0);
            for(std::size_t i=0;i<n;++i) {sum^=current[permutations[r][i]];next[i]=sum;}
            current.swap(next);
        }
        return current;
    }
};
static unsigned dot(const std::vector<block>& a,const std::vector<block>& b) {
    unsigned v=0;
    for(std::size_t i=0;i<a.size();++i) {
        const auto x=a[i].get<u64>(),y=b[i].get<u64>();
        v^=std::popcount(x[0]&y[0])^std::popcount(x[1]&y[1]);
    }
    return v&1;
}
static void fill(std::vector<block>& v,u64 seed) {
    for(auto& x:v) {const auto lo=randomWord(seed),hi=randomWord(seed);x=block(hi,lo);}
}
static void test(unsigned r,bool bucketed) {
    Cascade code(1024,r,bucketed);
    std::vector<block> input(1024),message(256),output(256);fill(input,71);fill(message,79);
    code.encode(input.data(),output.data());
    if(dot(code.forward(message),input)!=dot(message,output)) throw std::runtime_error("cascade adjoint failed");
    // All message basis coordinates, each in one bit lane, against a dense input.
    for(unsigned i=0;i<256;++i) {
        std::fill(message.begin(),message.end(),block(0,0));message[i]=block(0,1);
        if(dot(code.forward(message),input)!=(output[i].get<u64>()[0]&1)) throw std::runtime_error("cascade basis adjoint failed");
    }
    auto inplace=input;code.encode(inplace.data(),inplace.data());
    if(!std::equal(output.begin(),output.end(),inplace.begin()) || !std::equal(input.begin()+256,input.end(),inplace.begin()+256))
        throw std::runtime_error("cascade inplace/suffix failed");
    std::fill(input.begin(),input.end(),block(0,0));code.encode(input.data(),output.data());
    if(std::any_of(output.begin(),output.end(),[](block x){return x!=block(0,0);})) throw std::runtime_error("cascade zero failed");
}
int main(int argc,char** argv) {
    try {
        const unsigned rounds=argc>1?std::stoul(argv[1]):2;
        const bool bucketed=argc>2?std::stoul(argv[2])!=0:false;
        const int lock=open("/tmp/bare-spin-benchmark.lock",O_CREAT|O_RDWR,0600);
        if(lock<0 || flock(lock,LOCK_EX|LOCK_NB)) throw std::runtime_error("benchmark lock unavailable");
        cpu_set_t cpus;CPU_ZERO(&cpus);CPU_SET(15,&cpus);
        if(sched_setaffinity(0,sizeof(cpus),&cpus)) throw std::runtime_error("cannot pin CPU 15");
        guard();test(rounds,bucketed);
        constexpr std::size_t n=1<<22;Cascade code(n,rounds,bucketed);std::vector<block> input(n);fill(input,123);
        if(bucketed) for(auto& p:code.permutations) std::vector<u32>().swap(p);
        for(unsigned i=0;i<3;++i) code.encode(input.data(),input.data());
        std::vector<double> samples;samples.reserve(101);
        for(unsigned i=0;i<101;++i) {
            const auto begin=std::chrono::steady_clock::now();code.encode(input.data(),input.data());
            samples.push_back(std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-begin).count());
            sink=input[(997*i)%(n/4)].get<u64>()[0];
        }
        std::sort(samples.begin(),samples.end());u64 hash=0xcbf29ce484222325ULL;
        for(std::size_t i=0;i<n/4;++i) for(auto x:input[i].get<u64>()) {hash^=x;hash*=0x100000001b3ULL;}
        std::cout<<std::setprecision(10)<<"{\"configuration\":\"BCH128x32_accumulator_NOT_CERTIFIED\",\"rounds\":"<<rounds
            <<",\"m\":20,\"inplace\":true,\"trials\":101,\"median_ms\":"<<samples[50]
            <<",\"bucketed\":"<<(bucketed?"true":"false")
            <<",\"p10_ms\":"<<samples[10]<<",\"p90_ms\":"<<samples[90]
            <<",\"setup_bytes\":"<<(bucketed?2:1)*rounds*n*sizeof(u32)<<",\"workspace_bytes\":"<<2*n*sizeof(block)
            <<",\"correctness\":\"PASS\",\"output_hash\":\""<<std::hex<<hash<<std::dec<<"\"}\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
