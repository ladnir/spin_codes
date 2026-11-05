// Same code and same sampled route as the retained near-half implementation.
// No feedback work is timed here. Research diagnostics, not a library API.
#define SPIN_OR_WINDOW_NINE 1
#define SPIN_OR_WINDOW_COORDINATES 1
#define SPIN_OR_WINDOW_RAW_SCALE 1
#define main retained_near_half_main
#include "near_half_native.cpp"
#undef main
#include "RsOuterOptimized.h"
#include <sys/mman.h>
#include <fstream>

namespace ro {
using namespace native;
namespace opt=spin::research::one_route::rs_outer;
using Fn=void(*)(const Plan&,const E*,E*);
template<unsigned R,unsigned M,bool T=false>
NI void outer(const Plan& p,const E* in,E* out){opt::outer<R,M,T>(p.blocks,in,out);}
template<unsigned R,bool T=false>
NI void route(const Plan& p,const E* in,E* out){opt::route<R,T>(p.blocks,p.route.data(),in,out);}
template<unsigned R>
NI void streamed(const Plan& p,const E* in,E* out){opt::outer<R,2,false,true>(p.blocks,in,out);}
template<unsigned R,bool S=false>
NI void single(const Plan& p,const E* in,E* out){opt::outerOne<R,S>(p.blocks,in,out);}
template<unsigned P>
NI void batchedRoute(const Plan& p,const E* in,E* out){opt::routeBatched<P>(p.blocks,p.route.data(),in,out);}
template<unsigned R>
NI void fused(const Plan& p,const E* in,E* out){opt::outerFused<R>(p.blocks,in,out);}
struct Buffer {
    void* mapping;std::size_t bytes;E* pointer;
    explicit Buffer(std::size_t count,bool huge):bytes(count*sizeof(E)+(1U<<21)){
        mapping=mmap(nullptr,bytes,PROT_READ|PROT_WRITE,MAP_PRIVATE|MAP_ANONYMOUS,-1,0);
        require(mapping!=MAP_FAILED,"mmap failed");
        pointer=reinterpret_cast<E*>((reinterpret_cast<std::uintptr_t>(mapping)+(1U<<21)-1)&~std::uintptr_t((1U<<21)-1));
        const auto advised=(count*sizeof(E)+4095)&~std::size_t(4095);
        require(madvise(pointer,advised,huge?MADV_HUGEPAGE:MADV_NOHUGEPAGE)==0,"madvise failed");
        std::memset(pointer,0,count*sizeof(E));
        for(unsigned j=0;j<4;++j)pointer[count+j]={0x7d571ca9,0xa92f1803};
    }
    ~Buffer(){munmap(mapping,bytes);}
    E* data(){return pointer;}E& operator[](std::size_t i){return pointer[i];}
    void checkGuard(std::size_t count) const {for(unsigned j=0;j<4;++j)require(pointer[count+j]==E{0x7d571ca9,0xa92f1803},"output guard changed");}
    std::size_t hugeKiB() const {
        std::ifstream smaps("/proc/self/smaps");std::string line;bool inside=false;
        while(std::getline(smaps,line)){
            unsigned long lo,hi;
            if(std::sscanf(line.c_str(),"%lx-%lx",&lo,&hi)==2)
                inside=reinterpret_cast<std::uintptr_t>(pointer)>=lo&&reinterpret_cast<std::uintptr_t>(pointer)<hi;
            if(inside&&line.rfind("AnonHugePages:",0)==0)return std::stoull(line.substr(14));
        }
        return 0;
    }
};
struct Case {const char* name;Fn encode;Fn permute;unsigned rows;bool tiled;};
const Case cases[]={
    {"frozen",native::outer,native::routeOnly,4,false},
    {"copy_control",outer<4,0>,route<4>,4,false},
    {"leaf8",outer<4,1>,route<4>,4,false},
    {"fused_input",outer<4,2>,route<4>,4,false},
    {"blocked16",outer<16,2>,route<16>,16,false},
    {"blocked32",outer<32,2>,route<32>,32,false},
    {"blocked64",outer<64,2>,route<64>,64,false},
    {"blocked128",outer<128,2>,route<128>,128,false},
    {"blocked256",outer<256,2>,route<256>,256,false},
    {"tiled4",outer<4,2,true>,route<4,true>,4,true},
    {"tiled32",outer<32,2,true>,route<32,true>,32,true},
    {"tiled128",outer<128,2,true>,route<128,true>,128,true},
    {"stream16",streamed<16>,route<16>,16,false},
    {"stream32",streamed<32>,route<32>,32,false},
    {"stream64",streamed<64>,route<64>,64,false},
    {"stream128",streamed<128>,route<128>,128,false},
    {"single16",single<16>,route<16>,16,false},
    {"single32",single<32>,route<32>,32,false},
    {"single64",single<64>,route<64>,64,false},
    {"single128",single<128>,route<128>,128,false},
    {"single_stream32",single<32,true>,route<32>,32,false},
    {"single_stream64",single<64,true>,route<64>,64,false},
    {"single64_route16",single<64>,batchedRoute<0>,64,false},
    {"single64_prefill",single<64>,batchedRoute<1>,64,false},
    {"single64_lookahead",single<64>,batchedRoute<2>,64,false},
    {"frozen_route16",native::outer,batchedRoute<0>,4,false},
    {"fused16",fused<16>,batchedRoute<1>,16,false},
    {"fused32",fused<32>,batchedRoute<1>,32,false},
    {"fused64",fused<64>,batchedRoute<1>,64,false},
    {"fused128",fused<128>,batchedRoute<1>,128,false},
    {"fused128_route16",fused<128>,batchedRoute<0>,128,false},
    {"fused128_warm_loads",fused<128>,batchedRoute<3>,128,false},
    {"fused128_far_prefetch",fused<128>,batchedRoute<4>,128,false}
};
void kernelChecks(){
    alignas(64) E input[512],ref[1024],out[1024];
    std::mt19937_64 rng(293);
    for(unsigned t=0;t<684;++t){
        std::fill(input,input+512,E{});
        if(t<512)input[t]={1,0};
        else if(t<640){const auto bit=t-512;input[(t%4)*128+(t%128)]=bit<64?E{1ULL<<bit,0}:E{0,1ULL<<(bit-64)};}
        else for(auto& x:input)x=random(rng);
        opt::encodeFour<0>(input,ref);
        opt::encodeFour<1>(input,out);require(std::memcmp(ref,out,sizeof(ref))==0,"leaf8 differs");
        opt::encodeFour<2>(input,out);require(std::memcmp(ref,out,sizeof(ref))==0,"fused input differs");
        alignas(64) E one[256];
        for(unsigned lane=0;lane<4;++lane){
            opt::encodeOne(input+128*lane,one);
            for(unsigned j=0;j<256;++j)require(one[j]==ref[4*j+lane],"single-block differs");
        }
    }
    // Independent root-product generator matrix, not an FFT-based oracle.
    const auto& g=oracleTables().g;
    opt::encodeFour<2>(input,out);
    for(unsigned lane=0;lane<4;++lane)for(unsigned j=0;j<256;++j){
        E sum{};for(unsigned k=0;k<128;++k)sum^=byteScale(input[lane*128+k],g[k][j]);
        require(out[4*j+lane]==sum,"independent RS oracle failed");
    }
}
void tailChecks(){
    std::mt19937_64 rng(774);
    for(std::size_t blocks:{4,12,36,132,260}){
        Aligned<E> input(blocks*128),ref(blocks*256),out(blocks*256),routed(blocks*256);
        Aligned<std::uint32_t> labels(blocks*256);
        for(std::size_t i=0;i<blocks*128;++i)input[i]=random(rng);
        alignas(64) E tmp[1024];
        for(std::size_t r=0;r<blocks;r+=4){
            opt::encodeFour<0>(input.data()+r*128,tmp);
            for(unsigned j=0;j<256;++j)for(unsigned lane=0;lane<4;++lane)ref[j*blocks+r+lane]=tmp[4*j+lane];
        }
        opt::encodeRegional(blocks,input.data(),out.data());
        require(std::memcmp(out.data(),ref.data(),blocks*256*sizeof(E))==0,"partial outer tile failed");
        for(unsigned j=0;j<256;++j){
            auto* ix=labels.data()+j*blocks;std::iota(ix,ix+blocks,0U);std::shuffle(ix,ix+blocks,rng);
        }
        opt::routeRegional(blocks,labels.data(),out.data(),routed.data());
        for(std::size_t i=0;i<blocks*256;++i)require(routed[i]==ref[(i/blocks)*blocks+labels[i]],"partial route failed");
    }
}
double med(std::vector<double> v){std::sort(v.begin(),v.end());return v[v.size()/2];}
// Integrate only the exact outer replacement. The original fused routing and
// recursive inner remain untouched; no standalone-routing gain is assumed.
NI void fullOptimized(const Plan& p,const E* input,E* workspace,E* output){
    opt::encodeRegional(p.blocks,input,workspace);
    native::inner<true,false,true>(p,workspace,output);
}
void run(unsigned logK,unsigned reps,bool huge,bool finalOnly,bool fullIntegration){
    require(logK>=12&&logK<=20&&reps>=3&&reps<=101,"bounded benchmark arguments");
    kernelChecks();tailChecks();Plan p(std::size_t(1)<<logK);
    Buffer in(p.K,huge),work(p.N,huge),out(p.N,huge),refOuter(p.N,huge),refRoute(p.N,huge);
    std::mt19937_64 rng(931);for(std::size_t i=0;i<p.K;++i)in[i]=random(rng);
    const auto inputSum=checksum(in.data(),p.K);
    native::outer(p,in.data(),refOuter.data());native::routeOnly(p,refOuter.data(),refRoute.data());
    for(auto& c:cases){
        c.encode(p,in.data(),work.data());
        for(unsigned j=0;j<256;++j)for(std::size_t r=0;r<p.blocks;++r){
            const auto first=(r/c.rows)*c.rows,width=std::min<std::size_t>(c.rows,p.blocks-first);
            const auto at=c.tiled?first*256+j*width+r-first:j*p.blocks+r;
            require(work[at]==refOuter[j*p.blocks+r],"outer map changed");
        }
        c.permute(p,work.data(),out.data());
        require(std::memcmp(out.data(),refRoute.data(),p.N*sizeof(E))==0,"route map changed");
    }
    std::cout<<"K="<<p.K<<" reps="<<reps<<" huge_requested="<<huge
        <<" input_work_output_huge_KiB="<<in.hugeKiB()<<","<<work.hugeKiB()<<","<<out.hugeKiB()
        <<" kernel_and_full_checks=PASS checksum="<<checksum(refRoute.data(),p.N)<<"\n"<<std::flush;
    constexpr auto n=std::size(cases);
    std::vector<double> solo[n],combined[n];
    std::uint64_t sink=0;
    for(unsigned trial=0;trial<reps+3;++trial){
        for(std::size_t k=0;k<n;++k){
            const auto i=(trial&1)?n-1-k:k;const auto& c=cases[i];
            if(finalOnly&&std::string(c.name)!="frozen"&&std::string(c.name)!="fused128"&&std::string(c.name)!="fused128_warm_loads")continue;
            auto a=Clock::now();c.encode(p,in.data(),work.data());auto b=Clock::now();
            const auto alone=std::chrono::duration<double,std::milli>(b-a).count();
            sink^=checksum(work.data(),p.N);
            a=Clock::now();c.encode(p,in.data(),work.data());c.permute(p,work.data(),out.data());b=Clock::now();
            const auto together=std::chrono::duration<double,std::milli>(b-a).count();
            require(std::memcmp(out.data(),refRoute.data(),p.N*sizeof(E))==0,"timed output mismatch");
            if(trial>=3){solo[i].push_back(alone);combined[i].push_back(together);}
        }
    }
    std::cout<<std::fixed<<std::setprecision(6);
    for(std::size_t i=0;i<n;++i)if(!solo[i].empty())std::cout<<cases[i].name<<" outer_ms="<<med(solo[i])<<" outer_route_ms="<<med(combined[i])<<"\n";
    if(fullIntegration){
        native::full<false,true>(p,in.data(),work.data(),refRoute.data());
        fullOptimized(p,in.data(),work.data(),out.data());
        require(std::memcmp(out.data(),refRoute.data(),p.N*sizeof(E))==0,"full encoder output changed");
        std::vector<double> oldFull,newFull;
        for(unsigned trial=0;trial<12;++trial){
            for(unsigned k=0;k<2;++k){
                const bool updated=(k^(trial&1))!=0;
                const auto start=Clock::now();
                if(updated)fullOptimized(p,in.data(),work.data(),out.data());
                else native::full<false,true>(p,in.data(),work.data(),out.data());
                const auto ms=std::chrono::duration<double,std::milli>(Clock::now()-start).count();
                require(std::memcmp(out.data(),refRoute.data(),p.N*sizeof(E))==0,"timed full output changed");
                if(trial>=3)(updated?newFull:oldFull).push_back(ms);
            }
        }
        std::cout<<"full_control_ms="<<med(oldFull)<<" full_optimized_outer_ms="<<med(newFull)
            <<" full_samples=9 full_checksum="<<checksum(out.data(),p.N)<<"\n";
    }
    require(checksum(in.data(),p.K)==inputSum,"input changed");
    in.checkGuard(p.K);work.checkGuard(p.N);out.checkGuard(p.N);refOuter.checkGuard(p.N);refRoute.checkGuard(p.N);
    std::cout<<"input_and_guards=PASS sink="<<sink<<"\n";
}
}
int main(int argc,char**argv){try{if(argc!=4&&argc!=5)return 2;ro::run(unsigned(std::stoul(argv[1])),unsigned(std::stoul(argv[2])),std::stoul(argv[3])!=0,argc==5,argc==5&&std::string(argv[4])=="full");return 0;}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
