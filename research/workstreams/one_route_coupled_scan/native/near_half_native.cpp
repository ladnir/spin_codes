// Research-only implementation of NEAR_HALF_WINDOW.md, not a library profile.
// F128, fixed RS[256,128], one stratified route, exact Toeplitz L1024/b256.
#include "ToeplitzWindow1024.h"
#if SPIN_OR_WINDOW_NINE
#include "ToeplitzWindow1024Nine.h"
#endif
#ifndef SPIN_OR_WINDOW_RAW_SCALE
#define SPIN_OR_WINDOW_RAW_SCALE 0
#endif
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace native {
namespace sb=spin::research::char2::field128_lift::sparse_byte;
using E=sb::Element;
#if SPIN_OR_WINDOW_NINE
using Kernel=spin::research::one_route::ToeplitzWindow1024Nine;
#else
using Kernel=spin::research::one_route::ToeplitzWindow1024;
#endif
using FFT=spin::research::one_route::Transform256Cantor;
using Clock=std::chrono::steady_clock;
#define NI __attribute__((noinline))
template<class T> struct Aligned {
    struct Release {void operator()(T* p)const noexcept{std::free(p);}};
    std::unique_ptr<T,Release> pointer;
    std::size_t count;
    explicit Aligned(std::size_t n):count(n){
        pointer.reset(static_cast<T*>(std::aligned_alloc(64,(n*sizeof(T)+63)&~std::size_t(63))));
        if(!pointer)throw std::bad_alloc();
        std::memset(pointer.get(),0,n*sizeof(T));
    }
    T* data(){return pointer.get();} const T* data()const{return pointer.get();}
    T& operator[](std::size_t i){return data()[i];}
    const T& operator[](std::size_t i)const{return data()[i];}
};
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
E random(std::mt19937_64& rng){return {rng(),rng()};}
E byteScale(E x,std::uint8_t a){
    E result;
    _mm_storeu_si128(reinterpret_cast<__m128i*>(&result),
        _mm_gf2p8mul_epi8(_mm_loadu_si128(reinterpret_cast<const __m128i*>(&x)),_mm_set1_epi8(char(a))));
    return result;
}
std::uint8_t power(std::uint8_t x,unsigned n){
    std::uint8_t y=1;
    for(;n;n>>=1){if(n&1)y=sb::byteMultiply(y,x);x=sb::byteMultiply(x,x);}
    return y;
}

struct Plan {
    std::size_t K,N,blocks,batches;
    Aligned<std::uint32_t> route;
    Aligned<Kernel::PreparedBatch> keys;
    std::vector<E> raw;
    double routeMs=0,keyMs=0;
    explicit Plan(std::size_t k,bool retainRaw=false):K(k),N(2*k),blocks(k/128),batches(N/256),
        route(N),keys(batches){
        require(K>=512&&K%512==0,"reference K must be a multiple of512");
        require(blocks<256||blocks%256==0,"large regions must contain whole batches");
        std::mt19937_64 rng(113);
        // Reproducible fixtures only. The theorem assumes ideal uniform setup;
        // this benchmark generator is not a production cryptographic sampler.
        auto start=Clock::now();
        for(unsigned column=0;column<256;++column){
            auto* p=route.data()+column*blocks;
            std::iota(p,p+blocks,0U);std::shuffle(p,p+blocks,rng);
        }
        routeMs=std::chrono::duration<double,std::milli>(Clock::now()-start).count();
        if(retainRaw)raw.resize(batches*1279);
        start=Clock::now();
        alignas(64) E diagonals[1279];
        for(std::size_t t=0;t<batches;++t){
            for(auto& x:diagonals)x=random(rng);
            if(retainRaw)std::memcpy(raw.data()+1279*t,diagonals,sizeof(diagonals));
            Kernel::prepare(diagonals,keys[t]);
        }
        keyMs=std::chrono::duration<double,std::milli>(Clock::now()-start).count();
    }
};

// Novel-basis coefficients parameterize the same degree<128 RS space.
// Four independent blocks share a schedule, not message values or randomness.
NI void outer(const Plan& p,const E* input,E* planes){
    alignas(64) E scratch[256*4];
    for(std::size_t row=0;row<p.blocks;row+=4){
        for(unsigned j=0;j<128;++j){
            scratch[4*j]=input[(row+0)*128+j];
            scratch[4*j+1]=input[(row+1)*128+j];
            scratch[4*j+2]=input[(row+2)*128+j];
            scratch[4*j+3]=input[(row+3)*128+j];
        }
        std::memset(scratch+512,0,512*sizeof(E));
        FFT::evaluateNovelFour(scratch);
        for(unsigned j=0;j<256;++j)
            _mm512_store_si512(planes+j*p.blocks+row,_mm512_load_si512(scratch+4*j));
    }
}

template<bool Routed,bool Gfni,bool Deferred> NI void inner(const Plan& p,const E* input,E* output){
    Kernel::State state;state.reset();
    Kernel::Scratch scratch;
    alignas(64) E feedback[256];
    std::size_t base=0,nextRegion=p.blocks;
    for(std::size_t t=0;t<p.batches;++t){
        state.template feedback<Gfni,Deferred>(p.keys[t],feedback,scratch);
        const auto offset=256*t;
        if constexpr(Routed){
            if(p.blocks>=256){
                if(offset==nextRegion){base=nextRegion;nextRegion+=p.blocks;}
                const auto* source=input+base;const auto* labels=p.route.data()+offset;
                for(unsigned j=0;j<256;j+=4){
                    auto x=_mm512_castsi128_si512(_mm_loadu_si128(reinterpret_cast<const __m128i*>(source+labels[j])));
                    x=_mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(source+labels[j+1])),1);
                    x=_mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(source+labels[j+2])),2);
                    x=_mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(source+labels[j+3])),3);
                    _mm512_store_si512(output+offset+j,_mm512_xor_si512(x,_mm512_load_si512(feedback+j)));
                }
            }else{
                for(unsigned j=0;j<256;++j){
                    const auto index=offset+j,region=(index/p.blocks)*p.blocks;
                    output[index]=input[region+p.route[index]]^feedback[j];
                }
            }
        }else{
            for(unsigned j=0;j<256;j+=4)
                _mm512_store_si512(output+offset+j,_mm512_xor_si512(
                    _mm512_load_si512(input+offset+j),_mm512_load_si512(feedback+j)));
        }
        state.append(output+offset);
    }
}
NI void routeOnly(const Plan& p,const E* input,E* output){
    for(unsigned j=0;j<256;++j){
        const auto base=j*p.blocks;
        for(std::size_t i=0;i<p.blocks;++i)output[base+i]=input[base+p.route[base+i]];
    }
}
template<bool Gfni,bool Deferred> NI void full(const Plan& p,const E* input,E* workspace,E* output){
    outer(p,input,workspace);inner<true,Gfni,Deferred>(p,workspace,output);
}
void directOracle(const E* state,const E* diagonals,E* output){
    std::fill(output,output+256,E{});
    for(unsigned j=0;j<1024;++j)if(state[j]!=E{})
        for(unsigned i=0;i<256;++i)output[i]^=sb::multiplyOracle(state[j],diagonals[1023+i-j]);
}

// Independent root-product construction of the novel polynomial basis.
// No FFT implementation or its precomputed twiddles occurs in this oracle.
struct OuterOracle {
    std::uint8_t g[128][256]{};
    std::uint8_t basis[256][256]{};
    OuterOracle(){
        std::uint8_t s[8][129]{};
        for(unsigned level=0;level<8;++level){
            const auto h=1U<<level;s[level][0]=1;
            for(unsigned root=0;root<h;++root){
                auto a=FFT::evaluationPoint(root);
                for(unsigned d=root+1;d>0;--d)s[level][d]=s[level][d-1]^sb::byteMultiply(s[level][d],a);
                s[level][0]=sb::byteMultiply(s[level][0],a);
            }
            auto a=FFT::evaluationPoint(h);std::uint8_t value=0;
            for(unsigned d=h+1;d-- >0;)value=sb::byteMultiply(value,a)^s[level][d];
            require(value!=0,"zero subspace normalization");
            auto inv=power(value,254);
            for(unsigned d=0;d<=h;++d)s[level][d]=sb::byteMultiply(s[level][d],inv);
        }
        basis[0][0]=1;
        for(unsigned degree=1;degree<256;++degree){
            unsigned level=0;while((2U<<level)<=degree)++level;
            unsigned h=1U<<level,low=degree-h;
            for(unsigned i=0;i<=low;++i)for(unsigned j=0;j<=h;++j)
                basis[degree][i+j]^=sb::byteMultiply(basis[low][i],s[level][j]);
        }
        for(unsigned degree=0;degree<128;++degree)for(unsigned point=0;point<256;++point){
            std::uint8_t value=0;
            for(unsigned d=degree+1;d-- >0;)value=sb::byteMultiply(value,FFT::evaluationPoint(point))^basis[degree][d];
            g[degree][point]=value;
        }
    }
};
const OuterOracle& oracleTables(){static const OuterOracle t;return t;}
// Independent polynomial-coordinate oracle, using root products and triangular
// division. It does not call any FFT conversion/evaluation/interpolation code.
void inputCoordinates(const E* input,E* output){
    const auto& b=oracleTables().basis;
    std::fill(output,output+1024,E{});
    for(unsigned half=0;half<8;++half)for(unsigned j=0;j<128;++j){
        if(input[128*half+j]==E{})continue;
        for(unsigned d=0;d<=j;++d)if(b[j][d])
            output[128*half+d]^=byteScale(input[128*half+j],b[j][d]);
    }
}
void outputCoordinates(E* output){
    const auto& b=oracleTables().basis;
    for(unsigned half=0;half<2;++half)for(unsigned j=128;j-- >0;){
        const auto value=byteScale(output[128*half+j],power(b[128+j][128+j],254));
        output[128*half+j]=value;
        for(unsigned d=0;d<j;++d)if(b[128+j][128+d])
            output[128*half+d]^=byteScale(value,b[128+j][128+d]);
    }
}
void oracle(const E* state,const E* diagonals,E* output){
    if constexpr(SPIN_OR_WINDOW_COORDINATES){
        alignas(64) E converted[1024];
        inputCoordinates(state,converted);directOracle(converted,diagonals,output);
        outputCoordinates(output);
    }else directOracle(state,diagonals,output);
}
#if SPIN_OR_WINDOW_NINE
template<unsigned C> void checkRawScale(){
    namespace nine=spin::research::one_route::window_nine_detail;
    const auto zero=_mm512_setzero_si512();
    alignas(64) E actual[4],plain[4];
    for(unsigned word=0;word<65536;++word){
        const auto x=_mm512_set1_epi16(static_cast<short>(word));
        sb::deferred::Product raw{x,zero,zero,zero};
        _mm512_store_si512(actual,sb::deferred::reduce(nine::scaleRaw<C>(raw)));
        _mm512_store_si512(plain,sb::deferred::reduce(raw));
        for(unsigned lane=0;lane<4;++lane)
            require(actual[lane]==byteScale(plain[lane],C),"raw coefficient-slot scaling mismatch");
    }
    // All512 binary coefficient directions, including the extension tail.
    for(unsigned bit=0;bit<512;++bit){
        alignas(64) std::uint64_t words[32]{};
        words[(bit/128)*8+(bit%128)/64]=std::uint64_t(1)<<(bit%64);
        const sb::deferred::Product raw{_mm512_load_si512(words),_mm512_load_si512(words+8),
            _mm512_load_si512(words+16),_mm512_load_si512(words+24)};
        _mm512_store_si512(actual,sb::deferred::reduce(nine::scaleRaw<C>(raw)));
        _mm512_store_si512(plain,sb::deferred::reduce(raw));
        for(unsigned lane=0;lane<4;++lane)
            require(actual[lane]==byteScale(plain[lane],C),"raw extension-slot scaling mismatch");
    }
}
#endif
NI void check(){
#if SPIN_OR_WINDOW_NINE
    checkRawScale<FFT::evaluationPoint(2)>();checkRawScale<FFT::evaluationPoint(4)>();
#endif
    std::mt19937_64 rng(227);
    Kernel::PreparedBatch keys;Kernel::State state;Kernel::Scratch scratch;
    alignas(64) E raw[1279],s[1024]{},actual[256],expected[256],other[256];
    for(auto& x:raw)x=random(rng);
    Kernel::prepare(raw,keys);
    for(unsigned j=0;j<1024;++j){
        std::fill(s,s+1024,E{});s[j]={1,0};state.load(s);
        state.feedback<false,false>(keys,actual,scratch);
        if constexpr(SPIN_OR_WINDOW_COORDINATES)oracle(s,raw,expected);
        else for(unsigned i=0;i<256;++i)expected[i]=raw[1023+i-j];
        require(std::memcmp(actual,expected,sizeof(actual))==0,"unit-state map mismatch");
        state.feedback<false,true>(keys,other,scratch);
        require(std::memcmp(actual,other,sizeof(actual))==0,"deferred unit-state mismatch");
    }
    for(unsigned fixture=0;fixture<2;++fixture){
        for(auto& x:s)x=random(rng);
        state.load(s);
        oracle(s,raw,expected);
        state.feedback<false,false>(keys,actual,scratch);
        require(std::memcmp(actual,expected,sizeof(actual))==0,"dense independent Toeplitz mismatch");
        state.feedback<false,true>(keys,other,scratch);
        require(std::memcmp(actual,other,sizeof(actual))==0,"deferred dense mismatch");
        state.feedback<true,false>(keys,other,scratch);
        require(std::memcmp(actual,other,sizeof(actual))==0,"GFNI dense mismatch");
    }
    for(unsigned diagonal:{0U,1U,126U,127U,128U,255U,256U,1022U,1023U,1277U,1278U}){
        std::fill(raw,raw+1279,E{});raw[diagonal]={1,0};Kernel::prepare(raw,keys);
        state.load(s);state.feedback<false,true>(keys,actual,scratch);
        if constexpr(SPIN_OR_WINDOW_COORDINATES)oracle(s,raw,expected);
        else for(unsigned i=0;i<256;++i){int j=1023+int(i)-int(diagonal);expected[i]=j>=0&&j<1024?s[j]:E{};}
        require(std::memcmp(actual,expected,sizeof(actual))==0,"diagonal boundary mismatch");
    }
    Plan p(2048,true); //16 steps cross the four-batch ring repeatedly.
    Aligned<E> input(p.K),planes(p.N),out(p.N),reference(p.N),control(p.N);
    for(std::size_t i=0;i<p.K;++i)input[i]=random(rng);
    outer(p,input.data(),planes.data());
    const auto& table=oracleTables();
    for(unsigned j=0;j<256;++j){
        std::vector<bool> seen(p.blocks);
        for(std::size_t row=0;row<p.blocks;++row){
            E value{};for(unsigned d=0;d<128;++d)value^=byteScale(input[row*128+d],table.g[d][j]);
            require(value==planes[j*p.blocks+row],"RS polynomial oracle mismatch");
            auto label=p.route[j*p.blocks+row];
            require(label<p.blocks&&!seen[label],"regional route is not bijective");seen[label]=true;
        }
    }
    inner<true,false,true>(p,planes.data(),out.data());
    std::fill(s,s+1024,E{});
    for(std::size_t t=0;t<p.batches;++t){
        oracle(s,p.raw.data()+1279*t,expected);
        for(unsigned i=0;i<256;++i){
            auto index=256*t+i,base=(index/p.blocks)*p.blocks;
            reference[index]=planes[base+p.route[index]]^expected[i];
        }
        std::memmove(s,s+256,768*sizeof(E));
        std::memcpy(s+768,reference.data()+256*t,256*sizeof(E));
    }
    require(std::memcmp(out.data(),reference.data(),p.N*sizeof(E))==0,"full recursive oracle mismatch");
    inner<true,false,false>(p,planes.data(),control.data());
    require(std::memcmp(out.data(),control.data(),p.N*sizeof(E))==0,"full deferred/control mismatch");
    // Check the different, contiguous in-place API separately from routing.
    std::memcpy(control.data(),planes.data(),p.N*sizeof(E));
    inner<false,false,true>(p,control.data(),control.data());
    inner<false,false,false>(p,planes.data(),reference.data());
    require(std::memcmp(control.data(),reference.data(),p.N*sizeof(E))==0,"in-place contiguous mismatch");
    {
        // The large-region vector branch uses a different address schedule.
        Plan large(32768);
        Aligned<E> x(large.N),routed(large.N),a(large.N),b(large.N);
        for(std::size_t i=0;i<large.N;++i)x[i]=random(rng);
        routeOnly(large,x.data(),routed.data());
        inner<true,false,true>(large,x.data(),a.data());
        inner<false,false,true>(large,routed.data(),b.data());
        require(std::memcmp(a.data(),b.data(),large.N*sizeof(E))==0,"large-region routing mismatch");
    }
    std::cout<<"PASS: all1024 state-basis inputs; dense independent field oracle; diagonal endpoints; "
        <<"GFNI/deferred/control; RS root-product oracle; regional bijections; 16-step full recursion; in-place contiguous; fast regional routing\n";
}

std::uint64_t checksum(const E* x,std::size_t n){
    std::uint64_t s=0x243f6a8885a308d3ULL;
    for(std::size_t i=0;i<n;++i)s=(s^(x[i].lo+x[i].hi))*0x9e3779b97f4a7c15ULL;
    return s;
}
template<class F> double measure(unsigned repeats,F&& f){
    for(unsigned j=0;j<3;++j)f();
    std::vector<double> times;times.reserve(repeats);
    for(unsigned j=0;j<repeats;++j){auto start=Clock::now();f();
        times.push_back(std::chrono::duration<double,std::milli>(Clock::now()-start).count());}
    std::sort(times.begin(),times.end());return times[times.size()/2];
}
NI void transformStream(const E* input,E* output,std::size_t count){
    for(std::size_t i=0;i<count;i+=256){
        alignas(64) E a[256];std::memcpy(a,input+i,sizeof(a));
        if constexpr(SPIN_OR_WINDOW_COORDINATES){FFT::evaluateNovel(a);FFT::interpolateNovel(a);}
        else{FFT::evaluate(a);FFT::interpolate(a);}
        std::memcpy(output+i,a,sizeof(a));
    }
}
template<bool Gfni> NI void productStream(const E* a,const E* b,E* output,std::size_t count){
    for(std::size_t i=0;i<count;i+=4){
        auto x=_mm512_load_si512(a+i),y=_mm512_load_si512(b+i);
        _mm512_store_si512(output+i,Gfni?sb::multiply4Gfni(x,y):sb::multiply4(x,y));
    }
}
NI void readTape(const Kernel::PreparedBatch* tape,std::size_t batches,E* output){
    const auto* p=reinterpret_cast<const __m512i*>(tape);
    const auto n=batches*sizeof(Kernel::PreparedBatch)/sizeof(__m512i);
    auto a=_mm512_setzero_si512(),b=a,c=a,d=a;
    for(std::size_t i=0;i<n;i+=4){
        a=_mm512_xor_si512(a,_mm512_load_si512(p+i));
        b=_mm512_xor_si512(b,_mm512_load_si512(p+i+1));
        c=_mm512_xor_si512(c,_mm512_load_si512(p+i+2));
        d=_mm512_xor_si512(d,_mm512_load_si512(p+i+3));
    }
    _mm512_store_si512(output,_mm512_xor_si512(_mm512_xor_si512(a,b),_mm512_xor_si512(c,d)));
}
void micro(unsigned logK,unsigned repeats){
    require(logK>=12&&logK<=20&&repeats>=3&&repeats<=101,"bounded micro args");
    constexpr std::size_t count=1<<18;
    Aligned<E> a(count),b(count),c(count);
    std::mt19937_64 rng(331);for(std::size_t i=0;i<count;++i){a[i]=random(rng);b[i]=random(rng);}
    const double batches=double((std::size_t(2)<<logK)/256);
    auto show=[&](const char* label,double ms){std::cout<<std::fixed<<std::setprecision(6)<<label<<"="<<ms
        <<" checksum="<<checksum(c.data(),count)<<"\n"<<std::flush;};
    // Independent hot streams, normalized only by operation counts. These are
    // not additive phase timings of the full encoder (which streams384MiB keys).
    show("two_transform_pairs_normalized_ms",measure(repeats,[&]{transformStream(a.data(),c.data(),count);})*batches*2/(count/256));
    show("twelve_products_per_output_clmul_normalized_ms",measure(repeats,[&]{productStream<false>(a.data(),b.data(),c.data(),count);})*batches*3072/count);
    show("twelve_products_per_output_gfni_normalized_ms",measure(repeats,[&]{productStream<true>(a.data(),b.data(),c.data(),count);})*batches*3072/count);
}
void bench(unsigned logK,unsigned repeats,const std::string& phase){
    require(logK>=12&&logK<=20&&repeats>=3&&repeats<=101,"bounded bench: logK12..20, repeats3..101");
    require(phase=="all"||phase=="full"||phase=="deferred"||phase=="gfni"||phase=="keys","unknown phase");
    auto start=Clock::now();Plan p(std::size_t(1)<<logK);
    std::cout<<std::fixed<<std::setprecision(6)<<"K="<<p.K<<" N="<<p.N<<" B=256 L=1024 b=256 setup_ms="
        <<std::chrono::duration<double,std::milli>(Clock::now()-start).count()
        <<" route_setup_ms="<<p.routeMs<<" key_setup_ms="<<p.keyMs
        <<" prepared_keys_MiB="<<p.keys.count*sizeof(Kernel::PreparedBatch)/1048576.0
        <<" repetitions="<<repeats<<" products_per_output="<<Kernel::productsPerBatch/256
        <<" coordinate_variant="<<SPIN_OR_WINDOW_COORDINATES<<" raw_scale="<<SPIN_OR_WINDOW_RAW_SCALE<<"\n"<<std::flush;
    Aligned<E> input(p.K),workspace(p.N),output(p.N);
    std::mt19937_64 rng(931);for(std::size_t i=0;i<p.K;++i)input[i]=random(rng);
    outer(p,input.data(),workspace.data());
    auto show=[&](const char* name,double ms){std::cout<<name<<"_ms="<<ms<<" checksum="<<checksum(output.data(),p.N)<<"\n"<<std::flush;};
    if(phase=="all"){
        show("outer",measure(repeats,[&]{outer(p,input.data(),output.data());}));
        show("route",measure(repeats,[&]{routeOnly(p,workspace.data(),output.data());}));
        show("inner_contiguous",measure(repeats,[&]{inner<false,false,false>(p,workspace.data(),output.data());}));
        show("inner_routed",measure(repeats,[&]{inner<true,false,false>(p,workspace.data(),output.data());}));
        show("inner_deferred",measure(repeats,[&]{inner<true,false,true>(p,workspace.data(),output.data());}));
    }
    if(phase=="all"||phase=="full")show("full",measure(repeats,[&]{full<false,false>(p,input.data(),workspace.data(),output.data());}));
    if(phase=="all"||phase=="deferred")show("full_deferred",measure(repeats,[&]{full<false,true>(p,input.data(),workspace.data(),output.data());}));
    if(phase=="gfni")show("full_gfni",measure(repeats,[&]{full<true,false>(p,input.data(),workspace.data(),output.data());}));
    if(phase=="keys")show("key_tape_read",measure(repeats,[&]{readTape(p.keys.data(),p.batches,output.data());}));
}
} //namespace native
int main(int argc,char** argv){
    try{
        if(argc==2&&std::string(argv[1])=="--check"){native::check();return 0;}
        if(argc==5&&std::string(argv[1])=="--bench"){
            if(std::string(argv[4])=="micro"){
                native::micro(unsigned(std::stoul(argv[2])),unsigned(std::stoul(argv[3])));return 0;}
            native::bench(unsigned(std::stoul(argv[2])),unsigned(std::stoul(argv[3])),argv[4]);return 0;}
        std::cerr<<"usage: near_half_native --check | --bench logK repeats all|full|deferred|gfni|micro|keys\n";return 2;
    }catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<"\n";return 1;}
}
