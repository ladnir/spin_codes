// Research-only exact implementations of STRATIFIED_CLOSED_CERTIFICATE.md
// and SPECTRAL_BATCH.md: F128 feedback, L=128,b=64, fixed RS[64,32].
// Compile-time controls select Toeplitz or independently sampled spectral maps.
// No production profile or small-coefficient substitution is made here.
#include "Transform256.h"
#include "SparseByteField128Gfni.h"
#include <algorithm>
#include <array>
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

#ifndef SPIN_OR_SPECTRAL_POINTS
#define SPIN_OR_SPECTRAL_POINTS 0
#endif
#ifndef SPIN_OR_SPLIT192
#define SPIN_OR_SPLIT192 0
#endif
#ifndef SPIN_OR_NOVEL_OUTER
#define SPIN_OR_NOVEL_OUTER 0
#endif
static_assert(SPIN_OR_SPECTRAL_POINTS==0||SPIN_OR_SPECTRAL_POINTS==192||SPIN_OR_SPECTRAL_POINTS==256);
static_assert(!SPIN_OR_SPLIT192||SPIN_OR_SPECTRAL_POINTS==192);

namespace native {
namespace sb=spin::research::char2::field128_lift::sparse_byte;
using E=sb::Element;
using spin::research::one_route::Transform;
using Clock=std::chrono::steady_clock;
#define NI __attribute__((noinline))
constexpr unsigned keyCount=SPIN_OR_SPECTRAL_POINTS?SPIN_OR_SPECTRAL_POINTS:256;
constexpr unsigned feedbackOffset=SPIN_OR_SPECTRAL_POINTS?192:127;

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
std::uint8_t power(std::uint8_t x,unsigned e){
    std::uint8_t y=1;
    for(;e;e>>=1){if(e&1)y=sb::byteMultiply(y,x);x=sb::byteMultiply(x,x);}
    return y;
}

// Four independent extension-field products, no hidden field-basis conversion.
template<bool Gfni> NI void pointwise(E* values,const E* keys){
    for(unsigned i=0;i<keyCount;i+=4){
        const auto a=_mm512_load_si512(values+i),b=_mm512_load_si512(keys+i);
        const auto c=Gfni?sb::multiply4Gfni(a,b):sb::multiply4(a,b);
        _mm512_store_si512(values+i,c);
    }
}

template<bool Gfni> inline __attribute__((always_inline))
void feedback(const E* state,const E* keys,E* scratch){
#if SPIN_OR_SPLIT192
    constexpr auto c=Transform<256>::novelBasisValue(64,128);
    for(unsigned j=0;j<64;j+=4){
        const auto a=_mm512_load_si512(state+j),b=_mm512_load_si512(state+64+j);
        _mm512_store_si512(scratch+j,a);
        _mm512_store_si512(scratch+64+j,_mm512_xor_si512(a,b));
        _mm512_store_si512(scratch+128+j,_mm512_xor_si512(a,_mm512_gf2p8mul_epi8(b,_mm512_set1_epi8(char(c)))));
    }
    Transform<64>::evaluateNovelCoset<0>(scratch);
    Transform<64>::evaluateNovelCoset<64>(scratch+64);
    Transform<64>::evaluateNovelCoset<128>(scratch+128);
    pointwise<Gfni>(scratch,keys);
    Transform<64>::interpolateNovelCoset<0>(scratch);
    Transform<64>::interpolateNovelCoset<64>(scratch+64);
    Transform<64>::interpolateNovelCoset<128>(scratch+128);
    for(unsigned j=0;j<64;j+=4)
        _mm512_store_si512(scratch+192+j,_mm512_ternarylogic_epi64(
            _mm512_load_si512(scratch+j),_mm512_load_si512(scratch+64+j),_mm512_load_si512(scratch+128+j),0x96));
#else
    std::memcpy(scratch,state,128*sizeof(E));
    std::memset(scratch+128,0,128*sizeof(E));
    if constexpr(SPIN_OR_SPECTRAL_POINTS)Transform<256>::evaluateNovel(scratch);
    else Transform<256>::evaluate(scratch);
    pointwise<Gfni>(scratch,keys);
    if constexpr(SPIN_OR_SPECTRAL_POINTS==192)std::memset(scratch+192,0,64*sizeof(E));
    if constexpr(SPIN_OR_SPECTRAL_POINTS)Transform<256>::interpolateNovel(scratch);
    else Transform<256>::interpolate(scratch);
#endif
}

template<bool Gfni> NI void middle(const E* state,const E* transformed,E* output){
    alignas(64) E scratch[256];
    feedback<Gfni>(state,transformed,scratch);
    std::memcpy(output,scratch+feedbackOffset,64*sizeof(E));
}

struct Plan {
    std::size_t K,N,blocks,batches;
    Aligned<std::uint32_t> route;
    Aligned<E> keys;
    std::vector<E> raw;
    explicit Plan(std::size_t k,bool retainRaw=false):K(k),N(2*k),blocks(k/32),batches(N/64),
        route(N),keys(batches*keyCount){
        require(K>=128&&K%128==0,"K must be a positive multiple of128");
        require(blocks<64||blocks%64==0,"large regions must contain whole64-symbol batches");
        std::mt19937_64 rng(113);
        // Reproducible experiment only. Certificates concern ideal uniform setup.
        for(unsigned column=0;column<64;++column){
            auto* p=route.data()+column*blocks;
            std::iota(p,p+blocks,0U);
            std::shuffle(p,p+blocks,rng);
        }
        if(retainRaw&&!SPIN_OR_SPECTRAL_POINTS)raw.resize(batches*191);
        for(std::size_t t=0;t<batches;++t){
            auto* p=keys.data()+t*keyCount;
            if constexpr(SPIN_OR_SPECTRAL_POINTS){
                for(unsigned j=0;j<keyCount;++j)p[j]=random(rng);
            }else{
                for(unsigned j=0;j<191;++j){p[j]=random(rng);if(retainRaw)raw[t*191+j]=p[j];}
                // D coefficient r+127 is diagonal r. The high65 coefficients are zero.
                Transform<256>::evaluate(p);
            }
        }
    }
};

// Four outer blocks share the instruction schedule, not their message values.
// Every output cache line contains one RS coordinate from four adjacent blocks.
NI void outer(const Plan& p,const E* input,E* planes){
    alignas(64) E scratch[64*4];
    for(std::size_t row=0;row<p.blocks;row+=4){
        for(unsigned j=0;j<32;++j){
            scratch[4*j]=input[(row+0)*32+j];
            scratch[4*j+1]=input[(row+1)*32+j];
            scratch[4*j+2]=input[(row+2)*32+j];
            scratch[4*j+3]=input[(row+3)*32+j];
        }
        std::memset(scratch+128,0,128*sizeof(E));
        if constexpr(SPIN_OR_NOVEL_OUTER)Transform<64>::evaluateNovelFour(scratch);
        else Transform<64>::evaluateFour(scratch);
        for(unsigned j=0;j<64;++j)
            _mm512_store_si512(planes+j*p.blocks+row,_mm512_load_si512(scratch+4*j));
    }
}

template<bool Routed,bool Gfni> NI void inner(const Plan& p,const E* input,E* output){
    alignas(64) E state[128]{};
    alignas(64) E transformed[256];
    std::size_t regionBase=0,nextRegion=p.blocks;
    for(std::size_t t=0;t<p.batches;++t){
        feedback<Gfni>(state,p.keys.data()+keyCount*t,transformed);
        std::memcpy(state,state+64,64*sizeof(E));
        const auto offset=64*t;
        if constexpr(Routed){
            if(p.blocks>=64){
                if(offset==nextRegion){regionBase=nextRegion;nextRegion+=p.blocks;}
                const auto* source=input+regionBase;
                const auto* labels=p.route.data()+offset;
                for(unsigned j=0;j<64;j+=4){
                    auto x=_mm512_castsi128_si512(_mm_loadu_si128(reinterpret_cast<const __m128i*>(source+labels[j])));
                    x=_mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(source+labels[j+1])),1);
                    x=_mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(source+labels[j+2])),2);
                    x=_mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(source+labels[j+3])),3);
                    const auto value=_mm512_xor_si512(x,_mm512_loadu_si512(transformed+feedbackOffset+j));
                    _mm512_store_si512(output+offset+j,value);
                    _mm512_store_si512(state+64+j,value);
                }
            }else{
                // Small oracle fixture: several region boundaries occur within a batch.
                for(unsigned j=0;j<64;++j){
                    const auto index=offset+j,base=(index/p.blocks)*p.blocks;
                    const auto value=input[base+p.route[index]]^transformed[feedbackOffset+j];
                    output[index]=value;state[64+j]=value;
                }
            }
        }else{
            for(unsigned j=0;j<64;j+=4){
                const auto value=_mm512_xor_si512(_mm512_load_si512(input+offset+j),
                                                _mm512_loadu_si512(transformed+feedbackOffset+j));
                _mm512_store_si512(output+offset+j,value);
                _mm512_store_si512(state+64+j,value);
            }
        }
    }
}
NI void routeOnly(const Plan& p,const E* input,E* output){
    for(unsigned column=0;column<64;++column){
        const auto base=column*p.blocks;
        for(std::size_t j=0;j<p.blocks;++j)output[base+j]=input[base+p.route[base+j]];
    }
}
template<bool Gfni> NI void full(const Plan& p,const E* input,E* workspace,E* output){
    outer(p,input,workspace);inner<true,Gfni>(p,workspace,output);
}

template<unsigned N> void checkTransform(std::mt19937_64& rng){
    alignas(64) E values[N],original[N];
    for(unsigned degree=0;degree<N;++degree){
        std::memset(values,0,sizeof(values));values[degree]=random(rng);
        std::memcpy(original,values,sizeof(values));
        Transform<N>::evaluate(values);
        for(unsigned point=0;point<N;++point)
            require(values[point]==byteScale(original[degree],power(Transform<N>::evaluationPoint(point),degree)),
                    "transform monomial basis mismatch");
        Transform<N>::interpolate(values);
        require(std::memcmp(values,original,sizeof(values))==0,"inverse transform mismatch");
    }
    alignas(64) E src[4*N],expected[4*N];
    for(auto& x:src)x=random(rng);
    for(unsigned lane=0;lane<4;++lane){
        for(unsigned j=0;j<N;++j)values[j]=src[4*j+lane];
        Transform<N>::evaluate(values);
        for(unsigned j=0;j<N;++j)expected[4*j+lane]=values[j];
    }
    Transform<N>::evaluateFour(src);
    require(std::memcmp(src,expected,sizeof(src))==0,"four-stream transform mismatch");
}

// Slow, independently constructed polynomial basis for the oracle only.
// This uses root products and triangular polynomial division, never FFTs.
struct OracleTables {
    std::array<std::array<std::uint8_t,256>,256> basis{};
    std::array<std::array<std::uint8_t,256>,128> g{};
    std::array<std::array<std::uint8_t,256>,64> h{};
    OracleTables(){
        std::array<std::array<std::uint8_t,129>,8> subspace{};
        for(unsigned level=0;level<8;++level){
            const unsigned half=1U<<level;auto& s=subspace[level];s[0]=1;
            for(unsigned root=0;root<half;++root){
                const auto a=Transform<256>::evaluationPoint(root);
                for(unsigned d=root+1;d>0;--d)s[d]=s[d-1]^sb::byteMultiply(s[d],a);
                s[0]=sb::byteMultiply(s[0],a);
            }
            const auto a=Transform<256>::evaluationPoint(half);
            std::uint8_t value=0;
            for(unsigned d=half+1;d-- >0;)value=sb::byteMultiply(value,a)^s[d];
            require(value!=0,"oracle subspace normalization zero");
            const auto inv=power(value,254);
            for(unsigned d=0;d<=half;++d)s[d]=sb::byteMultiply(s[d],inv);
        }
        basis[0][0]=1;
        for(unsigned degree=1;degree<256;++degree){
            unsigned level=0;while((2U<<level)<=degree)++level;
            const auto half=1U<<level,low=degree-half;
            for(unsigned i=0;i<=low;++i)
                for(unsigned j=0;j<=half;++j)
                    basis[degree][i+j]^=sb::byteMultiply(basis[low][i],subspace[level][j]);
        }
        for(unsigned point=0;point<256;++point){
            const auto a=Transform<256>::evaluationPoint(point);
            for(unsigned degree=0;degree<128;++degree){
                std::uint8_t value=0;
                for(unsigned d=degree+1;d-- >0;)value=sb::byteMultiply(value,a)^basis[degree][d];
                g[degree][point]=value;
            }
            std::array<std::uint8_t,64> top{};
            // Top64 coefficients of (Z^256-Z)/(Z-a), since derivative=1.
            for(unsigned d=0;d<64;++d)top[d]=power(a,63-d);
            for(unsigned d=64;d-- >0;){
                const auto value=sb::byteMultiply(top[d],power(basis[192+d][192+d],254));
                h[d][point]=value;
                for(unsigned j=0;j<d;++j)top[j]^=sb::byteMultiply(value,basis[192+d][192+j]);
            }
        }
    }
};

void spectralOracle(const OracleTables& tables,const E* state,const E* keys,E* output){
    std::fill(output,output+64,E{});
    for(unsigned a=0;a<keyCount;++a){
        E value{};
        for(unsigned j=0;j<128;++j)value^=byteScale(state[j],tables.g[j][a]);
        value=sb::multiplyOracle(value,keys[a]);
        for(unsigned i=0;i<64;++i)output[i]^=byteScale(value,tables.h[i][a]);
    }
}

NI void check(){
    std::mt19937_64 rng(227);
    checkTransform<64>(rng);checkTransform<256>(rng);
    const OracleTables tables;
    alignas(64) E state[128],keys[256],diagonal[191],out[64],other[64];
    if constexpr(SPIN_OR_SPECTRAL_POINTS){
        for(auto& x:keys)x=random(rng);
        alignas(64) E polynomial[256];
        for(unsigned degree=0;degree<128;++degree){
            std::memset(polynomial,0,sizeof(polynomial));polynomial[degree]={1,0};
            Transform<256>::evaluateNovel(polynomial);
            for(unsigned point=0;point<256;++point)
                require(polynomial[point]==E{tables.g[degree][point],0},"novel evaluation polynomial oracle mismatch");
        }
        for(unsigned point=0;point<256;++point){
            std::memset(polynomial,0,sizeof(polynomial));polynomial[point]={1,0};
            Transform<256>::interpolateNovel(polynomial);
            for(unsigned j=0;j<64;++j)
                require(polynomial[192+j]==E{tables.h[j][point],0},"novel extraction polynomial oracle mismatch");
        }
    }else{
        for(auto& x:diagonal)x=random(rng);
        std::memcpy(keys,diagonal,sizeof(diagonal));std::memset(keys+191,0,65*sizeof(E));
        Transform<256>::evaluate(keys);
        for(unsigned j=0;j<128;++j){
            std::memset(state,0,sizeof(state));state[j]={1,0};
            middle<false>(state,keys,out);
            for(unsigned i=0;i<64;++i)require(out[i]==diagonal[127+i-j],"middle product state-basis mismatch");
        }
    }
    for(unsigned fixture=0;fixture<4;++fixture){
        for(auto& x:state)x=random(rng);
        middle<false>(state,keys,out);middle<true>(state,keys,other);
        require(std::memcmp(out,other,sizeof(out))==0,"CLMUL/GFNI middle product mismatch");
        if constexpr(SPIN_OR_SPECTRAL_POINTS){
            spectralOracle(tables,state,keys,other);
            require(std::memcmp(out,other,sizeof(out))==0,"spectral dense polynomial oracle mismatch");
        }else{
            for(unsigned i=0;i<64;++i){
                E expected{};
                for(unsigned j=0;j<128;++j)expected^=sb::multiplyOracle(state[j],diagonal[127+i-j]);
                require(out[i]==expected,"middle product dense oracle mismatch");
            }
        }
    }
    Plan p(256,true);
    Aligned<E> message(p.K),planes(p.N),encoded(p.N),reference(p.N);
    for(std::size_t j=0;j<p.K;++j)message[j]=random(rng);
    outer(p,message.data(),planes.data());
    for(unsigned column=0;column<64;++column){
        std::vector<bool> used(p.blocks);
        for(std::size_t block=0;block<p.blocks;++block){
            E value{};
            if constexpr(SPIN_OR_NOVEL_OUTER){
                for(unsigned degree=0;degree<32;++degree)value^=byteScale(message[32*block+degree],tables.g[degree][column]);
            }else{
                for(unsigned degree=32;degree-- >0;)value=byteScale(value,Transform<64>::evaluationPoint(column))^message[32*block+degree];
            }
            require(planes[column*p.blocks+block]==value,"outer Horner oracle mismatch");
            auto label=p.route[column*p.blocks+block];
            require(label<p.blocks&&!used[label],"route not a regional bijection");used[label]=true;
        }
    }
    inner<true,false>(p,planes.data(),encoded.data());
    std::memset(state,0,sizeof(state));
    for(std::size_t t=0;t<p.batches;++t){
        if constexpr(SPIN_OR_SPECTRAL_POINTS)spectralOracle(tables,state,p.keys.data()+keyCount*t,other);
        for(unsigned i=0;i<64;++i){
            const auto index=64*t+i,base=(index/p.blocks)*p.blocks;
            E value=planes[base+p.route[index]];
            if constexpr(SPIN_OR_SPECTRAL_POINTS)value^=other[i];
            else for(unsigned j=0;j<128;++j)value^=sb::multiplyOracle(state[j],p.raw[t*191+127+i-j]);
            reference[index]=value;
        }
        std::memcpy(state,state+64,64*sizeof(E));
        std::memcpy(state+64,reference.data()+64*t,64*sizeof(E));
    }
    require(std::memcmp(encoded.data(),reference.data(),p.N*sizeof(E))==0,"full stream oracle mismatch");
    std::cout<<"PASS: transforms64/256, inverse, four-stream outer, "
        <<(SPIN_OR_SPECTRAL_POINTS?"novel evaluation/extraction polynomial oracles, ":"128 state basis middle products, ")
        <<"4 dense feedback oracles, GFNI/CLMUL equality, regional bijections, full stream/free terminal oracle\n";
}

std::uint64_t checksum(const E* data,std::size_t n){
    std::uint64_t s=0x243f6a8885a308d3ULL;
    for(std::size_t j=0;j<n;++j)s=(s^(data[j].lo+data[j].hi))*0x9e3779b97f4a7c15ULL;
    return s;
}
template<class F> double measure(unsigned repeats,F&& f){
    for(unsigned j=0;j<3;++j)f();
    std::vector<double> times;times.reserve(repeats);
    for(unsigned j=0;j<repeats;++j){auto start=Clock::now();f();
        times.push_back(std::chrono::duration<double,std::milli>(Clock::now()-start).count());}
    std::sort(times.begin(),times.end());return times[times.size()/2];
}
NI void transformStream(const E* input,E* output,std::size_t batches){
    for(std::size_t t=0;t<batches;++t){
        alignas(64) E values[256];std::memcpy(values,input+256*t,sizeof(values));
        if constexpr(SPIN_OR_SPECTRAL_POINTS){Transform<256>::evaluateNovel(values);Transform<256>::interpolateNovel(values);}
        else {Transform<256>::evaluate(values);Transform<256>::interpolate(values);}
        std::memcpy(output+256*t,values,sizeof(values));
    }
}
template<bool Gfni> NI void productStream(const E* input,const E* keys,E* output,std::size_t count){
    for(std::size_t i=0;i<count;i+=4){
        const auto a=_mm512_load_si512(input+i),b=_mm512_load_si512(keys+i);
        _mm512_store_si512(output+i,Gfni?sb::multiply4Gfni(a,b):sb::multiply4(a,b));
    }
}
void bench(unsigned logK,unsigned repetitions,const std::string& phase){
    require(logK>=12&&logK<=20&&repetitions>=3&&repetitions<=501,"bounded bench expects logK12..20,repetitions3..501");
    require(phase=="all"||phase=="outer"||phase=="route"||phase=="inner"||phase=="full"||phase=="gfni"||phase=="micro","unknown benchmark phase");
    const auto start=Clock::now();Plan p(std::size_t(1)<<logK);
    std::cout<<std::fixed<<std::setprecision(6)<<"K="<<p.K<<" N="<<p.N<<" L=128 b=64 B=64 setup_ms="
        <<std::chrono::duration<double,std::milli>(Clock::now()-start).count()
        <<" keys_MiB="<<p.keys.count*sizeof(E)/1048576.0<<" repetitions="<<repetitions
        <<" spectral_points="<<SPIN_OR_SPECTRAL_POINTS<<" split192="<<SPIN_OR_SPLIT192
        <<" novel_outer="<<SPIN_OR_NOVEL_OUTER<<" cantor="<<SPIN_OR_TRANSFORM_CANTOR<<"\n"<<std::flush;
    Aligned<E> input(p.K),workspace(p.N),output(p.N);
    std::mt19937_64 rng(931);for(std::size_t j=0;j<p.K;++j)input[j]=random(rng);
    outer(p,input.data(),workspace.data());
    const auto show=[&](const char* label,double time){
        std::cout<<label<<"_ms="<<time<<" checksum="<<checksum(output.data(),p.N)<<"\n"<<std::flush;};
    if(phase=="all"||phase=="outer"){
        show("outer",measure(repetitions,[&]{outer(p,input.data(),output.data());}));
    }
    if(phase=="all"||phase=="route")show("route",measure(repetitions,[&]{routeOnly(p,workspace.data(),output.data());}));
    if(phase=="all"||phase=="inner"){
        show("inner_contiguous_clmul",measure(repetitions,[&]{inner<false,false>(p,workspace.data(),output.data());}));
        show("inner_routed_clmul",measure(repetitions,[&]{inner<true,false>(p,workspace.data(),output.data());}));
    }
    if(phase=="all"||phase=="full")show("full_clmul",measure(repetitions,[&]{full<false>(p,input.data(),workspace.data(),output.data());}));
    if(phase=="gfni")show("full_gfni",measure(repetitions,[&]{full<true>(p,input.data(),workspace.data(),output.data());}));
    if(phase=="micro"){
        // Fixed4-MiB streams: normalize to the K20 number of256-point batches.
        constexpr std::size_t micro=1<<18;
        Aligned<E> a(micro),b(micro),c(micro);
        for(std::size_t i=0;i<micro;++i){a[i]=random(rng);b[i]=random(rng);}
        const double scale=double(p.batches)/(micro/256);
        auto ms=measure(repetitions,[&]{transformStream(a.data(),c.data(),micro/256);});
        std::cout<<"unpruned_transform_pair_normalized_ms="<<ms*scale<<" checksum="<<checksum(c.data(),micro)<<"\n";
        ms=measure(repetitions,[&]{productStream<false>(a.data(),b.data(),c.data(),micro);});
        std::cout<<"pointwise_normalized_ms="<<ms*scale*keyCount/256.0<<" checksum="<<checksum(c.data(),micro)<<"\n";
    }
}
} // namespace native
int main(int argc,char** argv){
    try{
        if(argc==2&&std::string(argv[1])=="--check"){native::check();return 0;}
        if(argc==5&&std::string(argv[1])=="--bench"){
            native::bench(unsigned(std::stoul(argv[2])),unsigned(std::stoul(argv[3])),argv[4]);return 0;}
        std::cerr<<"usage: toeplitz_native --check | --bench logK repetitions all|outer|route|inner|full|gfni|micro\n";return 2;
    }catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<"\n";return 1;}
}
