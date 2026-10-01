// Minimal retained-driver transplant. New research code; no production changes.
#define SPIN_JOINT_NO_MAIN
#include "joint.cpp"
#undef SPIN_JOINT_NO_MAIN

namespace spin::detail::kernel {
void bchPackedFullGl32(const block*,block*,const std::uint64_t*);
}
namespace pd {
struct Gl32 {
    std::vector<std::uint64_t> coeff;
    explicit Gl32(std::size_t tiles,std::uint64_t seed):coeff(tiles*1024) {
        k::setup::Words rng(seed^0x75a1dc09ULL);
        for(std::size_t t=0;t<tiles;++t)for(unsigned g=0;g<32;++g) {
            std::array<std::uint32_t,32> a{},e{};unsigned rank;
            do {
                for(auto& v:a)v=std::uint32_t(rng());e=a;rank=0;
                for(unsigned c=0;c<32;++c) {
                    unsigned p=rank;while(p<32 && !(e[p]&(1U<<c)))++p;if(p==32)continue;
                    std::swap(e[rank],e[p]);for(unsigned j=rank+1;j<32;++j)if(e[j]&(1U<<c))e[j]^=e[rank];++rank;
                }
            }while(rank!=32);
            for(unsigned d=0;d<4;++d)for(unsigned lane=0;lane<4;++lane) {
                std::uint64_t m=0;for(unsigned j=0;j<8;++j)m|=std::uint64_t((a[8*lane+j]>>(8*((lane+d)%4)))&255)<<(8*(7-j));
                coeff[t*1024+32*g+8*d+2*lane]=m;coeff[t*1024+32*g+8*d+2*lane+1]=m;
            }
        }
    }
};
// The original lambda was sometimes outlined by GCC, depending on the other
// template instantiations in the translation unit. Force this tiny fixed-width
// emitter into each statically expanded inner point, retaining its SIMD shape.
template<bool Field> struct Emitter {
    k::block* scratch;
    const unsigned* bases;
    const std::uint64_t* controls;
    __m128i v1{},v2{},v3{};
    SPIN_FORCEINLINE void operator()(std::size_t i,k::block v) {
        switch(i&3) {
            case 3:v3=v.mData;break;case 2:v2=v.mData;break;case 1:v1=v.mData;break;
            case 0:{
                auto x=_mm512_castsi128_si512(v.mData);
                x=_mm512_inserti32x4(x,v1,1);x=_mm512_inserti32x4(x,v2,2);x=_mm512_inserti32x4(x,v3,3);
                if constexpr(Field)x=spin::research::gf16::transpose(x,std::uint32_t(controls[i/4]));
                _mm512_stream_si512(reinterpret_cast<__m512i*>(scratch+bases[i/4]),x);break;
            }
        }
    }
};
template<bool Mix,unsigned R=2> static SPIN_NOINLINE void encode(k::block* input,k::block* scratch,std::size_t n,
    const unsigned* masks,const Packets<4>& packets,const Gl32& setup,double* phases) {
    static_assert(R>=2 && R<=4);
    const auto start=phases?Clock::now():Clock::time_point{};
    Emitter<!Mix> emit{scratch,packets.bases.data(),packets.controls.data()};
    reverseInner<R>(input,n,masks,emit);_mm_sfence();
    const auto middle=phases?Clock::now():Clock::time_point{};
    for(std::size_t tile=0;tile<n/1024;++tile) {
        if constexpr(Mix)k::bchPackedFullGl32(scratch+tile*tileStride,input+tile*512,setup.coeff.data()+tile*1024);
        else k::bchTranspose4GfniBlend(scratch+tile*tileStride,input+tile*512);
    }
    if(phases){phases[0]+=std::chrono::duration<double,std::milli>(middle-start).count();phases[1]+=std::chrono::duration<double,std::milli>(Clock::now()-middle).count();}
}
static void scalarMix(const k::block* in,k::block* out,const std::uint64_t* coeff) {
    for(unsigned g=0;g<32;++g)for(unsigned row=0;row<32;++row) {
        const unsigned lane=row/8,j=row%8;std::uint32_t mask=0;
        for(unsigned d=0;d<4;++d)mask|=std::uint32_t((coeff[32*g+8*d+2*lane]>>(8*(7-j)))&255)<<(8*((lane+d)%4));
        auto sum=_mm_setzero_si128();while(mask){const auto c=std::countr_zero(mask);sum=_mm_xor_si128(sum,in[4*(8*g+c%8)+c/8].mData);mask&=mask-1;}
        out[4*(8*g+j)+lane]=k::block(sum);
    }
}
template<unsigned Mode,unsigned R=2> static void experiment(unsigned exponent,std::uint64_t seed,unsigned calls,bool profile) {
    static_assert(Mode<=2 && R>=2 && R<=4);
    const unsigned n=2U<<exponent;
    const Packets<4> packets(n/256,seed,true,true);
    std::vector<unsigned> masks(std::size_t(2)*R*(n/128));k::setup::Words rng(seed^0x3f625a92ULL);
    for(std::size_t i=0;i<masks.size();i+=2) {
        unsigned u;do{u=unsigned(rng())&((1U<<19)-1);}while(!u);unsigned v=unsigned(rng())&((1U<<19)-1);
        if(std::popcount(u&v)&1)v^=u&-u;masks[i]=u;masks[i+1]=v;
    }
    // Match the retained driver's vector allocation order and input alignment.
    // Only routing scratch requests the workspace advice, as in joint.cpp.
    std::vector<k::block> input(n),inner(n),actual(n),other(n),forward(n),expected(n),materialized(n);
    const auto scratchBlocks=(n/1024)*tileStride;std::vector<k::block> storage(scratchBlocks+4);
    auto* scratch=reinterpret_cast<k::block*>((reinterpret_cast<std::uintptr_t>(storage.data())+63)&~std::uintptr_t(63));
    k::workspace_routing::adviseOwned(scratch,scratchBlocks*16);
    const Gl32 setup(n/1024,seed);
    auto selected=[&](double* phases) {
        if constexpr(Mode==0)::encode<4,R,true,true,true>(input.data(),scratch,n,masks.data(),packets,phases);
        else encode<Mode==2,R>(input.data(),scratch,n,masks.data(),packets,setup,phases);
    };
    fill(other,971);forwardReference<R>(other.data(),forward.data(),n,masks.data());
    for(unsigned pattern=0;pattern<(calls?1U:3U);++pattern) {
        fill(input,pattern==2?997:913);if(pattern==1){std::fill(input.begin(),input.end(),k::block{});input[n-3]=k::block(1,3);}expected=input;
        reverseReference<R>(input.data(),inner.data(),n,masks.data());
        reverseInner<R>(input.data(),n,masks.data(),[&](std::size_t i,k::block v){actual[i]=v;});
        if(std::memcmp(actual.data(),inner.data(),n*16))throw std::runtime_error("dense inner mismatch");
        auto lhs=_mm_setzero_si128(),rhs=lhs;
        for(std::size_t i=0;i<n;++i){lhs=_mm_xor_si128(lhs,_mm_and_si128(input[i].mData,forward[i].mData));rhs=_mm_xor_si128(rhs,_mm_and_si128(other[i].mData,inner[i].mData));}
        if(_mm_movemask_epi8(_mm_cmpeq_epi8(lhs,rhs))!=65535)throw std::runtime_error("inner adjoint mismatch");
        for(std::size_t i=0;i<n;++i) {
            auto value=inner[i].mData;
            if constexpr(Mode!=2) {
                value=_mm_setzero_si128();const auto col=spin::research::gf16::multiply(packets.multipliers[i/4],1U<<(i&3));
                for(unsigned bit=0;bit<4;++bit)if((col>>bit)&1)value=_mm_xor_si128(value,inner[(i&~std::size_t(3))+bit].mData);
            }
            const auto x=packets.inverse[i];materialized[(x/1024)*1024+4*(x%256)+(x/256)%4]=k::block(value);
        }
        alignas(64)k::block mixed[1024];
        for(std::size_t tile=0;tile<n/1024;++tile) {
            if constexpr(Mode==2){scalarMix(materialized.data()+1024*tile,mixed,setup.coeff.data()+1024*tile);k::bchTranspose4(mixed,expected.data()+512*tile);}
            else k::bchTranspose4(materialized.data()+1024*tile,expected.data()+512*tile);
        }
        if constexpr(Mode==1){actual=input;::encode<4,R,true,true,true>(actual.data(),scratch,n,masks.data(),packets,nullptr);}
        selected(nullptr);
        if(std::memcmp(input.data(),expected.data(),n*16))throw std::runtime_error("independent full/suffix reference mismatch");
        if constexpr(Mode==1)if(std::memcmp(input.data(),actual.data(),n*16))throw std::runtime_error("forced-inline differs from legacy");
    }
    std::cout<<"checks: independent dense inner, inner adjoint, scalar route/mixer, production BCH, full buffer/suffix PASS; updates="<<R<<'\n';
    if(!calls)return;
    fill(input,913);double phases[2]{};auto run=[&]{selected(profile?phases:nullptr);};
    for(unsigned i=0;i<3;++i)run();std::vector<double> times;times.reserve(calls);
    for(unsigned i=0;i<calls;++i){const auto start=Clock::now();run();times.push_back(std::chrono::duration<double,std::milli>(Clock::now()-start).count());}
    std::sort(times.begin(),times.end());
    std::cout<<"exponent,mode,seed,calls,median_ms,p10_ms,p90_ms,checksum,updates\n"<<exponent<<','<<Mode<<','<<seed<<','<<calls<<','<<std::fixed<<std::setprecision(6)<<times[calls/2]<<','<<times[calls/10]<<','<<times[9*calls/10]<<','<<std::hex<<hash(input)<<std::dec<<','<<R<<'\n';
    if(profile)std::cout<<"phase_means_including_warmups,inner_route,"<<phases[0]/(calls+3)<<",bch,"<<phases[1]/(calls+3)<<'\n';
}
// Select the template once, outside setup, checks, and all timed loops.
template<unsigned R> static void dispatch(unsigned exponent,unsigned mode,std::uint64_t seed,unsigned calls,bool profile) {
    if(mode==0)experiment<0,R>(exponent,seed,calls,profile);
    else if(mode==1)experiment<1,R>(exponent,seed,calls,profile);
    else experiment<2,R>(exponent,seed,calls,profile);
}
}
int main(int argc,char** argv){try {
    if(argc<5||argc>7)throw std::invalid_argument("usage: packed_driver exponent mode seed calls [profile [updates=2]]");
    const unsigned exponent=std::stoul(argv[1]),mode=std::stoul(argv[2]),calls=std::stoul(argv[4]);const auto seed=std::stoull(argv[3]);const bool profile=argc>=6&&std::stoul(argv[5]);
    const auto updates=argc==7?std::stoul(argv[6]):2UL;
    if(exponent<14||exponent>20||mode>2||updates<2||updates>4)throw std::invalid_argument("unsupported case (updates must be 2, 3, or 4)");
    if(!__builtin_cpu_supports("avx512f")||!__builtin_cpu_supports("avx512vl")||!__builtin_cpu_supports("avx512bw")||!__builtin_cpu_supports("gfni"))throw std::runtime_error("AVX512/GFNI required");
    checkGf16();if(updates==2)pd::dispatch<2>(exponent,mode,seed,calls,profile);else if(updates==3)pd::dispatch<3>(exponent,mode,seed,calls,profile);else pd::dispatch<4>(exponent,mode,seed,calls,profile);
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
