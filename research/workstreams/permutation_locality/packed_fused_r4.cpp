// Separate exact-map optimization probe; the retained R4 source is unchanged.
#define main packed_fused_r4_retained_driver_main
#include "packed_driver.cpp"
#undef main
#include "FusedR4.h"

namespace fr=spin::research::fused_r4;
namespace fused_probe {
struct Setup {
    std::vector<fr::ScalarRow> scalar;
    std::vector<fr::WideRow> wide;
    Setup(const unsigned* masks,std::size_t epochs):scalar(epochs),wide(epochs) {
        for(std::size_t epoch=0;epoch<epochs;++epoch) {
            scalar[epoch]=fr::makeScalar(masks+8*epoch);
            wide[epoch]=fr::makeWide(scalar[epoch]);
        }
    }
};
template<class Row,class Emit> SPIN_FORCEINLINE void reverseInner(
    const k::block* in,std::size_t n,const Row* rows,Emit&& emit) {
    // Identical epoch order, fixed emission circuit, zeta, and feedback order.
    alignas(32) __m128i state[19]{},values[128],syndrome[19];
    for(std::size_t epoch=n/128;epoch-->0;) {
        const auto base=128*epoch;
        if(epoch+1==n/128) {
            for(unsigned p=128;p-->0;){values[p]=in[base+p].mData;emit(base+p,in[base+p]);}
        } else k::imtReversePoints<Map>(in+base,values,state,base,emit,std::make_index_sequence<128>{});
        if(!epoch)break;
        k::zeta<128>(values);Map::finish(values,syndrome);
        if(epoch+1==n/128)std::memcpy(state,syndrome,sizeof(state));
        else fr::step<true>(state,rows[epoch],syndrome);
    }
}
template<unsigned Mode> static SPIN_NOINLINE void encode(
    k::block* input,k::block* scratch,std::size_t n,const unsigned*,
    const Packets<4>& packets,const pd::Gl32& gl,const Setup& setup,double* phases) {
    static_assert(Mode==1 || Mode==2);
    const auto start=phases?Clock::now():Clock::time_point{};
    pd::Emitter<false> emit{scratch,packets.bases.data(),packets.controls.data()};
    if constexpr(Mode==1)reverseInner(input,n,setup.scalar.data(),emit);
    else reverseInner(input,n,setup.wide.data(),emit);
    _mm_sfence();
    const auto middle=phases?Clock::now():Clock::time_point{};
    for(std::size_t tile=0;tile<n/1024;++tile)
        k::bchPackedFullGl32(scratch+tile*tileStride,input+tile*512,gl.coeff.data()+tile*1024);
    if(phases) {
        phases[0]+=std::chrono::duration<double,std::milli>(middle-start).count();
        phases[1]+=std::chrono::duration<double,std::milli>(Clock::now()-middle).count();
    }
}
static void equal(const void* a,const void* b,std::size_t bytes,const char* why) {
    if(std::memcmp(a,b,bytes))throw std::runtime_error(why);
}
static void checkAdjoint(const k::block* x,const k::block* y,
    const k::block* forward,const k::block* reverse,std::size_t n) {
    auto lhs=_mm_setzero_si128(),rhs=lhs;
    for(std::size_t i=0;i<n;++i) {
        lhs=_mm_xor_si128(lhs,_mm_and_si128(x[i].mData,forward[i].mData));
        rhs=_mm_xor_si128(rhs,_mm_and_si128(y[i].mData,reverse[i].mData));
    }
    if(_mm_movemask_epi8(_mm_cmpeq_epi8(lhs,rhs))!=65535)
        throw std::runtime_error("original forward / fused reverse adjoint mismatch");
}
static void checkMaps(const unsigned* masks,const Setup& setup,std::size_t epochs) {
    // Encode all 19 basis vectors simultaneously in the low 19 bits of each
    // SIMD word; its high half supplies an additional independent dense input.
    // Every epoch and both orientations are checked against ORIGINAL masks.
    k::setup::Words rng(0xcaa4635U);
    for(std::size_t epoch=0;epoch<epochs;++epoch) {
        alignas(32) __m128i input[19],expected[19],scalar[19],wide[19];
        for(unsigned j=0;j<19;++j)input[j]=k::block(rng(),1ULL<<j).mData;
        for(unsigned transpose=0;transpose<2;++transpose) {
            std::memcpy(expected,input,sizeof(input));
            std::memcpy(scalar,input,sizeof(input));std::memcpy(wide,input,sizeof(input));
            const auto* m=masks+8*epoch;
            for(unsigned r=0;r<4;++r) {
                const auto i=transpose?3-r:r;
                const unsigned dot=m[2*i+(transpose?0:1)],out=m[2*i+(transpose?1:0)];
                auto d=_mm_setzero_si128();
                for(unsigned j=0;j<19;++j)if((dot>>j)&1)d=_mm_xor_si128(d,expected[j]);
                for(unsigned j=0;j<19;++j)if((out>>j)&1)expected[j]=_mm_xor_si128(expected[j],d);
            }
            const auto row=transpose?setup.scalar[epoch]:fr::makeScalar(m,false);
            const auto packed=transpose?setup.wide[epoch]:fr::makeWide(row);
            fr::step(scalar,row);fr::step(wide,packed);
            equal(expected,scalar,sizeof(input),"scalar fused basis-map mismatch");
            equal(expected,wide,sizeof(input),"wide fused basis-map mismatch");
        }
    }
}
template<unsigned Mode> static void experiment(unsigned exponent,std::uint64_t seed,unsigned calls,bool profile) {
    static_assert(Mode<=2);
    const unsigned n=2U<<exponent;
    const Packets<4> packets(n/256,seed,true,true);
    std::vector<unsigned> masks(std::size_t(8)*(n/128));k::setup::Words rng(seed^0x3f625a92ULL);
    for(std::size_t i=0;i<masks.size();i+=2) {
        unsigned u;do{u=unsigned(rng())&((1U<<19)-1);}while(!u);
        unsigned v=unsigned(rng())&((1U<<19)-1);
        if(std::popcount(u&v)&1)v^=u&-u;masks[i]=u;masks[i+1]=v;
    }
    // Preserve the baseline's main vector allocation order and scratch alignment.
    std::vector<k::block> input(n),inner(n),actual(n),other(n),forward(n),expected(n),materialized(n);
    const auto scratchBlocks=(n/1024)*tileStride;std::vector<k::block> storage(scratchBlocks+4);
    auto* scratch=reinterpret_cast<k::block*>((reinterpret_cast<std::uintptr_t>(storage.data())+63)&~std::uintptr_t(63));
    k::workspace_routing::adviseOwned(scratch,scratchBlocks*16);
    const pd::Gl32 gl(n/1024,seed);
    const Setup setup(masks.data(),n/128);
    checkMaps(masks.data(),setup,n/128);
    auto selected=[&](double* phases) {
        if constexpr(Mode==0)pd::encode<true,4>(input.data(),scratch,n,masks.data(),packets,gl,phases);
        else encode<Mode>(input.data(),scratch,n,masks.data(),packets,gl,setup,phases);
    };
    fill(other,971);forwardReference<4>(other.data(),forward.data(),n,masks.data());
    for(unsigned pattern=0;pattern<(calls?1U:3U);++pattern) {
        fill(input,pattern==2?997:913);
        if(pattern==1){std::fill(input.begin(),input.end(),k::block{});input[n-3]=k::block(1,3);}
        expected=input;
        reverseReference<4>(input.data(),inner.data(),n,masks.data());
        ::reverseInner<4>(input.data(),n,masks.data(),[&](std::size_t i,k::block v){actual[i]=v;});
        equal(actual.data(),inner.data(),n*16,"original dense/sequential inner mismatch");
        reverseInner(input.data(),n,setup.scalar.data(),[&](std::size_t i,k::block v){actual[i]=v;});
        equal(actual.data(),inner.data(),n*16,"scalar fused / original dense inner mismatch");
        checkAdjoint(input.data(),other.data(),forward.data(),actual.data(),n);
        reverseInner(input.data(),n,setup.wide.data(),[&](std::size_t i,k::block v){actual[i]=v;});
        equal(actual.data(),inner.data(),n*16,"wide fused / original dense inner mismatch");
        checkAdjoint(input.data(),other.data(),forward.data(),actual.data(),n);
        for(std::size_t i=0;i<n;++i) {
            const auto x=packets.inverse[i];
            materialized[(x/1024)*1024+4*(x%256)+(x/256)%4]=inner[i];
        }
        alignas(64) k::block mixed[1024];
        for(std::size_t tile=0;tile<n/1024;++tile) {
            pd::scalarMix(materialized.data()+1024*tile,mixed,gl.coeff.data()+1024*tile);
            k::bchTranspose4(mixed,expected.data()+512*tile);
        }
        // Full encoders share ONLY the route/GL32/BCH, never transformed masks
        // on the oracle side. All n words are compared, including untouched suffix.
        actual=input;pd::encode<true,4>(actual.data(),scratch,n,masks.data(),packets,gl,nullptr);
        equal(actual.data(),expected.data(),n*16,"original full/dense reference mismatch");
        actual=input;encode<1>(actual.data(),scratch,n,masks.data(),packets,gl,setup,nullptr);
        equal(actual.data(),expected.data(),n*16,"scalar fused full/reference/suffix mismatch");
        actual=input;encode<2>(actual.data(),scratch,n,masks.data(),packets,gl,setup,nullptr);
        equal(actual.data(),expected.data(),n*16,"wide fused full/reference/suffix mismatch");
        selected(nullptr);
        equal(input.data(),expected.data(),n*16,"selected full/reference/suffix mismatch");
    }
    std::cout<<"checks: every-epoch forward+transpose basis maps, original sequential+dense inner, "
        "both fused adjoints, original+both fused full encoder, scalar GL32/production BCH, suffix PASS; "
        "updates=4; epochs="<<n/128<<"; scalar_row_bytes="<<sizeof(fr::ScalarRow)
        <<"; wide_row_bytes="<<sizeof(fr::WideRow)<<'\n';
    if(!calls)return;
    fill(input,913);double phases[2]{};auto run=[&]{selected(profile?phases:nullptr);};
    for(unsigned i=0;i<3;++i)run();std::vector<double> times;times.reserve(calls);
    for(unsigned i=0;i<calls;++i) {
        const auto start=Clock::now();run();
        times.push_back(std::chrono::duration<double,std::milli>(Clock::now()-start).count());
    }
    std::sort(times.begin(),times.end());
    std::cout<<"exponent,mode,seed,calls,median_ms,p10_ms,p90_ms,checksum,updates\n"
        <<exponent<<','<<Mode<<','<<seed<<','<<calls<<','<<std::fixed<<std::setprecision(6)
        <<times[calls/2]<<','<<times[calls/10]<<','<<times[9*calls/10]<<','<<std::hex
        <<hash(input)<<std::dec<<",4\n";
    if(profile)std::cout<<"phase_means_including_warmups,inner_route,"<<phases[0]/(calls+3)
        <<",bch,"<<phases[1]/(calls+3)<<'\n';
}
}
int main(int argc,char** argv) {try {
    if(argc<4 || argc>6)throw std::invalid_argument(
        "usage: packed_fused_r4 exponent seed calls [mode=0 [profile=0]]; modes 0=original R4,1=fused scalar table,2=fused AVX512 table; calls=0 checks only");
    const unsigned exponent=std::stoul(argv[1]),calls=std::stoul(argv[3]);
    const auto seed=std::stoull(argv[2]);const unsigned mode=argc>=5?std::stoul(argv[4]):0;
    const bool profile=argc==6 && std::stoul(argv[5]);
    if(exponent<14 || exponent>20 || mode>2)throw std::invalid_argument("unsupported case");
    if(!__builtin_cpu_supports("avx512f") || !__builtin_cpu_supports("avx512vl") ||
       !__builtin_cpu_supports("avx512bw") || !__builtin_cpu_supports("gfni"))
        throw std::runtime_error("AVX512/GFNI required");
    if(mode==0)fused_probe::experiment<0>(exponent,seed,calls,profile);
    else if(mode==1)fused_probe::experiment<1>(exponent,seed,calls,profile);
    else fused_probe::experiment<2>(exponent,seed,calls,profile);
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
