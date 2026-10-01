// Exact-map optimization experiment. Neither certified sources nor the prior
// winning implementation are edited. All dispatch is compile-time in hot code.
#include "PacketReference.h"
#include "InnerWideState.h"
#include "InnerPacketStream.h"
#include "InnerDirectFeedback.h"
#include "InnerVbmiState.h"
namespace fusionprobe {
using namespace packetprobe;

template<int Mode> struct Options {
    static constexpr bool pack=Mode==1 || Mode==2 || Mode==9 || Mode==10;
    static constexpr bool unpack=Mode==1 || Mode==3 || Mode==5 || Mode==9 || Mode==10 || (Mode>=12 && Mode<=14);
    static constexpr bool direct=Mode==4 || Mode==5 || (Mode>=11 && Mode<=14);
    static constexpr bool stream=(Mode>=6 && Mode<=15) || (Mode>=18 && Mode<=21) || Mode>=24;
    static constexpr bool incremental=Mode!=7 && Mode!=10 && Mode!=13 && Mode!=19 && Mode<24;
    static constexpr bool pruned=Mode==8 || Mode==14 || Mode==21 || Mode==24 || Mode==26;
    static constexpr bool early=Mode==15 || Mode==16 || Mode==20 || Mode==21 || Mode>=25;
    static constexpr bool vbmiPack=(Mode>=17 && Mode<=20) || Mode==22 || Mode>=24;
    static constexpr bool vbmiUnpack=(Mode>=17 && Mode<=20) || Mode==23 || Mode>=24;
};

template<int Mode,class Emit>
SPIN_FORCEINLINE void reverseNew(const k::block* input,std::size_t n,
    const gf::Dense16Row* rows,Emit&& emit) {
    using O=Options<Mode>;
    alignas(64) __m128i words[16],moments[64],syndrome[16];
    alignas(64) __m512i packets[16],high[16];
    gf::Packed<16> state{},feedback;
    for(std::size_t epoch=n/64;epoch-->0;) {
        const auto* raw=input+64*epoch;
        const bool first=epoch+1==n/64;
        if(first) {
            loadPackets(raw,packets,std::make_index_sequence<16>{});
            for(unsigned h=16;h-->0;)emit(16*epoch+h,packets[h]);
        } else {
            if constexpr(O::vbmiUnpack)ip::wideUnpackVbmi(state,words);
            else if constexpr(O::unpack)ip::wideUnpack(state,words);
            else gf::unpack(state,words);
            // Emission consumes old coordinate words, whereas the next packed
            // state can already compute its independent matrix product.
            if constexpr(O::early)if(epoch)gf::denseStep16<false>(state,rows[epoch]);
            if constexpr(O::stream) {
                if constexpr(O::direct)
                    ip::streamStepHigh<O::pruned,O::incremental>(words,raw,16*epoch,emit,high);
                else ip::streamStep<O::pruned,O::incremental>(words,raw,16*epoch,emit,moments);
            } else {
                ip::Maps<true>::template emission<false>(words,packets);
                outputPackets<false>(raw,packets,moments,64*epoch,emit,std::make_index_sequence<16>{});
            }
        }
        if(!epoch)break;
        if constexpr(O::direct) {
            if constexpr(O::stream) {
                if(first)ip::directFeedback(packets,feedback);
                else ip::directFeedbackHigh(high,feedback);
            } else ip::directFeedback(packets,feedback);
        } else {
            if constexpr(O::stream) {if(first)ip::packetMoments(packets,moments);}
            else ip::packetMoments(packets,moments);
            ip::Maps<true>::finish(moments,syndrome);
            if constexpr(O::vbmiPack)ip::widePackVbmi(syndrome,feedback);
            else if constexpr(O::pack)ip::widePack(syndrome,feedback);
            else gf::pack<16>(syndrome,feedback);
        }
        if(first)state=feedback;
        else if constexpr(O::early) {
            state.v[0]=_mm512_xor_si512(state.v[0],feedback.v[0]);
            state.v[1]=_mm512_xor_si512(state.v[1],feedback.v[1]);
            state.v[2]=_mm512_xor_si512(state.v[2],feedback.v[2]);
            state.v[3]=_mm512_xor_si512(state.v[3],feedback.v[3]);
        } else gf::denseStep16<true>(state,rows[epoch],&feedback);
    }
}

template<int Mode> static SPIN_NOINLINE void encodeNew(k::block* input,k::block* scratch,std::size_t n,
    const Packets<4>& route,const Coefficients& coeff,const gf::Dense16Row* rows,double* phases) {
    const auto start=phases?Clock::now():Clock::time_point{};
    Route emit{scratch,route.bases.data()};
    reverseNew<Mode>(input,n,rows,emit);_mm_sfence();
    const auto middle=phases?Clock::now():Clock::time_point{};
    for(std::size_t tile=0;tile<n/1024;++tile)
        k::bchPackedCoeffCompact(scratch+tile*tileStride,input+tile*512,coeff.compact+512*tile);
    if(phases) {
        phases[0]+=std::chrono::duration<double,std::milli>(middle-start).count();
        phases[1]+=std::chrono::duration<double,std::milli>(Clock::now()-middle).count();
    }
}

template<int Mode> static void checkInnerNew(const Updates& setup,const Setup& transformed) {
    constexpr unsigned n=256;std::vector<k::block> input(n),expected(n),actual(n);
    for(unsigned bit=0;bit<n;++bit) {
        std::fill(input.begin(),input.end(),k::block{});input[bit]=k::block(0x729aca378bc42fedULL,0x942163aca892317bULL);
        innerOracle<true>(input.data(),expected.data(),n,setup);
        auto emit=[&](std::size_t i,__m512i x){_mm512_storeu_si512(actual.data()+4*i,x);};
        reverseNew<Mode>(input.data(),n,transformed.packed.data(),emit);
        equal(actual.data(),expected.data(),n*16,"four-epoch exhaustive input basis");
    }
}

static void checkFeedback() {
    for(unsigned coordinate=0;coordinate<64;++coordinate)for(unsigned bit=0;bit<128;++bit) {
        alignas(64) __m512i packets[16]{};
        auto* bytes=reinterpret_cast<unsigned char*>(packets);
        bytes[16*coordinate+bit/8]=static_cast<unsigned char>(1U<<(bit%8));
        alignas(64) __m128i moments[64]{},words[16];
        gf::Packed<16> actual,expected;
        ip::packetMoments(packets,moments);ip::Maps<true>::finish(moments,words);
        gf::pack<16>(words,expected);ip::directFeedback(packets,actual);
        equal(&actual,&expected,sizeof(actual),"direct packed feedback complete bit basis");
    }
}

template<int Mode> static void experiment(unsigned exponent,std::uint64_t seed,unsigned calls,bool profile) {
    const unsigned n=2U<<exponent;const Packets<4> route(n/256,seed,true,true);
    if(!calls) {
        if(!ip::wideStateSelfCheck())throw std::runtime_error("wide state full bit basis");
        if(!ip::vbmiStateSelfCheck())throw std::runtime_error("VBMI state full bit basis");
        checkFeedback();
    }
    const Updates setup(n/64,seed);const Setup transformed(setup);
    const pd::Gl32 gl(n/1024,seed);const Coefficients coeff(gl);
    std::vector<k::block> input(n),expected(n),initial(n),guarded((n/1024)*tileStride+12);
    auto* scratch=reinterpret_cast<k::block*>((reinterpret_cast<std::uintptr_t>(guarded.data()+4)+63)&~std::uintptr_t(63));
    const auto scratchBlocks=(n/1024)*tileStride;
    k::workspace_routing::adviseOwned(scratch,scratchBlocks*16);
    if constexpr(Mode!=0)checkInnerNew<Mode>(setup,transformed);
    auto selected=[&](double* phases) {
        if constexpr(Mode==0)packetprobe::encode<true,2,false,false,false>(input.data(),scratch,n,route,coeff,transformed.packed.data(),phases);
        else encodeNew<Mode>(input.data(),scratch,n,route,coeff,transformed.packed.data(),phases);
    };
    const k::block canary(0x781aca29312357bfULL,0x928fafcac352697bULL);
    for(unsigned pattern=0;pattern<(calls?1U:4U);++pattern) {
        fill(input,pattern==2?997:913);
        if(pattern==1){std::fill(input.begin(),input.end(),k::block{});input[n-3]=k::block(1,3);}
        if(pattern==3){std::fill(input.begin(),input.end(),k::block{});input[63]=k::block(1,2);input[64]=k::block(4,8);input[n-65]=k::block(16,32);}
        expected=input;initial=input;
        t64probe::encode(expected.data(),scratch,n,route,coeff,setup.packed.data(),nullptr);
        scratch[-1]=canary;scratch[scratchBlocks]=canary;selected(nullptr);
        equal(input.data(),expected.data(),n*16,"complete retained encoder/output suffix");
        equal(scratch-1,&canary,16,"leading scratch canary");equal(scratch+scratchBlocks,&canary,16,"trailing scratch canary");
        if(!calls) {
            innerOracle<true>(initial.data(),expected.data(),n,setup);
            for(std::size_t i=0;i<n;++i)initial[canonical(route.inverse[i])]=expected[i];
            alignas(64) k::block mixed[1024];
            for(std::size_t tile=0;tile<n/1024;++tile) {
                pd::scalarMix(initial.data()+1024*tile,mixed,gl.coeff.data()+1024*tile);
                k::bchTranspose4(mixed,expected.data()+512*tile);
            }
            equal(input.data(),expected.data(),n*8,"independent complete scalar reference");
        }
    }
    std::cout<<"checks PASS mode="<<Mode<<" seed="<<seed<<" exponent="<<exponent
             <<(Mode?" four-epoch basis;":"")<<" complete retained output/suffix; scratch canaries"
             <<(calls?"":"; exhaustive helpers and independent scalar complete encoder")<<'\n';
    if(!calls)return;
    fill(input,913);double phases[2]{};
    for(unsigned i=0;i<3;++i)selected(profile?phases:nullptr);
    std::vector<double> times;times.reserve(calls);
    for(unsigned i=0;i<calls;++i){const auto start=Clock::now();selected(profile?phases:nullptr);
        times.push_back(std::chrono::duration<double,std::milli>(Clock::now()-start).count());}
    std::sort(times.begin(),times.end());
    std::cout<<"mode,exponent,seed,calls,median_ms,p10_ms,p90_ms,checksum\n"
             <<Mode<<','<<exponent<<','<<seed<<','<<calls<<','<<std::fixed<<std::setprecision(6)
             <<times[calls/2]<<','<<times[calls/10]<<','<<times[9*calls/10]<<','<<std::hex<<hash(input)<<std::dec<<'\n';
    if(profile)std::cout<<"phase_means_including_warmups,inner_route,"<<phases[0]/(calls+3)
                      <<",packed_mixer_bch,"<<phases[1]/(calls+3)<<'\n';
}
}
int main(int argc,char** argv) {try {
    if(argc<5 || argc>6)throw std::invalid_argument("usage: inner-fusion mode exponent seed calls [profile]");
    if(!__builtin_cpu_supports("avx512f") || !__builtin_cpu_supports("avx512vl") ||
       !__builtin_cpu_supports("avx512bw") || !__builtin_cpu_supports("gfni") ||
       !__builtin_cpu_supports("avx512vbmi"))throw std::runtime_error("AVX512/GFNI/VBMI required");
    const unsigned mode=std::stoul(argv[1]),exponent=std::stoul(argv[2]),calls=std::stoul(argv[4]);
    const auto seed=std::stoull(argv[3]);const bool profile=argc==6 && std::stoul(argv[5]);
    if(exponent<14 || exponent>20)throw std::invalid_argument("exponent range 14..20");
    switch(mode) {
#define FUSION_CASE(M) case M:fusionprobe::experiment<M>(exponent,seed,calls,profile);break
        FUSION_CASE(0);FUSION_CASE(1);FUSION_CASE(2);FUSION_CASE(3);FUSION_CASE(4);
        FUSION_CASE(5);FUSION_CASE(6);FUSION_CASE(7);FUSION_CASE(8);FUSION_CASE(9);
        FUSION_CASE(10);FUSION_CASE(11);FUSION_CASE(12);FUSION_CASE(13);FUSION_CASE(14);
        FUSION_CASE(15);FUSION_CASE(16);FUSION_CASE(17);FUSION_CASE(18);FUSION_CASE(19);
        FUSION_CASE(20);FUSION_CASE(21);FUSION_CASE(22);FUSION_CASE(23);
        FUSION_CASE(24);FUSION_CASE(25);FUSION_CASE(26);
#undef FUSION_CASE
        default:throw std::invalid_argument("mode range 0..26");
    }
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
