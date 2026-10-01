// Exact-map small-size tuning. The setup distribution and encoded map are fixed.
#include "PacketReference.h"
#include "PacketInnerKernel.h"
#include "PacketInnerBoundaryTests.h"
namespace spin::detail::kernel {
void bchPackedCoeffAligned(const block*,block*,const std::uint64_t*);
}
namespace sizeprobe {
using namespace packetprobe;

struct ExpandedCoefficients {
    std::vector<std::uint64_t> storage;
    std::uint64_t* data;
    explicit ExpandedCoefficients(const pd::Gl32& original):storage(original.coeff.size()+8),
        data(reinterpret_cast<std::uint64_t*>((reinterpret_cast<std::uintptr_t>(storage.data())+63)&~std::uintptr_t(63))) {
        std::memcpy(data,original.coeff.data(),original.coeff.size()*sizeof(std::uint64_t));
    }
};

template<bool Cached,unsigned Lead=0> struct Route {
    static_assert(!Lead || Cached,"do not bring NT destinations into cache");
    k::block* scratch;const unsigned* bases;
    SPIN_FORCEINLINE void operator()(std::size_t i,__m512i value) {
        if constexpr(Lead)if(i>=Lead)__builtin_prefetch(scratch+bases[i-Lead],1,3);
        if constexpr(Cached)_mm512_store_si512(scratch+bases[i],value);
        else _mm512_stream_si512(reinterpret_cast<__m512i*>(scratch+bases[i]),value);
    }
};

template<unsigned Mode> struct Options {
    static constexpr unsigned recipe=Mode<64?Mode%4:4+(Mode-64)%3;
    static constexpr bool cached=Mode>=64 || (Mode&4)!=0,expanded=Mode<64?(Mode&8)!=0:Mode>=67;
    static constexpr unsigned lead=Mode>=64?0:(Mode/16==0?0:(Mode/16==1?8:(Mode/16==2?16:32)));
    static constexpr ip::InnerRecipe inner=recipe==0?ip::InnerRecipe::Packet:
        recipe==1?ip::InnerRecipe::PacketVbmi:recipe==2?ip::InnerRecipe::StreamVbmi:
        recipe==3?ip::InnerRecipe::StreamLegacy:recipe==4?ip::InnerRecipe::StreamVbmiEarly:
        recipe==5?ip::InnerRecipe::StreamVbmiPruned:ip::InnerRecipe::StreamVbmiPrunedEarly;
};

template<unsigned Mode> static SPIN_NOINLINE void encodeSize(k::block* input,k::block* scratch,std::size_t n,
    const Packets<4>& route,const Coefficients& coeff,const ExpandedCoefficients& expanded,
    const ip::PreparedPacketUpdates16& updates,double* phases) {
    using O=Options<Mode>;
    const auto start=phases?Clock::now():Clock::time_point{};
    Route<O::cached,O::lead> emit{scratch,route.bases.data()};
    // These compile-time recipes are the same kernels at every supported size.
    ip::transposeInner<O::inner>(input,n,updates.conjugatedRows(),emit);
    if constexpr(!O::cached)_mm_sfence();
    const auto middle=phases?Clock::now():Clock::time_point{};
    for(std::size_t tile=0;tile<n/1024;++tile) {
        if constexpr(O::expanded)k::bchPackedCoeffAligned(scratch+tile*tileStride,input+tile*512,expanded.data+1024*tile);
        else k::bchPackedCoeffCompact(scratch+tile*tileStride,input+tile*512,coeff.compact+512*tile);
    }
    if(phases) {
        phases[0]+=std::chrono::duration<double,std::milli>(middle-start).count();
        phases[1]+=std::chrono::duration<double,std::milli>(Clock::now()-middle).count();
    }
}

template<unsigned Mode> static void checkInnerSize(const Updates& setup,const ip::PreparedPacketUpdates16& transformed) {
    constexpr unsigned n=256;std::vector<k::block> input(n),expected(n),actual(n);
    ip::PreparedPacketUpdates16 shortUpdates(std::span<const Matrix>(setup.reverse.data(),n/64));
    for(unsigned coordinate=0;coordinate<n;++coordinate) {
        std::fill(input.begin(),input.end(),k::block{});input[coordinate]=k::block(0x729aca378bc42fedULL,0x942163aca892317bULL);
        innerOracle<true>(input.data(),expected.data(),n,setup);
        auto emit=[&](std::size_t i,__m512i value){_mm512_storeu_si512(actual.data()+4*i,value);};
        using O=Options<Mode>;
        ip::transposeInner<O::inner>(input.data(),n,shortUpdates.conjugatedRows(),emit);
        equal(actual.data(),expected.data(),n*16,"four-epoch basis");
    }
}

template<unsigned Mode> static void experiment(unsigned exponent,std::uint64_t seed,unsigned calls,bool profile) {
    using O=Options<Mode>;
    const unsigned n=2U<<exponent;const Packets<4> route(n/256,seed,true,true);
    const Updates setup(n/64,seed);const Setup retained(setup);
    const ip::PreparedPacketUpdates16 transformed(setup.reverse);
    equal(transformed.conjugatedRows().data(),retained.packed.data(),retained.packed.size()*sizeof(gf::Dense16Row),"prepared basis tables");
    const pd::Gl32 gl(n/1024,seed);const Coefficients coeff(gl);const ExpandedCoefficients expanded(gl);
    std::vector<k::block> input(n),expected(n),initial(n),guarded((n/1024)*tileStride+12);
    auto* scratch=reinterpret_cast<k::block*>((reinterpret_cast<std::uintptr_t>(guarded.data()+4)+63)&~std::uintptr_t(63));
    const auto scratchBlocks=(n/1024)*tileStride;
    k::workspace_routing::adviseOwned(scratch,scratchBlocks*16);
    checkInnerSize<Mode>(setup,transformed);
#if defined(__AVX512VBMI__)
    if(!calls && !ip::vbmiStateSelfCheck())throw std::runtime_error("VBMI exhaustive basis");
#endif
    auto selected=[&](double* phases){encodeSize<Mode>(input.data(),scratch,n,route,coeff,expanded,transformed,phases);};
    const k::block canary(0x781aca29312357bfULL,0x928fafcac352697bULL);
    for(unsigned pattern=0;pattern<(calls?1U:4U);++pattern) {
        fill(input,pattern==2?997:913);
        if(pattern==1){std::fill(input.begin(),input.end(),k::block{});input[n-3]=k::block(1,3);}
        if(pattern==3){std::fill(input.begin(),input.end(),k::block{});input[63]=k::block(1,2);input[64]=k::block(4,8);input[n-65]=k::block(16,32);}
        expected=input;initial=input;
        t64probe::encode(expected.data(),scratch,n,route,coeff,setup.packed.data(),nullptr);
        scratch[-1]=canary;scratch[scratchBlocks]=canary;selected(nullptr);
        equal(input.data(),expected.data(),n*16,"full retained output and suffix");
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
             <<" prepared matrices; four-epoch basis; whole output/suffix/canaries"
             <<(calls?"":"; independent scalar full encoder; boundary/API checks")<<'\n';
    if(!calls)return;
    fill(input,913);double phases[2]{};
    for(unsigned i=0;i<5;++i)selected(profile?phases:nullptr);
    std::vector<double> times;times.reserve(calls);
    for(unsigned i=0;i<calls;++i){const auto start=Clock::now();selected(profile?phases:nullptr);
        times.push_back(std::chrono::duration<double,std::milli>(Clock::now()-start).count());}
    std::sort(times.begin(),times.end());
    std::cout<<"mode,exponent,seed,calls,median_ms,p10_ms,p90_ms,checksum\n"
             <<Mode<<','<<exponent<<','<<seed<<','<<calls<<','<<std::fixed<<std::setprecision(6)
             <<times[calls/2]<<','<<times[calls/10]<<','<<times[9*calls/10]<<','<<std::hex<<hash(input)<<std::dec<<'\n';
    if(profile)std::cout<<"phase_means_including_warmups,inner_route,"<<phases[0]/(calls+5)
                      <<",packed_mixer_bch,"<<phases[1]/(calls+5)<<'\n';
}
}
int main(int argc,char** argv) {try {
    if(argc<5 || argc>6)throw std::invalid_argument("usage: inner-size mode exponent seed calls [profile]");
    if(!__builtin_cpu_supports("avx512f") || !__builtin_cpu_supports("avx512vl") ||
       !__builtin_cpu_supports("avx512bw") || !__builtin_cpu_supports("gfni"))
        throw std::runtime_error("AVX512/GFNI required");
#if defined(__AVX512VBMI__)
    if(!__builtin_cpu_supports("avx512vbmi"))throw std::runtime_error("VBMI build requires VBMI CPU");
#endif
    const unsigned mode=std::stoul(argv[1]),exponent=std::stoul(argv[2]),calls=std::stoul(argv[4]);
    const auto seed=std::stoull(argv[3]);const bool profile=argc==6 && std::stoul(argv[5]);
    if(exponent<14 || exponent>20)throw std::invalid_argument("exponent range 14..20");
    if(!calls)sizeprobe::checkInnerBoundaries();
    switch(mode) {
#define SIZE_CASE(M) case M:sizeprobe::experiment<M>(exponent,seed,calls,profile);break
        SIZE_CASE(0);SIZE_CASE(3);SIZE_CASE(4);SIZE_CASE(7);SIZE_CASE(8);SIZE_CASE(11);SIZE_CASE(12);SIZE_CASE(15);
        SIZE_CASE(20);SIZE_CASE(23);SIZE_CASE(28);SIZE_CASE(31);SIZE_CASE(36);SIZE_CASE(39);SIZE_CASE(44);SIZE_CASE(47);
        SIZE_CASE(52);SIZE_CASE(55);SIZE_CASE(60);SIZE_CASE(63);
#if defined(__AVX512VBMI__)
        SIZE_CASE(1);SIZE_CASE(2);SIZE_CASE(5);SIZE_CASE(6);SIZE_CASE(9);SIZE_CASE(10);SIZE_CASE(13);SIZE_CASE(14);
        SIZE_CASE(21);SIZE_CASE(22);SIZE_CASE(29);SIZE_CASE(30);SIZE_CASE(37);SIZE_CASE(38);SIZE_CASE(45);SIZE_CASE(46);
        SIZE_CASE(53);SIZE_CASE(54);SIZE_CASE(61);SIZE_CASE(62);
        SIZE_CASE(64);SIZE_CASE(65);SIZE_CASE(66);SIZE_CASE(67);SIZE_CASE(68);SIZE_CASE(69);
#endif
#undef SIZE_CASE
        default:throw std::invalid_argument("invalid mode/prefetch with NT store");
    }
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
