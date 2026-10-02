// Generated isolated harness: creative mode4/6 uses wide emission.
// Isolated exact-map comparisons. All modes share one setup/allocation path.
// Include the references first: their nested main macros must not erase ours.
#include "PacketReference.h"
#define main creative_unused_size_driver_main
#include "inner_size_probe.cpp"
#undef main
#include "further_mapWide.h"
#include "RoutePacked.h"
#include <limits>

namespace spin::detail::kernel {
void bchGatherCoeffCompact(const block*,block*,const std::uint64_t*,const std::uint32_t*);
}

namespace creativeprobe {
using namespace sizeprobe;
namespace rp = spin::research::route_packed;

// 0: retained size mode 6; 1/2: composed packet/stream;
// 3/4: composed-wide packet/stream; 5: retained stream + routed packing;
// 6: composed-wide stream + routed packing; 7/8: retained stream writing
// contiguous packets, then inverse-route outer loads (cached/NT stores).
template<bool NonTemporal> struct SequentialEmitter {
    k::block* scratch;
    SPIN_FORCEINLINE void operator()(std::size_t i,__m512i value) const {
        if constexpr(NonTemporal)
            _mm512_stream_si512(reinterpret_cast<__m512i*>(scratch+4*i),value);
        else _mm512_store_si512(scratch+4*i,value);
    }
};

template<unsigned Mode, class Emit>
static SPIN_FORCEINLINE void inner(const k::block* input, std::size_t n,
    const ip::PreparedPacketUpdates16& updates,
    const ip::PreparedComposedUpdates16& composed, Emit& emit) {
    static_assert(Mode <= 8);
    if constexpr(Mode == 0 || Mode == 5 || Mode >= 7)
        ip::transposeInner<ip::InnerRecipe::StreamVbmi>(input,n,updates.conjugatedRows(),emit);
    else if constexpr(Mode == 1 || Mode == 2)
        ip::transposeInnerComposed<Mode == 1 ? ip::ComposedRecipe::Packet : ip::ComposedRecipe::Stream>(
            input,n,composed.rows(),emit);
    else
        ip::transposeInnerFurtherWideRecipe<Mode == 3 ? ip::ComposedRecipe::Packet : ip::ComposedRecipe::Stream>(
            input,n,composed.rows(),emit);
}

template<unsigned Mode> static SPIN_NOINLINE void encodeCreative(
    k::block* input,k::block* scratch,std::size_t n,const Packets<4>& route,
    const Coefficients& coeff,const ExpandedCoefficients& expanded,
    const ip::PreparedPacketUpdates16& updates,
    const ip::PreparedComposedUpdates16& composed,const std::uint32_t* byteBases,
    const std::uint32_t* gatherOffsets,
    double* phases) {
    if constexpr(Mode == 0) {
        // Call the retained implementation itself, not a rewritten control.
        encodeSize<6>(input,scratch,n,route,coeff,expanded,updates,phases);
    } else {
        const auto start=phases?Clock::now():Clock::time_point{};
        if constexpr(Mode >= 7) {
            SequentialEmitter<Mode == 8> emit{scratch};
            inner<Mode>(input,n,updates,composed,emit);
            if constexpr(Mode == 8)_mm_sfence();
        } else if constexpr(Mode >= 5) {
            rp::Emitter emit{reinterpret_cast<unsigned char*>(scratch),byteBases};
            inner<Mode>(input,n,updates,composed,emit);
        } else {
            sizeprobe::Route<true> emit{scratch,route.bases.data()};
            inner<Mode>(input,n,updates,composed,emit);
        }
        const auto middle=phases?Clock::now():Clock::time_point{};
        for(std::size_t tile=0;tile<n/1024;++tile) {
            if constexpr(Mode >= 7)
                k::bchGatherCoeffCompact(scratch,input+tile*512,coeff.compact+512*tile,gatherOffsets+256*tile);
            else if constexpr(Mode >= 5)
                k::bchRoutePackedCoeffCompact(scratch+tile*tileStride,input+tile*512,coeff.compact+512*tile);
            else
                k::bchPackedCoeffCompact(scratch+tile*tileStride,input+tile*512,coeff.compact+512*tile);
        }
        if(phases) {
            phases[0]+=std::chrono::duration<double,std::milli>(middle-start).count();
            phases[1]+=std::chrono::duration<double,std::milli>(Clock::now()-middle).count();
        }
    }
}

template<unsigned Mode> static void checkShortInner(const Updates& setup) {
    // Covers the one-epoch case, first/last physical steps, and non-power-of-two
    // epoch counts. The oracle reads original-coordinate matrices, not our tables.
    for(unsigned n : {64U,128U,192U,256U}) {
        ip::PreparedPacketUpdates16 updates(std::span<const Matrix>(setup.reverse.data(),n/64));
        ip::PreparedComposedUpdates16 composed(updates.conjugatedRows());
        std::vector<k::block> input(n),expected(n),actual(n);
        std::vector<unsigned> writes(n/4);
        for(unsigned coordinate=0;coordinate<n;++coordinate) {
            std::fill(input.begin(),input.end(),k::block{});
            std::fill(actual.begin(),actual.end(),k::block(1,1));
            std::fill(writes.begin(),writes.end(),0U);
            input[coordinate]=k::block(0x729aca378bc42fedULL,0x942163aca892317bULL);
            innerOracle<true>(input.data(),expected.data(),n,setup);
            std::size_t next=n/4;
            auto emit=[&](std::size_t i,__m512i value) {
                if(i>=n/4 || next==0 || i!=--next || writes[i]++)
                    throw std::runtime_error("short inner packet index/order/duplication");
                _mm512_storeu_si512(actual.data()+4*i,value);
            };
            inner<Mode>(input.data(),n,updates,composed,emit);
            if(next)throw std::runtime_error("short inner missing packets");
            equal(actual.data(),expected.data(),n*16,"short original-coordinate inner basis");
        }
    }
}

static void scalarBchTranspose(const k::block* mixed,k::block* output) {
    // Direct, literal BCH matrix action: no optimized BCH/GFNI routine.
    for(unsigned lane=0;lane<4;++lane)
        for(unsigned row=0;row<128;++row) {
            k::block value{};
            for(unsigned column=0;column<256;++column)
                if((k::BchRows[row][column/64]>>(column%64))&1U)
                    value^=mixed[4*column+lane];
            output[128*lane+row]=value;
        }
}

template<bool Sequential>
static void setScratchCanaries(k::block* scratch,std::size_t tiles,const k::block& canary) {
    for(unsigned i=0;i<4;++i) {
        scratch[-4+int(i)]=canary;
        scratch[tiles*tileStride+i]=canary;
    }
    if constexpr(Sequential) {
        for(unsigned i=0;i<4;++i)scratch[tiles*1024+i]=canary;
    } else {
        for(std::size_t tile=0;tile<tiles;++tile)
            for(unsigned i=0;i<4;++i)scratch[tile*tileStride+1024+i]=canary;
    }
}

template<bool Sequential>
static void checkScratchCanaries(const k::block* scratch,std::size_t tiles,const k::block& canary) {
    for(unsigned i=0;i<4;++i) {
        equal(scratch-4+i,&canary,16,"leading scratch canary");
        equal(scratch+tiles*tileStride+i,&canary,16,"trailing scratch canary");
    }
    if constexpr(Sequential) {
        for(unsigned i=0;i<4;++i)
            equal(scratch+tiles*1024+i,&canary,16,"contiguous scratch end");
    } else {
        for(std::size_t tile=0;tile<tiles;++tile)
            for(unsigned i=0;i<4;++i)
                equal(scratch+tile*tileStride+1024+i,&canary,16,"scratch tile padding");
    }
}

template<unsigned Mode> static void experiment(unsigned exponent,std::uint64_t seed,unsigned calls,bool profile) {
    const unsigned n=2U<<exponent;
    // This allocation/setup sequence is identical for EVERY mode, including
    // control: neither composed rows nor the converted route are timed.
    const Packets<4> route(n/256,seed,true,true);
    const Updates setup(n/64,seed);const Setup retained(setup);
    const ip::PreparedPacketUpdates16 transformed(setup.reverse);
    const ip::PreparedComposedUpdates16 composed(transformed.conjugatedRows());
    std::vector<std::uint32_t> byteBases(route.bases.size());
    std::vector<std::uint32_t> gatherOffsets(route.bases.size(),std::numeric_limits<std::uint32_t>::max());
    for(std::size_t i=0;i<byteBases.size();++i) {
        const auto base=route.bases[i];
        byteBases[i]=rp::packedBaseBytes(base);
        const auto destination=(base/tileStride)*256+(base%tileStride)/4;
        if(destination>=gatherOffsets.size() || gatherOffsets[destination]!=std::numeric_limits<std::uint32_t>::max())
            throw std::runtime_error("inverse gather route is not a permutation");
        gatherOffsets[destination]=static_cast<std::uint32_t>(4*i);
    }
    equal(transformed.conjugatedRows().data(),retained.packed.data(),
          retained.packed.size()*sizeof(gf::Dense16Row),"prepared basis tables");
    const pd::Gl32 gl(n/1024,seed);const Coefficients coeff(gl);const ExpandedCoefficients expanded(gl);
    std::vector<k::block> input(n),expected(n),initial(n),guarded((n/1024)*tileStride+12);
    auto* scratch=reinterpret_cast<k::block*>((reinterpret_cast<std::uintptr_t>(guarded.data()+4)+63)&~std::uintptr_t(63));
    const auto scratchBlocks=(n/1024)*tileStride;
    k::workspace_routing::adviseOwned(scratch,scratchBlocks*16);
    checkShortInner<Mode>(setup);
    auto selected=[&](double* phases) {
        encodeCreative<Mode>(input.data(),scratch,n,route,coeff,expanded,transformed,composed,
                             byteBases.data(),gatherOffsets.data(),phases);
    };
    const k::block canary(0x781aca29312357bfULL,0x928fafcac352697bULL);
    for(unsigned pattern=0;pattern<(calls?1U:4U);++pattern) {
        fill(input,pattern==2?997:913);
        if(pattern==1){std::fill(input.begin(),input.end(),k::block{});input[n-3]=k::block(1,3);}
        if(pattern==3){std::fill(input.begin(),input.end(),k::block{});input[63]=k::block(1,2);input[64]=k::block(4,8);input[n-65]=k::block(16,32);}
        expected=input;initial=input;
        t64probe::encode(expected.data(),scratch,n,route,coeff,setup.packed.data(),nullptr);
        setScratchCanaries<(Mode >= 7)>(scratch,n/1024,canary);
        selected(nullptr);
        equal(input.data(),expected.data(),n*16,"full retained output and suffix");
        equal(input.data()+n/2,initial.data()+n/2,n*8,"untouched input suffix");
        checkScratchCanaries<(Mode >= 7)>(scratch,n/1024,canary);
        if(!calls) {
            innerOracle<true>(initial.data(),expected.data(),n,setup);
            for(std::size_t i=0;i<n;++i)initial[canonical(route.inverse[i])]=expected[i];
            alignas(64) k::block mixed[1024];
            for(std::size_t tile=0;tile<n/1024;++tile) {
                pd::scalarMix(initial.data()+1024*tile,mixed,gl.coeff.data()+1024*tile);
                scalarBchTranspose(mixed,expected.data()+512*tile);
            }
            equal(input.data(),expected.data(),n*8,"independent scalar full encoder");
        }
    }
    std::cout<<"checks PASS mode="<<Mode<<" seed="<<seed<<" exponent="<<exponent
             <<" prepared matrices; 1/2/3/4-epoch basis; whole output/suffix/canaries/padding"
             <<(calls?"":"; scalar full encoder and native helper bases")<<'\n';
    if(!calls)return;
    fill(input,913);double phases[2]{};
    for(unsigned i=0;i<5;++i)selected(profile?phases:nullptr);
    std::vector<double> times;times.reserve(calls);
    for(unsigned i=0;i<calls;++i) {
        const auto start=Clock::now();selected(profile?phases:nullptr);
        times.push_back(std::chrono::duration<double,std::milli>(Clock::now()-start).count());
    }
    std::sort(times.begin(),times.end());
    std::cout<<"mode,exponent,seed,calls,median_ms,p10_ms,p90_ms,checksum\n"
             <<Mode<<','<<exponent<<','<<seed<<','<<calls<<','<<std::fixed<<std::setprecision(6)
             <<times[calls/2]<<','<<times[calls/10]<<','<<times[9*calls/10]<<','<<std::hex<<hash(input)<<std::dec<<'\n';
    if(profile)std::cout<<"phase_means_including_warmups,inner_route,"<<phases[0]/(calls+5)
                      <<",packed_mixer_bch,"<<phases[1]/(calls+5)<<'\n';
}

static std::uint64_t unsignedArgument(const char* text) {
    if(!*text || *text=='-' || *text=='+')throw std::invalid_argument("expected an unsigned decimal integer");
    for(const char* p=text;*p;++p)if(*p<'0' || *p>'9')throw std::invalid_argument("invalid decimal argument");
    return std::stoull(text);
}
} // namespace creativeprobe

int main(int argc,char** argv) {try {
    if(argc<5 || argc>6)throw std::invalid_argument("usage: creative mode exponent seed calls [profile]");
    if(!__builtin_cpu_supports("avx512f") || !__builtin_cpu_supports("avx512vl") ||
       !__builtin_cpu_supports("avx512bw") || !__builtin_cpu_supports("avx512vbmi") || !__builtin_cpu_supports("gfni"))
        throw std::runtime_error("AVX512 F/VL/BW/VBMI and GFNI required");
    const auto mode=creativeprobe::unsignedArgument(argv[1]);
    const auto exponent=creativeprobe::unsignedArgument(argv[2]);
    const auto seed=creativeprobe::unsignedArgument(argv[3]);
    const auto calls=creativeprobe::unsignedArgument(argv[4]);
    const auto profile=argc==6?creativeprobe::unsignedArgument(argv[5]):0;
    if(mode>8 || exponent<14 || exponent>20 || calls>std::numeric_limits<unsigned>::max()-5U || profile>1)
        throw std::invalid_argument("mode 0..8, exponent 14..20, calls <= UINT32_MAX-5, profile 0/1 required");
    if(!calls) {
        sizeprobe::checkInnerBoundaries();
        if(!ip::vbmiStateSelfCheck())throw std::runtime_error("VBMI state basis");
        if(!ip::composedStateSelfCheck() || !ip::furtherMapWideSelfCheck())throw std::runtime_error("composed state/feedback basis");
        if(!ip::packetMomentsSelfCheck())throw std::runtime_error("packet feedback basis");
        if(!creativeprobe::rp::routePackedSelfCheck())throw std::runtime_error("routed packing basis/canaries");
        packetprobe::checkMaps<true,false>(true);
    }
    switch(mode) {
#define CREATIVE_CASE(M) case M:creativeprobe::experiment<M>(unsigned(exponent),seed,unsigned(calls),profile!=0);break
        CREATIVE_CASE(0);CREATIVE_CASE(1);CREATIVE_CASE(2);CREATIVE_CASE(3);
        CREATIVE_CASE(4);CREATIVE_CASE(5);CREATIVE_CASE(6);CREATIVE_CASE(7);CREATIVE_CASE(8);
#undef CREATIVE_CASE
    }
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
