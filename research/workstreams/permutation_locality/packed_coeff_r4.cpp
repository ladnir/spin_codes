// Research-only coefficient storage probe. Retained sources stay unchanged.
// Reuse original sampling, scalar GL32 oracle, and dense inner references.
#define main packed_coeff_retained_driver_main
#include "packed_driver.cpp"
#undef main
#include "FusedR4.h"

namespace spin::detail::kernel {
void bchPackedCoeffOriginal(const block*,block*,const std::uint64_t*);
void bchPackedCoeffAligned(const block*,block*,const std::uint64_t*);
void bchPackedCoeffCompact(const block*,block*,const std::uint64_t*);
unsigned packedCoeffTileMode();
}
namespace coeff_probe {
namespace fr=spin::research::fused_r4;
struct AlignedWords {
    std::vector<std::uint64_t> storage;
    std::uint64_t* data;
    explicit AlignedWords(std::size_t words):storage(words+8),data(reinterpret_cast<std::uint64_t*>(
        (reinterpret_cast<std::uintptr_t>(storage.data())+63)&~std::uintptr_t(63))) {}
};
struct Coefficients {
    AlignedWords expanded,compact;
    explicit Coefficients(const pd::Gl32& original):expanded(original.coeff.size()),compact(original.coeff.size()/2) {
        std::memcpy(expanded.data,original.coeff.data(),original.coeff.size()*8);
        for(std::size_t i=0;i<original.coeff.size()/2;++i) {
            if(original.coeff[2*i]!=original.coeff[2*i+1])
                throw std::runtime_error("original GFNI coefficient pair differs");
            compact.data[i]=original.coeff[2*i];
        }
    }
};
// Match the retained fused-R4 metadata preparation and allocation order.
struct FusedSetup {
    std::vector<fr::ScalarRow> scalar;
    std::vector<fr::WideRow> wide;
    FusedSetup(const unsigned* masks,std::size_t epochs):scalar(epochs),wide(epochs) {
        for(std::size_t epoch=0;epoch<epochs;++epoch) {
            scalar[epoch]=fr::makeScalar(masks+8*epoch);
            wide[epoch]=fr::makeWide(scalar[epoch]);
        }
    }
};
template<class Emit> SPIN_FORCEINLINE void reverseInner(
    const k::block* in,std::size_t n,const fr::ScalarRow* rows,Emit&& emit) {
    // Exact scalar-fused recurrence from packed_fused_r4.cpp: raw feedback is
    // added after the ORIGINAL four-transvection product, never before it.
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
template<unsigned Mode> SPIN_FORCEINLINE void kernel(const k::block* in,k::block* out,const std::uint64_t* coeff) {
    static_assert(Mode<4);
    if constexpr(Mode==0)k::bchPackedCoeffOriginal(in,out,coeff);
    else if constexpr(Mode==1)k::bchPackedCoeffAligned(in,out,coeff);
    else if constexpr(Mode==2)k::bchPackedCoeffCompact(in,out,coeff);
    else k::bchPackedFullGl32(in,out,coeff);
}
template<unsigned Mode> static SPIN_NOINLINE void encode(
    k::block* input,k::block* scratch,std::size_t n,const Packets<4>& packets,
    const pd::Gl32& original,const Coefficients& coeff,const FusedSetup& setup,double* phases) {
    const auto start=phases?Clock::now():Clock::time_point{};
    pd::Emitter<false> emit{scratch,packets.bases.data(),packets.controls.data()};
    reverseInner(input,n,setup.scalar.data(),emit);_mm_sfence();
    const auto middle=phases?Clock::now():Clock::time_point{};
    const auto* matrices=(Mode==0 || Mode==3)?original.coeff.data():(Mode==1?coeff.expanded.data:coeff.compact.data);
    constexpr std::size_t stride=Mode==2?512:1024;
    for(std::size_t tile=0;tile<n/1024;++tile)
        kernel<Mode>(scratch+tile*tileStride,input+tile*512,matrices+tile*stride);
    if(phases) {
        phases[0]+=std::chrono::duration<double,std::milli>(middle-start).count();
        phases[1]+=std::chrono::duration<double,std::milli>(Clock::now()-middle).count();
    }
}
static void equal(const void* a,const void* b,std::size_t bytes,const char* why) {
    if(std::memcmp(a,b,bytes))throw std::runtime_error(why);
}
static void basisCheck(std::uint64_t seed) {
    const pd::Gl32 original(1,seed);const Coefficients coeff(original);
    alignas(64) k::block input[1024]{},mixed[1024],expected[512],actual[512];
    std::vector<k::block> images(1024*512);
    for(unsigned coordinate=0;coordinate<1024;++coordinate) {
        input[coordinate]=k::block(~0ULL,~0ULL);
        pd::scalarMix(input,mixed,original.coeff.data());
        k::bchTranspose4(mixed,images.data()+512*coordinate);
        input[coordinate]=k::block{};
    }
    for(unsigned coordinate=0;coordinate<1024;++coordinate)for(unsigned bit=0;bit<128;++bit) {
        const auto unit=bit<64?k::block(0,1ULL<<bit):k::block(1ULL<<(bit-64),0);
        input[coordinate]=unit;
        for(unsigned row=0;row<512;++row)
            expected[row]=k::block(_mm_and_si128(images[512*coordinate+row].mData,unit.mData));
        kernel<0>(input,actual,original.coeff.data());equal(actual,expected,sizeof(actual),"original coefficient physical basis mismatch");
        kernel<1>(input,actual,coeff.expanded.data);equal(actual,expected,sizeof(actual),"aligned coefficient physical basis mismatch");
        kernel<2>(input,actual,coeff.compact.data);equal(actual,expected,sizeof(actual),"compact coefficient physical basis mismatch");
        input[coordinate]=k::block{};
    }
    std::cout<<"checks: 131072 physical basis vectors, all three coefficient layouts, original scalar GL32/production BCH PASS; seed="
        <<seed<<"; tile_mode="<<k::packedCoeffTileMode()<<'\n';
    // The API permits 16-byte alignment, not just the 64-byte alignment of
    // the basis arrays above. Guard both ends of a deliberately offset output.
    alignas(64) k::block guarded[514];
    const k::block guard(0x49a5c13bb409528dULL,0x952ad786204fb321ULL);
    k::setup::Words rng(seed^0x55d410cc);
    for(auto& word:input)word=k::block(rng(),rng());
    pd::scalarMix(input,mixed,original.coeff.data());k::bchTranspose4(mixed,expected);
    auto checkGuard=[&] {
        equal(guarded+1,expected,sizeof(expected),"16-byte-offset output mismatch");
        equal(guarded,&guard,sizeof(guard),"output prefix guard changed");
        equal(guarded+513,&guard,sizeof(guard),"output suffix guard changed");
    };
    guarded[0]=guarded[513]=guard;
    kernel<0>(input,guarded+1,original.coeff.data());checkGuard();
    kernel<1>(input,guarded+1,coeff.expanded.data);checkGuard();
    kernel<2>(input,guarded+1,coeff.compact.data);checkGuard();
    kernel<3>(input,guarded+1,original.coeff.data());checkGuard();
    std::cout<<"checks: 16-byte-offset output and prefix/suffix guards PASS\n";
}
template<unsigned Mode> static void experiment(unsigned exponent,std::uint64_t seed,unsigned calls,bool profile) {
    const unsigned n=2U<<exponent;
    const Packets<4> packets(n/256,seed,true,true);
    std::vector<unsigned> masks(std::size_t(8)*(n/128));k::setup::Words rng(seed^0x3f625a92ULL);
    for(std::size_t i=0;i<masks.size();i+=2) {
        unsigned u;do{u=unsigned(rng())&((1U<<19)-1);}while(!u);
        unsigned v=unsigned(rng())&((1U<<19)-1);
        if(std::popcount(u&v)&1)v^=u&-u;masks[i]=u;masks[i+1]=v;
    }
    // Every mode creates every setup form in this exact order, outside timing.
    std::vector<k::block> input(n),inner(n),actual(n),other(n),forward(n),expected(n),materialized(n);
    const auto scratchBlocks=(n/1024)*tileStride;std::vector<k::block> storage(scratchBlocks+4);
    auto* scratch=reinterpret_cast<k::block*>((reinterpret_cast<std::uintptr_t>(storage.data())+63)&~std::uintptr_t(63));
    k::workspace_routing::adviseOwned(scratch,scratchBlocks*16);
    const pd::Gl32 original(n/1024,seed);
    const FusedSetup setup(masks.data(),n/128);
    const Coefficients coeff(original);
    auto selected=[&](double* phases){encode<Mode>(input.data(),scratch,n,packets,original,coeff,setup,phases);};
    fill(other,971);forwardReference<4>(other.data(),forward.data(),n,masks.data());
    for(unsigned pattern=0;pattern<(calls?1U:3U);++pattern) {
        fill(input,pattern==2?997:913);
        if(pattern==1){std::fill(input.begin(),input.end(),k::block{});input[n-3]=k::block(1,3);}
        expected=input;
        reverseReference<4>(input.data(),inner.data(),n,masks.data());
        reverseInner(input.data(),n,setup.scalar.data(),[&](std::size_t i,k::block v){actual[i]=v;});
        equal(actual.data(),inner.data(),n*16,"fused inner / original dense reference mismatch");
        auto lhs=_mm_setzero_si128(),rhs=lhs;
        for(std::size_t i=0;i<n;++i) {
            lhs=_mm_xor_si128(lhs,_mm_and_si128(input[i].mData,forward[i].mData));
            rhs=_mm_xor_si128(rhs,_mm_and_si128(other[i].mData,actual[i].mData));
        }
        if(_mm_movemask_epi8(_mm_cmpeq_epi8(lhs,rhs))!=65535)
            throw std::runtime_error("original forward/fused reverse adjoint mismatch");
        for(std::size_t i=0;i<n;++i) {
            const auto x=packets.inverse[i];
            materialized[(x/1024)*1024+4*(x%256)+(x/256)%4]=inner[i];
        }
        alignas(64) k::block mixed[1024],retained[512],generated[512];
        for(std::size_t tile=0;tile<n/1024;++tile) {
            const auto* source=materialized.data()+1024*tile;
            const auto* matrices=original.coeff.data()+1024*tile;
            pd::scalarMix(source,mixed,matrices);k::bchTranspose4(mixed,expected.data()+512*tile);
            k::bchPackedFullGl32(source,retained,matrices);
            kernel<0>(source,generated,matrices);
            equal(retained,generated,sizeof(retained),"common tile schedule / retained kernel mismatch");
        }
        actual=input;pd::encode<true,4>(actual.data(),scratch,n,masks.data(),packets,original,nullptr);
        equal(actual.data(),expected.data(),n*16,"original sequential R4 full reference mismatch");
        actual=input;encode<0>(actual.data(),scratch,n,packets,original,coeff,setup,nullptr);
        equal(actual.data(),expected.data(),n*16,"original coefficient full/suffix mismatch");
        actual=input;encode<1>(actual.data(),scratch,n,packets,original,coeff,setup,nullptr);
        equal(actual.data(),expected.data(),n*16,"aligned coefficient full/suffix mismatch");
        actual=input;encode<2>(actual.data(),scratch,n,packets,original,coeff,setup,nullptr);
        equal(actual.data(),expected.data(),n*16,"compact coefficient full/suffix mismatch");
        actual=input;encode<3>(actual.data(),scratch,n,packets,original,coeff,setup,nullptr);
        equal(actual.data(),expected.data(),n*16,"retained BCH / fused R4 full/suffix mismatch");
        selected(nullptr);equal(input.data(),expected.data(),n*16,"selected full/suffix mismatch");
    }
    std::cout<<"checks: original dense inner/adjoint, original scalar GL32/production BCH, retained kernel, "
        "all coefficient layouts, full buffer/suffix PASS; updates=4; tile_mode="<<k::packedCoeffTileMode()
        <<"; original_mod64="<<(reinterpret_cast<std::uintptr_t>(original.coeff.data())&63)
        <<"; expanded_bytes="<<original.coeff.size()*8<<"; compact_bytes="<<original.coeff.size()*4<<'\n';
    if(!calls)return;
    fill(input,913);double phases[2]{};auto run=[&]{selected(profile?phases:nullptr);};
    for(unsigned i=0;i<3;++i)run();std::vector<double> times;times.reserve(calls);
    for(unsigned i=0;i<calls;++i) {
        const auto start=Clock::now();run();
        times.push_back(std::chrono::duration<double,std::milli>(Clock::now()-start).count());
    }
    std::sort(times.begin(),times.end());
    std::cout<<"exponent,mode,seed,calls,median_ms,p10_ms,p90_ms,checksum,updates,tile_mode\n"
        <<exponent<<','<<Mode<<','<<seed<<','<<calls<<','<<std::fixed<<std::setprecision(6)
        <<times[calls/2]<<','<<times[calls/10]<<','<<times[9*calls/10]<<','<<std::hex
        <<hash(input)<<std::dec<<",4,"<<k::packedCoeffTileMode()<<'\n';
    if(profile)std::cout<<"phase_means_including_warmups,inner_route,"<<phases[0]/(calls+3)
        <<",packed_mixer_bch,"<<phases[1]/(calls+3)<<'\n';
}
}
int main(int argc,char** argv) {try {
    if(!__builtin_cpu_supports("avx512f") || !__builtin_cpu_supports("avx512vl") ||
       !__builtin_cpu_supports("avx512bw") || !__builtin_cpu_supports("gfni"))
        throw std::runtime_error("AVX512/GFNI required");
    if(argc==3 && std::string(argv[1])=="basis") {
        coeff_probe::basisCheck(std::stoull(argv[2]));return 0;
    }
    if(argc<4 || argc>6)throw std::invalid_argument(
        "usage: packed_coeff_r4 exponent seed calls [mode=0 [profile=0]]; or basis seed; "
        "modes 0=original expanded,1=aligned expanded,2=aligned compact,3=retained BCH; calls=0 checks all layouts");
    const unsigned exponent=std::stoul(argv[1]),calls=std::stoul(argv[3]);
    const auto seed=std::stoull(argv[2]);const unsigned mode=argc>=5?std::stoul(argv[4]):0;
    const bool profile=argc==6 && std::stoul(argv[5]);
    if(exponent<14 || exponent>20 || mode>3)throw std::invalid_argument("unsupported coefficient probe");
    if(mode==0)coeff_probe::experiment<0>(exponent,seed,calls,profile);
    else if(mode==1)coeff_probe::experiment<1>(exponent,seed,calls,profile);
    else if(mode==2)coeff_probe::experiment<2>(exponent,seed,calls,profile);
    else coeff_probe::experiment<3>(exponent,seed,calls,profile);
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
