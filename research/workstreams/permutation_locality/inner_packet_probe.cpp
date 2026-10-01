// Isolated exact-map inner experiment. T64Reference.h is mechanically extracted
// by inner_packet_codegen.py; certified sources and production code are untouched.
#include "T64Reference.h"
#include "InnerPacketMaps.h"
#include "InnerPacketFeedback.h"
#include "InnerWordUpdate.h"
namespace ip = spin::research::packet_inner;
namespace packetprobe {
using namespace t64probe;

static Matrix multiply(const Matrix& a,const Matrix& b) {
    Matrix result{};
    for(unsigned i=0;i<16;++i)for(unsigned j=0;j<16;++j)
        if((a[i]>>j)&1)result[i]^=b[j];
    return result;
}
struct Setup {
    std::vector<gf::Dense16Row> packed;
    std::vector<std::uint16_t> wordOriginal,wordBasis;
    explicit Setup(const Updates& original):packed(original.packed.size()),
        wordOriginal(16*packed.size()),wordBasis(16*packed.size()) {
        Matrix p{},inverse{};
        for(unsigned j=0;j<16;++j) {
            p[j]=ip::Maps<true>::basisRows[j];inverse[j]=ip::Maps<true>::inverseRows[j];
        }
        const auto id=multiply(p,inverse);
        for(unsigned j=0;j<16;++j)if(id[j]!=(1U<<j))throw std::runtime_error("basis inverse");
        for(std::size_t e=0;e<packed.size();++e) {
            const auto matrix=multiply(multiply(p,original.reverse[e]),inverse);
            packed[e]=gf::makeDense16(matrix.data());
            for(unsigned j=0;j<16;++j) {
                wordOriginal[16*e+j]=original.reverse[e][j];wordBasis[16*e+j]=matrix[j];
            }
        }
    }
};
template<bool Basis> static void basisState(const __m128i* original,__m128i* result) {
    for(unsigned j=0;j<16;++j) {
        result[j]=_mm_setzero_si128();
        const auto mask=ip::Maps<Basis>::basisRows[j];
        for(unsigned i=0;i<16;++i)if((mask>>i)&1)result[j]=_mm_xor_si128(result[j],original[i]);
    }
}
template<std::size_t... I> SPIN_FORCEINLINE void loadPackets(
    const k::block* input,__m512i* out,std::index_sequence<I...>) {
    ((out[I]=_mm512_loadu_si512(input+4*I)),...);
}
template<bool KeepRaw,unsigned H,class Emit> SPIN_FORCEINLINE void outputPacket(
    const k::block* input,__m512i* expansion,__m128i* rawWords,std::size_t base,Emit& emit) {
    const auto raw=_mm512_loadu_si512(input+4*H);
    if constexpr(KeepRaw)_mm512_store_si512(rawWords+4*H,raw);
    expansion[H]=_mm512_xor_si512(expansion[H],raw);
    emit(base/4+H,expansion[H]);
}
template<bool KeepRaw,class Emit,std::size_t... I> SPIN_FORCEINLINE void outputPackets(
    const k::block* input,__m512i* expansion,__m128i* rawWords,std::size_t base,Emit& emit,std::index_sequence<I...>) {
    // Reverse issue order matches the retained routing stream. Each store is
    // already a full cache line; there is no sub-cache-line scatter to combine.
    (outputPacket<KeepRaw,15-I>(input,expansion,rawWords,base,emit),...);
}
struct Route {
    k::block* scratch;const unsigned* bases;
    SPIN_FORCEINLINE void operator()(std::size_t i,__m512i x) {
        _mm512_stream_si512(reinterpret_cast<__m512i*>(scratch+bases[i]),x);
    }
};
// Feedback=0 retains the scalar-block zeta, 1 uses packets from raw input,
// 2 uses emitted packets (valid because A^T A=0). Pruned changes only expansion.
template<bool Basis,unsigned Feedback,bool Pruned,bool Word=false,bool Group=false,class Emit>
SPIN_FORCEINLINE void reverse(const k::block* input,std::size_t n,
    const gf::Dense16Row* rows,Emit&& emit,const std::uint16_t* wordRows=nullptr) {
    alignas(64) __m128i words[16],moments[64],syndrome[16];
    alignas(64) __m512i packets[16];
    gf::Packed<16> state{},feedback;
    for(std::size_t epoch=n/64;epoch-->0;) {
        const auto* raw=input+64*epoch;
        if(epoch+1==n/64) {
            loadPackets(raw,packets,std::make_index_sequence<16>{});
            for(unsigned h=16;h-->0;) {
                if constexpr(Feedback==0)_mm512_store_si512(moments+4*h,packets[h]);
                emit(16*epoch+h,packets[h]);
            }
        } else {
            if constexpr(!Word)gf::unpack(state,words);
            ip::Maps<Basis>::template emission<Pruned>(words,packets);
            outputPackets<Feedback==0>(raw,packets,moments,64*epoch,emit,std::make_index_sequence<16>{});
        }
        if(!epoch)break;
        if constexpr(Feedback==0) {
            k::zeta<64>(moments);
        } else {
            if constexpr(Feedback==1)loadPackets(raw,packets,std::make_index_sequence<16>{});
            ip::packetMoments(packets,moments);
        }
        if constexpr(Group) {
            static_assert(!Word);
            if(epoch+1==n/64)ip::Maps<Basis>::finishPacked(moments,state);
            else {ip::Maps<Basis>::finishPacked(moments,feedback);gf::denseStep16<true>(state,rows[epoch],&feedback);}
        } else if constexpr(Word) {
            ip::Maps<Basis>::finish(moments,syndrome);
            if(epoch+1==n/64)std::memcpy(words,syndrome,sizeof(words));
            else ip::wordUpdate(words,wordRows+16*epoch,syndrome);
        } else {
            ip::Maps<Basis>::finish(moments,syndrome);
            if(epoch+1==n/64)gf::pack<16>(syndrome,state);
            else {gf::pack<16>(syndrome,feedback);gf::denseStep16<true>(state,rows[epoch],&feedback);}
        }
    }
}
template<bool Basis,unsigned Feedback,bool Pruned,bool Word=false,bool Group=false>
static SPIN_NOINLINE void encode(k::block* input,k::block* scratch,std::size_t n,
    const Packets<4>& route,const Coefficients& coeff,const gf::Dense16Row* rows,double* phases,
    const std::uint16_t* wordRows=nullptr) {
    const auto start=phases?Clock::now():Clock::time_point{};
    Route emit{scratch,route.bases.data()};
    reverse<Basis,Feedback,Pruned,Word,Group>(input,n,rows,emit,wordRows);_mm_sfence();
    const auto middle=phases?Clock::now():Clock::time_point{};
    for(std::size_t tile=0;tile<n/1024;++tile)
        k::bchPackedCoeffCompact(scratch+tile*tileStride,input+tile*512,coeff.compact+512*tile);
    if(phases) {
        phases[0]+=std::chrono::duration<double,std::milli>(middle-start).count();
        phases[1]+=std::chrono::duration<double,std::milli>(Clock::now()-middle).count();
    }
}
template<bool Basis,bool Pruned> static void checkMaps(bool fullBits) {
    for(unsigned bit=0;bit<16;++bit)for(unsigned payload=0;payload<(fullBits?128U:1U);++payload) {
        alignas(64) __m128i original[16]{},state[16];alignas(64) __m512i out[16];
        original[bit]=fullBits?_mm_set_epi64x(payload>=64?(1ULL<<(payload&63)):0,payload<64?(1ULL<<payload):0)
            :k::block(0x123456789abcdef0ULL,0xfedcba9876543210ULL).mData;
        basisState<Basis>(original,state);ip::Maps<Basis>::template emission<Pruned>(state,out);
        alignas(64) k::block words[64];std::memcpy(words,out,sizeof(out));
        for(unsigned p=0;p<64;++p) {
            const auto expected=((Fixed::expansionRows[bit]>>p)&1)?original[bit]:_mm_setzero_si128();
            equal(words+p,&expected,16,"packet expansion physical basis");
        }
    }
    for(unsigned p=0;p<64;++p)for(unsigned payload=0;payload<(fullBits?128U:1U);++payload) {
        alignas(64) k::block raw[64]{};alignas(64) __m512i packets[16];
        alignas(64) __m128i gotMoments[64]{},reference[64],got[16],orig[16],want[16];
        raw[p]=fullBits?k::block(_mm_set_epi64x(payload>=64?(1ULL<<(payload&63)):0,payload<64?(1ULL<<payload):0))
            :k::block(0xa15c89e630725bf4ULL,0x5afe897654321fedULL);
        loadPackets(raw,packets,std::make_index_sequence<16>{});
        ip::packetMoments(packets,gotMoments);std::memcpy(reference,raw,sizeof(raw));k::zeta<64>(reference);
        for(unsigned m=0;m<64;++m)if(__builtin_popcount(m)<=2)
            equal(gotMoments+m,reference+m,16,"packet moment physical basis");
        Fixed::finish(reference,orig);basisState<Basis>(orig,want);
        ip::Maps<Basis>::finish(gotMoments,got);equal(got,want,sizeof(got),"packet feedback basis");
        gf::Packed<16> packed;ip::Maps<Basis>::finishPacked(gotMoments,packed);
        gf::unpack(packed,got);equal(got,want,sizeof(got),"grouped packed feedback basis");
    }
}
template<bool Basis,unsigned Feedback,bool Pruned,bool Word=false,bool Group=false> static void checkInner(const Updates& setup,const Setup& transformed) {
    // Every input coordinate at four epochs exercises initialization, TWO
    // conjugated updates, final emission, and the omitted final feedback.
    constexpr unsigned n=256;std::vector<k::block> input(n),expected(n),actual(n);
    for(unsigned bit=0;bit<n;++bit) {
        std::fill(input.begin(),input.end(),k::block{});input[bit]=k::block(0x729aca378bc42fedULL,0x942163aca892317bULL);
        innerOracle<true>(input.data(),expected.data(),n,setup);
        auto emit=[&](std::size_t i,__m512i x){_mm512_storeu_si512(actual.data()+4*i,x);};
        reverse<Basis,Feedback,Pruned,Word,Group>(input.data(),n,Basis?transformed.packed.data():setup.packed.data(),emit,
            Basis?transformed.wordBasis.data():transformed.wordOriginal.data());
        equal(actual.data(),expected.data(),n*16,"four-epoch exhaustive input basis");
    }
}
template<int Mode> static void experiment(unsigned exponent,std::uint64_t seed,unsigned calls,bool profile) {
    constexpr bool basis=Mode==2 || Mode==4 || Mode==5 || Mode==7 || Mode==8 || Mode==9 || Mode==13 || Mode==15 || Mode>=16;
    constexpr unsigned feedback=(Mode==3 || Mode==4 || Mode==8)?1:((Mode==5 || Mode>=9)?2:0);
    constexpr bool pruned=(Mode>=6 && Mode<=9) || Mode==11 || Mode==14 || Mode==15 || Mode==17;
    constexpr bool word=Mode>=12 && Mode<=15;
    constexpr bool group=Mode>=16;
    const unsigned n=2U<<exponent;const Packets<4> route(n/256,seed,true,true);
    if(!calls && !ip::packetMomentsSelfCheck())throw std::runtime_error("packet feedback full bit basis");
    if(!calls && word && !ip::wordUpdateSelfCheck())throw std::runtime_error("word update full basis/masks");
    const Updates setup(n/64,seed);const Setup transformed(setup);
    const pd::Gl32 gl(n/1024,seed);const Coefficients coeff(gl);
    std::vector<k::block> input(n),expected(n),initial(n),guarded((n/1024)*tileStride+12);
    auto* scratch=reinterpret_cast<k::block*>((reinterpret_cast<std::uintptr_t>(guarded.data()+4)+63)&~std::uintptr_t(63));
    const auto scratchBlocks=(n/1024)*tileStride;
    k::workspace_routing::adviseOwned(scratch,scratchBlocks*16);
    if constexpr(Mode!=0) {checkMaps<basis,pruned>(!calls);checkInner<basis,feedback,pruned,word,group>(setup,transformed);}
    auto selected=[&](double* phases) {
        if constexpr(Mode==0)t64probe::encode(input.data(),scratch,n,route,coeff,setup.packed.data(),phases);
        else encode<basis,feedback,pruned,word,group>(input.data(),scratch,n,route,coeff,
            basis?transformed.packed.data():setup.packed.data(),phases,
            basis?transformed.wordBasis.data():transformed.wordOriginal.data());
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
            // Independent scalar inner oracle plus explicit permutation, scalar
            // GL32 and retained BCH: no candidate basis or GFNI state tables.
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
             <<(Mode?" maps, four-epoch basis,":"")<<" complete retained output/suffix, scratch canaries"
             <<(calls?"":"; independent scalar complete encoder")<<'\n';
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
    if(argc<5 || argc>6)throw std::invalid_argument("usage: inner-packet mode exponent seed calls [profile]");
    if(!__builtin_cpu_supports("avx512f") || !__builtin_cpu_supports("avx512vl") ||
       !__builtin_cpu_supports("avx512bw") || !__builtin_cpu_supports("gfni"))throw std::runtime_error("AVX512/GFNI required");
    const unsigned mode=std::stoul(argv[1]),exponent=std::stoul(argv[2]),calls=std::stoul(argv[4]);
    const auto seed=std::stoull(argv[3]);
    const bool profile=argc==6 && std::stoul(argv[5]);
    if(exponent<14 || exponent>20)throw std::invalid_argument("exponent range 14..20");
    switch(mode) {
#define PACKET_CASE(M) case M:packetprobe::experiment<M>(exponent,seed,calls,profile);break
        PACKET_CASE(0);PACKET_CASE(1);PACKET_CASE(2);PACKET_CASE(3);PACKET_CASE(4);
        PACKET_CASE(5);PACKET_CASE(6);PACKET_CASE(7);PACKET_CASE(8);PACKET_CASE(9);PACKET_CASE(10);PACKET_CASE(11);
        PACKET_CASE(12);PACKET_CASE(13);PACKET_CASE(14);PACKET_CASE(15);
        PACKET_CASE(16);PACKET_CASE(17);
#undef PACKET_CASE
        default:throw std::invalid_argument("mode range 0..17");
    }
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
