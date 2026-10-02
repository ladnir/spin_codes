// Isolated exact-map GFNI R4 probe; no retained source or production edits.
#define main gfni_r4_retained_driver_main
#include "packed_driver.cpp"
#undef main
#include "FusedR4.h"
#include "FusedR4Gfni.h"

namespace spin::detail::kernel {
void bchPackedCoeffCompact(const block*,block*,const std::uint64_t*);
unsigned packedCoeffTileMode();
}
namespace gp {
namespace fr=spin::research::fused_r4;
namespace gf=spin::research::gfni_r4;
struct AlignedWords {
    std::vector<std::uint64_t> storage;std::uint64_t* data;
    explicit AlignedWords(std::size_t words):storage(words+8),data(reinterpret_cast<std::uint64_t*>(
        (reinterpret_cast<std::uintptr_t>(storage.data())+63)&~std::uintptr_t(63))) {}
};
struct Coefficients {
    AlignedWords expanded,compact;
    explicit Coefficients(const pd::Gl32& original):expanded(original.coeff.size()),compact(original.coeff.size()/2) {
        std::memcpy(expanded.data,original.coeff.data(),original.coeff.size()*8);
        for(std::size_t i=0;i<original.coeff.size()/2;++i) {
            if(original.coeff[2*i]!=original.coeff[2*i+1])throw std::runtime_error("coefficient pair differs");
            compact.data[i]=original.coeff[2*i];
        }
    }
};
struct ScalarSetup {
    std::vector<fr::ScalarRow> scalar;std::vector<fr::WideRow> wide;
    ScalarSetup(const unsigned* masks,std::size_t epochs):scalar(epochs),wide(epochs) {
        for(std::size_t i=0;i<epochs;++i){scalar[i]=fr::makeScalar(masks+8*i);wide[i]=fr::makeWide(scalar[i]);}
    }
};
static void equal(const void* a,const void* b,std::size_t size,const char* why) {
    if(std::memcmp(a,b,size))throw std::runtime_error(why);
}
template<unsigned S> static void sample(unsigned* masks,k::setup::Words& rng) {
    for(unsigned i=0;i<4;++i) {
        unsigned u;do{u=unsigned(rng())&((1U<<S)-1);}while(!u);
        unsigned v=unsigned(rng())&((1U<<S)-1);
        if(std::popcount(u&v)&1)v^=u&-u;masks[2*i]=u;masks[2*i+1]=v;
    }
}
template<unsigned S,unsigned R=4> static void scalarOracle(__m128i* state,const unsigned* masks,bool transpose) {
    for(unsigned r=0;r<R;++r) {
        const unsigned i=transpose?R-1-r:r;
        const auto dot=masks[2*i+(transpose?0:1)],out=masks[2*i+(transpose?1:0)];
        auto value=_mm_setzero_si128();
        for(unsigned j=0;j<S;++j)if((dot>>j)&1)value=_mm_xor_si128(value,state[j]);
        for(unsigned j=0;j<S;++j)if((out>>j)&1)state[j]=_mm_xor_si128(state[j],value);
    }
}
template<unsigned S> static void checkPadding(const gf::Packed<S>& state) {
    if constexpr(S==19) {
        const auto mask=_mm512_set1_epi8(char(0xf8));
        if(_mm512_test_epi64_mask(state.v[4],mask) || _mm512_test_epi64_mask(state.v[5],mask))
            throw std::runtime_error("GFNI nonzero coordinate padding");
    }
}
template<unsigned S> static void checkCore(std::uint64_t seed) {
    k::setup::Words rng(seed);unsigned masks[8];sample<S>(masks,rng);
    alignas(64) __m128i input[S]{},syndrome[S],oracle[S],actual[S+2];
    const auto guard=k::block(0x530d7740ca1a5491ULL,0x73a236c2aba07cddULL).mData;
    actual[0]=actual[S+1]=guard;
    gf::Packed<S> packed{},feedback{};
    for(unsigned coordinate=0;coordinate<S;++coordinate)for(unsigned bit=0;bit<128;++bit) {
        input[coordinate]=(bit<64?k::block(0,1ULL<<bit):k::block(1ULL<<(bit-64),0)).mData;
        gf::pack<S>(input,packed);checkPadding(packed);gf::unpack(packed,actual+1);
        equal(input,actual+1,sizeof(input),"GFNI pack/unpack physical basis mismatch");
        for(unsigned transpose=0;transpose<2;++transpose) {
            std::memcpy(oracle,input,sizeof(input));scalarOracle<S>(oracle,masks,transpose);
            gf::pack<S>(input,packed);gf::step(packed,gf::makeRow<S>(masks,transpose));
            checkPadding(packed);gf::unpack(packed,actual+1);
            equal(oracle,actual+1,sizeof(input),"GFNI physical basis map mismatch");
        }
        input[coordinate]=_mm_setzero_si128();
    }
    // Fresh masks, both orientations, arbitrary payloads, post-update feedback,
    // and input states covering the complete coordinate basis simultaneously.
    for(unsigned trial=0;trial<256;++trial) {
        sample<S>(masks,rng);
        for(unsigned j=0;j<S;++j) {
            input[j]=k::block(rng(),1ULL<<j).mData;
            syndrome[j]=k::block(rng(),rng()).mData;
        }
        gf::pack<S>(syndrome,feedback);
        for(unsigned transpose=0;transpose<2;++transpose) {
            std::memcpy(oracle,input,sizeof(input));scalarOracle<S>(oracle,masks,transpose);
            for(unsigned j=0;j<S;++j)oracle[j]=_mm_xor_si128(oracle[j],syndrome[j]);
            gf::pack<S>(input,packed);gf::step<S,true>(packed,gf::makeRow<S>(masks,transpose),&feedback);
            checkPadding(packed);gf::unpack(packed,actual+1);
            equal(oracle,actual+1,sizeof(input),"GFNI original-mask/feedback oracle mismatch");
        }
    }
    equal(actual,&guard,sizeof(guard),"GFNI unpack prefix guard changed");
    equal(actual+S+1,&guard,sizeof(guard),"GFNI unpack suffix guard changed");
    std::cout<<"core PASS: S="<<S<<",physical_basis="<<128*S
        <<",forward+transpose,256 fresh masks,post-update syndrome,zero padding,unpack guards\n";
}
template<unsigned R> static void checkDense16(std::uint64_t seed) {
    static_assert(R==4 || R==8);
    k::setup::Words rng(seed);unsigned masks[2*R];
    for(unsigned trial=0;trial<64;++trial) {
        for(unsigned r=0;r<R;r+=4)sample<16>(masks+2*r,rng);
        unsigned forwardRows[16]{},transposeRows[16]{};
        for(unsigned j=0;j<16;++j) {
            unsigned column=1U<<j;
            for(unsigned r=0;r<R;++r)
                if(std::popcount(column&masks[2*r+1])&1)column^=masks[2*r];
            transposeRows[j]=column;
            for(unsigned i=0;i<16;++i)forwardRows[i]|=((column>>i)&1)<<j;
        }
        const auto forward=gf::makeDense16(forwardRows),transpose=gf::makeDense16(transposeRows);
        alignas(32) __m128i input[16],syndrome[16],expected[16],actual[16];
        for(unsigned j=0;j<16;++j)syndrome[j]=k::block(rng(),rng()).mData;
        gf::Packed<16> state,feedback;gf::pack<16>(syndrome,feedback);
        const unsigned patterns=trial?1:2048;
        for(unsigned pattern=0;pattern<patterns;++pattern) {
            for(unsigned j=0;j<16;++j)input[j]=trial?k::block(rng(),1ULL<<j).mData:_mm_setzero_si128();
            if(!trial) {
                const auto bit=pattern%128;
                input[pattern/128]=(bit<64?k::block(0,1ULL<<bit):k::block(1ULL<<(bit-64),0)).mData;
            }
            for(unsigned t=0;t<2;++t) {
                std::memcpy(expected,input,sizeof(input));scalarOracle<16,R>(expected,masks,t);
                gf::pack<16>(input,state);gf::denseStep16(state,t?transpose:forward);gf::unpack(state,actual);
                equal(actual,expected,sizeof(actual),"dense16/original transvection product mismatch");
                for(unsigned j=0;j<16;++j)expected[j]=_mm_xor_si128(expected[j],syndrome[j]);
                gf::pack<16>(input,state);gf::denseStep16<true>(state,t?transpose:forward,&feedback);gf::unpack(state,actual);
                equal(actual,expected,sizeof(actual),"dense16 post-update feedback mismatch");
            }
        }
    }
    std::cout<<"dense16 core PASS: R="<<R<<",2048 physical basis,64 fresh products,forward+transpose,post-update feedback\n";
}
static void checkEpochs(const unsigned* masks,const std::vector<gf::Row<19>>& rows) {
    for(std::size_t i=0;i<rows.size();++i) {
        alignas(32) __m128i source[19],expected[19],actual[19];
        for(unsigned j=0;j<19;++j)source[j]=k::block(1ULL<<j,1ULL<<j).mData;
        std::memcpy(expected,source,sizeof(source));scalarOracle<19>(expected,masks+8*i,true);
        gf::Packed<19> state;gf::pack<19>(source,state);gf::step(state,rows[i]);
        checkPadding(state);gf::unpack(state,actual);
        equal(expected,actual,sizeof(expected),"per-epoch GFNI original-mask basis mismatch");
    }
}
template<unsigned Mode,class Emit> SPIN_FORCEINLINE void reverse(
    const k::block* input,std::size_t n,const ScalarSetup& setup,
    const gf::Row<19>* rows,Emit&& emit) {
    static_assert(Mode<=2);
    alignas(32) __m128i words[19]{},values[128],syndrome[19];
    gf::Packed<19> state{},feedback{};
    for(std::size_t epoch=n/128;epoch-->0;) {
        const auto base=128*epoch;
        if(epoch+1==n/128) {
            for(unsigned p=128;p-->0;){values[p]=input[base+p].mData;emit(base+p,input[base+p]);}
        } else {
            if constexpr(Mode==1)gf::unpack(state,words);
            k::imtReversePoints<Map>(input+base,values,words,base,emit,std::make_index_sequence<128>{});
        }
        if(!epoch)break;
        k::zeta<128>(values);Map::finish(values,syndrome);
        if(epoch+1==n/128) {
            if constexpr(Mode==1)gf::pack<19>(syndrome,state);
            else std::memcpy(words,syndrome,sizeof(words));
        } else if constexpr(Mode==0)fr::step<true>(words,setup.scalar[epoch],syndrome);
        else {
            if constexpr(Mode==2)gf::pack<19>(words,state);
            gf::pack<19>(syndrome,feedback);gf::step<19,true>(state,rows[epoch],&feedback);
            if constexpr(Mode==2)gf::unpack(state,words);
        }
    }
}
template<unsigned Mode> static SPIN_NOINLINE void encode(
    k::block* input,k::block* scratch,std::size_t n,const Packets<4>& packets,
    const Coefficients& coeff,const ScalarSetup& setup,const gf::Row<19>* rows,double* phases) {
    const auto start=phases?Clock::now():Clock::time_point{};
    pd::Emitter<false> emit{scratch,packets.bases.data(),packets.controls.data()};
    reverse<Mode>(input,n,setup,rows,emit);_mm_sfence();
    const auto middle=phases?Clock::now():Clock::time_point{};
    for(std::size_t tile=0;tile<n/1024;++tile)
        k::bchPackedCoeffCompact(scratch+tile*tileStride,input+tile*512,coeff.compact.data+tile*512);
    if(phases) {
        phases[0]+=std::chrono::duration<double,std::milli>(middle-start).count();
        phases[1]+=std::chrono::duration<double,std::milli>(Clock::now()-middle).count();
    }
}
template<unsigned Mode> static void experiment(unsigned exponent,std::uint64_t seed,unsigned calls,bool profile) {
    const unsigned n=2U<<exponent;
    const Packets<4> packets(n/256,seed,true,true);
    std::vector<unsigned> masks(std::size_t(8)*(n/128));k::setup::Words rng(seed^0x3f625a92ULL);
    for(std::size_t i=0;i<n/128;++i)sample<19>(masks.data()+8*i,rng);
    // Same principal allocation order as packed_coeff_r4.cpp. Every mode also
    // allocates every GFNI row AFTER the unchanged original setup allocations.
    std::vector<k::block> input(n),inner(n),actual(n),other(n),forward(n),expected(n),materialized(n);
    const auto scratchBlocks=(n/1024)*tileStride;std::vector<k::block> storage(scratchBlocks+4);
    auto* scratch=reinterpret_cast<k::block*>((reinterpret_cast<std::uintptr_t>(storage.data())+63)&~std::uintptr_t(63));
    k::workspace_routing::adviseOwned(scratch,scratchBlocks*16);
    const pd::Gl32 original(n/1024,seed);const ScalarSetup setup(masks.data(),n/128);
    const Coefficients coeff(original);std::vector<gf::Row<19>> rows(n/128);
    for(std::size_t i=0;i<rows.size();++i)rows[i]=gf::makeRow<19>(masks.data()+8*i);
    checkEpochs(masks.data(),rows);
    auto selected=[&](double* phases){encode<Mode>(input.data(),scratch,n,packets,coeff,setup,rows.data(),phases);};
    fill(other,971);forwardReference<4>(other.data(),forward.data(),n,masks.data());
    for(unsigned pattern=0;pattern<(calls?1U:3U);++pattern) {
        fill(input,pattern==2?997:913);
        if(pattern==1){std::fill(input.begin(),input.end(),k::block{});input[n-3]=k::block(1,3);}
        expected=input;reverseReference<4>(input.data(),inner.data(),n,masks.data());
        auto checkInner=[&]<unsigned M>() {
            reverse<M>(input.data(),n,setup,rows.data(),[&](std::size_t i,k::block v){actual[i]=v;});
            equal(actual.data(),inner.data(),n*16,"GFNI inner/original dense reference mismatch");
            auto lhs=_mm_setzero_si128(),rhs=lhs;
            for(std::size_t i=0;i<n;++i) {
                lhs=_mm_xor_si128(lhs,_mm_and_si128(input[i].mData,forward[i].mData));
                rhs=_mm_xor_si128(rhs,_mm_and_si128(other[i].mData,actual[i].mData));
            }
            equal(&lhs,&rhs,16,"GFNI inner/original forward adjoint mismatch");
        };
        checkInner.template operator()<0>();checkInner.template operator()<1>();checkInner.template operator()<2>();
        for(std::size_t i=0;i<n;++i) {
            const auto x=packets.inverse[i];
            materialized[(x/1024)*1024+4*(x%256)+(x/256)%4]=inner[i];
        }
        alignas(64) k::block mixed[1024];
        for(std::size_t tile=0;tile<n/1024;++tile) {
            pd::scalarMix(materialized.data()+1024*tile,mixed,original.coeff.data()+1024*tile);
            k::bchTranspose4(mixed,expected.data()+512*tile);
        }
        actual=input;pd::encode<true,4>(actual.data(),scratch,n,masks.data(),packets,original,nullptr);
        equal(actual.data(),expected.data(),n*16,"original sequential full/reference mismatch");
        auto checkFull=[&]<unsigned M>() {
            actual=input;encode<M>(actual.data(),scratch,n,packets,coeff,setup,rows.data(),nullptr);
            equal(actual.data(),expected.data(),n*16,"GFNI full/original reference/suffix mismatch");
        };
        checkFull.template operator()<0>();checkFull.template operator()<1>();checkFull.template operator()<2>();
        selected(nullptr);equal(input.data(),expected.data(),n*16,"selected full/reference mismatch");
    }
    std::cout<<"checks PASS: every-epoch basis, all3 inner modes/original dense+forward adjoint, "
        "original sequential/full scalar GL32+production BCH reference,suffix; gfni_row_bytes="<<sizeof(gf::Row<19>)
        <<"; packed_state_bytes="<<sizeof(gf::Packed<19>)<<"; gfni_setup_bytes="<<rows.size()*sizeof(rows[0])
        <<"; tile_mode="<<k::packedCoeffTileMode()<<'\n';
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
        <<times[calls/2]<<','<<times[calls/10]<<','<<times[9*calls/10]<<','<<std::hex<<hash(input)
        <<std::dec<<",4,"<<k::packedCoeffTileMode()<<'\n';
    if(profile)std::cout<<"phase_means_including_warmups,inner_route,"<<phases[0]/(calls+3)
        <<",packed_mixer_bch,"<<phases[1]/(calls+3)<<'\n';
}
}
int main(int argc,char** argv) {try {
    if(!__builtin_cpu_supports("avx512f") || !__builtin_cpu_supports("avx512vl") ||
       !__builtin_cpu_supports("avx512bw") || !__builtin_cpu_supports("gfni"))
        throw std::runtime_error("AVX512/GFNI required");
    if(argc==3 && (std::string(argv[1])=="core" || std::string(argv[1])=="basis")) {
        const auto seed=std::stoull(argv[2]);gp::checkCore<16>(seed);gp::checkCore<19>(seed);
        gp::checkDense16<4>(seed);gp::checkDense16<8>(seed);return 0;
    }
    if(argc<5 || argc>6)throw std::invalid_argument(
        "usage: gfni_r4_probe exponent seed calls mode [profile]; or core seed; "
        "0=scalar fused,1=persistent packed GFNI,2=roundtrip GFNI; calls=0 checks all modes");
    const unsigned exponent=std::stoul(argv[1]),calls=std::stoul(argv[3]),mode=std::stoul(argv[4]);
    const auto seed=std::stoull(argv[2]);const bool profile=argc==6 && std::stoul(argv[5]);
    if(exponent<14 || exponent>20 || mode>2)throw std::invalid_argument("unsupported GFNI probe");
    if(mode==0)gp::experiment<0>(exponent,seed,calls,profile);
    else if(mode==1)gp::experiment<1>(exponent,seed,calls,profile);
    else gp::experiment<2>(exponent,seed,calls,profile);
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
