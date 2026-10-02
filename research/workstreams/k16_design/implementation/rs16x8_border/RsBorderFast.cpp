#include "RsBorder.h"
#include "BorderEvaluator.h"
namespace spin::research::rsborder {
namespace {
using namespace detail::packet::large;
struct PackedExtra {__m512i lo,hi;};
template<unsigned Extra> static SPIN_FORCEINLINE PackedExtra packExtra(const __m128i* words) {
    const auto zero=_mm_setzero_si128();
    const auto a=join(words[0],Extra>=2?words[1]:zero,Extra>=3?words[2]:zero,Extra>=4?words[3]:zero);
    const auto b=_mm512_setzero_si512();
    const auto basis=_mm512_set1_epi64(0x0102040810204080ULL);
    return {_mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(a,packIndex<0>(),b),0),
            _mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(a,packIndex<8>(),b),0)};
}
template<unsigned Extra> static SPIN_FORCEINLINE void unpackExtra(const PackedExtra& packed,__m128i* words) {
    const auto basis=_mm512_set1_epi64(0x8040201008040201ULL);
    const auto lo=_mm512_gf2p8affine_epi64_epi8(basis,packed.lo,0);
    const auto hi=_mm512_gf2p8affine_epi64_epi8(basis,packed.hi,0);
    const auto joint=_mm512_permutex2var_epi8(lo,unpackIndex<0>(),hi);
    words[0]=_mm512_castsi512_si128(joint);
    if constexpr(Extra>=2)words[1]=_mm512_extracti32x4_epi32(joint,1);
    if constexpr(Extra>=3)words[2]=_mm512_extracti32x4_epi32(joint,2);
    if constexpr(Extra>=4)words[3]=_mm512_extracti32x4_epi32(joint,3);
}
static SPIN_FORCEINLINE void update(PackedState& state,PackedExtra& extra,const PackedState& feedback,
    const PackedExtra& extraFeedback,const BorderRow& row) {
    const auto lower0=_mm512_set1_epi64(row.lower[0]),lower1=_mm512_set1_epi64(row.lower[1]);
    const auto corner=_mm512_set1_epi64(row.corner);
    PackedExtra next;
    next.lo=_mm512_xor_si512(_mm512_ternarylogic_epi64(
        _mm512_gf2p8affine_epi64_epi8(state.v[0],lower0,0),
        _mm512_gf2p8affine_epi64_epi8(state.v[2],lower1,0),
        _mm512_gf2p8affine_epi64_epi8(extra.lo,corner,0),0x96),extraFeedback.lo);
    next.hi=_mm512_xor_si512(_mm512_ternarylogic_epi64(
        _mm512_gf2p8affine_epi64_epi8(state.v[1],lower0,0),
        _mm512_gf2p8affine_epi64_epi8(state.v[3],lower1,0),
        _mm512_gf2p8affine_epi64_epi8(extra.hi,corner,0),0x96),extraFeedback.hi);
    const auto upper0=_mm512_set1_epi64(row.upper[0]),upper1=_mm512_set1_epi64(row.upper[1]);
    const auto u0=_mm512_gf2p8affine_epi64_epi8(extra.lo,upper0,0);
    const auto u1=_mm512_gf2p8affine_epi64_epi8(extra.hi,upper0,0);
    const auto u2=_mm512_gf2p8affine_epi64_epi8(extra.lo,upper1,0);
    const auto u3=_mm512_gf2p8affine_epi64_epi8(extra.hi,upper1,0);
    denseStep16(state,row.base,&feedback);
    state.v[0]=_mm512_xor_si512(state.v[0],u0);state.v[1]=_mm512_xor_si512(state.v[1],u1);
    state.v[2]=_mm512_xor_si512(state.v[2],u2);state.v[3]=_mm512_xor_si512(state.v[3],u3);
    extra=next;
}
struct StreamingRoute {
    Block* scratch;const std::uint32_t* bases;
    SPIN_FORCEINLINE void operator()(std::size_t packet,__m512i value) {
        _mm512_stream_si512(reinterpret_cast<__m512i*>(scratch+bases[packet]),value);
    }
};
template<unsigned Extra,class Emit>
static SPIN_NOINLINE void reverseBorder(const Block* input,std::size_t n,const BorderRow* rows,Emit& emit) {
    alignas(64) __m128i words[16],extraWords[4],moments[64],syndrome[16],extraSyndrome[4];
    alignas(64) __m512i packets[16];
    PackedState state{},feedback;PackedExtra extra{};
    for(std::size_t epoch=n/64;epoch-->0;) {
        const auto* raw=input+64*epoch;const bool first=epoch+1==n/64;
        if(first) {
            loadPackets(raw,packets,std::make_index_sequence<16>{});
            for(unsigned h=16;h-->0;)emit(16*epoch+h,packets[h]);
        } else {
            wideUnpackVbmi(state,words);unpackExtra<Extra>(extra,extraWords);
            fast::streamStepBorder<Extra>(words,extraWords,raw,16*epoch,emit,moments);
        }
        if(!epoch)break;
        if(first)packetMoments(packets,moments);
        finish(moments,syndrome);widePackVbmi(syndrome,feedback);
        extraSyndrome[0]=moments[3];
        if constexpr(Extra>=2)extraSyndrome[1]=moments[5];
        if constexpr(Extra>=3)extraSyndrome[2]=moments[9];
        if constexpr(Extra>=4)extraSyndrome[3]=moments[17];
        const auto extraFeedback=packExtra<Extra>(extraSyndrome);
        if(first){state=feedback;extra=extraFeedback;}
        else update(state,extra,feedback,extraFeedback,rows[epoch]);
    }
}
}
void reverseRoute(const Block* input,Block* scratch,const Plan& plan) {
    if(plan.stateBits==16){rswide::reverseRoute(input,scratch,plan.outer);return;}
    StreamingRoute emit{scratch,plan.outer.route.data()};
    if(plan.stateBits==17)reverseBorder<1>(input,plan.outer.n,plan.updates.data(),emit);
    else if(plan.stateBits==18)reverseBorder<2>(input,plan.outer.n,plan.updates.data(),emit);
    else if(plan.stateBits==19)reverseBorder<3>(input,plan.outer.n,plan.updates.data(),emit);
    else reverseBorder<4>(input,plan.outer.n,plan.updates.data(),emit);
    _mm_sfence();
}
void transposeFast(const Block* input,Block* output,Block* scratch,const Plan& plan) {
    reverseRoute(input,scratch,plan);rswide::outerFast(scratch,output,plan.outer);
}
}
