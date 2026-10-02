#include "K16Paired15NoAlias.h"
#include "T64Paired15ShuffleMap.h"
#include "T64PairedOptimizedMap.h"
#include <stdexcept>
namespace spin::research::k16codesign {
namespace {
using namespace detail::packet::fast;
struct Four {__m512i a,b,c,d;};
static SPIN_FORCEINLINE Four transpose4(__m512i a,__m512i b,__m512i c,__m512i d) {
    const auto ab0=_mm512_shuffle_i32x4(a,b,0x44),ab1=_mm512_shuffle_i32x4(a,b,0xee);
    const auto cd0=_mm512_shuffle_i32x4(c,d,0x44),cd1=_mm512_shuffle_i32x4(c,d,0xee);
    return {_mm512_shuffle_i32x4(ab0,cd0,0x88),_mm512_shuffle_i32x4(ab0,cd0,0xdd),
            _mm512_shuffle_i32x4(ab1,cd1,0x88),_mm512_shuffle_i32x4(ab1,cd1,0xdd)};
}
static SPIN_FORCEINLINE __m512i totals4(__m512i a,__m512i b,__m512i c,__m512i d) {
    const auto t=transpose4(a,b,c,d);
    return _mm512_ternarylogic_epi64(t.a,t.b,_mm512_xor_si512(t.c,t.d),0x96);
}
static SPIN_FORCEINLINE void finishWide(const __m512i* high,WordPackets16& feedback) {
    const auto l0=_mm512_mask_xor_epi64(high[0],0x33,high[0],_mm512_shuffle_i32x4(high[0],high[0],0xf5));
    const auto low=_mm512_mask_xor_epi64(l0,0x0f,l0,_mm512_shuffle_i32x4(l0,l0,0xee));
    const auto middle=transpose4(high[1],high[2],high[4],high[8]);
    const auto odd=_mm512_xor_si512(middle.b,middle.d);
    const auto upper=_mm512_xor_si512(middle.c,middle.d);
    const auto total=_mm512_ternarylogic_epi64(middle.a,middle.c,odd,0x96);
    const auto quadratic=totals4(high[3],high[5],high[6],high[9]);
    const auto pair=_mm512_xor_si512(_mm512_shuffle_i32x4(high[10],high[12],0x44),
                                    _mm512_shuffle_i32x4(high[10],high[12],0xee));
    const auto both=_mm512_xor_si512(pair,_mm512_shuffle_i32x4(pair,pair,0xb1));
    // low=(z0,z1,z2,z3); total=(z4,z8,z16,z32)
    feedback.v[0]=low;
    feedback.v[1]=total;
    // odd=(z5,z9,z17,z33), upper=(z6,z10,z18,z34)
    // quadratic=(z12,z20,z24,z36), both=(z40,z40,z48,z48)
    // v2=(z9,0,z33+z12,z17+z36), v3=(z10+z20,z34+z24,z6+z40,z48).
    const auto left2=_mm512_maskz_shuffle_i32x4(0xff0f,odd,odd,0xb1);
    const auto right2=_mm512_shuffle_i32x4(quadratic,quadratic,0xc0);
    feedback.v[2]=_mm512_mask_xor_epi64(left2,0xf0,left2,right2);
    const auto left3=_mm512_shuffle_i32x4(upper,both,0x8d);
    const auto right3=_mm512_shuffle_i32x4(quadratic,upper,0x09);
    feedback.v[3]=_mm512_mask_xor_epi64(left3,0x3f,left3,right3);
}
static SPIN_FORCEINLINE void highFromPackets(const __m512i* packets,__m512i* high) {
    const auto a=firstStage<0>(packets),b=firstStage<4>(packets),c=firstStage<8>(packets),d=firstStage<12>(packets);
    const auto z4=_mm512_xor_si512(b.z0,d.z0),z8=_mm512_xor_si512(c.z0,d.z0);
    high[0]=_mm512_ternarylogic_epi64(a.z0,c.z0,z4,0x96);high[4]=z4;high[8]=z8;high[12]=d.z0;
    const auto z5=_mm512_xor_si512(b.z1,d.z1),z9=_mm512_xor_si512(c.z1,d.z1);
    high[1]=_mm512_ternarylogic_epi64(a.z1,c.z1,z5,0x96);high[5]=z5;high[9]=z9;
    const auto z6=_mm512_xor_si512(b.z2,d.z2),z10=_mm512_xor_si512(c.z2,d.z2);
    high[2]=_mm512_ternarylogic_epi64(a.z2,c.z2,z6,0x96);high[6]=z6;high[10]=z10;
    high[3]=_mm512_ternarylogic_epi64(a.z3,b.z3,_mm512_xor_si512(c.z3,d.z3),0x96);
}
template<bool Field>
static SPIN_FORCEINLINE void updateWide(const PackedState& before,const std::uint64_t* row,
    const WordPackets16& feedback,WordPackets16& words) {
    __m512i lo0,hi0,lo1,hi1;
    if constexpr(Field) {
        const auto m0=_mm512_set1_epi64(row[0]),m1=_mm512_set1_epi64(row[1]),m2=_mm512_set1_epi64(row[2]);
        const auto p00=_mm512_gf2p8affine_epi64_epi8(m0,before.v[0],0),p01=_mm512_gf2p8affine_epi64_epi8(m0,before.v[1],0);
        const auto p10=_mm512_gf2p8affine_epi64_epi8(m1,before.v[2],0),p11=_mm512_gf2p8affine_epi64_epi8(m1,before.v[3],0);
        const auto p20=_mm512_gf2p8affine_epi64_epi8(m2,_mm512_xor_si512(before.v[0],before.v[2]),0);
        const auto p21=_mm512_gf2p8affine_epi64_epi8(m2,_mm512_xor_si512(before.v[1],before.v[3]),0);
        lo0=_mm512_xor_si512(p00,p10);hi0=_mm512_xor_si512(p01,p11);
        lo1=_mm512_xor_si512(p00,p20);hi1=_mm512_xor_si512(p01,p21);
    } else {
        const auto m00=_mm512_set1_epi64(row[0]),m01=_mm512_set1_epi64(row[1]);
        const auto m10=_mm512_set1_epi64(row[2]),m11=_mm512_set1_epi64(row[3]);
        lo0=_mm512_xor_si512(_mm512_gf2p8affine_epi64_epi8(m00,before.v[0],0),_mm512_gf2p8affine_epi64_epi8(m01,before.v[2],0));
        hi0=_mm512_xor_si512(_mm512_gf2p8affine_epi64_epi8(m00,before.v[1],0),_mm512_gf2p8affine_epi64_epi8(m01,before.v[3],0));
        lo1=_mm512_xor_si512(_mm512_gf2p8affine_epi64_epi8(m10,before.v[0],0),_mm512_gf2p8affine_epi64_epi8(m11,before.v[2],0));
        hi1=_mm512_xor_si512(_mm512_gf2p8affine_epi64_epi8(m10,before.v[1],0),_mm512_gf2p8affine_epi64_epi8(m11,before.v[3],0));
    }
    words.v[0]=_mm512_xor_si512(_mm512_permutex2var_epi8(lo0,unpackIndex<0>(),hi0),feedback.v[0]);
    words.v[1]=_mm512_xor_si512(_mm512_permutex2var_epi8(lo0,unpackIndex<4>(),hi0),feedback.v[1]);
    words.v[2]=_mm512_xor_si512(_mm512_permutex2var_epi8(lo1,unpackIndex<0>(),hi1),feedback.v[2]);
    words.v[3]=_mm512_xor_si512(_mm512_permutex2var_epi8(lo1,unpackIndex<4>(),hi1),feedback.v[3]);
}
struct ByteRouteWide {
    std::byte* scratch;const std::uint32_t* offsets;
    SPIN_FORCEINLINE void operator()(std::size_t packet,__m512i value) {
        _mm512_store_si512(reinterpret_cast<__m512i*>(scratch+offsets[packet]),value);
    }
};
template<unsigned Variant>
static SPIN_NOINLINE void reverseWide(const Block* __restrict input,Block* __restrict scratch,const rs::Plan& plan,
    const PairedTables& tables,const PairedOptimizedTables& optimized) {
    static_assert(Variant<2);
    constexpr bool field=false,rawFeedback=false;
    ByteRouteWide emit{reinterpret_cast<std::byte*>(scratch),optimized.routeBytes.data()};
    alignas(64) __m512i packets[16],high[16];
    PackedState before;WordPackets16 state,feedback;
    for(std::size_t epoch=plan.n/64;epoch-->0;) {
        const auto* raw=input+64*epoch;const bool first=epoch+1==plan.n/64;
        if(first) {
            loadPackets(raw,packets,std::make_index_sequence<16>{});
            for(unsigned h=16;h-->0;)emit(16*epoch+h,packets[h]);
            if(!epoch)break;
            highFromPackets(packets,high);finishWide(high,state);continue;
        }
        if constexpr(rawFeedback) {
            if(epoch) {
                highFromPackets(reinterpret_cast<const __m512i*>(raw),high);finishWide(high,feedback);
                packWordGroup<0>(state,before);packWordGroup<1>(state,before);
                WordPackets16 next;
                const auto* row=field?optimized.fieldUpdates.data()+3*epoch:tables.updates.data()+4*epoch;
                updateWide<field>(before,row,feedback,next);
                pairedopt::emitOnly(state,raw,16*epoch,emit);state=next;
            } else pairedopt::emitOnly(state,raw,0,emit);
        } else {
            if(epoch){packWordGroup<0>(state,before);packWordGroup<1>(state,before);}
            paired15shuffle::streamHigh<(Variant!=0)>(state,raw,16*epoch,emit,high);
            if(!epoch)break;
            finishWide(high,feedback);
            const auto* row=field?optimized.fieldUpdates.data()+3*epoch:tables.updates.data()+4*epoch;
            updateWide<field>(before,row,feedback,state);
        }
    }
}
static SPIN_NOINLINE void reversePeeled(const Block* __restrict input,Block* __restrict scratch,
    const rs::Plan& plan,const PairedTables& tables,const PairedOptimizedTables& optimized) {
    const auto epochs=plan.n/64;
    if(!epochs)return;
    ByteRouteWide emit{reinterpret_cast<std::byte*>(scratch),optimized.routeBytes.data()};
    alignas(64) __m512i packets[16],high[16];
    PackedState before;WordPackets16 state,feedback;
    const auto first=epochs-1;
    loadPackets(input+64*first,packets,std::make_index_sequence<16>{});
    for(unsigned h=16;h-->0;)emit(16*first+h,packets[h]);
    if(!first)return;
    highFromPackets(packets,high);finishWide(high,state);
    for(std::size_t epoch=first;--epoch;) {
        packWordGroup<0>(state,before);packWordGroup<1>(state,before);
        paired15shuffle::streamHigh<true>(state,input+64*epoch,16*epoch,emit,high);
        finishWide(high,feedback);
        updateWide<false>(before,tables.updates.data()+4*epoch,feedback,state);
    }
    paired15shuffle::streamHigh<true>(state,input,0,emit,high);
}
}
void reverseRoutePaired15NoAlias(const Block* in,Block* scratch,const rs::Plan& plan,const PairedTables& tables,
    const PairedOptimizedTables& optimized,unsigned variant) {
    switch(variant) {
        case 0:reverseWide<1>(in,scratch,plan,tables,optimized);break;
        case 1:reversePeeled(in,scratch,plan,tables,optimized);break;
        default:throw std::invalid_argument("unknown no-alias paired15 variant");
    }
}
}
