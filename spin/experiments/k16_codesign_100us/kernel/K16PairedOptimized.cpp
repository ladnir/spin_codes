#include "K16PairedOptimized.h"
#include "T64PairedOptimizedMap.h"
#include <stdexcept>
namespace spin::research::k16codesign {
namespace {
using namespace detail::packet::fast;
struct ByteRoute {
    std::byte* scratch;const std::uint32_t* offsets;
    SPIN_FORCEINLINE void operator()(std::size_t packet,__m512i value) {
        _mm512_store_si512(reinterpret_cast<__m512i*>(scratch+offsets[packet]),value);
    }
};
// Karatsuba over GF256, performed with the GFNI operand roles swapped. Each
// changing packed state qword is the matrix operand; multiplication matrices
// are the data operands. Thus six affine instructions produce ordinary-word
// bytes directly, preserving the fused inverse transpose of the dense path.
static SPIN_FORCEINLINE void fieldUpdate(const PackedState& before,const std::uint64_t* row,
    const __m128i* feedback,WordPackets16& words) {
    const auto m0=_mm512_set1_epi64(row[0]),m1=_mm512_set1_epi64(row[1]),m2=_mm512_set1_epi64(row[2]);
    const auto p00=_mm512_gf2p8affine_epi64_epi8(m0,before.v[0],0);
    const auto p01=_mm512_gf2p8affine_epi64_epi8(m0,before.v[1],0);
    const auto p10=_mm512_gf2p8affine_epi64_epi8(m1,before.v[2],0);
    const auto p11=_mm512_gf2p8affine_epi64_epi8(m1,before.v[3],0);
    const auto p20=_mm512_gf2p8affine_epi64_epi8(m2,_mm512_xor_si512(before.v[0],before.v[2]),0);
    const auto p21=_mm512_gf2p8affine_epi64_epi8(m2,_mm512_xor_si512(before.v[1],before.v[3]),0);
    const auto lo0=_mm512_xor_si512(p00,p10),hi0=_mm512_xor_si512(p01,p11);
    const auto lo1=_mm512_xor_si512(p00,p20),hi1=_mm512_xor_si512(p01,p21);
    words.v[0]=_mm512_xor_si512(_mm512_permutex2var_epi8(lo0,unpackIndex<0>(),hi0),join(feedback[0],feedback[1],feedback[2],feedback[3]));
    words.v[1]=_mm512_xor_si512(_mm512_permutex2var_epi8(lo0,unpackIndex<4>(),hi0),join(feedback[4],feedback[5],feedback[6],feedback[7]));
    words.v[2]=_mm512_xor_si512(_mm512_permutex2var_epi8(lo1,unpackIndex<0>(),hi1),join(feedback[8],feedback[9],feedback[10],feedback[11]));
    words.v[3]=_mm512_xor_si512(_mm512_permutex2var_epi8(lo1,unpackIndex<4>(),hi1),join(feedback[12],feedback[13],feedback[14],feedback[15]));
}
template<unsigned Variant>
static SPIN_NOINLINE void reverseOptimized(const Block* input,Block* scratch,const rs::Plan& plan,
    const PairedTables& tables,const PairedOptimizedTables& optimized) {
    constexpr bool rawFeedback=Variant==1 || Variant==4;
    constexpr bool field=Variant==3 || Variant==4;
    ByteRoute emit{reinterpret_cast<std::byte*>(scratch),optimized.routeBytes.data()};
    alignas(64) __m128i moments[64],syndrome[16];
    alignas(64) __m512i packets[16];
    PackedState before;WordPackets16 state;
    for(std::size_t epoch=plan.n/64;epoch-->0;) {
        const auto* raw=input+64*epoch;
        const bool first=epoch+1==plan.n/64;
        if(first) {
            loadPackets(raw,packets,std::make_index_sequence<16>{});
            for(unsigned h=16;h-->0;)emit(16*epoch+h,packets[h]);
            if(!epoch)break;
            packetMoments(packets,moments);paired::finishMono(moments,syndrome);
            wordPackets(syndrome,state);
            continue;
        }
        if constexpr(rawFeedback) {
            // The source loads are safe without a temporary packet array:
            // firstStage uses unaligned loads and reads exactly 64 records.
            if(epoch) {
                packetMoments(reinterpret_cast<const __m512i*>(raw),moments);
                paired::finishMono(moments,syndrome);
                packWordGroup<0>(state,before);packWordGroup<1>(state,before);
                WordPackets16 next;
                if constexpr(field)fieldUpdate(before,optimized.fieldUpdates.data()+3*epoch,syndrome,next);
                else applyPackedToPackets(before,tables.updates.data()+4*epoch,syndrome,next);
                pairedopt::emitOnly(state,raw,16*epoch,emit);
                state=next;
            } else pairedopt::emitOnly(state,raw,0,emit);
        } else {
            if(epoch){packWordGroup<0>(state,before);packWordGroup<1>(state,before);}
            if constexpr(Variant==2)pairedopt::streamTernary(state,raw,16*epoch,emit,moments);
            else paired::stream(state,raw,16*epoch,emit,moments);
            if(!epoch)break;
            paired::finishMono(moments,syndrome);
            if constexpr(field)fieldUpdate(before,optimized.fieldUpdates.data()+3*epoch,syndrome,state);
            else applyPackedToPackets(before,tables.updates.data()+4*epoch,syndrome,state);
        }
    }
}
}
void reverseRoutePairedOptimized(const Block* in,Block* scratch,const rs::Plan& plan,
    const PairedTables& tables,const PairedOptimizedTables& optimized,unsigned variant) {
    switch(variant) {
        case 0:reverseOptimized<0>(in,scratch,plan,tables,optimized);break;
        case 1:reverseOptimized<1>(in,scratch,plan,tables,optimized);break;
        case 2:reverseOptimized<2>(in,scratch,plan,tables,optimized);break;
        case 3:reverseOptimized<3>(in,scratch,plan,tables,optimized);break;
        case 4:reverseOptimized<4>(in,scratch,plan,tables,optimized);break;
        default:throw std::invalid_argument("unknown paired optimized variant");
    }
}
}
