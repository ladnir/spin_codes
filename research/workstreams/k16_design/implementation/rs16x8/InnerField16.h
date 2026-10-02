// Exact transitive state-map family. Only the dense update differs from
// PacketLargeInner::reverse; evaluation, feedback and routing stay retained.
#pragma once
#include "Tower32ByteRandomizer.h"
#include "../../../../../spin/src/packet/PacketLargeInner.h"
namespace spin::research::rswide {
template<class Emit>
SPIN_FORCEINLINE void reverseField16(const Block* input,std::size_t n,
    const std::array<std::uint8_t,3>* rows,Emit& emit) {
    using namespace detail::packet::large;
    alignas(64) __m128i words[16],moments[64],syndrome[16];
    alignas(64) __m512i packets[16];
    PackedState state{},feedback;
    for(std::size_t epoch=n/64;epoch-->0;) {
        const auto* raw=input+64*epoch;
        const bool first=epoch+1==n/64;
        if(first) {
            loadPackets(raw,packets,std::make_index_sequence<16>{});
            for(unsigned h=16;h-->0;)emit(16*epoch+h,packets[h]);
        } else {
            wideUnpackVbmi(state,words);
            streamStep(words,raw,16*epoch,emit,moments);
        }
        if(!epoch)break;
        if(first)packetMoments(packets,moments);
        finish(moments,syndrome);
        widePackVbmi(syndrome,feedback);
        if(first)state=feedback;
        else {
            const auto lo=rs::tower32byte::quadraticMultiply(state.v[0],state.v[2],rows[epoch].data());
            const auto hi=rs::tower32byte::quadraticMultiply(state.v[1],state.v[3],rows[epoch].data());
            state.v[0]=_mm512_xor_si512(lo.lo,feedback.v[0]);
            state.v[1]=_mm512_xor_si512(hi.lo,feedback.v[1]);
            state.v[2]=_mm512_xor_si512(lo.hi,feedback.v[2]);
            state.v[3]=_mm512_xor_si512(hi.hi,feedback.v[3]);
        }
    }
}
}
