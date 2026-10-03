#include "Packet8.h"
#include <cstring>

namespace spin::research::k16codesign {
void outerScalar(const rs::Block*,rs::Block*,const rs::Plan&);
void outerForwardScalar(const rs::Block*,rs::Block*,const rs::Plan&);
}
namespace spin::research::packet8 {
namespace {
// Build each product from literal field columns. No GFNI tables or packed
// kernel helpers are used by these independent binary-coordinate oracles.
void addMultiply(Block* output,const Block* input,std::uint8_t scalar,bool adjoint) {
    for(unsigned column=0;column<8;++column) {
        const auto image=multiply(scalar,std::uint8_t(1U<<column));
        for(unsigned row=0;row<8;++row)if((image>>row)&1U)
            output[adjoint?column:row]^=input[adjoint?row:column];
    }
}
template<bool Forward>
void inner(const Block* input,Block* output,const Plan& plan) {
    Block state[16]{};
    for(std::size_t step=0;step<plan.n()/64;++step) {
        const auto epoch=Forward?step:plan.n()/64-1-step;
        Block feedback[16]{};
        for(unsigned packet=0;packet<8;++packet) {
            const auto linear=64*epoch+8*packet;
            const auto routed=plan.route[8*epoch+packet];
            Block raw[8],value[8];
            std::memcpy(raw,input+(Forward?routed:linear),sizeof(raw));
            for(unsigned b=0;b<8;++b) {
                value[b]=raw[b]; value[b]^=state[b]; feedback[b]^=raw[b];
            }
            addMultiply(value,state+8,std::uint8_t(packet),!Forward);
            addMultiply(feedback+8,raw,std::uint8_t(packet),!Forward);
            std::memcpy(output+(Forward?linear:routed),value,sizeof(value));
        }
        Block next[16];std::memcpy(next,feedback,sizeof(next));
        const auto& m=plan.updates[epoch].forward;
        if constexpr(Forward) {
            addMultiply(next,state,m[0],false);addMultiply(next,state+8,m[1],false);
            addMultiply(next+8,state,m[2],false);addMultiply(next+8,state+8,m[3],false);
        } else {
            addMultiply(next,state,m[0],true);addMultiply(next,state+8,m[2],true);
            addMultiply(next+8,state,m[1],true);addMultiply(next+8,state+8,m[3],true);
        }
        std::memcpy(state,next,sizeof(state));
    }
}
}
void reverseRouteScalar(const Block* in,Block* raw,const Plan& plan){inner<false>(in,raw,plan);}
void transposeScalar(const Block* in,Block* out,Block* raw,const Plan& plan) {
    reverseRouteScalar(in,raw,plan);k16codesign::outerScalar(raw,out,plan.outer);
}
void forwardScalar(const Block* message,Block* out,const Plan& plan) {
    std::vector<Block> raw(plan.scratchBlocks());
    k16codesign::outerForwardScalar(message,raw.data(),plan.outer);
    inner<true>(raw.data(),out,plan);
}
void packScratchScalar(const Block* raw,Block* packed,const Plan& plan) {
    for(std::size_t g=0;g<plan.outer.groups;++g)for(unsigned packet=0;packet<32;++packet) {
        const auto index=g*rs::groupStride+8*packet;
        auto* bytes=reinterpret_cast<std::uint8_t*>(packed+index);
        const auto* input=reinterpret_cast<const std::uint8_t*>(raw+index);
        for(unsigned half=0;half<2;++half)for(unsigned byte=0;byte<8;++byte)
            for(unsigned bit=0;bit<8;++bit) {
                unsigned value=0;
                for(unsigned c=0;c<8;++c)value|=((input[16*c+8*half+byte]>>bit)&1U)<<c;
                // The retained GFNI pack uses descending payload-bit lanes.
                bytes[64*half+8*byte+7-bit]=std::uint8_t(value);
            }
    }
}
void unpackScratchScalar(const Block* packed,Block* raw,const Plan& plan) {
    for(std::size_t g=0;g<plan.outer.groups;++g)for(unsigned packet=0;packet<32;++packet) {
        const auto index=g*rs::groupStride+8*packet;
        const auto* bytes=reinterpret_cast<const std::uint8_t*>(packed+index);
        auto* out=reinterpret_cast<std::uint8_t*>(raw+index);
        for(unsigned c=0;c<8;++c)for(unsigned byte=0;byte<16;++byte) {
            unsigned value=0;
            for(unsigned bit=0;bit<8;++bit)value|=((bytes[64*(byte/8)+8*(byte%8)+7-bit]>>c)&1U)<<bit;
            out[16*c+byte]=std::uint8_t(value);
        }
    }
}
}
