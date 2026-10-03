#include "Packet8Wide24.h"
#include <bit>
#include <cstring>

namespace spin::research::packet8wide24 {
namespace {
// Literal field columns, independent of the SIMD multiplication circuits.
void addMultiply8(Block* output,const Block* input,std::uint8_t scalar,bool adjoint) {
    for(unsigned column=0; column<8; ++column) {
        const auto image=multiply8(scalar,std::uint8_t(1U<<column));
        for(unsigned row=0; row<8; ++row) if((image>>row)&1U)
            output[adjoint?column:row]^=input[adjoint?row:column];
    }
}
void addMultiply24(Block* output,const Block* input,std::uint32_t scalar,bool adjoint) {
    for(unsigned column=0; column<24; ++column) {
        auto image=multiply24(scalar,std::uint32_t(1U<<column));
        for(; image; image&=image-1) {
            const auto row=std::countr_zero(image);
            output[adjoint?column:row]^=input[adjoint?row:column];
        }
    }
}

constexpr unsigned multiply4(unsigned a,unsigned b) {
    unsigned result=0;
    for(unsigned i=0; i<4; ++i) {
        if(b&1U) result^=a;
        b>>=1; a<<=1; if(a&16U) a^=0x13U;
    }
    return result;
}
constexpr unsigned inverse4(unsigned value) {
    unsigned result=1;
    for(unsigned i=0; i<14; ++i) result=multiply4(result,value);
    return result;
}
// Eight independent RS[16,8] rows use this plain Lagrange generator.
// A message basis bit generates a 64-bit word in one RS row.
constexpr auto rsRows=[] {
    std::array<std::uint64_t,32> result{};
    for(unsigned messageBit=0; messageBit<32; ++messageBit)
        for(unsigned symbol=0; symbol<16; ++symbol) {
            unsigned numerator=1,denominator=1;
            for(unsigned other=0; other<8; ++other) if(other!=messageBit/4) {
                numerator=multiply4(numerator,symbol^other);
                denominator=multiply4(denominator,(messageBit/4)^other);
            }
            const auto coefficient=multiply4(numerator,inverse4(denominator));
            result[messageBit]|=std::uint64_t(multiply4(coefficient,1U<<(messageBit%4)))<<(4*symbol);
        }
    return result;
}();

void outerForward(const Block* message,Block* raw,const Plan& plan) {
    for(std::size_t group=0; group<plan.groups; ++group) {
        Block coded[512];
        for(unsigned lane=0; lane<8; ++lane)
            for(unsigned bit=0; bit<64; ++bit) {
                Block value{};
                for(unsigned messageBit=0; messageBit<32; ++messageBit)
                    if((rsRows[messageBit]>>bit)&1U)
                        value^=message[256*group+32*lane+messageBit];
                coded[32*(bit/4)+4*lane+bit%4]=value;
            }
        const auto* rows=plan.outerRows.data()+16*group;
        auto* output=raw+groupStride*group;
        for(unsigned symbol=0; symbol<16; ++symbol)
            for(unsigned physical=0; physical<32; ++physical) {
                Block value{};
                for(unsigned logical=0; logical<32; ++logical)
                    if((rows[symbol][logical]>>physical)&1U)
                        value^=coded[32*symbol+logical];
                output[32*symbol+physical]=value;
            }
    }
}

template<bool Forward>
void inner(const Block* input,Block* output,const Plan& plan) {
    Block state[24]{};
    for(std::size_t step=0; step<plan.n/64; ++step) {
        const auto epoch=Forward?step:plan.n/64-1-step;
        Block feedback[24]{};
        for(unsigned packet=0; packet<8; ++packet) {
            const auto linear=64*epoch+8*packet;
            const auto routed=plan.route[8*epoch+packet];
            Block raw[8],value[8];
            std::memcpy(raw,input+(Forward?routed:linear),sizeof(raw));
            for(unsigned bit=0; bit<8; ++bit) {
                value[bit]=raw[bit]; value[bit]^=state[bit]; feedback[bit]^=raw[bit];
            }
            const auto h=std::uint8_t(packet),square=multiply8(h,h);
            addMultiply8(value,state+8,h,!Forward);
            addMultiply8(value,state+16,square,!Forward);
            addMultiply8(feedback+8,raw,h,!Forward);
            addMultiply8(feedback+16,raw,square,!Forward);
            std::memcpy(output+(Forward?linear:routed),value,sizeof(value));
        }
        Block next[24]; std::memcpy(next,feedback,sizeof(next));
        // The forward sampled update is the binary adjoint of the ordinary
        // field multiplication stored for reverse evaluation.
        addMultiply24(next,state,plan.updates[epoch].scalar,Forward);
        std::memcpy(state,next,sizeof(state));
    }
}
}

void reverseRouteScalar(const Block* input,Block* raw,const Plan& plan) {
    inner<false>(input,raw,plan);
}
void outerScalar(const Block* raw,Block* output,const Plan& plan) {
    for(std::size_t group=0; group<plan.groups; ++group) {
        const auto* input=raw+groupStride*group;
        const auto* rows=plan.outerRows.data()+16*group;
        Block mixed[512];
        for(unsigned symbol=0; symbol<16; ++symbol)
            for(unsigned logical=0; logical<32; ++logical) {
                Block value{};
                for(auto mask=rows[symbol][logical]; mask; mask&=mask-1)
                    value^=input[32*symbol+std::countr_zero(mask)];
                mixed[32*symbol+logical]=value;
            }
        for(unsigned lane=0; lane<8; ++lane)
            for(unsigned messageBit=0; messageBit<32; ++messageBit) {
                Block value{};
                for(auto mask=rsRows[messageBit]; mask; mask&=mask-1) {
                    const auto bit=std::countr_zero(mask);
                    value^=mixed[32*(bit/4)+4*lane+bit%4];
                }
                output[256*group+32*lane+messageBit]=value;
            }
    }
}
void transposeScalar(const Block* input,Block* output,Block* raw,const Plan& plan) {
    reverseRouteScalar(input,raw,plan);
    outerScalar(raw,output,plan);
}
void forwardScalar(const Block* message,Block* output,const Plan& plan) {
    std::vector<Block> raw(plan.scratchBlocks());
    outerForward(message,raw.data(),plan);
    inner<true>(raw.data(),output,plan);
}
void packScratchScalar(const Block* raw,Block* packed,const Plan& plan) {
    for(std::size_t group=0; group<plan.groups; ++group)
        for(unsigned packet=0; packet<64; ++packet) {
            const auto index=group*groupStride+8*packet;
            auto* output=reinterpret_cast<std::uint8_t*>(packed+index);
            const auto* input=reinterpret_cast<const std::uint8_t*>(raw+index);
            for(unsigned half=0; half<2; ++half)
                for(unsigned byte=0; byte<8; ++byte)
                    for(unsigned bit=0; bit<8; ++bit) {
                        unsigned value=0;
                        for(unsigned coordinate=0; coordinate<8; ++coordinate)
                            value|=((input[16*coordinate+8*half+byte]>>bit)&1U)<<coordinate;
                        output[64*half+8*byte+7-bit]=std::uint8_t(value);
                    }
        }
}
void unpackScratchScalar(const Block* packed,Block* raw,const Plan& plan) {
    for(std::size_t group=0; group<plan.groups; ++group)
        for(unsigned packet=0; packet<64; ++packet) {
            const auto index=group*groupStride+8*packet;
            const auto* input=reinterpret_cast<const std::uint8_t*>(packed+index);
            auto* output=reinterpret_cast<std::uint8_t*>(raw+index);
            for(unsigned coordinate=0; coordinate<8; ++coordinate)
                for(unsigned byte=0; byte<16; ++byte) {
                    unsigned value=0;
                    for(unsigned bit=0; bit<8; ++bit)
                        value|=((input[64*(byte/8)+8*(byte%8)+7-bit]>>coordinate)&1U)<<bit;
                    output[16*coordinate+byte]=std::uint8_t(value);
                }
        }
}
}
