#include "PacketPlan.h"
#include <bit>
#include <cstring>

namespace spin::detail::packet {
namespace {
// Literal expansion map of the measured t64/s20 inner. The last four
// coordinates are x0*x1, x0*x2, x0*x3, x0*x4; feedback is its transpose.
constexpr std::uint16_t columns[64]={
    0x1,0x3,0x5,0x3507,0x9,0x118b,0x940d,0xb08f,
    0x11,0x5993,0xe415,0x8897,0xd999,0x919b,0xa99d,0xd49f,
    0x21,0x2ca3,0xb8a5,0xa127,0x1a29,0x272b,0x36ad,0x3eaf,
    0x6631,0x1333,0x3ab5,0x7ab7,0xa5b9,0xc13b,0x6d3d,0x3cbf,
    0x41,0xcfc3,0xe45,0xf4c7,0xee49,0x304b,0x744d,0x9f4f,
    0xb6d1,0x20d3,0x5cd5,0xffd7,0x8159,0x6db,0xff5d,0x4ddf,
    0x78e1,0x9be3,0xce65,0x1867,0x8ce9,0x7e6b,0xae6d,0x69ef,
    0xa871,0x12f3,0xfaf5,0x7577,0x85f9,0x2efb,0x437d,0xdd7f};
constexpr unsigned column(unsigned p) {
    return unsigned(columns[p])|(((p&1U)*((p>>1)&1U))<<16)
        |(((p&1U)*((p>>2)&1U))<<17)|(((p&1U)*((p>>3)&1U))<<18)
        |(((p&1U)*((p>>4)&1U))<<19);
}
constexpr unsigned multiply(unsigned a,unsigned b) {
    unsigned result=0;
    for(unsigned i=0;i<4;++i) {
        if(b&1)result^=a;
        b>>=1;a<<=1;if(a&16)a^=0x13;
    }
    return result;
}
constexpr unsigned inverse(unsigned a) {
    unsigned result=1;
    for(unsigned i=0;i<14;++i)result=multiply(result,a);
    return result;
}
// Independent Lagrange generator, not the optimized parity factorization.
constexpr auto rsRows=[] {
    std::array<std::uint64_t,32> rows{};
    for(unsigned messageBit=0;messageBit<32;++messageBit)
        for(unsigned symbol=0;symbol<16;++symbol) {
            unsigned numerator=1,denominator=1;
            for(unsigned other=0;other<8;++other)
                if(other!=messageBit/4) {
                    numerator=multiply(numerator,symbol^other);
                    denominator=multiply(denominator,(messageBit/4)^other);
                }
            rows[messageBit]|=std::uint64_t(multiply(multiply(numerator,inverse(denominator)),1U<<(messageBit%4)))<<(4*symbol);
        }
    return rows;
}();
constexpr unsigned index(unsigned symbol,unsigned coordinate) {
    return 32*symbol+16*(coordinate/16)+4*(coordinate%4)+(coordinate%16)/4;
}

void outerForwardGroup(const Block* message,Block* routed,const std::array<std::uint32_t,32>* matrices) {
    Block coded[512];
    for(unsigned lane=0;lane<8;++lane)
        for(unsigned bit=0;bit<64;++bit) {
            auto value=_mm_setzero_si128();
            for(unsigned messageBit=0;messageBit<32;++messageBit)
                if((rsRows[messageBit]>>bit)&1)value=_mm_xor_si128(value,message[32*lane+messageBit].mData);
            coded[index(bit/4,4*lane+bit%4)]=Block(value);
        }
    for(unsigned symbol=0;symbol<16;++symbol)
        for(unsigned row=0;row<32;++row) {
            auto value=_mm_setzero_si128();
            for(unsigned c=0;c<32;++c)
                if((matrices[symbol][c]>>row)&1)value=_mm_xor_si128(value,coded[index(symbol,c)].mData);
            routed[index(symbol,row)]=Block(value);
        }
}

void outerScalar(const Block* scratch,Block* output,const Plan& plan) {
    for(std::size_t group=0;group<plan.groups;++group) {
        const auto* routed=scratch+groupStride*group;const auto* matrices=plan.outerMatrices.data()+16*group;
        Block mixed[512];
        for(unsigned symbol=0;symbol<16;++symbol)
            for(unsigned row=0;row<32;++row) {
                auto value=_mm_setzero_si128();
                for(auto mask=matrices[symbol][row];mask;mask&=mask-1)
                    value=_mm_xor_si128(value,routed[index(symbol,std::countr_zero(mask))].mData);
                mixed[index(symbol,row)]=Block(value);
            }
        for(unsigned lane=0;lane<8;++lane)
            for(unsigned messageBit=0;messageBit<32;++messageBit) {
                auto value=_mm_setzero_si128();
                for(auto mask=rsRows[messageBit];mask;mask&=mask-1) {
                    const auto bit=std::countr_zero(mask);
                    value=_mm_xor_si128(value,mixed[index(bit/4,4*lane+bit%4)].mData);
                }
                output[256*group+32*lane+messageBit]=Block(value);
            }
    }
}
}

void transposeScalar(const Block* input,Block* output,Block* scratch,const Plan& plan) {
    __m128i state[20]{},next[20]{},feedback[20];
    for(std::size_t epoch=plan.n/64;epoch-->0;) {
        for(auto& value:feedback)value=_mm_setzero_si128();
        for(unsigned p=0;p<64;++p) {
            const auto i=64*epoch+p;const auto raw=input[i].mData;auto value=raw;
            for(unsigned mask=column(p);mask;mask&=mask-1) {
                const auto j=std::countr_zero(mask);value=_mm_xor_si128(value,state[j]);
                feedback[j]=_mm_xor_si128(feedback[j],raw);
            }
            scratch[plan.route[i/4]+(i&3)]=Block(value);
        }
        if(!epoch)break;
        for(unsigned j=0;j<20;++j) {
            auto value=feedback[j];
            for(unsigned mask=plan.reverseMatrices[epoch][j];mask;mask&=mask-1)value=_mm_xor_si128(value,state[std::countr_zero(mask)]);
            next[j]=value;
        }
        std::memcpy(state,next,sizeof(state));
    }
    outerScalar(scratch,output,plan);
}

void forwardScalar(const Block* message,Block* encoded,const Plan& plan) {
    std::vector<Block> routed(plan.scratchBlocks());
    forwardScalar(message,encoded,routed.data(),plan);
}

void forwardScalar(const Block* message,Block* encoded,Block* routed,const Plan& plan) {
    for(std::size_t group=0;group<plan.groups;++group)
        outerForwardGroup(message+256*group,routed+groupStride*group,plan.outerMatrices.data()+16*group);
    for(std::size_t i=0;i<plan.n;++i)encoded[i]=routed[plan.route[i/4]+(i&3)];
    __m128i state[20]{},next[20]{},feedback[20];
    for(std::size_t epoch=0;epoch<plan.n/64;++epoch) {
        for(auto& value:feedback)value=_mm_setzero_si128();
        for(unsigned p=0;p<64;++p) {
            const auto i=64*epoch+p;const auto raw=encoded[i].mData;auto value=raw;
            for(unsigned mask=column(p);mask;mask&=mask-1) {
                const auto j=std::countr_zero(mask);value=_mm_xor_si128(value,state[j]);
                feedback[j]=_mm_xor_si128(feedback[j],raw);
            }
            encoded[i]=Block(value);
        }
        if(epoch+1==plan.n/64)break;
        for(unsigned j=0;j<20;++j) {
            auto value=feedback[j];
            for(unsigned c=0;c<20;++c)if((plan.reverseMatrices[epoch][c]>>j)&1)value=_mm_xor_si128(value,state[c]);
            next[j]=value;
        }
        std::memcpy(state,next,sizeof(state));
    }
}
template<unsigned L> static void scalarWide(const Block* message,Block* encoded,Block* routed,const Plan& plan) {
    alignas(16) Block local[256];
    for(unsigned lane=0;lane<L;++lane) {
        for(std::size_t group=0;group<plan.groups;++group) {
            for(unsigned i=0;i<256;++i)local[i]=message[(256*group+i)*L+lane];
            outerForwardGroup(local,routed+groupStride*group,plan.outerMatrices.data()+16*group);
        }
        __m128i state[20]{},next[20]{},feedback[20];
        for(std::size_t epoch=0;epoch<plan.n/64;++epoch) {
            for(auto& value:feedback)value=_mm_setzero_si128();
            for(unsigned p=0;p<64;++p) {
                const auto i=64*epoch+p;const auto raw=routed[plan.route[i/4]+(i&3)].mData;auto value=raw;
                for(unsigned mask=column(p);mask;mask&=mask-1) {
                    const auto j=std::countr_zero(mask);value=_mm_xor_si128(value,state[j]);
                    feedback[j]=_mm_xor_si128(feedback[j],raw);
                }
                encoded[i*L+lane]=Block(value);
            }
            if(epoch+1==plan.n/64)break;
            for(unsigned j=0;j<20;++j) {
                auto value=feedback[j];
                for(unsigned c=0;c<20;++c)if((plan.reverseMatrices[epoch][c]>>j)&1)value=_mm_xor_si128(value,state[c]);
                next[j]=value;
            }
            std::memcpy(state,next,sizeof(state));
        }
    }
}
void forwardScalarWide(const Block* in,Block* out,Block* scratch,const Plan& p,unsigned lanes) {
    if(lanes==2)scalarWide<2>(in,out,scratch,p);else scalarWide<4>(in,out,scratch,p);
}
}
