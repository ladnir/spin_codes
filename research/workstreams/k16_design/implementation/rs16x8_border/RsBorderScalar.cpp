#include "RsBorder.h"
#include <bit>
#include <cstring>
namespace spin::research::rsborder {
namespace {
constexpr std::uint16_t columns[64]={
    0x1,0x3,0x5,0x3507,0x9,0x118b,0x940d,0xb08f,
    0x11,0x5993,0xe415,0x8897,0xd999,0x919b,0xa99d,0xd49f,
    0x21,0x2ca3,0xb8a5,0xa127,0x1a29,0x272b,0x36ad,0x3eaf,
    0x6631,0x1333,0x3ab5,0x7ab7,0xa5b9,0xc13b,0x6d3d,0x3cbf,
    0x41,0xcfc3,0xe45,0xf4c7,0xee49,0x304b,0x744d,0x9f4f,
    0xb6d1,0x20d3,0x5cd5,0xffd7,0x8159,0x6db,0xff5d,0x4ddf,
    0x78e1,0x9be3,0xce65,0x1867,0x8ce9,0x7e6b,0xae6d,0x69ef,
    0xa871,0x12f3,0xfaf5,0x7577,0x85f9,0x2efb,0x437d,0xdd7f};
unsigned column(unsigned p,unsigned bits) {
    auto result=unsigned(columns[p]);
    if(bits>=17)result|=((p&1U)*((p>>1)&1U))<<16;
    if(bits>=18)result|=((p&1U)*((p>>2)&1U))<<17;
    if(bits>=19)result|=((p&1U)*((p>>3)&1U))<<18;
    if(bits>=20)result|=((p&1U)*((p>>4)&1U))<<19;
    return result;
}
// Undo the validated s16 forward inner to recover the outer's routed word.
// This lets the extended oracle share only the already independently checked
// outer implementation; the new inner below is literal binary matrix algebra.
void undoBaseInner(Block* encoded,const rswide::Plan& base) {
    __m128i state[16]{},next[16],feedback[16];
    for(std::size_t epoch=0;epoch<base.n/64;++epoch) {
        for(auto& v:feedback)v=_mm_setzero_si128();
        for(unsigned p=0;p<64;++p) {
            auto raw=encoded[64*epoch+p].mData;
            for(unsigned mask=columns[p];mask;mask&=mask-1)raw=_mm_xor_si128(raw,state[std::countr_zero(mask)]);
            encoded[64*epoch+p]=Block(raw);
            for(unsigned mask=columns[p];mask;mask&=mask-1)feedback[std::countr_zero(mask)]=_mm_xor_si128(feedback[std::countr_zero(mask)],raw);
        }
        if(epoch+1==base.n/64)break;
        for(unsigned j=0;j<16;++j) {
            auto value=feedback[j];
            for(unsigned c=0;c<16;++c)if((base.reverseMatrices[epoch][c]>>j)&1)value=_mm_xor_si128(value,state[c]);
            next[j]=value;
        }
        std::memcpy(state,next,sizeof(state));
    }
}
}
void transposeScalar(const Block* input,Block* output,Block* scratch,const Plan& plan) {
    if(plan.stateBits==16){rswide::transposeScalar(input,output,scratch,plan.outer);return;}
    const auto bits=plan.stateBits;__m128i state[20]{},next[20]{},feedback[20];
    for(std::size_t epoch=plan.outer.n/64;epoch-->0;) {
        for(auto& value:feedback)value=_mm_setzero_si128();
        for(unsigned p=0;p<64;++p) {
            const auto i=64*epoch+p;const auto raw=input[i].mData;auto value=raw;
            for(unsigned mask=column(p,bits);mask;mask&=mask-1) {
                const auto j=std::countr_zero(mask);value=_mm_xor_si128(value,state[j]);
                feedback[j]=_mm_xor_si128(feedback[j],raw);
            }
            scratch[plan.outer.route[i/4]+(i&3)]=Block(value);
        }
        if(!epoch)break;
        for(unsigned j=0;j<bits;++j) {
            auto value=feedback[j];
            for(unsigned mask=plan.reverseMatrices[epoch][j];mask;mask&=mask-1)value=_mm_xor_si128(value,state[std::countr_zero(mask)]);
            next[j]=value;
        }
        std::memcpy(state,next,sizeof(state));
    }
    rswide::outerScalar(scratch,output,plan.outer);
}
void forwardScalar(const Block* message,Block* encoded,const Plan& plan) {
    rswide::forwardScalar(message,encoded,plan.outer);
    if(plan.stateBits==16)return;
    undoBaseInner(encoded,plan.outer);
    const auto bits=plan.stateBits;__m128i state[20]{},next[20]{},feedback[20];
    for(std::size_t epoch=0;epoch<plan.outer.n/64;++epoch) {
        for(auto& value:feedback)value=_mm_setzero_si128();
        for(unsigned p=0;p<64;++p) {
            const auto i=64*epoch+p;const auto raw=encoded[i].mData;auto value=raw;
            for(unsigned mask=column(p,bits);mask;mask&=mask-1) {
                const auto j=std::countr_zero(mask);value=_mm_xor_si128(value,state[j]);
                feedback[j]=_mm_xor_si128(feedback[j],raw);
            }
            encoded[i]=Block(value);
        }
        if(epoch+1==plan.outer.n/64)break;
        for(unsigned j=0;j<bits;++j) {
            auto value=feedback[j];
            for(unsigned c=0;c<bits;++c)if((plan.reverseMatrices[epoch][c]>>j)&1)value=_mm_xor_si128(value,state[c]);
            next[j]=value;
        }
        std::memcpy(state,next,sizeof(state));
    }
}
}
