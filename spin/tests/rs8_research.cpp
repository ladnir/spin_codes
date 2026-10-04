// Optional frozen research oracle. This file is not part of the installed
// library. It deliberately retains literal binary-matrix/Lagrange operations.
#include <spin/Code.h>
#include <spin/PacketCode.h>
#include "../experiments/rs8_sizes/Kernels.h"
#include "../experiments/k20_two_ms/rs8/Rs8.h"
#include "../src/certified22/productionMap.h"
#include <algorithm>
#include <array>
#include <bit>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>

namespace code=spin::detail::certified22;
namespace p15=spin::detail::paired15;
namespace sizes=spin::research::rs8sizes;
namespace frozen=spin::research::twoms::rs8;
using Block=code::Block;
namespace {
struct Case {std::size_t k;std::uint64_t route,inner;};
Case context{};
Block* blocks(spin::Buffer& b) {return reinterpret_cast<Block*>(b.bytes().data());}
void require(bool value,const char* message) {
    if(value)return;
    std::cerr<<"RS8 research K="<<context.k<<" route="<<context.route<<" inner="<<context.inner<<'\n';
    throw std::runtime_error(message);
}
void fill(spin::Buffer& buffer,std::uint64_t state) {
    for(std::size_t i=0;i<buffer.bytes().size();i+=8) {
        state+=0x9e3779b97f4a7c15ULL;auto word=state;
        word=(word^(word>>30))*0xbf58476d1ce4e5b9ULL;
        word=(word^(word>>27))*0x94d049bb133111ebULL;word^=word>>31;
        std::memcpy(buffer.bytes().data()+i,&word,8);
    }
}
void compare(std::span<const std::byte> expected,std::span<const std::byte> actual,const char* what) {
    require(expected.size()==actual.size(),"oracle result shape");
    if(!std::memcmp(expected.data(),actual.data(),expected.size()))return;
    std::size_t first=0;while(expected[first]==actual[first])++first;
    std::cerr<<"direction="<<what<<" first_byte="<<first<<" record="<<first/16<<'\n';
    require(false,"public RS8 differs from the frozen research map");
}
std::uint64_t checksum(std::span<const std::byte> bytes) {
    std::uint64_t result=0;
    for(std::size_t i=0;i<bytes.size();i+=8) {
        std::uint64_t word;std::memcpy(&word,bytes.data()+i,8);
        result=(result^word)*0x100000001b3ULL;
    }
    return result;
}

// The inner oracle uses literal columns, not the GFNI basis/slot layout.
constexpr std::uint16_t columns15[64]={
    0x1,0x3,0x5,0x87,0x9,0xb,0x400d,0x408f,0x11,0x113,0x815,0x997,0x419,0x51b,0x4c1d,0x4d9f,
    0x21,0x1023,0x25,0x10a7,0x829,0x182b,0x482d,0x58af,0x2031,0x3133,0x2835,0x39b7,0x2c39,0x3d3b,0x643d,0x75bf,
    0x41,0x443,0x2045,0x24c7,0x1049,0x144b,0x704d,0x74cf,0x4051,0x4553,0x6855,0x6dd7,0x5459,0x515b,0x3c5d,0x39df,
    0x261,0x1663,0x2265,0x36e7,0x1a69,0xe6b,0x7a6d,0x6eef,0x6271,0x7773,0x4a75,0x5ff7,0x7e79,0x6b7b,0x167d,0x3ff};

template<bool S15,class Plan>
static void forwardInnerScalar(const Block* input,Block* output,const Plan& p) {
    constexpr unsigned stateSize=S15?16:22;
    alignas(16) __m128i state[stateSize]{},next[stateSize],feedback[stateSize];
    for(std::size_t epoch=0;epoch<p.n/64;++epoch) {
        for(auto& value:feedback) value=_mm_setzero_si128();
        for(unsigned j=0;j<64;++j) {
            const auto i=64*epoch+j,routed=p.route[i/4]+(i&3);
            const auto raw=input[routed].mData;auto value=raw;
            const auto column=[&] {
                if constexpr(S15) return std::uint32_t(columns15[j]);
                else return code::map::oldColumns[j];
            }();
            for(auto mask=column;mask;mask&=mask-1) {
                const auto c=std::countr_zero(mask);
                value=_mm_xor_si128(value,state[c]);
                feedback[c]=_mm_xor_si128(feedback[c],raw);
            }
            output[i]=Block(value);
        }
        if(epoch+1==p.n/64) break;
        for(unsigned j=0;j<stateSize;++j) {
            auto value=feedback[j];
            for(unsigned c=0;c<stateSize;++c)
                if((p.reverseMatrices[epoch][c]>>j)&1U)
                    value=_mm_xor_si128(value,state[c]);
            next[j]=value;
        }
        std::memcpy(state,next,sizeof(state));
    }
}


void run(Case selected) {
    context=selected;const auto [k,routeSeed,innerSeed]=selected;
    code::Plan geometry(k,routeSeed,innerSeed);
    frozen::Tables32 tables;frozen::setupField32(geometry,routeSeed,tables);
    std::optional<p15::Plan> paired;
    if(k<=65536) {
        paired.emplace(65536,routeSeed,innerSeed);
        sizes::resizePaired15(*paired,geometry,innerSeed);
    }
    spin::Buffer message(k*16),input(2*k*16),scratch(geometry.scratchBlocks()*16);
    spin::Buffer forward(2*k*16),transpose(k*16),gotForward(2*k*16),gotTranspose(k*16);
    fill(message,1949);fill(input,913);
    frozen::outerForwardScalar(blocks(message),blocks(scratch),geometry,tables);
    if(paired)forwardInnerScalar<true>(blocks(scratch),blocks(forward),*paired);
    else forwardInnerScalar<false>(blocks(scratch),blocks(forward),geometry);

    // Calling the retained literal full transpose also computes its old outer;
    // ignore that output and apply the frozen RS8 outer to the routed scratch.
    std::fill(scratch.bytes().begin(),scratch.bytes().end(),std::byte{0x3c});
    if(paired)p15::transposeScalar(blocks(input),blocks(transpose),blocks(scratch),*paired);
    else code::transposeScalar(blocks(input),blocks(transpose),blocks(scratch),geometry);
    frozen::outerScalar(blocks(scratch),blocks(transpose),geometry,tables);

    // The frozen SIMD route and outer must agree with this literal reference.
    spin::Buffer fastScratch(geometry.scratchBlocks()*16);
    std::fill(fastScratch.bytes().begin(),fastScratch.bytes().end(),std::byte{0x3c});
    if(paired)sizes::reverseRoute15(blocks(input),blocks(fastScratch),*paired);
    else code::reverseRouteFast(blocks(input),blocks(fastScratch),geometry);
    compare(scratch.bytes(),fastScratch.bytes(),"frozen routed scratch and padding");
    frozen::outerFast(blocks(fastScratch),blocks(gotTranspose),geometry,tables);
    compare(transpose.bytes(),gotTranspose.bytes(),"frozen SIMD transpose");

    for(auto backend:{spin::Backend::Portable,spin::Backend::Automatic}) {
        spin::Code publicCode({k,spin::Parameters::PacketRs8T64,routeSeed,innerSeed},{backend});
        auto work=publicCode.make_workspace();
        publicCode.forward_bytes(message.bytes(),gotForward.bytes(),work);
        publicCode.transpose_bytes(input.bytes(),gotTranspose.bytes(),work);
        compare(forward.bytes(),gotForward.bytes(),"public forward");
        compare(transpose.bytes(),gotTranspose.bytes(),"public transpose");
        // Reusing output as input exercises the exact public in-place wrapper.
        std::memcpy(gotForward.bytes().data(),input.bytes().data(),input.bytes().size());
        publicCode.transpose_inplace_bytes(gotForward.bytes(),work);
        compare(transpose.bytes(),gotForward.bytes().first(k*16),"public in-place prefix");
        compare(input.bytes().subspan(k*16),gotForward.bytes().subspan(k*16),"public in-place suffix");
    }
    std::cout<<"RS8 frozen map match K="<<k<<" route="<<routeSeed<<" inner="<<innerSeed
        <<" forward=0x"<<std::hex<<checksum(forward.bytes())<<" transpose=0x"
        <<checksum(transpose.bytes())<<std::dec<<'\n';
}
}
int main() {try {
    // The optional oracle objects themselves are compiled for AVX512/GFNI.
    if(!spin::packet_fast_available())return 77;
    require(frozen::selfTest(),"frozen RS8 flattened field arithmetic failed");
    for(auto selected:std::array<Case,9>{{
        {2176,17,43},{6144,913,1123},{65536,17,43},{65536,113,113},
        {65664,17,43},{262144,17,43},{262144,113,113},
        {1048576,41,41},{1048576,113,113}}})run(selected);
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}

