// Optional integration test against the independently frozen research sources.
// Keep this outside the installed library and do not copy library setup into it.
#include <spin/Code.h>
#include <spin/PacketCode.h>
#include "../experiments/k16_codesign_100us/kernel/K16Paired15Shuffle.h"
#include "../experiments/k16_codesign_100us/kernel/K16NativeOuter.h"
#include <array>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <utility>

namespace rs=spin::research::rs;
namespace oracle=spin::research::k16codesign;
namespace spin::research::k16codesign {
void outerScalar(const rs::Block*,rs::Block*,const rs::Plan&);
void outerForwardScalar(const rs::Block*,rs::Block*,const rs::Plan&);
}
static rs::Block* blocks(spin::Buffer& buffer) {return reinterpret_cast<rs::Block*>(buffer.bytes().data());}
static void fill(spin::Buffer& buffer,std::uint64_t seed) {
    for(std::size_t i=0;i<buffer.bytes().size();i+=8) {
        seed+=0x9e3779b97f4a7c15ULL;auto z=seed;
        z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;z=(z^(z>>27))*0x94d049bb133111ebULL;z^=z>>31;
        std::memcpy(buffer.bytes().data()+i,&z,8);
    }
}
static std::uint64_t hash(const spin::Buffer& buffer) {
    std::uint64_t result=0;
    for(std::size_t i=0;i<buffer.bytes().size();i+=8) {
        std::uint64_t word;std::memcpy(&word,buffer.bytes().data()+i,8);
        result=(result^word)*0x100000001b3ULL;
    }
    return result;
}
static void compare(const spin::Buffer& expected,const spin::Buffer& actual,const char* direction,
                    spin::Backend backend,std::uint64_t route,std::uint64_t inner) {
    if(!std::memcmp(expected.bytes().data(),actual.bytes().data(),expected.bytes().size()))return;
    std::size_t first=0;while(expected.bytes()[first]==actual.bytes()[first])++first;
    std::cerr<<"frozen research mismatch: direction="<<direction
        <<" backend="<<(backend==spin::Backend::Portable?"Portable":"Avx512")
        <<" route="<<route<<" inner="<<inner<<" width_bits=128 first_byte="<<first
        <<" record="<<first/16<<" expected="<<std::to_integer<unsigned>(expected.bytes()[first])
        <<" actual="<<std::to_integer<unsigned>(actual.bytes()[first])<<'\n';
    throw std::runtime_error("paired15 public encoding differs from frozen mode52 research map");
}
int main() {try {
    // Original research objects contain unrestricted AVX512 intrinsics.
    if(!spin::packet_fast_available() || !spin::capabilities().forward512)return 77;
    constexpr std::size_t k=65536;
    for(auto seeds:std::array<std::pair<std::uint64_t,std::uint64_t>,3>{{{1,1},{17,43},{913,1123}}}) {
        const auto [routeSeed,innerSeed]=seeds;
        rs::Plan plan(k,routeSeed,innerSeed,rs::Variant::Rs16Gf16,rs::InnerKernel::RetainedStreaming);
        oracle::PairedTables tables;oracle::PairedOptimizedTables optimized;oracle::NativeOuterTables native;
        oracle::customizePaired15Shuffle(plan,innerSeed,tables,optimized);
        oracle::customizeNativeField16(plan,routeSeed,native);
        spin::Buffer message(k*16),input(2*k*16),scratch(plan.scratchBlocks()*16);
        spin::Buffer forward(2*k*16),transpose(k*16),gotForward(2*k*16),gotTranspose(k*16);
        fill(message,1949);fill(input,913);
        oracle::outerForwardScalar(blocks(message),blocks(scratch),plan);
        oracle::forwardInnerPaired15Scalar(blocks(scratch),blocks(forward),plan);
        oracle::reverseRoutePaired15Scalar(blocks(input),blocks(scratch),plan);
        oracle::outerScalar(blocks(scratch),blocks(transpose),plan);
        for(auto backend:{spin::Backend::Portable,spin::Backend::Automatic}) {
            spin::Code code({k,spin::Parameters::PacketRsT64S15K16,routeSeed,innerSeed},{backend});
            auto work=code.make_workspace();
            code.forward_bytes(message.bytes(),gotForward.bytes(),work);
            code.transpose_bytes(input.bytes(),gotTranspose.bytes(),work);
            compare(forward,gotForward,"forward",code.backend(),routeSeed,innerSeed);
            compare(transpose,gotTranspose,"transpose",code.backend(),routeSeed,innerSeed);
        }
        std::cout<<"paired15 research match route="<<routeSeed<<" inner="<<innerSeed
            <<" forward=0x"<<std::hex<<hash(forward)<<" transpose=0x"<<hash(transpose)<<std::dec<<'\n';
    }
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
