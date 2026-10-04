// Optional validation against the frozen literal research map. This is not a
// benchmark and is never installed with the library.
#include <spin/Code.h>
#include <spin/PacketCode.h>
#include "../src/certified22/Plan.h"
#include "../experiments/k18_certified22/K18Certified22.h"
#include "../experiments/k16_codesign_100us/kernel/K16NativeOuter.h"
#include <algorithm>
#include <array>
#include <cstring>
#include <iomanip>
#include <iostream>

namespace rs=spin::research::rs;
namespace outer=spin::research::k16codesign;
namespace oracle=spin::research::k18cert22;
namespace production=spin::detail::certified22;
namespace spin::research::k16codesign {
void outerScalar(const rs::Block*,rs::Block*,const rs::Plan&);
void outerForwardScalar(const rs::Block*,rs::Block*,const rs::Plan&);
}
namespace {
static_assert(std::is_same_v<rs::Block,production::Block>);
struct Case {std::size_t k;std::uint64_t route,inner;};
Case context{};
rs::Block* blocks(spin::Buffer& buffer) {return reinterpret_cast<rs::Block*>(buffer.bytes().data());}
void require(bool condition,const char* message) {
    if(condition)return;
    std::cerr<<"K="<<context.k<<" route="<<context.route<<" inner="<<context.inner<<'\n';
    throw std::runtime_error(message);
}
void fill(spin::Buffer& buffer,std::uint64_t seed) {
    for(std::size_t i=0;i<buffer.bytes().size();i+=8) {
        seed+=0x9e3779b97f4a7c15ULL;auto z=seed;
        z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;z=(z^(z>>27))*0x94d049bb133111ebULL;z^=z>>31;
        std::memcpy(buffer.bytes().data()+i,&z,8);
    }
}
std::uint64_t hash(const spin::Buffer& buffer) {
    std::uint64_t result=0;
    for(std::size_t i=0;i<buffer.bytes().size();i+=8) {
        std::uint64_t word;std::memcpy(&word,buffer.bytes().data()+i,8);
        result=(result^word)*0x100000001b3ULL;
    }
    return result;
}
void compare(const spin::Buffer& expected,const spin::Buffer& actual,const char* direction) {
    require(expected.bytes().size()==actual.bytes().size(),"research comparison shape mismatch");
    if(!std::memcmp(expected.bytes().data(),actual.bytes().data(),expected.bytes().size()))return;
    std::size_t first=0;while(expected.bytes()[first]==actual.bytes()[first])++first;
    std::cerr<<"direction="<<direction<<" first_byte="<<first<<" record="<<first/16
        <<" expected="<<std::to_integer<unsigned>(expected.bytes()[first])
        <<" actual="<<std::to_integer<unsigned>(actual.bytes()[first])<<'\n';
    require(false,"promoted certified22 map differs from frozen literal research oracle");
}
void run(Case selected) {
    context=selected;const auto [k,routeSeed,innerSeed]=selected;
    rs::Plan plan(k,routeSeed,innerSeed,rs::Variant::Rs16Gf16,rs::InnerKernel::RetainedStreaming);
    oracle::Tables tables;oracle::prepare(plan,innerSeed,tables);
    outer::NativeOuterTables native;outer::customizeNativeField16(plan,routeSeed,native);
    production::Plan promoted(k,routeSeed,innerSeed);
    require(promoted.k==plan.k && promoted.n==plan.n && promoted.groups==plan.groups,"setup geometry changed");
    require(promoted.scratchBlocks()==plan.scratchBlocks(),"scratch layout changed");
    require(promoted.route==plan.route,"sampled route changed");
    require(promoted.routeBytes==tables.routeBytes,"route byte offsets changed");
    require(promoted.reverseMatrices==tables.reverseMatrices,"literal GL22 distribution/seed changed");
    require(promoted.updates==tables.updates,"basis-conjugated GL22 tables changed");
    require(promoted.outerMatrices16==plan.outerMatrices16,"literal native outer matrices changed");
    require(promoted.outerField==native.field,"native outer field coefficients changed");
    spin::Buffer message(k*16),input(2*k*16),scratch(plan.scratchBlocks()*16),fastScratch(plan.scratchBlocks()*16);
    spin::Buffer forward(2*k*16),transpose(k*16),gotForward(2*k*16),gotTranspose(k*16);
    fill(message,1949);fill(input,913);
    outer::outerForwardScalar(blocks(message),blocks(scratch),plan);
    oracle::forwardInnerScalar(blocks(scratch),blocks(forward),plan,tables);
    // Check the routed scratch, including all unused padding: the outer has a
    // kernel and can otherwise hide an incorrect inner or route implementation.
    std::fill(scratch.bytes().begin(),scratch.bytes().end(),std::byte{0x3c});
    std::fill(fastScratch.bytes().begin(),fastScratch.bytes().end(),std::byte{0x3c});
    oracle::reverseRouteScalar(blocks(input),blocks(scratch),plan,tables);
    production::reverseRouteFast(blocks(input),blocks(fastScratch),promoted);
    compare(scratch,fastScratch,"routed scratch and padding");
    outer::outerScalar(blocks(scratch),blocks(transpose),plan);
    for(auto backend:{spin::Backend::Portable,spin::Backend::Automatic}) {
        spin::Code code({k,spin::Parameters::PacketRsT64S22,routeSeed,innerSeed},{backend});auto work=code.make_workspace();
        code.forward_bytes(message.bytes(),gotForward.bytes(),work);
        code.transpose_bytes(input.bytes(),gotTranspose.bytes(),work);
        compare(forward,gotForward,"public forward");compare(transpose,gotTranspose,"public transpose");
    }
    std::cout<<"certified22 research match K="<<k<<" route="<<routeSeed<<" inner="<<innerSeed
        <<" forward=0x"<<std::hex<<hash(forward)<<" transpose=0x"<<hash(transpose)<<std::dec<<'\n';
}
}
int main() {try {
    // Research objects are compiled for AVX512/GFNI, unlike the installed fallback.
    if(!spin::packet_fast_available() || !spin::capabilities().forward512)return 77;
    for(auto selected:std::array<Case,7>{{{2048,1,1},{2048,17,43},{2048,913,1123},
                                         {262144,17,43},{262144,113,113},
                                         {1048576,41,41},{1048576,113,113}}})run(selected);
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
