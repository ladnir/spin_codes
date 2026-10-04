#include "K20Variants.h"
#include <spin/Code.h>
#include <spin/PacketCode.h>
#include <algorithm>
#include <bit>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace rs = spin::research::rs;
namespace inner = spin::research::k18cert22;
namespace outer = spin::research::k16codesign;
namespace tuned = spin::research::k20small;
namespace spin::research::k16codesign {
void outerScalar(const rs::Block*, rs::Block*, const rs::Plan&);
void outerForwardScalar(const rs::Block*, rs::Block*, const rs::Plan&);
}
using Block = rs::Block;
static Block* blocks(spin::Buffer& b) { return reinterpret_cast<Block*>(b.bytes().data()); }
static void fill(spin::Buffer& b, std::uint64_t state) {
    for (std::size_t i = 0; i < b.bytes().size(); i += 8) {
        state += 0x9e3779b97f4a7c15ULL; auto v = state;
        v = (v ^ (v >> 30)) * 0xbf58476d1ce4e5b9ULL;
        v = (v ^ (v >> 27)) * 0x94d049bb133111ebULL; v ^= v >> 31;
        std::memcpy(b.bytes().data()+i, &v, 8);
    }
}
static unsigned dot(const Block* a, const Block* b, std::size_t n) {
    unsigned p = 0;
    for (std::size_t i = 0; i < n; ++i) {
        std::uint64_t x[2], y[2]; std::memcpy(x,a+i,16); std::memcpy(y,b+i,16);
        p ^= std::popcount(x[0]&y[0]) ^ std::popcount(x[1]&y[1]);
    }
    return p&1;
}
static void require(bool value, const char* why) {
    if (!value) throw std::runtime_error(why);
}
template<unsigned RouteNt, unsigned OuterNt>
static void checkVariant(const rs::Plan& plan, const inner::Tables& tables,
                         const outer::NativeOuterTables& native, spin::Buffer& input,
                         spin::Buffer& scratch, spin::Buffer& output,
                         spin::Buffer& scalarScratch, spin::Buffer& expected,
                         const std::vector<std::byte>& saved) {
    const auto encode = [&](const Block* x, Block* y) {
        tuned::route<RouteNt>(x,blocks(scratch),plan,tables);
        tuned::finish<OuterNt>(blocks(scratch),y,plan,native);
    };
    std::memcpy(input.bytes().data(),saved.data(),saved.size());
    std::fill(scratch.bytes().begin(),scratch.bytes().end(),std::byte{0x3c});
    tuned::route<RouteNt>(blocks(input)+4,blocks(scratch),plan,tables);
    require(!std::memcmp(scratch.bytes().data(),scalarScratch.bytes().data(),scratch.bytes().size()),
            "routed scalar mismatch or padding overwrite");
    for (unsigned offset = 0; offset < 4; ++offset) {
        std::fill(output.bytes().begin(),output.bytes().end(),std::byte{0xa5});
        std::memcpy(blocks(input)+4+offset,saved.data()+64,plan.n*16);
        auto* out = blocks(output)+4+offset;
        encode(blocks(input)+4+offset,out);
        require(!std::memcmp(out,blocks(expected),plan.k*16),"scalar output mismatch");
        const auto guard = [](std::byte x) { return x == std::byte{0xa5}; };
        require(std::all_of(output.bytes().begin(),output.bytes().begin()+(4+offset)*16,guard)
            && std::all_of(output.bytes().begin()+(4+offset+plan.k)*16,output.bytes().end(),guard),
            "output guard mismatch");
        require(!std::memcmp(blocks(input)+4+offset,saved.data()+64,plan.n*16),"input modified");
        encode(blocks(input)+4+offset,blocks(input)+4+offset);
        require(!std::memcmp(blocks(input)+4+offset,blocks(expected),plan.k*16),"in-place output mismatch");
        require(!std::memcmp(blocks(input)+4+offset+plan.k,saved.data()+64+plan.k*16,
                            (plan.n-plan.k)*16),"in-place tail modified");
    }
}
static void check(std::size_t k, std::uint64_t seed) {
    rs::Plan plan(k,seed,seed,rs::Variant::Rs16Gf16,rs::InnerKernel::RetainedStreaming);
    inner::Tables tables; inner::prepare(plan,seed,tables);
    outer::NativeOuterTables native; outer::customizeNativeField16(plan,seed,native);
    spin::Buffer input((plan.n+12)*16), scratch(plan.scratchBlocks()*16),
        output((plan.k+12)*16), scalarScratch(plan.scratchBlocks()*16), expected(plan.k*16),
        message(plan.k*16), encoded(plan.n*16), forwardScratch(plan.scratchBlocks()*16);
    fill(input,913); fill(message,1949);
    const std::vector<std::byte> saved(input.bytes().begin(),input.bytes().end());
    std::fill(scalarScratch.bytes().begin(),scalarScratch.bytes().end(),std::byte{0x3c});
    inner::reverseRouteScalar(blocks(input)+4,blocks(scalarScratch),plan,tables);
    outer::outerScalar(blocks(scalarScratch),blocks(expected),plan);
    outer::outerForwardScalar(blocks(message),blocks(forwardScratch),plan);
    inner::forwardInnerScalar(blocks(forwardScratch),blocks(encoded),plan,tables);
    require(dot(blocks(encoded),blocks(input)+4,plan.n)==dot(blocks(message),blocks(expected),plan.k),
            "literal-map adjoint mismatch");
    checkVariant<false,false>(plan,tables,native,input,scratch,output,scalarScratch,expected,saved);
    checkVariant<true,false>(plan,tables,native,input,scratch,output,scalarScratch,expected,saved);
    checkVariant<false,true>(plan,tables,native,input,scratch,output,scalarScratch,expected,saved);
    checkVariant<true,true>(plan,tables,native,input,scratch,output,scalarScratch,expected,saved);
    checkVariant<1,2>(plan,tables,native,input,scratch,output,scalarScratch,expected,saved);
    checkVariant<1,3>(plan,tables,native,input,scratch,output,scalarScratch,expected,saved);
    checkVariant<2,1>(plan,tables,native,input,scratch,output,scalarScratch,expected,saved);
    checkVariant<2,3>(plan,tables,native,input,scratch,output,scalarScratch,expected,saved);
    checkVariant<3,0>(plan,tables,native,input,scratch,output,scalarScratch,expected,saved);
    checkVariant<3,4>(plan,tables,native,input,scratch,output,scalarScratch,expected,saved);
    checkVariant<0,4>(plan,tables,native,input,scratch,output,scalarScratch,expected,saved);
    std::cout << "PASS eleven variants, scalar route/output, guards, alignments, adjoint, in-place K="
        << k << " seed=" << seed << '\n';
}
int main(int argc, char** argv) { try {
    if (!spin::packet_fast_available() || !spin::capabilities().forward512) return 77;
    if (argc != 3) throw std::invalid_argument("K seed");
    require(tuned::packedCarrySelfTest(), "packed-carry basis identity failed");
    check(std::stoull(argv[1]),std::stoull(argv[2]));
    return 0;
} catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; } }
