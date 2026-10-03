#include "Packet8.h"
#include "../../src/kernels/SetupRandom.h"
#include <algorithm>
#include <numeric>

namespace spin::research::packet8 {
namespace {
using detail::kernel::setup::Words;
using detail::kernel::setup::Divisor;
template<class Container>
void shuffle(Container& values,Words& random,const std::vector<Divisor>& divisors) {
    std::iota(values.begin(),values.end(),0U);
    for(std::size_t i=values.size();i>1;--i)
        std::swap(values[i-1],values[divisors[i].sample(random,i)]);
}
}
Plan::Plan(std::size_t k,std::uint64_t seed)
    :outer(k,seed,seed,rs::Variant::Rs16Gf16,rs::InnerKernel::RetainedStreaming) {
    k16codesign::customizeNativeField16(outer,seed,native);
    Words routing(seed^0xa13df7412678bb08ULL);
    std::vector<Divisor> divisors(std::max(outer.groups,std::size_t(32))+1);
    for(std::size_t i=2;i<divisors.size();++i)divisors[i]=Divisor(i);
    std::vector<std::array<std::uint32_t,32>> columns(outer.groups);
    for(auto& p:columns)shuffle(p,routing,divisors);
    std::vector<std::uint32_t> positions(outer.groups);
    route.resize(n()/8);
    for(unsigned region=0;region<32;++region) {
        shuffle(positions,routing,divisors);
        for(std::size_t g=0;g<outer.groups;++g)
            route[region*outer.groups+positions[g]]=std::uint32_t(g*rs::groupStride+8*columns[g][region]);
    }
    Words inner(seed^0xb862dec09a47c103ULL);
    updates.resize(n()/64);
    for(auto& update:updates) {
        do {
            const auto value=inner();
            for(unsigned i=0;i<4;++i)update.forward[i]=std::uint8_t(value>>(8*i));
        } while((multiply(update.forward[0],update.forward[3])^
                 multiply(update.forward[1],update.forward[2]))==0);
        for(unsigned i=0;i<4;++i)update.adjoint[i]=adjointMatrix(update.forward[i]);
    }
}
}
