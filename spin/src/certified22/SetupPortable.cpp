// Portable replica of the exact certified-s22 inner setup (Setup.cpp
// prepareInner), plain C++ for any 64-bit CPU. Same SplitMix word stream and
// acceptance, so every table is bit-identical (spin_turbo_setup_test).
// Candidate c is words 22c..22c+21; epoch e takes the e-th full-rank
// candidate. Eight candidates are rank-tested together by incremental
// reduction (pivot = lowest set bit), as plain lane loops the compiler can
// vectorize; the basis transforms are those of SetupTransform.h.
#include "Plan.h"
#include "SetupTransform.h"

namespace spin::detail::certified22 {
namespace {
using namespace transform;
constexpr unsigned lanes=8;
struct Batch {
    std::uint32_t rows[22][lanes];
    std::uint32_t alive[lanes];
};
void rankBatch(std::uint64_t seed,std::uint64_t first,Batch& out) {
    std::uint32_t basis[22][lanes],pivot[22][lanes];
    for(unsigned l=0;l<lanes;++l) out.alive[l]=1;
    for(unsigned r=0;r<22;++r) {
        std::uint32_t v[lanes];
        for(unsigned l=0;l<lanes;++l) {
            v[l]=std::uint32_t(mix(seed+step*(22*(first+l)+r+1)))&0x3fffffU;
            out.rows[r][l]=v[l];
        }
        for(unsigned i=0;i<r;++i)
            for(unsigned l=0;l<lanes;++l) v[l]^=basis[i][l]&(0U-std::uint32_t((v[l]&pivot[i][l])!=0));
        for(unsigned l=0;l<lanes;++l) {
            out.alive[l]&=std::uint32_t(v[l]!=0);
            basis[r][l]=v[l];
            pivot[r][l]=v[l]&(0U-v[l]);
        }
    }
}
}

void prepareInnerPortable(Plan& p,std::uint64_t seed) {
    const auto epochs=p.n/64;
    p.reverseMatrices.resize(epochs);
    p.updates.resize(9*epochs);
    p.forwardUpdates.resize(9*epochs);
    // Setup.cpp draws from Words(seed^0x3f625a92).
    const auto streamSeed=seed^0x3f625a92ULL;
    Batch batch;
    std::size_t epoch=0;
    for(std::uint64_t first=0;epoch<epochs;first+=lanes) {
        rankBatch(streamSeed,first,batch);
        for(unsigned lane=0;lane<lanes && epoch<epochs;++lane) {
            if(!batch.alive[lane]) continue;
            auto& literal=p.reverseMatrices[epoch];
            for(unsigned r=0;r<22;++r) literal[r]=batch.rows[r][lane];
            std::uint32_t u[24],g[24];
            transforms(literal.data(),u,g);
            pack(u,g,p.updates.data()+9*epoch,p.forwardUpdates.data()+9*epoch);
            ++epoch;
        }
    }
}
}
