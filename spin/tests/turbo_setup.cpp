// Turbo (PacketRs8T64) setup differential test: the portable, AVX2 and AVX-512 setups
// must reproduce the scalar setup's tables bit for bit, and the table-driven
// outer coefficients must equal the literal bit-serial definitions.
#include "../src/rs8/Rs8.h"
#include "../src/rs8/Field.h"
#include "../src/Cpu.h"
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace spin::detail;
template<class A,class B> void same(const A& a,const B& b,const char* what,std::size_t k) {
    if(a.size()!=b.size() || !std::equal(a.begin(),a.end(),b.begin())) {
        std::cerr<<"setup mismatch: "<<what<<" K="<<k<<'\n';
        throw std::runtime_error(what);
    }
}
void check(std::size_t k,std::uint64_t route,std::uint64_t inner,certified22::SetupPath path) {
    const rs8::Plan fast(k,route,inner,path),scalar(k,route,inner,certified22::SetupPath::Scalar);
    const auto& f=*fast.certified22Plan;
    const auto& s=*scalar.certified22Plan;
    same(f.route,s.route,"route",k);
    same(f.routeBytes,s.routeBytes,"routeBytes",k);
    same(f.reverseMatrices,s.reverseMatrices,"reverseMatrices",k);
    same(f.updates,s.updates,"updates",k);
    same(f.forwardUpdates,s.forwardUpdates,"forwardUpdates",k);
    same(fast.scalars,scalar.scalars,"scalars",k);
    same(fast.field,scalar.field,"field",k);
    same(fast.adjoint,scalar.adjoint,"adjoint",k);
    for(std::size_t symbol=0; symbol<fast.scalars.size(); ++symbol) {
        const auto x=fast.scalars[symbol];
        if(!x || fast.field[symbol]!=rs8::field::coefficients(x) || fast.adjoint[symbol]!=rs8::literalAdjoint(x)) {
            std::cerr<<"outer literal mismatch K="<<k<<" symbol="<<symbol<<'\n';
            throw std::runtime_error("outer literal");
        }
    }
}
}
int main() {
    try {
        // Every fast setup available here against the literal Scalar setup:
        // Portable always, Avx2 where present, Automatic (AVX-512) where present.
        std::vector<certified22::SetupPath> paths{certified22::SetupPath::Portable};
        std::string names="portable";
        if(cpu_avx2()) { paths.push_back(certified22::SetupPath::Avx2); names+=", AVX2"; }
        if(cpu_packet512()) { paths.push_back(certified22::SetupPath::Automatic); names+=", AVX-512"; }
        for(auto path:paths)
            for(std::size_t k:{std::size_t{128},std::size_t{384},std::size_t{1}<<14,std::size_t{65664},std::size_t{1}<<18})
                for(std::uint64_t seed:{1ULL,0x9e3779b97f4a7c15ULL,~0ULL}) check(k,seed,seed*3+7,path);
        std::cout<<"setup: "<<names<<" and scalar Turbo setup tables identical\n";
        return 0;
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
