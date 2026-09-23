#pragma once
#include "Feistel.h"
#include "../../src/kernels/Inner.h"
#include "../../src/kernels/generated/BchCircuit.h"
#include <numeric>
#include <stdexcept>
#ifndef SPIN_FEISTEL_BATCH_ROUTE
#define SPIN_FEISTEL_BATCH_ROUTE 1
#endif
namespace spin::detail::kernel {void bchTranspose4(const block*,block*);}
namespace spin::experimental::feistel {
namespace kernel=detail::kernel;
inline std::vector<std::uint32_t> makeMasksK18(std::uint64_t seed) {
    std::vector<std::uint32_t> masks(2*((1U<<19)/128));
    kernelWords words(seed);
    for(unsigned i=0;i<masks.size();i+=2) {
        unsigned u;do{u=unsigned(words())&((1U<<19)-1);}while(!u);
        unsigned v=unsigned(words())&((1U<<19)-1);
        if(std::popcount(u&v)&1)v^=u&-u;
        masks[i]=u;masks[i+1]=v;
    }
    return masks;
}
// K=2^18 experiment. Keep the exact-uniform generator as a lean control, so
// gains from removing execution tables are not attributed to Feistel mixing.
inline std::vector<std::uint32_t> uniformK18(std::uint64_t seed) {
    constexpr unsigned rows=2048,n=rows*256;
    kernelWords words(seed);
    std::vector<kernel::setup::Divisor> divisors(rows+1);
    for(unsigned d=2;d<=rows;++d)divisors[d]=kernel::setup::Divisor(d);
    auto shuffle=[&](auto& p) {
        std::iota(p.begin(),p.end(),0);
        for(unsigned i=unsigned(p.size());i>1;--i)std::swap(p[i-1],p[divisors[i].sample(words,i)]);
    };
    std::vector<std::uint8_t> coords(n);std::array<std::uint8_t,256> row;
    for(unsigned j=0;j<rows;++j){shuffle(row);for(unsigned c=0;c<256;++c)coords[row[c]*rows+j]=std::uint8_t(c);}
    std::vector<std::uint32_t> route(n),positions(rows);
    for(unsigned r=0;r<256;++r){shuffle(positions);for(unsigned j=0;j<rows;++j)route[r*rows+positions[j]]=256*j+coords[r*rows+j];}
    return route;
}
template<unsigned Rounds,bool Materialized> class Route;
template<unsigned Rounds> class Route<Rounds,false> {
    Routing<11,Rounds> routing_;
public:
    explicit Route(std::uint64_t seed):routing_(seed) {}
    unsigned outer(unsigned i) const {return routing_.outer(i);}
    void outerBatch(unsigned first,unsigned* out) const {routing_.outerBatch(first,out);}
    std::size_t bytes() const{return routing_.bytes();}
};
template<unsigned Rounds> class Route<Rounds,true> {
    std::vector<std::uint32_t> route_;
public:
    explicit Route(std::uint64_t seed):route_(Routing<11,Rounds>(seed).materialize()) {}
    unsigned outer(unsigned i) const{return route_[i];}
    std::size_t bytes() const{return 4*route_.size();}
};
template<> class Route<0,true> {
    std::vector<std::uint32_t> route_;
public:
    explicit Route(std::uint64_t seed):route_(uniformK18(seed)) {}
    unsigned outer(unsigned i) const{return route_[i];}
    std::size_t bytes() const{return 4*route_.size();}
};
template<unsigned Rounds,bool Materialized,class RouteType=Route<Rounds,Materialized>> class Direct {
    RouteType route_;
    std::vector<std::uint32_t> masks_;
    bool four_;
    template<bool Four> void run(kernel::block* buffer,kernel::block* scratch) const {
        std::array<unsigned,16> routeBatch;
        kernel::innerReverse<kernel::Map128S19>(buffer,N,masks_.data(),[&](std::size_t i,kernel::block v) {
            unsigned x;
            if constexpr(!Materialized && SPIN_FEISTEL_BATCH_ROUTE) {
                // The S19 reverse emitter visits coordinates in decreasing order.
                if((i&15)==15)route_.outerBatch(unsigned(i)&~15U,routeBatch.data());
                x=routeBatch[i&15];
            }else x=route_.outer(unsigned(i));
            const auto at=Four?((x&~1023U)|((x&255U)<<2)|((x>>8)&3U)):x;
            scratch[at]=v;
        });
        if constexpr(Four) {
            for(unsigned i=0;i<N;i+=1024)kernel::bchTranspose4(scratch+i,buffer+i/2);
        }else {
            for(unsigned i=0;i<N;i+=512)kernel::bchTranspose2(scratch+i,scratch+i+256,buffer+i/2,buffer+i/2+128);
        }
    }
public:
    static constexpr unsigned N=1U<<19;
    Direct(std::uint64_t routeSeed,std::uint64_t innerSeed):Direct(RouteType(routeSeed),innerSeed) {}
    Direct(RouteType route,std::uint64_t innerSeed):route_(std::move(route)),masks_(makeMasksK18(innerSeed)),four_(kernel::bchAvx512Available()) {}
    std::size_t bytes() const{return route_.bytes()+4*masks_.size();}
    void encode(kernel::block* buffer,std::vector<kernel::block>& scratch) const {
        if(scratch.size()!=N)throw std::invalid_argument("experimental scratch size");
#if SPIN_BCH_AVX512
        if(four_){run<true>(buffer,scratch.data());return;}
#endif
        run<false>(buffer,scratch.data());
    }
};
}
