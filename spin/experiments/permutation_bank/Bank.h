#pragma once
#include "../../src/kernels/SetupRandom.h"
#include <algorithm>
#include <bit>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace spin::experimental::bank {
using Words=detail::kernel::setup::Words;
enum class Mode { Xor, XorAdd, RotateAdd, InputRotateAdd };
inline unsigned rotate(unsigned x,unsigned r,unsigned bits) {
    const unsigned mask=(1U<<bits)-1;
    return r?((x<<r)|(x>>(bits-r)))&mask:x;
}
struct Parameters {unsigned entry,a,b,c,e,ri,ro;};
struct Bank {
    unsigned bits,n,mask,count;
    std::vector<unsigned> p,inverse;
    Bank(unsigned d,unsigned size,std::uint64_t seed):bits(d),n(1U<<d),mask(n-1),count(size),p(n*size),inverse(n*size) {
        if(d<5 || d>16 || !std::has_single_bit(size))throw std::invalid_argument("bank geometry");
        Words words(seed);
        std::vector<detail::kernel::setup::Divisor> div(n+1);
        for(unsigned i=2;i<=n;++i)div[i]=detail::kernel::setup::Divisor(i);
        for(unsigned j=0;j<size;++j) {
            auto first=p.begin()+j*n;std::iota(first,first+n,0);
            for(unsigned i=n;i>1;--i)std::swap(first[i-1],first[div[i].sample(words,i)]);
            for(unsigned x=0;x<n;++x)inverse[j*n+first[x]]=x;
        }
    }
    Parameters sample(Words& words) const {
        detail::kernel::setup::Divisor rd(bits);
        return {unsigned(words())&(count-1),unsigned(words())&mask,unsigned(words())&mask,
            unsigned(words())&mask,unsigned(words())&mask,unsigned(rd.sample(words,bits)),unsigned(rd.sample(words,bits))};
    }
    template<Mode M> unsigned in(unsigned x,const Parameters& k) const {
        x^=k.a;
        if constexpr(M==Mode::RotateAdd || M==Mode::InputRotateAdd)x=rotate(x,k.ri,bits);
        if constexpr(M!=Mode::Xor)x=(x+k.b)&mask;
        return x;
    }
    template<Mode M> unsigned out(unsigned x,const Parameters& k) const {
        if constexpr(M==Mode::InputRotateAdd)return x;
        x^=k.c;
        if constexpr(M==Mode::RotateAdd)x=rotate(x,k.ro,bits);
        if constexpr(M!=Mode::Xor)x=(x+k.e)&mask;
        return x;
    }
    template<Mode M> unsigned forward(unsigned x,const Parameters& k) const {return out<M>(p[k.entry*n+in<M>(x,k)],k);}
    template<Mode M> unsigned backward(unsigned x,const Parameters& k) const {
        if constexpr(M!=Mode::InputRotateAdd) {
            if constexpr(M!=Mode::Xor)x=(x-k.e)&mask;
            if constexpr(M==Mode::RotateAdd)x=rotate(x,k.ro?bits-k.ro:0,bits);
            x^=k.c;
        }
        x=inverse[k.entry*n+x];
        if constexpr(M!=Mode::Xor)x=(x-k.b)&mask;
        if constexpr(M==Mode::RotateAdd || M==Mode::InputRotateAdd)x=rotate(x,k.ri?bits-k.ri:0,bits);
        return x^k.a;
    }
};
}
