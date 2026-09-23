#pragma once
#include <array>
#include <cstdint>
#include <limits>
#ifdef _MSC_VER
#include <intrin.h>
#endif

namespace spin::detail::kernel::setup {
// Preserve descriptor-v1's exact SplitMix stream, including rejected words.
// Independent arithmetic lanes allow the compiler to overlap the multiplications.
class Words {
    std::uint64_t seed_;
    std::array<std::uint64_t,8> words_;
    unsigned next_=8;
    void refill() noexcept {
        constexpr std::uint64_t step=0x9e3779b97f4a7c15ULL;
        for(unsigned i=0;i<8;++i) {
            auto v=seed_+step*(i+1);
            v=(v^(v>>30))*0xbf58476d1ce4e5b9ULL;
            v=(v^(v>>27))*0x94d049bb133111ebULL;
            words_[i]=v^(v>>31);
        }
        seed_+=8*step;next_=0;
    }
public:
    explicit Words(std::uint64_t seed):seed_(seed) {}
    std::uint64_t operator()() noexcept {
        if(next_==8)refill();
        return words_[next_++];
    }
};

// Bounds are in [2,2^32-1]. With m=floor(2^64/d), high(x*m) is
// floor(x/d) or one less. Thus one subtraction corrects the remainder.
// Precompute once per bound, not once per shuffle draw. Unlike multiply-high
// range sampling, this preserves the existing modulo sampler's exact outputs.
struct Divisor {
    std::uint64_t reciprocal=0,threshold=0;
    Divisor()=default;
    explicit Divisor(std::uint64_t d) {
        constexpr auto max=std::numeric_limits<std::uint64_t>::max();
        reciprocal=max/d;
        threshold=max%d+1;
        if(threshold==d){++reciprocal;threshold=0;}
    }
    std::uint64_t remainder(std::uint64_t x,std::uint64_t d) const noexcept {
#ifdef _MSC_VER
        const auto q=__umulh(x,reciprocal);
#else
        const auto q=std::uint64_t((__uint128_t(x)*reciprocal)>>64);
#endif
        const auto r=x-q*d;
        return r>=d?r-d:r;
    }
    template<class Source> std::uint64_t sample(Source& words,std::uint64_t d) const noexcept {
        std::uint64_t x;
        do{x=words();}while(x<threshold);
        return remainder(x,d);
    }
};
}
