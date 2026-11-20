#include "Plan.h"
#include "productionMap.h"
#include "../kernels/SetupRandom.h"
#include "../Cpu.h"
#include <algorithm>
#include <bit>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace spin::detail::certified22 {
namespace {
using kernel::setup::Words;
using kernel::setup::Divisor;
using Matrix = std::array<std::uint32_t,22>;

bool fullRank(Matrix rows) {
    for(unsigned column=0; column<22; ++column) {
        unsigned pivot=column;
        while(pivot<22 && !((rows[pivot]>>column)&1U)) ++pivot;
        if(pivot==22) return false;
        std::swap(rows[pivot],rows[column]);
        for(unsigned row=column+1; row<22; ++row)
            if((rows[row]>>column)&1U) rows[row]^=rows[column];
    }
    return true;
}
Matrix multiply(const Matrix& left,const Matrix& right) {
    Matrix result{};
    for(unsigned row=0; row<22; ++row)
        for(auto mask=left[row]; mask; mask&=mask-1)
            result[row]^=right[std::countr_zero(mask)];
    return result;
}

#if SPIN_AVX512
// The same SplitMix stream as Words, generated eight words per vector.
class BufferedWords {
    std::uint64_t seed_, next_ = 0;
    unsigned position_ = 256;
    alignas(64) std::uint64_t words_[256];
public:
    explicit BufferedWords(std::uint64_t seed) : seed_(seed) {}
    std::uint64_t operator()() noexcept {
        if(position_ == 256) { fillWords(seed_, next_, words_, 256); next_ += 256; position_ = 0; }
        return words_[position_++];
    }
};
#endif
template<class Container,class Source>
void shuffle(Container& values, Source& words, const std::vector<Divisor>& divisors) {
    std::iota(values.begin(), values.end(), 0U);
    for(std::size_t i = values.size(); i > 1; --i)
        std::swap(values[i - 1], values[divisors[i].sample(words, i)]);
}
template<class Source>
void prepareRoute(Plan& p, Source& words) {
    std::vector<Divisor> divisors(std::max(p.groups, std::size_t{64}) + 1);
    for(std::size_t i = 2; i < divisors.size(); ++i) divisors[i] = Divisor(i);
    std::vector<std::array<std::uint32_t,64>> columns(p.groups);
    for(auto& column : columns) shuffle(column, words, divisors);
    std::vector<std::uint32_t> destinations(p.groups);
    p.route.resize(p.n / 4);
    p.routeBytes.resize(p.n / 4);
    for(unsigned region = 0; region < 64; ++region) {
        shuffle(destinations, words, divisors);
        for(std::size_t group = 0; group < p.groups; ++group) {
            const auto i = region * p.groups + destinations[group];
            const auto offset = std::uint32_t(group * groupStride + 4 * columns[group][region]);
            p.route[i] = offset;
            p.routeBytes[i] = 16 * offset;
        }
    }
}

void prepareInner(Plan& p,std::uint64_t seed) {
    const auto epochs=p.n/64;
    p.reverseMatrices.resize(epochs);
    p.updates.resize(9*epochs);
    p.forwardUpdates.resize(9*epochs);
    Words words(seed^0x3f625a92ULL);
    const auto packBlocks=[](const Matrix& transformed,std::uint64_t* out) {
        for(unsigned o=0; o<3; ++o)
            for(unsigned in=0; in<3; ++in) {
                std::uint64_t packed=0;
                for(unsigned j=0; j<8; ++j) {
                    const auto row=8*o+j;
                    if(row<22) packed|=std::uint64_t((transformed[row]>>(8*in))&255U)<<(8*j);
                }
                out[3*o+in]=packed;
            }
    };
    for(std::size_t epoch=0; epoch<epochs; ++epoch) {
        auto& literal=p.reverseMatrices[epoch];
        // Exact uniform-GL22 rejection sampler from the certified setup.
        do {
            for(auto& row:literal) row=static_cast<std::uint32_t>(words())&0x3fffffU;
        } while(!fullRank(literal));
        packBlocks(multiply(multiply(map::basis,literal),map::inverseBasis),p.updates.data()+9*epoch);
        // The forward recurrence uses R^T in the same circuit basis.
        Matrix transposed{};
        for(unsigned row=0; row<22; ++row)
            for(unsigned column=0; column<22; ++column)
                transposed[column]|=((literal[row]>>column)&1U)<<row;
        packBlocks(multiply(multiply(map::basis,transposed),map::inverseBasis),p.forwardUpdates.data()+9*epoch);
    }
}
}
bool validMessageSize(std::size_t k) noexcept {
    if(!k || k%128) return false;
    const auto groups=k/128;
    // Each route points to a complete 64-byte packet. Also reserve workspace
    // for the widest supported records without overflowing size_t.
    constexpr auto maxByteRoute=std::numeric_limits<std::uint32_t>::max();
    constexpr auto maxSize=std::numeric_limits<std::size_t>::max();
    return groups<=maxByteRoute/(groupStride*sizeof(Block))
        && groups<=maxSize/(groupStride*sizeof(Block)*4)
        && k<=maxSize/(2*sizeof(Block)*4);
}
Plan::Plan(std::size_t size,std::uint64_t routeSeed,std::uint64_t innerSeed,SetupPath path) {
    if(!validMessageSize(size))
        throw std::invalid_argument("Turbo s22 K must be a positive multiple of 128 within the byte-route and workspace limits");
    k=size;
    n=2*size;
    groups=size/128;
#if SPIN_AVX512
    const bool packet512=path==SetupPath::Automatic && cpu_packet512();
#else
    const bool packet512=false;
#endif
    const bool avx2=SPIN_X86 && !packet512 && (path==SetupPath::Avx2 || (path==SetupPath::Automatic && cpu_avx2()));
#if SPIN_AVX512
    if(packet512) {
        BufferedWords words(routeSeed);
        prepareRoute(*this,words);
    } else
#endif
    { Words words(routeSeed); prepareRoute(*this,words); }
    if(packet512) {
#if SPIN_AVX512
        prepareInnerFast(*this,innerSeed);
#endif
    } else if(avx2) {
#if SPIN_X86
        prepareInnerAvx2(*this,innerSeed);
#endif
    } else if(path==SetupPath::Scalar) prepareInner(*this,innerSeed);
    else prepareInnerPortable(*this,innerSeed);
}
std::size_t Plan::setupBytes() const noexcept {
    return route.capacity()*sizeof(route[0])+routeBytes.capacity()*sizeof(routeBytes[0])
        +reverseMatrices.capacity()*sizeof(reverseMatrices[0])
        +updates.capacity()*sizeof(updates[0])
        +forwardUpdates.capacity()*sizeof(forwardUpdates[0]);
}
}
