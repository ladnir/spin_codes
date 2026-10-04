#include "Plan.h"
#include "productionMap.h"
#include <bit>
#include <cstring>

namespace spin::detail::certified22 {
namespace {
constexpr unsigned multiply4(unsigned a, unsigned b) {
    unsigned result = 0;
    for(unsigned i = 0; i < 4; ++i) {
        if(b & 1U) result ^= a;
        b >>= 1;
        a <<= 1;
        if(a & 16U) a ^= 0x13U;
    }
    return result;
}
constexpr unsigned inverse4(unsigned a) {
    unsigned result = 1;
    for(unsigned i = 0; i < 14; ++i) result = multiply4(result, a);
    return result;
}
// Full binary generator from Lagrange interpolation; no fast RS factor.
constexpr auto rsRows = [] {
    std::array<std::uint64_t,32> rows{};
    for(unsigned messageBit = 0; messageBit < 32; ++messageBit)
        for(unsigned symbol = 0; symbol < 16; ++symbol) {
            unsigned numerator = 1, denominator = 1;
            for(unsigned other = 0; other < 8; ++other)
                if(other != messageBit/4) {
                    numerator = multiply4(numerator, symbol ^ other);
                    denominator = multiply4(denominator, (messageBit/4) ^ other);
                }
            const auto coefficient = multiply4(numerator, inverse4(denominator));
            const auto value = multiply4(coefficient, 1U << (messageBit%4));
            rows[messageBit] |= std::uint64_t(value) << (4*symbol);
        }
    return rows;
}();

template<bool Forward>
void inner(const Block* input, Block* output, const Plan& p, unsigned lanes, unsigned lane) {
    alignas(16) __m128i state[22]{}, next[22], feedback[22];
    for(std::size_t step = 0; step < p.n/64; ++step) {
        const auto epoch = Forward ? step : p.n/64 - 1 - step;
        for(auto& value : feedback) value = _mm_setzero_si128();
        for(unsigned j = 0; j < 64; ++j) {
            const auto i = 64*epoch+j, routed = p.route[i/4]+(i&3);
            const auto raw = input[Forward ? routed : i].mData;
            auto value = raw;
            for(auto mask = map::oldColumns[j]; mask; mask &= mask-1) {
                const auto c = std::countr_zero(mask);
                value = _mm_xor_si128(value, state[c]);
                feedback[c] = _mm_xor_si128(feedback[c], raw);
            }
            output[Forward ? i*lanes+lane : routed] = Block(value);
        }
        if(step+1 == p.n/64) break;
        for(unsigned j = 0; j < 22; ++j) {
            auto value = feedback[j];
            if constexpr(Forward) {
                for(unsigned c = 0; c < 22; ++c)
                    if((p.reverseMatrices[epoch][c] >> j) & 1U)
                        value = _mm_xor_si128(value, state[c]);
            } else {
                for(unsigned mask = p.reverseMatrices[epoch][j]; mask; mask &= mask-1)
                    value = _mm_xor_si128(value, state[std::countr_zero(mask)]);
            }
            next[j] = value;
        }
        std::memcpy(state, next, sizeof(state));
    }
}
void outerTranspose(const Block* input, Block* output, const std::array<std::uint16_t,16>* matrices) {
    alignas(16) Block mixed[256];
    for(unsigned symbol = 0; symbol < 16; ++symbol)
        for(unsigned row = 0; row < 16; ++row) {
            auto value = _mm_setzero_si128();
            for(unsigned mask = matrices[symbol][row]; mask; mask &= mask-1) {
                const auto c = std::countr_zero(mask);
                value = _mm_xor_si128(value, input[4*(4*symbol+c%4)+c/4].mData);
            }
            mixed[4*(4*symbol+row%4)+row/4] = Block(value);
        }
    for(unsigned row = 0; row < 4; ++row)
        for(unsigned bit = 0; bit < 32; ++bit) {
            auto value = _mm_setzero_si128();
            for(auto mask = rsRows[bit]; mask; mask &= mask-1)
                value = _mm_xor_si128(value, mixed[4*std::countr_zero(mask)+row].mData);
            output[32*row+bit] = Block(value);
        }
}
void outerForward(const Block* message, Block* scratch, const std::array<std::uint16_t,16>* matrices,
                  unsigned lanes, unsigned lane) {
    alignas(16) Block coded[256];
    for(unsigned row = 0; row < 4; ++row)
        for(unsigned bit = 0; bit < 64; ++bit) {
            auto value = _mm_setzero_si128();
            for(unsigned m = 0; m < 32; ++m)
                if((rsRows[m] >> bit) & 1U)
                    value = _mm_xor_si128(value, message[(32*row+m)*lanes+lane].mData);
            coded[4*bit+row] = Block(value);
        }
    for(unsigned symbol = 0; symbol < 16; ++symbol)
        for(unsigned row = 0; row < 16; ++row) {
            auto value = _mm_setzero_si128();
            for(unsigned column = 0; column < 16; ++column)
                if((matrices[symbol][column] >> row) & 1U)
                    value = _mm_xor_si128(value, coded[4*(4*symbol+column%4)+column/4].mData);
            scratch[4*(4*symbol+row%4)+row/4] = Block(value);
        }
}
}
void reverseRouteScalar(const Block* input, Block* scratch, const Plan& p) {
    inner<false>(input, scratch, p, 1, 0);
}
void forwardInnerScalar(const Block* scratch, Block* output, const Plan& p, unsigned lanes) {
    for(unsigned lane = 0; lane < lanes; ++lane)
        inner<true>(scratch+lane*p.scratchBlocks(), output, p, lanes, lane);
}
void transposeScalar(const Block* input, Block* output, Block* scratch, const Plan& p) {
    reverseRouteScalar(input, scratch, p);
    for(std::size_t group = 0; group < p.groups; ++group)
        outerTranspose(scratch+group*groupStride, output+group*128,
                       p.outerMatrices16.data()+group*16);
}
void forwardScalar(const Block* message, Block* output, Block* scratch, const Plan& p, unsigned lanes) {
    for(std::size_t group = 0; group < p.groups; ++group)
        for(unsigned lane = 0; lane < lanes; ++lane)
            outerForward(message+group*128*lanes, scratch+lane*p.scratchBlocks()+group*groupStride,
                         p.outerMatrices16.data()+group*16, lanes, lane);
    forwardInnerScalar(scratch, output, p, lanes);
}
}
