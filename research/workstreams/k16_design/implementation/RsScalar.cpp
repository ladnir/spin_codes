#include "RsPrototype.h"
#include <bit>
#include <cstring>

namespace spin::research::rs {
namespace {
// Literal selected expansion columns; feedback is their binary transpose.
// Selected-map SHA256:
// 4652916f85eb3484d9fef6ce86fffeb0d6bca0e0ec0508eea801009e5b212fe7.
constexpr std::uint16_t columns[64] = {
    0x1,0x3,0x5,0x3507,0x9,0x118b,0x940d,0xb08f,
    0x11,0x5993,0xe415,0x8897,0xd999,0x919b,0xa99d,0xd49f,
    0x21,0x2ca3,0xb8a5,0xa127,0x1a29,0x272b,0x36ad,0x3eaf,
    0x6631,0x1333,0x3ab5,0x7ab7,0xa5b9,0xc13b,0x6d3d,0x3cbf,
    0x41,0xcfc3,0xe45,0xf4c7,0xee49,0x304b,0x744d,0x9f4f,
    0xb6d1,0x20d3,0x5cd5,0xffd7,0x8159,0x6db,0xff5d,0x4ddf,
    0x78e1,0x9be3,0xce65,0x1867,0x8ce9,0x7e6b,0xae6d,0x69ef,
    0xa871,0x12f3,0xfaf5,0x7577,0x85f9,0x2efb,0x437d,0xdd7f};

constexpr unsigned multiply(unsigned a, unsigned b) {
    unsigned result = 0;
    while(b) {
        if(b & 1) result ^= a;
        a <<= 1;
        if(a & 256) a ^= 0x11b;
        b >>= 1;
    }
    return result;
}

// Full systematic RS matrix from interpolation at points0,1,2,3 and
// evaluation at0,...,7. This oracle does not use the five-map factor.
constexpr unsigned generator[8][4] = {
    {1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1},
    {27,28,18,20},{28,27,20,18},{18,20,27,28},{20,18,28,27}};
constexpr auto rsRows = [] {
    std::array<std::uint64_t,32> rows{};
    for(unsigned messageBit = 0; messageBit < 32; ++messageBit)
        for(unsigned symbol = 0; symbol < 8; ++symbol) {
            const auto value = multiply(generator[symbol][messageBit / 8], 1U << (messageBit % 8));
            rows[messageBit] |= std::uint64_t(value) << (8 * symbol);
        }
    return rows;
}();

constexpr unsigned multiply16(unsigned a, unsigned b) {
    unsigned result = 0;
    for(unsigned i = 0; i < 4; ++i) {
        if(b & 1) result ^= a;
        b >>= 1;
        a <<= 1;
        if(a & 16) a ^= 0x13;
    }
    return result;
}

constexpr unsigned inverse16(unsigned a) {
    unsigned result = 1;
    for(unsigned i = 0; i < 14; ++i) result = multiply16(result, a);
    return result;
}

// Independently form the full32-by64 binary RS16 matrix by Lagrange
// interpolation. Neither the fast zeta factor nor its coefficients are used.
constexpr auto rsRows16 = [] {
    std::array<std::uint64_t,32> rows{};
    for(unsigned messageBit = 0; messageBit < 32; ++messageBit)
        for(unsigned symbol = 0; symbol < 16; ++symbol) {
            unsigned numerator = 1, denominator = 1;
            for(unsigned other = 0; other < 8; ++other)
                if(other != messageBit / 4) {
                    numerator = multiply16(numerator, symbol ^ other);
                    denominator = multiply16(denominator, (messageBit / 4) ^ other);
                }
            const auto coefficient = multiply16(numerator, inverse16(denominator));
            const auto value = multiply16(coefficient, 1U << (messageBit % 4));
            rows[messageBit] |= std::uint64_t(value) << (4 * symbol);
        }
    return rows;
}();

void reverseRoute(const Block* input, Block* scratch, const Plan& plan) {
    alignas(16) __m128i state[16]{}, next[16], feedback[16];
    for(std::size_t epoch = plan.n / 64; epoch-- > 0;) {
        for(auto& value : feedback) value = _mm_setzero_si128();
        for(unsigned p = 0; p < 64; ++p) {
            const auto i = 64 * epoch + p;
            const auto raw = input[i].mData;
            auto value = raw;
            for(unsigned mask = columns[p]; mask; mask &= mask - 1) {
                const auto j = std::countr_zero(mask);
                value = _mm_xor_si128(value, state[j]);
                // Raw input makes this independent of the fused C^T C=0 trick.
                feedback[j] = _mm_xor_si128(feedback[j], raw);
            }
            scratch[plan.route[i / 4] + (i & 3)] = Block(value);
        }
        if(!epoch) break;
        for(unsigned j = 0; j < 16; ++j) {
            auto value = feedback[j];
            for(unsigned mask = plan.reverseMatrices[epoch][j]; mask; mask &= mask - 1)
                value = _mm_xor_si128(value, state[std::countr_zero(mask)]);
            next[j] = value;
        }
        std::memcpy(state, next, sizeof(state));
    }
}

void outerTranspose(const Block* routed, Block* output,
                    const std::array<std::uint32_t,32>* matrices) {
    alignas(16) Block mixed[256];
    for(unsigned symbol = 0; symbol < 8; ++symbol)
        for(unsigned row = 0; row < 32; ++row) {
            auto value = _mm_setzero_si128();
            for(auto mask = matrices[symbol][row]; mask; mask &= mask - 1) {
                const auto c = std::countr_zero(mask);
                value = _mm_xor_si128(value, routed[4 * (8 * symbol + c % 8) + c / 8].mData);
            }
            mixed[4 * (8 * symbol + row % 8) + row / 8] = Block(value);
        }
    for(unsigned lane = 0; lane < 4; ++lane)
        for(unsigned messageBit = 0; messageBit < 32; ++messageBit) {
            auto value = _mm_setzero_si128();
            for(auto mask = rsRows[messageBit]; mask; mask &= mask - 1)
                value = _mm_xor_si128(value, mixed[4 * std::countr_zero(mask) + lane].mData);
            output[32 * lane + messageBit] = Block(value);
        }
}

void outerForward(const Block* message, Block* routed,
                  const std::array<std::uint32_t,32>* reverseMatrices) {
    alignas(16) Block coded[256];
    for(unsigned lane = 0; lane < 4; ++lane)
        for(unsigned bit = 0; bit < 64; ++bit) {
            auto value = _mm_setzero_si128();
            for(unsigned messageBit = 0; messageBit < 32; ++messageBit)
                if((rsRows[messageBit] >> bit) & 1)
                    value = _mm_xor_si128(value, message[32 * lane + messageBit].mData);
            coded[4 * bit + lane] = Block(value);
        }
    for(unsigned symbol = 0; symbol < 8; ++symbol)
        for(unsigned row = 0; row < 32; ++row) {
            auto value = _mm_setzero_si128();
            for(unsigned column = 0; column < 32; ++column)
                if((reverseMatrices[symbol][column] >> row) & 1)
                    value = _mm_xor_si128(value, coded[4 * (8 * symbol + column % 8) + column / 8].mData);
            routed[4 * (8 * symbol + row % 8) + row / 8] = Block(value);
        }
}

void outerTranspose16(const Block* routed, Block* output,
                      const std::array<std::uint16_t,16>* matrices) {
    alignas(16) Block mixed[256];
    for(unsigned symbol = 0; symbol < 16; ++symbol)
        for(unsigned row = 0; row < 16; ++row) {
            auto value = _mm_setzero_si128();
            for(unsigned mask = matrices[symbol][row]; mask; mask &= mask - 1) {
                const auto c = std::countr_zero(mask);
                value = _mm_xor_si128(value, routed[4 * (4 * symbol + c % 4) + c / 4].mData);
            }
            mixed[4 * (4 * symbol + row % 4) + row / 4] = Block(value);
        }
    for(unsigned lane = 0; lane < 4; ++lane)
        for(unsigned messageBit = 0; messageBit < 32; ++messageBit) {
            auto value = _mm_setzero_si128();
            for(auto mask = rsRows16[messageBit]; mask; mask &= mask - 1)
                value = _mm_xor_si128(value, mixed[4 * std::countr_zero(mask) + lane].mData);
            output[32 * lane + messageBit] = Block(value);
        }
}

void outerForward16(const Block* message, Block* routed,
                    const std::array<std::uint16_t,16>* reverseMatrices) {
    alignas(16) Block coded[256];
    for(unsigned lane = 0; lane < 4; ++lane)
        for(unsigned bit = 0; bit < 64; ++bit) {
            auto value = _mm_setzero_si128();
            for(unsigned messageBit = 0; messageBit < 32; ++messageBit)
                if((rsRows16[messageBit] >> bit) & 1)
                    value = _mm_xor_si128(value, message[32 * lane + messageBit].mData);
            coded[4 * bit + lane] = Block(value);
        }
    for(unsigned symbol = 0; symbol < 16; ++symbol)
        for(unsigned row = 0; row < 16; ++row) {
            auto value = _mm_setzero_si128();
            for(unsigned column = 0; column < 16; ++column)
                if((reverseMatrices[symbol][column] >> row) & 1)
                    value = _mm_xor_si128(value, coded[4 * (4 * symbol + column % 4) + column / 4].mData);
            routed[4 * (4 * symbol + row % 4) + row / 4] = Block(value);
        }
}
}

void transposeScalar(const Block* input, Block* output, Block* scratch, const Plan& plan) {
    reverseRoute(input, scratch, plan);
    if(plan.variant == Variant::Rs8Gf256) {
        for(std::size_t group = 0; group < plan.groups; ++group)
            outerTranspose(scratch + group * groupStride, output + group * 128,
                           plan.outerMatrices.data() + group * 8);
    } else {
        for(std::size_t group = 0; group < plan.groups; ++group)
            outerTranspose16(scratch + group * groupStride, output + group * 128,
                             plan.outerMatrices16.data() + group * 16);
    }
}

void forwardScalar(const Block* message, Block* encoded, const Plan& plan) {
    std::vector<Block> routed(plan.scratchBlocks());
    if(plan.variant == Variant::Rs8Gf256) {
        for(std::size_t group = 0; group < plan.groups; ++group)
            outerForward(message + group * 128, routed.data() + group * groupStride,
                         plan.outerMatrices.data() + group * 8);
    } else {
        for(std::size_t group = 0; group < plan.groups; ++group)
            outerForward16(message + group * 128, routed.data() + group * groupStride,
                           plan.outerMatrices16.data() + group * 16);
    }
    for(std::size_t i = 0; i < plan.n; ++i)
        encoded[i] = routed[plan.route[i / 4] + (i & 3)];
    alignas(16) __m128i state[16]{}, next[16], feedback[16];
    for(std::size_t epoch = 0; epoch < plan.n / 64; ++epoch) {
        for(auto& value : feedback) value = _mm_setzero_si128();
        for(unsigned p = 0; p < 64; ++p) {
            const auto i = 64 * epoch + p;
            const auto raw = encoded[i].mData;
            auto value = raw;
            for(unsigned mask = columns[p]; mask; mask &= mask - 1) {
                const auto j = std::countr_zero(mask);
                value = _mm_xor_si128(value, state[j]);
                feedback[j] = _mm_xor_si128(feedback[j], raw);
            }
            encoded[i] = Block(value);
        }
        if(epoch + 1 == plan.n / 64) break;
        for(unsigned j = 0; j < 16; ++j) {
            auto value = feedback[j];
            for(unsigned c = 0; c < 16; ++c)
                if((plan.reverseMatrices[epoch][c] >> j) & 1)
                    value = _mm_xor_si128(value, state[c]);
            next[j] = value;
        }
        std::memcpy(state, next, sizeof(state));
    }
}
}
