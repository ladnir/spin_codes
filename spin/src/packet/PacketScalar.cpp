#include "PacketPlan.h"
#include "../kernels/generated/BchCircuit.h"
#include <bit>
#include <cstring>

namespace spin::detail::packet {
namespace {
using storage::block;
// The physical t64/s16 expansion map, selected-map SHA256
// 4652916f85eb3484d9fef6ce86fffeb0d6bca0e0ec0508eea801009e5b212fe7.
// Feedback is its transpose. These are original, not fused-basis coordinates.
constexpr std::uint16_t columns[64] = {
    0x1,0x3,0x5,0x3507,0x9,0x118b,0x940d,0xb08f,
    0x11,0x5993,0xe415,0x8897,0xd999,0x919b,0xa99d,0xd49f,
    0x21,0x2ca3,0xb8a5,0xa127,0x1a29,0x272b,0x36ad,0x3eaf,
    0x6631,0x1333,0x3ab5,0x7ab7,0xa5b9,0xc13b,0x6d3d,0x3cbf,
    0x41,0xcfc3,0xe45,0xf4c7,0xee49,0x304b,0x744d,0x9f4f,
    0xb6d1,0x20d3,0x5cd5,0xffd7,0x8159,0x6db,0xff5d,0x4ddf,
    0x78e1,0x9be3,0xce65,0x1867,0x8ce9,0x7e6b,0xae6d,0x69ef,
    0xa871,0x12f3,0xfaf5,0x7577,0x85f9,0x2efb,0x437d,0xdd7f};

void innerRoute(const block* input, block* scratch, const Plan& plan) {
    alignas(16) __m128i state[16]{}, next[16], feedback[16];
    const auto epochs = plan.n / 64;
    for(std::size_t epoch = epochs; epoch-- > 0;) {
        for(auto& x : feedback) x = _mm_setzero_si128();
        for(unsigned p = 0; p < 64; ++p) {
            const auto i = 64 * epoch + p;
            const auto raw = input[i].mData;
            auto value = raw;
            for(unsigned mask = columns[p]; mask; mask &= mask - 1) {
                const auto j = std::countr_zero(mask);
                value = _mm_xor_si128(value, state[j]);
                // Use raw input: this oracle does not rely on A^T A=0.
                feedback[j] = _mm_xor_si128(feedback[j], raw);
            }
            scratch[plan.route[i / 4] + (i & 3)] = block(value);
        }
        if(!epoch) break; // No flushing and no update after the final output.
        for(unsigned j = 0; j < 16; ++j) {
            auto value = feedback[j];
            for(unsigned mask = plan.reverseMatrices[epoch][j]; mask; mask &= mask - 1)
                value = _mm_xor_si128(value, state[std::countr_zero(mask)]);
            next[j] = value;
        }
        std::memcpy(state, next, sizeof(state));
    }
}

void outerTile(const block* routed, block* output,
               const std::array<std::uint32_t, 32>* matrices) {
    alignas(16) block mixed[1024];
    for(unsigned group = 0; group < 32; ++group)
        for(unsigned row = 0; row < 32; ++row) {
            auto value = _mm_setzero_si128();
            for(auto mask = matrices[group][row]; mask; mask &= mask - 1) {
                const auto c = std::countr_zero(mask);
                value = _mm_xor_si128(value, routed[4 * (8 * group + c % 8) + c / 8].mData);
            }
            mixed[4 * (8 * group + row % 8) + row / 8] = block(value);
        }
    // Literal BCH action gives an independent fallback for the fused circuit.
    for(unsigned lane = 0; lane < 4; ++lane)
        for(unsigned row = 0; row < 128; ++row) {
            auto value = _mm_setzero_si128();
            for(unsigned limb = 0; limb < 4; ++limb)
                for(auto mask = kernel::BchRows[row][limb]; mask; mask &= mask - 1) {
                    const auto c = 64 * limb + std::countr_zero(mask);
                    value = _mm_xor_si128(value, mixed[4 * c + lane].mData);
                }
            output[128 * lane + row] = block(value);
        }
}
}

void transposeScalar(const storage::block* input, storage::block* output,
                     storage::block* scratch, const Plan& plan) {
    innerRoute(input, scratch, plan);
    for(std::size_t tile = 0; tile < plan.n / 1024; ++tile)
        outerTile(scratch + tile * tileStride, output + tile * 512,
                  plan.outerMatrices.data() + tile * 32);
}
}
