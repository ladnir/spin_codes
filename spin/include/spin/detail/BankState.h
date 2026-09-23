#pragma once
#include <spin/Code.h>
#include <bit>
#include <vector>

namespace spin::detail {
// Versioned heuristic family: one fixed permutation of rows per region,
// fixed independent row permutations, and fresh wrappers. Never a full fresh route.
struct BankState {
    CodeSpec spec;
    std::uint64_t bank_seed;
    std::size_t rows;
    std::vector<std::uint32_t> bank, row_keys, masks;
    // K18 native fast layout. Scalar/generic routing retains the natural table.
    std::vector<std::uint32_t> indexed, packed_keys;
    std::array<std::uint32_t,256> region{}, shift{};
    BankState(CodeSpec, std::uint64_t);
    void refresh(std::uint64_t routeSeed,std::uint64_t innerSeed) noexcept;
    std::size_t bytes() const noexcept {
        return sizeof(*this)+4*(bank.capacity()+row_keys.capacity()+masks.capacity()+indexed.capacity()+packed_keys.capacity());
    }
    static std::uint32_t column(std::uint32_t c,std::uint32_t key) noexcept {
        c=std::rotr(std::uint8_t(c),int(key>>24));
        return ((c*(key&255)+(key>>8))^(key>>16))&255;
    }
    std::uint32_t destination(std::size_t g,std::size_t p) const noexcept {
        p+=shift[g];if(p>=rows)p-=rows;
        const auto x=bank[region[g]*rows+p];
        return (x&~255U)|column(x&255,row_keys[x>>8]);
    }
};
}
