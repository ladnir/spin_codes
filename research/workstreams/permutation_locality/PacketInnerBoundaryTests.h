#pragma once

// Setup/check-only tests for the reusable exact inner API. No timing, main,
// executable target, or optimized-coordinate oracle is part of this header.
#include "PacketInnerKernel.h"
#include "GfniT64.h"
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <span>
#include <stdexcept>
#include <vector>

namespace sizeprobe {
namespace packet_boundary_detail {
namespace inner = spin::research::packet_inner;
using Block = spin::detail::kernel::block;
using Matrix = inner::PreparedPacketUpdates16::OriginalReverseMatrix;
using Word = std::array<std::uint64_t,2>;
using Fixed = spin::research::t64::Fixed;
static_assert(sizeof(Block) == sizeof(Word));

inline void require(bool condition, const char* message) {
    if(!condition) throw std::runtime_error(message);
}

inline Word xored(Word a, const Word& b) noexcept {
    a[0] ^= b[0]; a[1] ^= b[1]; return a;
}

inline unsigned scalarRank(Matrix rows) noexcept {
    unsigned rank = 0;
    for(unsigned bit = 0; bit < 16; ++bit) {
        unsigned pivot = rank;
        while(pivot < 16 && !(rows[pivot] & (1U << bit))) ++pivot;
        if(pivot == 16) continue;
        std::swap(rows[pivot], rows[rank]);
        for(unsigned row = rank + 1; row < 16; ++row)
            if(rows[row] & (1U << bit)) rows[row] ^= rows[rank];
        ++rank;
    }
    return rank;
}

inline bool symmetric(const Matrix& rows) noexcept {
    for(unsigned i = 0; i < 16; ++i)
        for(unsigned j = 0; j < 16; ++j)
            if(((rows[i] >> j) & 1U) != ((rows[j] >> i) & 1U)) return false;
    return true;
}

inline std::array<Matrix,3> originalReverseMatrices() {
    std::array<Matrix,3> result{};
    for(unsigned epoch = 0; epoch < result.size(); ++epoch) {
        auto& rows = result[epoch];
        for(unsigned i = 0; i < 16; ++i) rows[i] = 1U << i;
        // Sequential elementary row additions and swaps preserve rank.
        // Each epoch has a different order and offsets. These fixtures are
        // already the original-coordinate REVERSE matrices, not forward M.
        for(unsigned pass = 0; pass < 3; ++pass)
            for(unsigned i = 0; i < 16; ++i)
                rows[i] ^= rows[(i + 1 + 2 * epoch + 2 * pass) % 16];
        std::swap(rows[epoch], rows[9 + epoch]);
        require(scalarRank(rows) == 16, "packet boundary GL16 fixture lost rank");
        require(!symmetric(rows), "packet boundary GL16 fixture is symmetric");
        for(unsigned previous = 0; previous < epoch; ++previous)
            require(rows != result[previous], "packet boundary GL16 fixtures repeat");
    }
    return result;
}

// Independent original-coordinate recurrence. All payload arithmetic is
// scalar uint64 XOR. The only shared mathematical data are the fixed
// expansion rows; no state basis, packed table, GFNI helper, zeta circuit,
// emission table, packet moment, or optimized conversion participates.
inline std::vector<Word> scalarReverse(std::span<const Word> input,
                                       std::span<const Matrix> originalReverse) {
    require(input.size() % 64 == 0 && originalReverse.size() == input.size() / 64,
            "packet boundary scalar oracle shape");
    std::vector<Word> output(input.size());
    std::array<Word,16> state{};
    for(std::size_t epoch = originalReverse.size(); epoch-- > 0;) {
        std::array<Word,16> feedback{};
        for(unsigned p = 0; p < 64; ++p) {
            const auto& raw = input[64 * epoch + p];
            auto value = raw;
            for(unsigned j = 0; j < 16; ++j) {
                if((Fixed::expansionRows[j] >> p) & 1U) {
                    value = xored(value, state[j]);
                    feedback[j] = xored(feedback[j], raw);
                }
            }
            output[64 * epoch + p] = value;
        }
        if(!epoch) break;
        std::array<Word,16> next = feedback;
        for(unsigned j = 0; j < 16; ++j)
            for(unsigned k = 0; k < 16; ++k)
                if((originalReverse[epoch][j] >> k) & 1U)
                    next[j] = xored(next[j], state[k]);
        state = next;
    }
    return output;
}

inline std::uint64_t mix(std::uint64_t value) noexcept {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

inline std::vector<Word> denseInput(std::size_t n) {
    std::vector<Word> result(n);
    for(std::size_t i = 0; i < n; ++i)
        result[i] = {mix(2 * i + 0x21345ULL), mix(2 * i + 0x91765ULL)};
    return result;
}

template<inner::InnerRecipe Recipe>
inline void compareCase(std::span<const Word> values,
                        std::span<const Matrix> originalReverse,
                        const inner::PreparedPacketUpdates16& prepared,
                        bool nullForEmpty = false) {
    constexpr std::size_t guard = 4;
    const Word sentinel{0xa359814072bb0de6ULL,0xd7ca816593ef024bULL};
    const auto n = values.size();
    require(n <= 192, "packet boundary fixture exceeds its tested size range");
    std::vector<Block> input(n + 2 * guard), output(n + 2 * guard);
    for(std::size_t i = 0; i < input.size(); ++i) {
        std::memcpy(input.data() + i, sentinel.data(), sizeof(Block));
        std::memcpy(output.data() + i, sentinel.data(), sizeof(Block));
    }
    for(std::size_t i = 0; i < n; ++i)
        std::memcpy(input.data() + guard + i, values[i].data(), sizeof(Block));
    const auto savedInput = input;
    const auto expected = scalarReverse(values, originalReverse);
    std::size_t remaining = n / 4;
    auto emit = [&](std::size_t packet, __m512i value) {
        require(remaining != 0, "packet boundary emitted an extra packet");
        require(packet == remaining - 1, "packet boundary emission order or index");
        --remaining;
        _mm512_storeu_si512(output.data() + guard + 4 * packet, value);
    };
    const Block* source = nullForEmpty ? nullptr : input.data() + guard;
    require(!nullForEmpty || !n, "packet boundary null fixture must be empty");
    inner::transposeInner<Recipe>(source, n, prepared.conjugatedRows(), emit);
    require(remaining == 0, "packet boundary omitted an output packet");
    require(std::memcmp(input.data(), savedInput.data(), input.size() * sizeof(Block)) == 0,
            "packet boundary changed input or input guards");
    for(std::size_t i = 0; i < guard; ++i) {
        require(std::memcmp(output.data() + i, sentinel.data(), sizeof(Block)) == 0,
                "packet boundary leading output guard changed");
        require(std::memcmp(output.data() + guard + n + i, sentinel.data(), sizeof(Block)) == 0,
                "packet boundary trailing output guard changed");
    }
    for(std::size_t i = 0; i < n; ++i)
        require(std::memcmp(output.data() + guard + i, expected[i].data(), sizeof(Block)) == 0,
                "packet boundary output differs from original-coordinate scalar recurrence");
}

template<class Function>
inline void expectInvalid(Function&& function, const char* message) {
    bool rejected = false;
    try { function(); }
    catch(const std::invalid_argument&) { rejected = true; }
    require(rejected, message);
}

template<inner::InnerRecipe Recipe>
inline void invalidCalls(const std::array<Matrix,3>& matrices) {
    const auto all = std::span<const Matrix>(matrices);
    inner::PreparedPacketUpdates16 empty(all.first(0)), one(all.first(1)), two(all.first(2)), three(all);
    std::array<Block,192> input{};
    unsigned emitted = 0;
    auto emit = [&](std::size_t, __m512i) { ++emitted; };
    auto invalid = [&](const Block* source, std::size_t n,
                       inner::PreparedPacketUpdates16::ConjugatedRows updates) {
        expectInvalid([&] { inner::transposeInner<Recipe>(source, n, updates, emit); },
                      "packet boundary malformed inner call was accepted");
        require(emitted == 0, "packet boundary rejected call emitted output");
    };
    invalid(input.data(), 1, empty.conjugatedRows());
    invalid(input.data(), 63, empty.conjugatedRows());
    invalid(input.data(), 65, one.conjugatedRows());
    invalid(input.data(), 64, empty.conjugatedRows());
    invalid(input.data(), 0, one.conjugatedRows());
    invalid(input.data(), 128, three.conjugatedRows());
    invalid(input.data(), 192, two.conjugatedRows());
    invalid(nullptr, 64, one.conjugatedRows());
    inner::transposeInner<Recipe>(nullptr, 0, empty.conjugatedRows(), emit);
    require(emitted == 0, "packet boundary empty null input emitted output");
}

template<inner::InnerRecipe Recipe>
inline void recipeCases(const std::array<Matrix,3>& matrices) {
    const auto all = std::span<const Matrix>(matrices);
    for(const std::size_t n : {std::size_t(0), std::size_t(64), std::size_t(128), std::size_t(192)}) {
        const auto reverse = all.first(n / 64);
        inner::PreparedPacketUpdates16 prepared(reverse);
        require(prepared.epochs() == n / 64 && prepared.conjugatedRows().size() == n / 64,
                "packet boundary prepared epoch count");
        auto values = denseInput(n);
        compareCase<Recipe>(values, reverse, prepared);
        std::fill(values.begin(), values.end(), Word{});
        compareCase<Recipe>(values, reverse, prepared, n == 0);
        for(std::size_t coordinate = 0; coordinate < n; ++coordinate) {
            values[coordinate] = {0x729aca378bc42fedULL,0x942163aca892317bULL};
            compareCase<Recipe>(values, reverse, prepared);
            values[coordinate] = {};
        }
        // Explicit payload-bit bases at every endpoint and both sides of
        // a physical-step boundary detect lane/half/bit conversion mistakes.
        std::vector<std::size_t> boundaries;
        if(n) boundaries = {0, n - 1};
        for(std::size_t step = 64; step < n; step += 64) {
            boundaries.push_back(step - 1); boundaries.push_back(step);
        }
        std::sort(boundaries.begin(), boundaries.end());
        boundaries.erase(std::unique(boundaries.begin(), boundaries.end()), boundaries.end());
        for(const auto coordinate : boundaries)
            for(unsigned bit = 0; bit < 128; ++bit) {
                values[coordinate][bit / 64] = std::uint64_t(1) << (bit % 64);
                compareCase<Recipe>(values, reverse, prepared);
                values[coordinate] = {};
            }
    }
    invalidCalls<Recipe>(matrices);

    // Explicitly opting out of rank validation permits singular linear
    // updates, but still must implement the exact same recurrence. At three
    // epochs, index1 is genuinely used on a nonzero entering state.
    auto singular = matrices;
    singular[1].fill(0);
    inner::PreparedPacketUpdates16 allowed(std::span<const Matrix>(singular), false);
    auto values = denseInput(192);
    compareCase<Recipe>(values, singular, allowed);
    values.assign(192, Word{}); values[191] = {0x9bd2157408e3fa61ULL,0x183a7ec496d205bfULL};
    compareCase<Recipe>(values, singular, allowed);

    // Preparation owns its material: subsequent changes to the source
    // matrices must not alter the already prepared encoder.
    auto mutableSource = matrices;
    inner::PreparedPacketUpdates16 owned{std::span<const Matrix>(mutableSource)};
    for(auto& matrix : mutableSource) matrix.fill(0);
    values = denseInput(192);
    compareCase<Recipe>(values, matrices, owned);
}

inline void invalidPreparation(const std::array<Matrix,3>& matrices) {
    auto invalid = matrices;
    invalid[1][3] |= 1U << 16;
    expectInvalid([&] { inner::PreparedPacketUpdates16 rejected{std::span<const Matrix>(invalid)}; },
                  "packet boundary preparation accepted mask bit16");
    expectInvalid([&] { inner::PreparedPacketUpdates16 rejected(std::span<const Matrix>(invalid), false); },
                  "packet boundary rank opt-out bypassed mask-width validation");
    invalid = matrices; invalid[1].fill(0);
    expectInvalid([&] { inner::PreparedPacketUpdates16 rejected{std::span<const Matrix>(invalid)}; },
                  "packet boundary default preparation accepted singular matrix");
    // Unused boundary entries are also validated by the preparation API.
    invalid = matrices; invalid[0] = invalid[2] = Matrix{};
    expectInvalid([&] { inner::PreparedPacketUpdates16 rejected{std::span<const Matrix>(invalid)}; },
                  "packet boundary preparation skipped singular boundary matrices");
    const inner::PreparedPacketUpdates16 checked{std::span<const Matrix>(matrices)};
    const inner::PreparedPacketUpdates16 unchecked(std::span<const Matrix>(matrices), false);
    require(std::memcmp(checked.conjugatedRows().data(), unchecked.conjugatedRows().data(),
                        3 * sizeof(inner::PreparedPacketUpdates16::PackedRow)) == 0,
            "packet boundary rank option changed prepared matrix values");
}

} // namespace packet_boundary_detail

// Call only after the frontend has verified the required ISA. The explicit
// SPIN_PACKET_ENABLE_VBMI override must be consistent across translation
// units, as required by PacketInnerKernel.h. A non-VBMI build checks the
// retained Packet and StreamLegacy recipes without naming VBMI intrinsics.
// This is an ordinary check function with SIMD scratch, never a coroutine.
inline void checkInnerBoundaries() {
    namespace tests = packet_boundary_detail;
    namespace inner = spin::research::packet_inner;
    const auto matrices = tests::originalReverseMatrices();
    tests::invalidPreparation(matrices);
    tests::recipeCases<inner::InnerRecipe::Packet>(matrices);
    tests::recipeCases<inner::InnerRecipe::StreamLegacy>(matrices);
#if defined(__AVX512VBMI__) || (defined(SPIN_PACKET_ENABLE_VBMI) && SPIN_PACKET_ENABLE_VBMI)
    tests::recipeCases<inner::InnerRecipe::PacketVbmi>(matrices);
    tests::recipeCases<inner::InnerRecipe::StreamVbmi>(matrices);
    tests::recipeCases<inner::InnerRecipe::StreamVbmiEarly>(matrices);
    tests::recipeCases<inner::InnerRecipe::StreamVbmiPruned>(matrices);
    tests::recipeCases<inner::InnerRecipe::StreamVbmiPrunedEarly>(matrices);
#endif
}

} // namespace sizeprobe
