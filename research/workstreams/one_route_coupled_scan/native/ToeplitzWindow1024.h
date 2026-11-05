#pragma once

#include "Transform256.h"
#include "SparseByteField128Gfni.h"
#include "SparseByteField128ToeplitzDeferred.h"
#include <array>
#include <cstddef>
#include <cstring>
#ifndef SPIN_OR_WINDOW_COORDINATES
#define SPIN_OR_WINDOW_COORDINATES 0
#endif
static_assert(SPIN_OR_WINDOW_COORDINATES == 0 || SPIN_OR_WINDOW_COORDINATES == 1);

namespace spin::research::one_route {

#if defined(_MSC_VER)
#define SPIN_OR_WINDOW_INLINE __forceinline
#else
#define SPIN_OR_WINDOW_INLINE inline __attribute__((always_inline))
#endif

// Exact fixed-setup Toeplitz feedback: T[i,j]=raw[1023+i-j].
// This header changes neither the sampled1279 diagonals nor their distribution.
// The state is1024 preceding output symbols in oldest-first order; feedback
// returns256 symbols and does not modify that state. The caller adds the routed
// input, emits the complete batch, and only then calls append on those outputs.
struct ToeplitzWindow1024 {
    using Field = char2::field128_lift::sparse_byte::Element;
    static constexpr unsigned memory = 1024;
    static constexpr unsigned batch = 256;
    static constexpr unsigned rawCoefficients = memory + batch - 1;
    static constexpr unsigned preparedCoefficients = 12 * 256;
    static constexpr unsigned productsPerBatch = preparedCoefficients;
    static constexpr bool coordinateVariant = SPIN_OR_WINDOW_COORDINATES;
    static constexpr unsigned extraction = coordinateVariant ? 128 : 127;

    // Explicit alternative distribution, NOT the identical Toeplitz encoder.
    // Multiply evaluations by their point to replace D(Z) by Z*D(Z). Runtime
    // then uses novel input/output coordinates: T' = diag(U,U) T diag(V x8).
    // U and V are fixed invertible maps, so the uniform-feedback proof applies.
    // Setup evaluate() remains MONOMIAL, including in this alternative.
    static void prepareSpectrum(Field* values) noexcept {
        if constexpr (coordinateVariant) {
            constexpr std::uint64_t repeat = 0x0101010101010101ULL;
            for(unsigned i=0;i<256;i+=4){
                const auto a=repeat*Transform256Cantor::evaluationPoint(i);
                const auto b=repeat*Transform256Cantor::evaluationPoint(i+1);
                const auto c=repeat*Transform256Cantor::evaluationPoint(i+2);
                const auto d=repeat*Transform256Cantor::evaluationPoint(i+3);
                _mm512_store_si512(values+i,_mm512_gf2p8mul_epi8(
                    _mm512_load_si512(values+i),_mm512_set_epi64(d,d,c,c,b,b,a,a)));
            }
        }
    }

    // Four input chunks, each a256-square Toeplitz product. For one chunk,
    // split the square as [[A,B],[C,A]], where each block is128-square.
    // u=A(lo+hi), v=(A+C)lo, w=(A+B)hi gives (u+w,u+v).
    // The twelve transformed leaves are packed in their consumption order:
    // four evaluation points, four chunks, then A,A+C,A+B.
    struct alignas(64) PreparedBatch {
        Field values[64][4][3][4];
    };

    // Runtime scratch is separate from the cached state and setup. No hot-path
    // allocation occurs. Only two inverse transforms are needed: sum the four
    // chunk contributions and combine u+w/u+v before interpolation.
    struct alignas(64) Scratch {
        Field spectra[2][256];
    };

    // Setup-only temporary. Precomputing the leaf transforms is deliberately
    // separate from the critical encoding path. Every transform is MONOMIAL;
    // replacing evaluate/interpolate by their Novel variants changes the map.
    // raw and destination must not overlap.
    static void prepare(const Field* raw, PreparedBatch& destination) noexcept {
        alignas(64) Field leaves[3][256];
        for (unsigned chunk = 0; chunk != 4; ++chunk) {
            const unsigned base = 256 * (3 - chunk);
            for (unsigned k = 0; k != 255; ++k) {
                const auto a = raw[base + 128 + k];
                leaves[0][k] = a;
                leaves[1][k] = a ^ raw[base + 256 + k];
                leaves[2][k] = a ^ raw[base + k];
            }
            // Degree254, not255: degree255 would alias into extracted index127.
            leaves[0][255] = leaves[1][255] = leaves[2][255] = Field{};
            Transform256Cantor::evaluate(leaves[0]);
            Transform256Cantor::evaluate(leaves[1]);
            Transform256Cantor::evaluate(leaves[2]);
            prepareSpectrum(leaves[0]);
            prepareSpectrum(leaves[1]);
            prepareSpectrum(leaves[2]);
            for (unsigned p = 0; p != 64; ++p) {
                _mm512_store_si512(destination.values[p][chunk][0], _mm512_load_si512(leaves[0] + 4 * p));
                _mm512_store_si512(destination.values[p][chunk][1], _mm512_load_si512(leaves[1] + 4 * p));
                _mm512_store_si512(destination.values[p][chunk][2], _mm512_load_si512(leaves[2] + 4 * p));
            }
        }
    }

private:
    template<bool Gfni>
    static SPIN_OR_WINDOW_INLINE __m512i product(__m512i value, __m512i key) noexcept {
        if constexpr (Gfni)
            return char2::field128_lift::sparse_byte::multiply4Gfni(value, key);
        else
            return char2::field128_lift::sparse_byte::multiply4(value, key);
    }

    template<unsigned Chunk, bool Gfni>
    static SPIN_OR_WINDOW_INLINE void accumulate(const Field* state, const PreparedBatch& keys,
                                  unsigned p, __m512i& low, __m512i& high) noexcept {
        const auto x0 = _mm512_load_si512(state + 4 * p);
        const auto x1 = _mm512_load_si512(state + 256 + 4 * p);
        const auto shared = product<Gfni>(_mm512_xor_si512(x0, x1),
            _mm512_load_si512(keys.values[p][Chunk][0]));
        low = _mm512_xor_si512(low, shared);
        high = _mm512_xor_si512(high, shared);
        high = _mm512_xor_si512(high,
            product<Gfni>(x0, _mm512_load_si512(keys.values[p][Chunk][1])));
        low = _mm512_xor_si512(low,
            product<Gfni>(x1, _mm512_load_si512(keys.values[p][Chunk][2])));
    }

    template<unsigned Chunk>
    static SPIN_OR_WINDOW_INLINE void accumulateDeferred(const Field* state, const PreparedBatch& keys,
        unsigned p, char2::field128_lift::sparse_byte::deferred::Product& low,
        char2::field128_lift::sparse_byte::deferred::Product& high) noexcept {
        namespace deferred = char2::field128_lift::sparse_byte::deferred;
        const auto x0 = _mm512_load_si512(state + 4 * p);
        const auto x1 = _mm512_load_si512(state + 256 + 4 * p);
        {
            const auto shared = deferred::multiplyRaw4(_mm512_xor_si512(x0, x1),
                _mm512_load_si512(keys.values[p][Chunk][0]));
            low = deferred::sum(low, shared);
            high = deferred::sum(high, shared);
        }
        high = deferred::sum(high,
            deferred::multiplyRaw4(x0, _mm512_load_si512(keys.values[p][Chunk][1])));
        low = deferred::sum(low,
            deferred::multiplyRaw4(x1, _mm512_load_si512(keys.values[p][Chunk][2])));
    }

public:
    struct alignas(64) State {
        // Physical batches form a four-slot ring. A128-symbol half becomes a
        // degree<=127 polynomial, padded to256, then evaluated. A complete
        // half therefore occupies4096 bytes; all eight cached halves use32KiB.
        Field transformed[4][512]{};
        unsigned oldest = 0;

        void reset() noexcept {
            std::memset(transformed, 0, sizeof(transformed));
            oldest = 0;
        }

    private:
        static void transformBatch(const Field* input, Field* output) noexcept {
            std::memcpy(output, input, 128 * sizeof(Field));
            std::memset(output + 128, 0, 128 * sizeof(Field));
            std::memcpy(output + 256, input + 128, 128 * sizeof(Field));
            std::memset(output + 384, 0, 128 * sizeof(Field));
            if constexpr(coordinateVariant){
                Transform256Cantor::evaluateNovel(output);
                Transform256Cantor::evaluateNovel(output + 256);
            }else{
                Transform256Cantor::evaluate(output);
                Transform256Cantor::evaluate(output + 256);
            }
        }

    public:
        // Oracle/initialization interface. Input is1024 raw field symbols in
        // oldest-first order and must not overlap this State's internal cache.
        void load(const Field* rawState) noexcept {
            transformBatch(rawState, transformed[0]);
            transformBatch(rawState + 256, transformed[1]);
            transformBatch(rawState + 512, transformed[2]);
            transformBatch(rawState + 768, transformed[3]);
            oldest = 0;
        }

        // Exactly two new forward transforms per encoded batch. The previous
        // six half-transforms remain valid and are reused without conversion.
        // Input must not overlap this State's cache; ordinary output buffers
        // may immediately be appended after feedback and routed-input addition.
        void append(const Field* outputBatch) noexcept {
            transformBatch(outputBatch, transformed[oldest]);
            oldest = (oldest + 1) & 3U;
        }

        template<bool Gfni = false, bool Deferred = false>
        void feedback(const PreparedBatch& keys, Field* output, Scratch& scratch) const noexcept {
            // keys, output, scratch and this State must be mutually disjoint.
            static_assert(!(Gfni && Deferred), "Deferred is the CLMUL arithmetic alternative");
            const auto* s0 = transformed[oldest];
            const auto* s1 = transformed[(oldest + 1) & 3U];
            const auto* s2 = transformed[(oldest + 2) & 3U];
            const auto* s3 = transformed[(oldest + 3) & 3U];
            for (unsigned p = 0; p != 64; ++p) {
                if constexpr (Deferred) {
                    namespace deferred = char2::field128_lift::sparse_byte::deferred;
                    const auto zero = _mm512_setzero_si512();
                    deferred::Product low{zero, zero, zero, zero};
                    deferred::Product high{zero, zero, zero, zero};
                    accumulateDeferred<0>(s0, keys, p, low, high);
                    accumulateDeferred<1>(s1, keys, p, low, high);
                    accumulateDeferred<2>(s2, keys, p, low, high);
                    accumulateDeferred<3>(s3, keys, p, low, high);
                    // Both field reductions are GF2-linear. Sum all raw
                    // products first; no carry propagation or degree growth
                    // occurs under XOR. Each output spectrum is reduced once.
                    _mm512_store_si512(scratch.spectra[0] + 4 * p, deferred::reduce(low));
                    _mm512_store_si512(scratch.spectra[1] + 4 * p, deferred::reduce(high));
                } else {
                    auto low = _mm512_setzero_si512();
                    auto high = _mm512_setzero_si512();
                    accumulate<0, Gfni>(s0, keys, p, low, high);
                    accumulate<1, Gfni>(s1, keys, p, low, high);
                    accumulate<2, Gfni>(s2, keys, p, low, high);
                    accumulate<3, Gfni>(s3, keys, p, low, high);
                    _mm512_store_si512(scratch.spectra[0] + 4 * p, low);
                    _mm512_store_si512(scratch.spectra[1] + 4 * p, high);
                }
            }
            if constexpr(coordinateVariant){
                Transform256Cantor::interpolateNovel(scratch.spectra[0]);
                Transform256Cantor::interpolateNovel(scratch.spectra[1]);
            }else{
                Transform256Cantor::interpolate(scratch.spectra[0]);
                Transform256Cantor::interpolate(scratch.spectra[1]);
            }
            // Each leaf product has degree<=381. Interpolation at allGF256
            // points computes modulo X^256-X. Wrapped terms land only at
            // degrees1..126, so the requested middle127..254 is unchanged.
            for (unsigned j = 0; j != 128; j += 4) {
                _mm512_storeu_si512(output + j,
                    _mm512_loadu_si512(scratch.spectra[0] + extraction + j));
                _mm512_storeu_si512(output + 128 + j,
                    _mm512_loadu_si512(scratch.spectra[1] + extraction + j));
            }
        }
    };
};

static_assert(sizeof(ToeplitzWindow1024::PreparedBatch) == 12 * 256 * sizeof(Element));
static_assert(alignof(ToeplitzWindow1024::PreparedBatch) == 64);
static_assert(alignof(ToeplitzWindow1024::Scratch) == 64);
static_assert(alignof(ToeplitzWindow1024::State) == 64);

#undef SPIN_OR_WINDOW_INLINE
} // namespace spin::research::one_route
