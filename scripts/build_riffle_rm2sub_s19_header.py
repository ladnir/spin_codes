#!/usr/bin/env python3
"""Generate the optimized transposed RM2Sub-S19 C++ header.

The generator binds the selected s=19 constituent to the five-group A
partition recorded by the performance search.  It also synthesizes and
verifies the small quadratic projection used after the pruned RM transform.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import probe_bch_forward_xor_circuit as paar  # noqa: E402


DEFAULT_SELECTION = ROOT / (
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_"
    "rm2sub_t128_s20/receipts/min_state/s19_rm2sub_selection.json"
)
DEFAULT_OUTPUT = ROOT / (
    "constructions/riffle_parityshear12_bchperm_transpose_bitshuffle_"
    "splitstate_preaddmul_rm2sub_t128_s19/benchmark/RiffleRm2SubS19.h"
)
DEFAULT_RECEIPT = ROOT / (
    "constructions/riffle_parityshear12_bchperm_transpose_bitshuffle_"
    "splitstate_preaddmul_rm2sub_t128_s19/benchmark/"
    "RiffleRm2SubS19.generated.json"
)

# The 30,000-candidate search minimized the grouped A XOR ledger.
A_GROUPS = (
    (0, 16, 6, 2),
    (14, 18, 1, 10),
    (9, 12, 3, 5),
    (17, 13, 15, 8),
    (11, 7, 4),
)
MONOMIAL_POSITIONS = (
    0x03, 0x05, 0x09, 0x11, 0x21, 0x41,
    0x06, 0x0A, 0x12, 0x22, 0x42,
    0x0C, 0x14, 0x24, 0x44,
    0x18, 0x28, 0x48,
    0x30, 0x50,
    0x60,
)


def cpp_xor(signals: list[int], base_names: list[str], prefix: str) -> str:
    def name(signal: int) -> str:
        return base_names[signal] if signal < len(base_names) else f"{prefix}{signal}"

    expression = name(signals[0])
    for signal in signals[1:]:
        expression = f"rm2s19Xor({expression}, {name(signal)})"
    return expression


def projection_code(masks: list[int]) -> tuple[str, paar.Circuit]:
    paar.DIMENSION = 21
    circuit = paar.synthesize(128, 2, masks)
    base = [f"monomials[{index}]" for index in range(21)]
    lines: list[str] = []
    for signal, (left, right) in enumerate(circuit.gates, start=21):
        lines.append(
            f"\t\t\tconst auto q{signal} = rm2s19Xor("
            f"{base[left] if left < 21 else f'q{left}'}, "
            f"{base[right] if right < 21 else f'q{right}'});"
        )
    for output, signals in enumerate(circuit.output_signals):
        lines.append(
            f"\t\t\toutput[{8 + output}] = "
            f"{cpp_xor(signals, base, 'q')};"
        )
    return "\n".join(lines), circuit


def format_array(values: list[int], width: int, per_line: int = 8) -> str:
    rows = []
    for start in range(0, len(values), per_line):
        row = ", ".join(f"0x{value:0{width}x}" for value in values[start:start + per_line])
        rows.append(f"\t\t\t{row},")
    return "\n".join(rows)


def make_header(columns: list[int], masks: list[int], projection: str) -> str:
    groups = "\n".join(
        f"\t\tinline constexpr std::array<unsigned, {len(group)}> "
        f"Rm2Sub19Group{index}{{{', '.join(map(str, group))}}};"
        for index, group in enumerate(A_GROUPS)
    )
    column_array = format_array(columns, 5)
    mask_array = format_array(masks, 6, 6)
    monomial_array = ", ".join(f"0x{x:02x}" for x in MONOMIAL_POSITIONS)

    return f'''#pragma once

#include <array>
#include <bit>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <utility>
#include <vector>

#include <immintrin.h>
#include <cryptoTools/Common/Defines.h>

namespace osuCrypto
{{
\tnamespace detail
\t{{
\t\tinline constexpr std::array<u32, 11> Rm2Sub19QuadraticMasks{{
{mask_array}
\t\t}};

\t\tinline constexpr std::array<unsigned, 21> Rm2Sub19MonomialPositions{{
\t\t\t{monomial_array},
\t\t}};

\t\tinline constexpr std::array<u32, 128> Rm2Sub19Columns{{
{column_array}
\t\t}};

{groups}

\t\tstruct Rm2Sub19FieldTransposeSchedule
\t\t{{
\t\t\talignas(32) std::array<u32, 19> rowMasks{{}};
\t\t}};

\t\tOC_FORCEINLINE __m128i rm2s19Xor(__m128i left, __m128i right) noexcept
\t\t{{
\t\t\treturn _mm_xor_si128(left, right);
\t\t}}

\t\tinline u32 rm2s19FieldMultiplyScalar(u32 left, u32 right) noexcept
\t\t{{
\t\t\tconstexpr u32 modulusLow = 0x27; // x^19+x^5+x^2+x+1.
\t\t\tu32 result = 0;
\t\t\tfor (unsigned bit = 0; bit < 19; ++bit)
\t\t\t{{
\t\t\t\tresult ^= left & (0U - (right & 1U));
\t\t\t\tright >>= 1;
\t\t\t\tconst u32 carry = left >> 18;
\t\t\t\tleft = (left << 1) & 0x7ffffU;
\t\t\t\tleft ^= modulusLow & (0U - carry);
\t\t\t}}
\t\t\treturn result;
\t\t}}

\t\tinline Rm2Sub19FieldTransposeSchedule rm2s19MakeFieldTransposeSchedule(
\t\t\tu32 coefficient) noexcept
\t\t{{
\t\t\tRm2Sub19FieldTransposeSchedule schedule{{}};
\t\t\tfor (unsigned row = 0; row < 19; ++row)
\t\t\t\tschedule.rowMasks[row] = rm2s19FieldMultiplyScalar(1U << row, coefficient);
\t\t\treturn schedule;
\t\t}}

\t\tOC_FORCEINLINE void rm2s19BuildNibbleTable(
\t\t\tconst __m128i* input, __m128i* table) noexcept
\t\t{{
\t\t\ttable[0] = _mm_setzero_si128();
\t\t\ttable[1] = input[0];
\t\t\ttable[2] = input[1];
\t\t\ttable[3] = rm2s19Xor(input[0], input[1]);
\t\t\ttable[4] = input[2];
\t\t\ttable[5] = rm2s19Xor(input[2], table[1]);
\t\t\ttable[6] = rm2s19Xor(input[2], table[2]);
\t\t\ttable[7] = rm2s19Xor(input[2], table[3]);
\t\t\ttable[8] = input[3];
\t\t\ttable[9] = rm2s19Xor(input[3], table[1]);
\t\t\ttable[10] = rm2s19Xor(input[3], table[2]);
\t\t\ttable[11] = rm2s19Xor(input[3], table[3]);
\t\t\ttable[12] = rm2s19Xor(input[3], table[4]);
\t\t\ttable[13] = rm2s19Xor(input[3], table[5]);
\t\t\ttable[14] = rm2s19Xor(input[3], table[6]);
\t\t\ttable[15] = rm2s19Xor(input[3], table[7]);
\t\t}}

\t\tOC_FORCEINLINE void rm2s19BuildThreeBitTable(
\t\t\tconst __m128i* input, __m128i* table) noexcept
\t\t{{
\t\t\ttable[0] = _mm_setzero_si128();
\t\t\ttable[1] = input[0];
\t\t\ttable[2] = input[1];
\t\t\ttable[3] = rm2s19Xor(input[0], input[1]);
\t\t\ttable[4] = input[2];
\t\t\ttable[5] = rm2s19Xor(input[2], table[1]);
\t\t\ttable[6] = rm2s19Xor(input[2], table[2]);
\t\t\ttable[7] = rm2s19Xor(input[2], table[3]);
\t\t}}

\t\tOC_FORCEINLINE void rm2s19FieldMultiplyTranspose(
\t\t\t__m128i* state, const Rm2Sub19FieldTransposeSchedule& schedule) noexcept
\t\t{{
\t\t\talignas(32) __m128i table[4][16];
\t\t\tfor (unsigned group = 0; group < 4; ++group)
\t\t\t\trm2s19BuildNibbleTable(state + 4 * group, table[group]);
\t\t\talignas(32) __m128i finalTable[8];
\t\t\trm2s19BuildThreeBitTable(state + 16, finalTable);
\t\t\talignas(32) __m128i output[19];
\t\t\tfor (unsigned row = 0; row < 19; ++row)
\t\t\t{{
\t\t\t\tconst u32 mask = schedule.rowMasks[row];
\t\t\t\tauto value = rm2s19Xor(table[0][mask & 15], table[1][(mask >> 4) & 15]);
\t\t\t\tvalue = rm2s19Xor(value, table[2][(mask >> 8) & 15]);
\t\t\t\tvalue = rm2s19Xor(value, table[3][(mask >> 12) & 15]);
\t\t\t\toutput[row] = rm2s19Xor(value, finalTable[(mask >> 16) & 7]);
\t\t\t}}
\t\t\tstd::memcpy(state, output, sizeof(output));
\t\t}}

\t\ttemplate<unsigned Distance>
\t\tOC_FORCEINLINE void rm2s19ZetaPair(__m128i* values, unsigned low) noexcept
\t\t{{
\t\t\tauto lhs = _mm256_load_si256(reinterpret_cast<const __m256i*>(values + low));
\t\t\tconst auto rhs = _mm256_load_si256(
\t\t\t\treinterpret_cast<const __m256i*>(values + low + Distance));
\t\t\t_mm256_store_si256(reinterpret_cast<__m256i*>(values + low), _mm256_xor_si256(lhs, rhs));
\t\t}}

#if defined(__AVX512F__)
\t\ttemplate<unsigned Distance>
\t\tOC_FORCEINLINE void rm2s19ZetaQuad(__m128i* values, unsigned low) noexcept
\t\t{{
\t\t\tauto lhs = _mm512_loadu_si512(values + low);
\t\t\tconst auto rhs = _mm512_loadu_si512(values + low + Distance);
\t\t\t_mm512_storeu_si512(values + low, _mm512_xor_si512(lhs, rhs));
\t\t}}
#endif

\t\ttemplate<unsigned Distance>
\t\tOC_FORCEINLINE void rm2s19ZetaStage(__m128i* values) noexcept
\t\t{{
\t\t\tfor (unsigned base = 0; base < 128; base += 2 * Distance)
\t\t\t{{
\t\t\t\tif constexpr (Distance == 1)
\t\t\t\t\tvalues[base] = rm2s19Xor(values[base], values[base + 1]);
\t\t\t\telse
\t\t\t\t\tfor (unsigned offset = 0; offset < Distance; offset += 2)
\t\t\t\t\t\trm2s19ZetaPair<Distance>(values, base + offset);
\t\t\t}}
\t\t}}

\t\tOC_FORCEINLINE void rm2s19ZetaTransposePruned(__m128i* values) noexcept
\t\t{{
#if defined(__AVX512F__)
\t\t\tfor (unsigned offset = 0; offset < 64; offset += 4) rm2s19ZetaQuad<64>(values, offset);
\t\t\tfor (unsigned base = 0; base < 128; base += 64)
\t\t\t\tfor (unsigned offset = 0; offset < 32; offset += 4) rm2s19ZetaQuad<32>(values, base + offset);
\t\t\tfor (unsigned base = 0; base < 128; base += 32)
\t\t\t\tfor (unsigned offset = 0; offset < 16; offset += 4) rm2s19ZetaQuad<16>(values, base + offset);
\t\t\tfor (unsigned base = 0; base <= 96; base += 16)
\t\t\t\tfor (unsigned offset = 0; offset < 8; offset += 4) rm2s19ZetaQuad<8>(values, base + offset);
\t\t\tfor (unsigned base = 0; base <= 48; base += 8) rm2s19ZetaQuad<4>(values, base);
\t\t\tfor (unsigned base = 64; base <= 80; base += 8) rm2s19ZetaQuad<4>(values, base);
\t\t\trm2s19ZetaQuad<4>(values, 96);
#else
\t\t\trm2s19ZetaStage<64>(values);
\t\t\trm2s19ZetaStage<32>(values);
\t\t\trm2s19ZetaStage<16>(values);
\t\t\tfor (unsigned base = 0; base <= 96; base += 16)
\t\t\t\tfor (unsigned offset = 0; offset < 8; offset += 2) rm2s19ZetaPair<8>(values, base + offset);
\t\t\tfor (unsigned base = 0; base <= 48; base += 8)
\t\t\t\tfor (unsigned offset = 0; offset < 4; offset += 2) rm2s19ZetaPair<4>(values, base + offset);
\t\t\tfor (unsigned base = 64; base <= 80; base += 8)
\t\t\t\tfor (unsigned offset = 0; offset < 4; offset += 2) rm2s19ZetaPair<4>(values, base + offset);
\t\t\trm2s19ZetaPair<4>(values, 96); rm2s19ZetaPair<4>(values, 98);
#endif
\t\t\tfor (unsigned base = 0; base <= 24; base += 4) rm2s19ZetaPair<2>(values, base);
\t\t\tfor (unsigned base = 32; base <= 40; base += 4) rm2s19ZetaPair<2>(values, base);
\t\t\trm2s19ZetaPair<2>(values, 48);
\t\t\tfor (unsigned base = 64; base <= 72; base += 4) rm2s19ZetaPair<2>(values, base);
\t\t\trm2s19ZetaPair<2>(values, 80); rm2s19ZetaPair<2>(values, 96);
#define RIFFLE_RM2S19_D1(Low) values[Low] = rm2s19Xor(values[Low], values[(Low) + 1])
\t\t\tRIFFLE_RM2S19_D1(0);  RIFFLE_RM2S19_D1(2);  RIFFLE_RM2S19_D1(4);  RIFFLE_RM2S19_D1(6);
\t\t\tRIFFLE_RM2S19_D1(8);  RIFFLE_RM2S19_D1(10); RIFFLE_RM2S19_D1(12); RIFFLE_RM2S19_D1(16);
\t\t\tRIFFLE_RM2S19_D1(18); RIFFLE_RM2S19_D1(20); RIFFLE_RM2S19_D1(24); RIFFLE_RM2S19_D1(32);
\t\t\tRIFFLE_RM2S19_D1(34); RIFFLE_RM2S19_D1(36); RIFFLE_RM2S19_D1(40); RIFFLE_RM2S19_D1(48);
\t\t\tRIFFLE_RM2S19_D1(64); RIFFLE_RM2S19_D1(66); RIFFLE_RM2S19_D1(68); RIFFLE_RM2S19_D1(72);
\t\t\tRIFFLE_RM2S19_D1(80); RIFFLE_RM2S19_D1(96);
#undef RIFFLE_RM2S19_D1
\t\t}}

\t\tOC_FORCEINLINE void rm2s19FinishB(__m128i* values, __m128i* output) noexcept
\t\t{{
\t\t\trm2s19ZetaTransposePruned(values);
\t\t\toutput[0] = values[0];
\t\t\tfor (unsigned linear = 0; linear < 7; ++linear) output[1 + linear] = values[1U << linear];
\t\t\talignas(32) __m128i monomials[21];
\t\t\tfor (unsigned monomial = 0; monomial < 21; ++monomial)
\t\t\t\tmonomials[monomial] = values[Rm2Sub19MonomialPositions[monomial]];
{projection}
\t\t}}

\t\ttemplate<std::size_t Point>
\t\tOC_FORCEINLINE __m128i rm2s19GroupedAAddend(
\t\t\tconst __m128i table[4][16], const __m128i finalTable[8]) noexcept
\t\t{{
\t\t\tconstexpr u32 column = Rm2Sub19Columns[Point];
\t\t\tconstexpr unsigned i0 = ((column >> 0) & 1U) | (((column >> 16) & 1U) << 1) |
\t\t\t\t(((column >> 6) & 1U) << 2) | (((column >> 2) & 1U) << 3);
\t\t\tconstexpr unsigned i1 = ((column >> 14) & 1U) | (((column >> 18) & 1U) << 1) |
\t\t\t\t(((column >> 1) & 1U) << 2) | (((column >> 10) & 1U) << 3);
\t\t\tconstexpr unsigned i2 = ((column >> 9) & 1U) | (((column >> 12) & 1U) << 1) |
\t\t\t\t(((column >> 3) & 1U) << 2) | (((column >> 5) & 1U) << 3);
\t\t\tconstexpr unsigned i3 = ((column >> 17) & 1U) | (((column >> 13) & 1U) << 1) |
\t\t\t\t(((column >> 15) & 1U) << 2) | (((column >> 8) & 1U) << 3);
\t\t\tconstexpr unsigned i4 = ((column >> 11) & 1U) | (((column >> 7) & 1U) << 1) |
\t\t\t\t(((column >> 4) & 1U) << 2);
\t\t\tstatic_assert(i0 != 0);
\t\t\tauto addend = table[0][i0];
\t\t\tif constexpr (i1 != 0) addend = rm2s19Xor(addend, table[1][i1]);
\t\t\tif constexpr (i2 != 0) addend = rm2s19Xor(addend, table[2][i2]);
\t\t\tif constexpr (i3 != 0) addend = rm2s19Xor(addend, table[3][i3]);
\t\t\tif constexpr (i4 != 0) addend = rm2s19Xor(addend, finalTable[i4]);
\t\t\treturn addend;
\t\t}}

\t\ttemplate<std::size_t ReversePoint, typename Emit>
\t\tOC_FORCEINLINE void rm2s19EmitPoint(
\t\t\tconst block* input, __m128i* values, const __m128i table[4][16],
\t\t\tconst __m128i finalTable[8], u64 innerBase, Emit& emit) noexcept
\t\t{{
\t\t\tconstexpr std::size_t point = 127 - ReversePoint;
\t\t\tconst auto value = rm2s19Xor(input[point].mData, rm2s19GroupedAAddend<point>(table, finalTable));
\t\t\tvalues[point] = value;
\t\t\temit(innerBase + point, block(value));
\t\t}}

\t\ttemplate<typename Emit, std::size_t... ReversePoints>
\t\tOC_FORCEINLINE void rm2s19EmitPoints(
\t\t\tconst block* input, __m128i* values, const __m128i table[4][16],
\t\t\tconst __m128i finalTable[8], u64 innerBase, Emit& emit,
\t\t\tstd::index_sequence<ReversePoints...>) noexcept
\t\t{{
\t\t\t(rm2s19EmitPoint<ReversePoints>(input, values, table, finalTable, innerBase, emit), ...);
\t\t}}

\t\ttemplate<typename Emit>
\t\tOC_FORCEINLINE void rm2s19AddAGroupedAndEmit(
\t\t\tconst block* input, const __m128i* state, __m128i* values,
\t\t\tu64 innerBase, Emit& emit) noexcept
\t\t{{
\t\t\talignas(32) __m128i groupedState[4][4]{{
\t\t\t\t{{state[0], state[16], state[6], state[2]}},
\t\t\t\t{{state[14], state[18], state[1], state[10]}},
\t\t\t\t{{state[9], state[12], state[3], state[5]}},
\t\t\t\t{{state[17], state[13], state[15], state[8]}},
\t\t\t}};
\t\t\talignas(32) __m128i table[4][16];
\t\t\tfor (unsigned group = 0; group < 4; ++group) rm2s19BuildNibbleTable(groupedState[group], table[group]);
\t\t\talignas(32) __m128i finalState[3]{{state[11], state[7], state[4]}};
\t\t\talignas(32) __m128i finalTable[8];
\t\t\trm2s19BuildThreeBitTable(finalState, finalTable);
\t\t\trm2s19EmitPoints(input, values, table, finalTable, innerBase, emit, std::make_index_sequence<128>{{}});
\t\t}}

\t\ttemplate<typename Emit>
\t\tOC_FORCEINLINE void rm2s19CopyAndEmit(
\t\t\tconst block* input, __m128i* values, u64 innerBase, Emit& emit) noexcept
\t\t{{
\t\t\tfor (u64 offset = 128; offset-- > 0;)
\t\t\t{{
\t\t\t\tconst auto value = input[offset].mData;
\t\t\t\tvalues[offset] = value;
\t\t\t\temit(innerBase + offset, block(value));
\t\t\t}}
\t\t}}
\t}}

\tclass RiffleRm2SubS19Transpose
\t{{
\tpublic:
\t\tstatic constexpr u64 stepBlocks = 128;
\t\tstatic constexpr u64 stateBlocks = 19;

\t\tvoid init(u64 epochs, u64 coefficientSeed)
\t\t{{
\t\t\tif (epochs == 0) throw std::invalid_argument("RM2Sub epoch count must be positive");
\t\t\tmSchedules.resize(epochs);
\t\t\tfor (u64 epoch = 0; epoch < epochs; ++epoch)
\t\t\t{{
\t\t\t\tu32 coefficient;
\t\t\t\tdo coefficient = static_cast<u32>(splitmix64(coefficientSeed)) & 0x7ffffU;
\t\t\t\twhile (coefficient == 0);
\t\t\t\tmSchedules[epoch] = detail::rm2s19MakeFieldTransposeSchedule(coefficient);
\t\t\t}}
\t\t}}

\t\ttemplate<typename Emit>
\t\tOC_FORCEINLINE void emitReverse(
\t\t\tconst block* __restrict input, u64 codeBlocks, Emit&& emit) const noexcept
\t\t{{
\t\t\tconst u64 epochs = codeBlocks / stepBlocks;
\t\t\talignas(32) __m128i state[19]{{}};
\t\t\talignas(32) __m128i syndrome[19];
\t\t\talignas(32) __m128i values[128];
\t\t\tfor (u64 epoch = epochs; epoch-- > 0;)
\t\t\t{{
\t\t\t\tconst block* node = input + epoch * stepBlocks;
\t\t\t\tconst u64 innerBase = epoch * stepBlocks;
\t\t\t\tif (epoch + 1 == epochs) detail::rm2s19CopyAndEmit(node, values, innerBase, emit);
\t\t\t\telse detail::rm2s19AddAGroupedAndEmit(node, state, values, innerBase, emit);
\t\t\t\tif (epoch == 0) break;
\t\t\t\tdetail::rm2s19FinishB(values, syndrome);
\t\t\t\tif (epoch + 1 == epochs) std::memcpy(state, syndrome, sizeof(state));
\t\t\t\telse
\t\t\t\t{{
\t\t\t\t\tdetail::rm2s19FieldMultiplyTranspose(state, mSchedules[epoch]);
\t\t\t\t\tfor (unsigned bit = 0; bit < 19; ++bit) state[bit] = detail::rm2s19Xor(state[bit], syndrome[bit]);
\t\t\t\t}}
\t\t\t}}
\t\t}}

\t\tconst std::vector<detail::Rm2Sub19FieldTransposeSchedule>& schedules() const noexcept
\t\t{{
\t\t\treturn mSchedules;
\t\t}}

\tprivate:
\t\tstd::vector<detail::Rm2Sub19FieldTransposeSchedule> mSchedules;

\t\tstatic u64 splitmix64(u64& state) noexcept
\t\t{{
\t\t\tu64 value = (state += 0x9e3779b97f4a7c15ULL);
\t\t\tvalue = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
\t\t\tvalue = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
\t\t\treturn value ^ (value >> 31);
\t\t}}
\t}};
}}
'''


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    args = parser.parse_args()

    selection = json.loads(args.selection.read_text(encoding="utf-8"))
    if int(selection["parameters"]["state_bits"]) != 19:
        raise ValueError("selection is not state size 19")
    selected = selection["selected"]
    columns = [int(value, 16) for value in selected["B_columns_hex"]]
    masks = [int(value, 16) for value in selected["quadratic_masks_hex"]]
    projection, circuit = projection_code(masks)
    header = make_header(columns, masks, projection)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(header, encoding="utf-8")

    receipt = {
        "schema": "riffle-rm2sub-s19-generated-header-v1",
        "selection": str(args.selection.resolve().relative_to(ROOT)),
        "selected_seed": selected["seed"],
        "a_groups": A_GROUPS,
        "quadratic_projection_seed": 128,
        "quadratic_projection_shared_gates": len(circuit.gates),
        "quadratic_projection_xors": circuit.xor_count,
        "quadratic_projection_transpose_xors": paar.transpose_xor_count(circuit),
        "field_polynomial": "x^19+x^5+x^2+x+1",
        "output": str(args.output.resolve().relative_to(ROOT)),
    }
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(f"quadratic_projection_xors,{circuit.xor_count}")
    print(f"quadratic_projection_shared_gates,{len(circuit.gates)}")
    print(f"wrote,{args.output}")
    print(f"wrote,{args.receipt}")


if __name__ == "__main__":
    main()
