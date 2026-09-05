#!/usr/bin/env python3
"""Generate the exact ExtendedBch256x128-Eq3 transposed XOR circuit.

The parent is the parity extension of the primitive narrow-sense BCH(255,37)
code.  It has parameters [256,131,>=38].  Intersecting it with
x_0+x_1=x_0+x_2=x_0+x_3=0 gives the deterministic [256,128,>=38] subcode
used by the first end-to-end FieldCheckpoint implementation.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import probe_bch_forward_xor_circuit as paar
from check_bch_transposed_inner import extended_bch_check_rows, row_rank
from exact_bch_spectrum_small import nullspace_basis


DIMENSION = 128
LENGTH = 256
M = 8
DESIGNED_DISTANCE = 37
EXTRA_CHECK_COORDINATES = (1, 2, 3)


def build_rows() -> list[int]:
    parent, length, dimension, rank = extended_bch_check_rows(
        M, DESIGNED_DISTANCE
    )
    if (length, dimension, rank) != (256, 131, 125):
        raise RuntimeError(
            f"unexpected extended BCH parent {(length, dimension, rank)}"
        )
    extra = [(1 << 0) | (1 << coordinate) for coordinate in EXTRA_CHECK_COORDINATES]
    ladder = [row_rank(parent + extra[:count]) for count in range(4)]
    if ladder != [125, 126, 127, 128]:
        raise RuntimeError(f"extra checks are not independent: {ladder}")
    rows = list(nullspace_basis(parent + extra, LENGTH))
    if len(rows) != DIMENSION or row_rank(rows) != DIMENSION:
        raise RuntimeError("generator basis has the wrong rank")
    return rows


def columns(rows: list[int]) -> list[int]:
    result = [0] * LENGTH
    for row_index, row in enumerate(rows):
        while row:
            bit = (row & -row).bit_length() - 1
            result[bit] |= 1 << row_index
            row &= row - 1
    return result


def rows_sha256(rows: list[int]) -> str:
    return hashlib.sha256(
        b"".join(row.to_bytes(LENGTH // 8, "little") for row in rows)
    ).hexdigest()


def write_wrapper(path: Path, rows: list[int], transpose_xors: int) -> None:
    row_lines = []
    for row in rows:
        words = [(row >> (64 * word)) & ((1 << 64) - 1) for word in range(4)]
        row_lines.append(
            "\t\t\t{{" + ", ".join(f"0x{word:016x}ULL" for word in words) + "}},"
        )
    text = "\n".join(
        [
            "#pragma once",
            "",
            '#include "ExtendedBch256x128Eq3TransposeCircuit.h"',
            '#include "ExtendedBch256x128Eq3TransposeCircuit2.h"',
            "#if defined(LIBOTE_RIFFLE_EXPERIMENTAL_AVX512_OUTER) && defined(__AVX512F__)",
            '#include "ExtendedBch256x128Eq3TransposeCircuit4.h"',
            "#endif",
            "",
            "#include <array>",
            "#include <bit>",
            "#include <cryptoTools/Common/Defines.h>",
            "",
            "namespace osuCrypto",
            "{",
            "\t// Extended BCH(255,37), parity extended, then intersected with",
            "\t// x_0+x_1=x_0+x_2=x_0+x_3=0. Exact [256,128,>=38].",
            "\tstruct ExtendedBch256x128Eq3",
            "\t{",
            "\t\tstatic constexpr u64 dimension = 128;",
            "\t\tstatic constexpr u64 length = 256;",
            f"\t\tstatic constexpr u64 transposedEncoderXorCount = {transpose_xors};",
            f'\t\tstatic constexpr const char* generatorRowsSha256 = "{rows_sha256(rows)}";',
            "\t\tstatic constexpr std::array<std::array<u64, 4>, dimension> generatorRows = {{",
            *row_lines,
            "\t\t}};",
            "",
            "\t\tstatic void transposeBlock(const block* __restrict word, block* __restrict message)",
            "\t\t{",
            "\t\t\tdetail::extendedBch256x128Eq3TransposeCircuit(word, message);",
            "\t\t}",
            "",
            "\t\tstatic void transposeBlock2(",
            "\t\t\tconst block* __restrict word0,",
            "\t\t\tconst block* __restrict word1,",
            "\t\t\tblock* __restrict message0,",
            "\t\t\tblock* __restrict message1)",
            "\t\t{",
            "\t\t\tdetail::extendedBch256x128Eq3TransposeCircuit2(",
            "\t\t\t\tword0, word1, message0, message1);",
            "\t\t}",
            "",
            "#if defined(LIBOTE_RIFFLE_EXPERIMENTAL_AVX512_OUTER) && defined(__AVX512F__)",
            "\t\tstatic void transposeBlock4(",
            "\t\t\tconst block* __restrict word0,",
            "\t\t\tconst block* __restrict word1,",
            "\t\t\tconst block* __restrict word2,",
            "\t\t\tconst block* __restrict word3,",
            "\t\t\tblock* __restrict message0,",
            "\t\t\tblock* __restrict message1,",
            "\t\t\tblock* __restrict message2,",
            "\t\t\tblock* __restrict message3)",
            "\t\t{",
            "\t\t\tdetail::extendedBch256x128Eq3TransposeCircuit4(",
            "\t\t\t\tword0, word1, word2, word3,",
            "\t\t\t\tmessage0, message1, message2, message3);",
            "\t\t}",
            "#endif",
            "",
            "\t\ttemplate<typename T>",
            "\t\tstatic void transposeReference(const T* __restrict word, T* __restrict message)",
            "\t\t{",
            "\t\t\tfor (u64 row = 0; row < dimension; ++row)",
            "\t\t\t{",
            "\t\t\t\tmessage[row] = T{};",
            "\t\t\t\tfor (u64 chunk = 0; chunk < 4; ++chunk)",
            "\t\t\t\t{",
            "\t\t\t\t\tauto mask = generatorRows[row][chunk];",
            "\t\t\t\t\twhile (mask)",
            "\t\t\t\t\t{",
            "\t\t\t\t\t\tconst auto bit = static_cast<u64>(std::countr_zero(mask));",
            "\t\t\t\t\t\tmessage[row] ^= word[64 * chunk + bit];",
            "\t\t\t\t\t\tmask &= mask - 1;",
            "\t\t\t\t\t}",
            "\t\t\t\t}",
            "\t\t\t}",
            "\t\t}",
            "\t};",
            "}",
            "",
        ]
    )
    path.write_text(text, encoding="utf-8")


def write_packed_transpose_header(
    path: Path,
    circuit: paar.Circuit,
    seed: int,
    schedule: list[int] | None = None,
    namespace: str = "osuCrypto::detail",
) -> None:
    terms = paar.transpose_terms(circuit)
    schedule = paar.descending_transpose_schedule(circuit) if schedule is None else schedule

    def term_name(term: tuple[str, int]) -> str:
        kind, index = term
        return f"word[{index}]" if kind == "word" else f"a{index}"

    lines = [
        "#pragma once",
        "",
        "// Generated by scripts/build_bch256_128_eq3_circuit.py.",
        f"// Two-way AVX2 packing of the exact seed-{seed} transposed circuit.",
        "",
        "#include <immintrin.h>",
        "#include <cryptoTools/Common/Defines.h>",
        "",
        f"namespace {namespace}",
        "{",
        "\tstruct ExtendedBch256x128Eq3PackedBlock",
        "\t{",
        "\t\t__m256i value;",
        "\t};",
        "",
        "\tOC_FORCEINLINE ExtendedBch256x128Eq3PackedBlock operator^(",
        "\t\tExtendedBch256x128Eq3PackedBlock left,",
        "\t\tExtendedBch256x128Eq3PackedBlock right)",
        "\t{",
        "\t\treturn { _mm256_xor_si256(left.value, right.value) };",
        "\t}",
        "",
        "#if defined(_MSC_VER)",
        "#define LIBOTE_RIFFLE_NOINLINE __declspec(noinline)",
        "#elif defined(__GNUC__) || defined(__clang__)",
        "#define LIBOTE_RIFFLE_NOINLINE __attribute__((noinline))",
        "#else",
        "#define LIBOTE_RIFFLE_NOINLINE",
        "#endif",
        "\tLIBOTE_RIFFLE_NOINLINE inline void extendedBch256x128Eq3TransposeCircuit2(",
        "\t\tconst block* __restrict word0,",
        "\t\tconst block* __restrict word1,",
        "\t\tblock* __restrict message0,",
        "\t\tblock* __restrict message1)",
        "\t{",
        "\t\tExtendedBch256x128Eq3PackedBlock word[256];",
        "\t\tfor (u64 index = 0; index < 256; ++index)",
        "\t\t\tword[index].value = _mm256_set_m128i(word1[index].mData, word0[index].mData);",
    ]
    for signal in schedule:
        expression = " ^ ".join(term_name(term) for term in terms[signal])
        if signal < DIMENSION:
            lines.extend(
                [
                    f"\t\tconst ExtendedBch256x128Eq3PackedBlock output{signal} = {expression};",
                    f"\t\tmessage0[{signal}].mData = _mm256_castsi256_si128(output{signal}.value);",
                    f"\t\tmessage1[{signal}].mData = _mm256_extracti128_si256(output{signal}.value, 1);",
                ]
            )
        else:
            if not expression:
                raise ValueError("packed transpose contains an unused internal signal")
            lines.append(
                f"\t\tconst ExtendedBch256x128Eq3PackedBlock a{signal} = {expression};"
            )
    lines.extend(["\t}", "#undef LIBOTE_RIFFLE_NOINLINE", "}", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def write_packed4_transpose_header(
    path: Path,
    circuit: paar.Circuit,
    seed: int,
    schedule: list[int] | None = None,
    namespace: str = "osuCrypto::detail",
) -> None:
    terms = paar.transpose_terms(circuit)
    schedule = paar.descending_transpose_schedule(circuit) if schedule is None else schedule

    def term_name(term: tuple[str, int]) -> str:
        kind, index = term
        return f"word[{index}]" if kind == "word" else f"a{index}"

    lines = [
        "#pragma once",
        "",
        "// Generated by scripts/build_bch256_128_eq3_circuit.py.",
        f"// Four-way AVX-512 packing of the exact seed-{seed} transposed circuit.",
        "",
        "#include <immintrin.h>",
        "#include <cryptoTools/Common/Defines.h>",
        "",
        f"namespace {namespace}",
        "{",
        "\tstruct ExtendedBch256x128Eq3PackedBlock4",
        "\t{",
        "\t\t__m512i value;",
        "\t};",
        "",
        "\tOC_FORCEINLINE ExtendedBch256x128Eq3PackedBlock4 operator^(",
        "\t\tExtendedBch256x128Eq3PackedBlock4 left,",
        "\t\tExtendedBch256x128Eq3PackedBlock4 right)",
        "\t{",
        "\t\treturn { _mm512_xor_si512(left.value, right.value) };",
        "\t}",
        "",
        "#if defined(_MSC_VER)",
        "#define LIBOTE_RIFFLE_NOINLINE4 __declspec(noinline)",
        "#elif defined(__GNUC__) || defined(__clang__)",
        "#define LIBOTE_RIFFLE_NOINLINE4 __attribute__((noinline))",
        "#else",
        "#define LIBOTE_RIFFLE_NOINLINE4",
        "#endif",
        "\tLIBOTE_RIFFLE_NOINLINE4 inline void extendedBch256x128Eq3TransposeCircuit4(",
        "\t\tconst block* __restrict word0,",
        "\t\tconst block* __restrict word1,",
        "\t\tconst block* __restrict word2,",
        "\t\tconst block* __restrict word3,",
        "\t\tblock* __restrict message0,",
        "\t\tblock* __restrict message1,",
        "\t\tblock* __restrict message2,",
        "\t\tblock* __restrict message3)",
        "\t{",
        "\t\tExtendedBch256x128Eq3PackedBlock4 word[256];",
        "\t\tfor (u64 index = 0; index < 256; ++index)",
        "\t\t{",
        "\t\t\tconst auto low = _mm256_set_m128i(word1[index].mData, word0[index].mData);",
        "\t\t\tconst auto high = _mm256_set_m128i(word3[index].mData, word2[index].mData);",
        "\t\t\tword[index].value = _mm512_inserti64x4(",
        "\t\t\t\t_mm512_castsi256_si512(low), high, 1);",
        "\t\t}",
    ]
    for signal in schedule:
        expression = " ^ ".join(term_name(term) for term in terms[signal])
        if signal < DIMENSION:
            lines.extend(
                [
                    f"\t\tconst ExtendedBch256x128Eq3PackedBlock4 output{signal} = {expression};",
                    f"\t\tmessage0[{signal}].mData = _mm512_extracti32x4_epi32(output{signal}.value, 0);",
                    f"\t\tmessage1[{signal}].mData = _mm512_extracti32x4_epi32(output{signal}.value, 1);",
                    f"\t\tmessage2[{signal}].mData = _mm512_extracti32x4_epi32(output{signal}.value, 2);",
                    f"\t\tmessage3[{signal}].mData = _mm512_extracti32x4_epi32(output{signal}.value, 3);",
                ]
            )
        else:
            if not expression:
                raise ValueError("packed transpose contains an unused internal signal")
            lines.append(
                f"\t\tconst ExtendedBch256x128Eq3PackedBlock4 a{signal} = {expression};"
            )
    lines.extend(["\t}", "#undef LIBOTE_RIFFLE_NOINLINE4", "}", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    # Seed 114 was the best point in the recorded 128-seed search.
    parser.add_argument("--trials", type=int, default=1)
    parser.add_argument("--seed", type=int, default=114)
    parser.add_argument("--circuit-header", type=Path, required=True)
    parser.add_argument("--packed-circuit-header", type=Path, required=True)
    parser.add_argument("--packed4-circuit-header", type=Path)
    parser.add_argument("--wrapper-header", type=Path, required=True)
    parser.add_argument(
        "--transpose-schedule",
        choices=("descending", "greedy"),
        default="descending",
    )
    parser.add_argument("--schedule-seed", type=int, default=0)
    parser.add_argument("--schedule-trials", type=int, default=1)
    parser.add_argument("--packed-namespace", default="osuCrypto::detail")
    args = parser.parse_args()

    rows = build_rows()
    target = columns(rows)
    paar.DIMENSION = DIMENSION
    paar.LENGTH = LENGTH
    best = None
    best_seed = None
    for seed in range(args.seed, args.seed + args.trials):
        circuit = paar.synthesize(seed, 2, target)
        if best is None or paar.transpose_xor_count(circuit) < paar.transpose_xor_count(best):
            best = circuit
            best_seed = seed
    assert best is not None and best_seed is not None
    transpose_xors = paar.transpose_xor_count(best)
    descending_schedule = paar.descending_transpose_schedule(best)
    descending_peak, descending_area = paar.transpose_live_stats(
        best, descending_schedule
    )
    packed_schedule = descending_schedule
    packed_schedule_seed = None
    if args.transpose_schedule == "greedy":
        schedule_candidates = (
            (schedule_seed, paar.greedy_transpose_schedule(best, schedule_seed))
            for schedule_seed in range(
                args.schedule_seed, args.schedule_seed + args.schedule_trials
            )
        )
        packed_schedule_seed, packed_schedule = min(
            schedule_candidates,
            key=lambda item: paar.transpose_live_stats(best, item[1]),
        )
    packed_peak, packed_area = paar.transpose_live_stats(best, packed_schedule)
    paar.write_production_transpose_header(
        args.circuit_header,
        best,
        best_seed,
        2,
        "extendedBch256x128Eq3TransposeCircuit",
        False,
    )
    write_packed_transpose_header(
        args.packed_circuit_header,
        best,
        best_seed,
        packed_schedule,
        args.packed_namespace,
    )
    if args.packed4_circuit_header is not None:
        write_packed4_transpose_header(
            args.packed4_circuit_header,
            best,
            best_seed,
            packed_schedule,
            args.packed_namespace,
        )
    write_wrapper(args.wrapper_header, rows, transpose_xors)
    print(f"generator_rows_sha256,{rows_sha256(rows)}")
    print(f"best_seed,{best_seed}")
    print(f"forward_xors,{best.xor_count}")
    print(f"transpose_xors,{transpose_xors}")
    print(f"descending_transpose_peak_live,{descending_peak}")
    print(f"descending_transpose_live_area,{descending_area}")
    print(f"packed_transpose_schedule,{args.transpose_schedule}")
    if packed_schedule_seed is not None:
        print(f"packed_transpose_schedule_seed,{packed_schedule_seed}")
    print(f"packed_transpose_peak_live,{packed_peak}")
    print(f"packed_transpose_live_area,{packed_area}")
    print(f"circuit_header,{args.circuit_header}")
    print(f"packed_circuit_header,{args.packed_circuit_header}")
    if args.packed4_circuit_header is not None:
        print(f"packed4_circuit_header,{args.packed4_circuit_header}")
    print(f"wrapper_header,{args.wrapper_header}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
