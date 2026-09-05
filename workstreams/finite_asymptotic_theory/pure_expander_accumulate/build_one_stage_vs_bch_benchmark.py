#!/usr/bin/env python3
"""Generate a matched AVX2 benchmark for sparse-EA and optimized BCH.

The sparse path applies the transpose of ``C = A E`` for one sampled
degree-33 map ``E : F_2^256 -> F_2^512``.  The BCH path applies two copies
of the exact ``ExtendedBch256x128-Eq3`` transpose.  Both packed kernels
process two independent constituent calls per AVX2 lane pair.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import build_bch256_128_eq3_circuit as bch  # noqa: E402
import probe_bch_forward_xor_circuit as paar  # noqa: E402


SPARSE_INPUTS = 256
SPARSE_OUTPUTS = 512
SPARSE_DEGREE = 33
SPARSE_MATRIX_SEED = 20260903
SPARSE_OPTIMIZER_SEEDS = range(1, 5)
EXPECTED_SPARSE_ROWS_SHA256 = (
    "42508e8b7d1d74f52ef3b4b9d884655d113fec81408bc095e0908da995ae2fac"
)
EXPECTED_SPARSE_CIRCUIT_SHA256 = (
    "929e57413cf86e25534aa1a817315f4278253567ef906190f32bbddd32bbf62e"
)


def sample_rows() -> list[int]:
    rng = random.Random(SPARSE_MATRIX_SEED)
    rows = []
    for _ in range(SPARSE_OUTPUTS):
        row = 0
        for coordinate in rng.sample(range(SPARSE_INPUTS), SPARSE_DEGREE):
            row |= 1 << coordinate
        rows.append(row)
    return rows


def rows_hash(rows: list[int], width: int) -> str:
    byte_width = (width + 7) // 8
    return hashlib.sha256(
        b"".join(row.to_bytes(byte_width, "little") for row in rows)
    ).hexdigest()


def circuit_hash(circuit: paar.Circuit) -> str:
    encoded = json.dumps(
        {"gates": circuit.gates, "outputs": circuit.output_signals},
        separators=(",", ":"),
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def best_circuit(
    inputs: int, outputs: int, target: list[int], seeds: range
) -> tuple[paar.Circuit, int]:
    paar.DIMENSION = inputs
    paar.LENGTH = outputs
    candidates = [(paar.synthesize(seed, 2, target), seed) for seed in seeds]
    return min(candidates, key=lambda item: (item[0].xor_count, item[1]))


def best_schedule(
    circuit: paar.Circuit, inputs: int, outputs: int, trials: int
) -> tuple[list[int], int, int, int]:
    paar.DIMENSION = inputs
    paar.LENGTH = outputs
    candidates = []
    for seed in range(trials):
        schedule = paar.greedy_transpose_schedule(circuit, seed)
        peak, area = paar.transpose_live_stats(circuit, schedule)
        candidates.append((peak, area, seed, schedule))
    peak, area, seed, schedule = min(candidates)
    return schedule, seed, peak, area


def emit_matrix(name: str, rows: list[int], width: int) -> list[str]:
    chunks = (width + 63) // 64
    lines = [
        f"alignas(64) constexpr std::uint64_t {name}[{len(rows)}][{chunks}] = {{"
    ]
    for row in rows:
        words = [
            (row >> (64 * chunk)) & ((1 << 64) - 1)
            for chunk in range(chunks)
        ]
        lines.append(
            "    {" + ", ".join(f"0x{word:016x}ULL" for word in words) + "},"
        )
    lines.append("};")
    return lines


def emit_packed_transpose(
    name: str,
    type_name: str,
    circuit: paar.Circuit,
    inputs: int,
    outputs: int,
    schedule: list[int],
    accumulate: bool,
) -> list[str]:
    paar.DIMENSION = inputs
    paar.LENGTH = outputs
    terms = paar.transpose_terms(circuit)

    def term_name(term: tuple[str, int]) -> str:
        kind, index = term
        return f"word[{index}]" if kind == "word" else f"a{index}"

    lines = [
        f"struct {type_name} {{ __m256i value; }};",
        f"inline {type_name} operator^({type_name} left, {type_name} right) noexcept",
        "{",
        "    return { _mm256_xor_si256(left.value, right.value) };",
        "}",
        f"NOINLINE void {name}(",
        "    const Block* RESTRICT word0, const Block* RESTRICT word1,",
        "    Block* RESTRICT message0, Block* RESTRICT message1) noexcept",
        "{",
        f"    alignas(64) std::array<{type_name}, {outputs}> word;",
    ]
    if accumulate:
        lines.extend(
            [
                "    __m256i suffix = _mm256_setzero_si256();",
                f"    for (std::size_t index = {outputs}; index-- > 0;)",
                "    {",
                "        const __m256i packed = _mm256_set_m128i(",
                "            word1[index].value, word0[index].value);",
                "        suffix = _mm256_xor_si256(suffix, packed);",
                "        word[index].value = suffix;",
                "    }",
            ]
        )
    else:
        lines.extend(
            [
                f"    for (std::size_t index = 0; index < {outputs}; ++index)",
                "        word[index].value = _mm256_set_m128i(",
                "            word1[index].value, word0[index].value);",
            ]
        )

    for signal in schedule:
        expression = " ^ ".join(term_name(term) for term in terms[signal])
        if signal < inputs:
            lines.extend(
                [
                    f"    const {type_name} output{signal} = {expression};",
                    f"    message0[{signal}].value = _mm256_castsi256_si128(output{signal}.value);",
                    f"    message1[{signal}].value = _mm256_extracti128_si256(output{signal}.value, 1);",
                ]
            )
        else:
            if not expression:
                raise RuntimeError("transpose schedule contains an unused signal")
            lines.append(f"    const {type_name} a{signal} = {expression};")
    lines.append("}")
    return lines


def generate(output: Path, schedule_trials: int) -> dict[str, object]:
    sparse_rows = sample_rows()
    sparse_row_hash = rows_hash(sparse_rows, SPARSE_INPUTS)
    if sparse_row_hash != EXPECTED_SPARSE_ROWS_SHA256:
        raise RuntimeError(f"sparse row hash changed: {sparse_row_hash}")
    sparse_circuit, sparse_seed = best_circuit(
        SPARSE_INPUTS, SPARSE_OUTPUTS, sparse_rows, SPARSE_OPTIMIZER_SEEDS
    )
    sparse_circuit_digest = circuit_hash(sparse_circuit)
    if sparse_circuit_digest != EXPECTED_SPARSE_CIRCUIT_SHA256:
        raise RuntimeError(f"sparse circuit hash changed: {sparse_circuit_digest}")
    sparse_schedule, sparse_schedule_seed, sparse_peak, sparse_area = best_schedule(
        sparse_circuit,
        SPARSE_INPUTS,
        SPARSE_OUTPUTS,
        schedule_trials,
    )
    paar.DIMENSION = SPARSE_INPUTS
    paar.LENGTH = SPARSE_OUTPUTS
    sparse_transpose_xors = paar.transpose_xor_count(sparse_circuit)

    bch_rows = bch.build_rows()
    bch_target = bch.columns(bch_rows)
    bch_circuit, bch_seed = best_circuit(
        bch.DIMENSION, bch.LENGTH, bch_target, range(114, 115)
    )
    bch_schedule, bch_schedule_seed, bch_peak, bch_area = best_schedule(
        bch_circuit, bch.DIMENSION, bch.LENGTH, schedule_trials
    )

    lines = [
        "// Generated by build_one_stage_vs_bch_benchmark.py.",
        "#include <algorithm>",
        "#include <array>",
        "#include <chrono>",
        "#include <cmath>",
        "#include <cstddef>",
        "#include <cstdint>",
        "#include <cstdlib>",
        "#include <cstring>",
        "#include <iomanip>",
        "#include <immintrin.h>",
        "#include <iostream>",
        "#include <numeric>",
        "#include <stdexcept>",
        "#include <string>",
        "#include <vector>",
        "#if defined(__linux__)",
        "#include <pthread.h>",
        "#include <sched.h>",
        "#endif",
        "",
        "#if defined(_MSC_VER)",
        "#define NOINLINE __declspec(noinline)",
        "#define RESTRICT __restrict",
        "#else",
        "#define NOINLINE __attribute__((noinline))",
        "#define RESTRICT __restrict__",
        "#endif",
        "",
        "struct alignas(16) Block { __m128i value; };",
        "",
        *emit_matrix("sparseRows", sparse_rows, SPARSE_INPUTS),
        "",
        *emit_matrix("bchRows", bch_rows, bch.LENGTH),
        "",
        *emit_packed_transpose(
            "sparseEaTranspose2",
            "SparsePackedBlock",
            sparse_circuit,
            SPARSE_INPUTS,
            SPARSE_OUTPUTS,
            sparse_schedule,
            True,
        ),
        "",
        *emit_packed_transpose(
            "bchTranspose2",
            "BchPackedBlock",
            bch_circuit,
            bch.DIMENSION,
            bch.LENGTH,
            bch_schedule,
            False,
        ),
        "",
        r'''std::uint64_t splitmix64(std::uint64_t& state) noexcept
{
    std::uint64_t value = (state += 0x9e3779b97f4a7c15ULL);
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

Block randomBlock(std::uint64_t& state) noexcept
{
    const auto low = splitmix64(state);
    const auto high = splitmix64(state);
    return { _mm_set_epi64x(static_cast<long long>(high), static_cast<long long>(low)) };
}

Block xorBlock(Block left, Block right) noexcept
{
    return { _mm_xor_si128(left.value, right.value) };
}

bool equalBlock(Block left, Block right) noexcept
{
    const auto different = _mm_xor_si128(left.value, right.value);
    return _mm_test_all_zeros(different, different) != 0;
}

void sparseReference(const Block* word, Block* message) noexcept
{
    std::array<Block, 512> suffix;
    Block running{ _mm_setzero_si128() };
    for (std::size_t index = 512; index-- > 0;)
    {
        running = xorBlock(running, word[index]);
        suffix[index] = running;
    }
    for (std::size_t column = 0; column < 256; ++column)
    {
        Block value{ _mm_setzero_si128() };
        for (std::size_t row = 0; row < 512; ++row)
            if ((sparseRows[row][column / 64] >> (column % 64)) & 1ULL)
                value = xorBlock(value, suffix[row]);
        message[column] = value;
    }
}

void bchReference(const Block* word, Block* message) noexcept
{
    for (std::size_t row = 0; row < 128; ++row)
    {
        Block value{ _mm_setzero_si128() };
        for (std::size_t coordinate = 0; coordinate < 256; ++coordinate)
            if ((bchRows[row][coordinate / 64] >> (coordinate % 64)) & 1ULL)
                value = xorBlock(value, word[coordinate]);
        message[row] = value;
    }
}

void selfTest()
{
    std::uint64_t state = 0x5350494e2d45412dULL;
    std::array<Block, 512> sparseWord0;
    std::array<Block, 512> sparseWord1;
    std::array<Block, 256> sparseActual0;
    std::array<Block, 256> sparseActual1;
    std::array<Block, 256> sparseExpected0;
    std::array<Block, 256> sparseExpected1;
    for (auto& value : sparseWord0) value = randomBlock(state);
    for (auto& value : sparseWord1) value = randomBlock(state);
    sparseEaTranspose2(
        sparseWord0.data(), sparseWord1.data(),
        sparseActual0.data(), sparseActual1.data());
    sparseReference(sparseWord0.data(), sparseExpected0.data());
    sparseReference(sparseWord1.data(), sparseExpected1.data());
    for (std::size_t index = 0; index < 256; ++index)
        if (!equalBlock(sparseActual0[index], sparseExpected0[index]) ||
            !equalBlock(sparseActual1[index], sparseExpected1[index]))
            throw std::runtime_error("sparse-EA self-test failed");

    std::array<Block, 256> bchWord0;
    std::array<Block, 256> bchWord1;
    std::array<Block, 128> bchActual0;
    std::array<Block, 128> bchActual1;
    std::array<Block, 128> bchExpected0;
    std::array<Block, 128> bchExpected1;
    for (auto& value : bchWord0) value = randomBlock(state);
    for (auto& value : bchWord1) value = randomBlock(state);
    bchTranspose2(
        bchWord0.data(), bchWord1.data(), bchActual0.data(), bchActual1.data());
    bchReference(bchWord0.data(), bchExpected0.data());
    bchReference(bchWord1.data(), bchExpected1.data());
    for (std::size_t index = 0; index < 128; ++index)
        if (!equalBlock(bchActual0[index], bchExpected0[index]) ||
            !equalBlock(bchActual1[index], bchExpected1[index]))
            throw std::runtime_error("BCH self-test failed");
}

NOINLINE void runSparse(const Block* RESTRICT word, Block* RESTRICT message) noexcept
{
    constexpr std::size_t constituents = (std::size_t{1} << 21) / 512;
    for (std::size_t outer = 0; outer < constituents; outer += 2)
        sparseEaTranspose2(
            word + outer * 512,
            word + (outer + 1) * 512,
            message + outer * 256,
            message + (outer + 1) * 256);
}

NOINLINE void runBch(const Block* RESTRICT word, Block* RESTRICT message) noexcept
{
    constexpr std::size_t constituents = (std::size_t{1} << 21) / 256;
    for (std::size_t outer = 0; outer < constituents; outer += 2)
        bchTranspose2(
            word + outer * 256,
            word + (outer + 1) * 256,
            message + outer * 128,
            message + (outer + 1) * 128);
}

std::uint64_t checksum(const std::vector<Block>& values) noexcept
{
    __m128i sum = _mm_setzero_si128();
    for (const auto value : values) sum = _mm_xor_si128(sum, value.value);
    alignas(16) std::uint64_t words[2];
    _mm_store_si128(reinterpret_cast<__m128i*>(words), sum);
    return words[0] ^ words[1];
}

void pinThread(unsigned cpu)
{
#if defined(__linux__)
    cpu_set_t set;
    CPU_ZERO(&set);
    CPU_SET(cpu, &set);
    if (pthread_setaffinity_np(pthread_self(), sizeof(set), &set) != 0)
        throw std::runtime_error("failed to pin benchmark thread");
#else
    static_cast<void>(cpu);
#endif
}

int main(int argc, char** argv)
{
    if (argc != 5)
    {
        std::cerr << "usage: bench sparse|bch trials warmups cpu\n";
        return 2;
    }
    const std::string variant = argv[1];
    const auto trials = static_cast<unsigned>(std::stoul(argv[2]));
    const auto warmups = static_cast<unsigned>(std::stoul(argv[3]));
    const auto cpu = static_cast<unsigned>(std::stoul(argv[4]));
    if (variant != "sparse" && variant != "bch")
        throw std::invalid_argument("variant must be sparse or bch");
    pinThread(cpu);
    selfTest();

    constexpr std::size_t inputBlocks = std::size_t{1} << 21;
    constexpr std::size_t outputBlocks = std::size_t{1} << 20;
    std::vector<Block> input(inputBlocks);
    std::vector<Block> output(outputBlocks);
    std::uint64_t state = 0x70656163682d3739ULL;
    for (auto& value : input) value = randomBlock(state);

    const auto run = [&] {
        if (variant == "sparse") runSparse(input.data(), output.data());
        else runBch(input.data(), output.data());
    };
    for (unsigned index = 0; index < warmups; ++index) run();

    std::vector<double> samples;
    samples.reserve(trials);
    std::uint64_t sink = 0;
    for (unsigned trial = 0; trial < trials; ++trial)
    {
        input[trial % input.size()].value = _mm_xor_si128(
            input[trial % input.size()].value,
            _mm_set_epi64x(0, static_cast<long long>(trial + 1)));
        std::atomic_signal_fence(std::memory_order_seq_cst);
        const auto begin = std::chrono::steady_clock::now();
        run();
        const auto end = std::chrono::steady_clock::now();
        std::atomic_signal_fence(std::memory_order_seq_cst);
        samples.push_back(std::chrono::duration<double, std::milli>(end - begin).count());
        sink ^= checksum(output);
    }
    std::sort(samples.begin(), samples.end());
    const double median = samples[samples.size() / 2];
    const double mean = std::accumulate(samples.begin(), samples.end(), 0.0) / samples.size();
    const double variance = std::accumulate(
        samples.begin(), samples.end(), 0.0,
        [mean](double sum, double value) {
            const double difference = value - mean;
            return sum + difference * difference;
        }) / samples.size();
    std::cout << std::setprecision(17)
        << "{\n"
        << "  \"schema\": \"one-stage-sparse-ea-vs-bch-benchmark-v1\",\n"
        << "  \"variant\": \"" << variant << "\",\n"
        << "  \"input_blocks\": " << inputBlocks << ",\n"
        << "  \"output_blocks\": " << outputBlocks << ",\n"
        << "  \"trials\": " << trials << ",\n"
        << "  \"warmups\": " << warmups << ",\n"
        << "  \"cpu\": " << cpu << ",\n"
        << "  \"median_ms\": " << median << ",\n"
        << "  \"mean_ms\": " << mean << ",\n"
        << "  \"stddev_ms\": " << std::sqrt(variance) << ",\n"
        << "  \"minimum_ms\": " << samples.front() << ",\n"
        << "  \"maximum_ms\": " << samples.back() << ",\n"
        << "  \"sink\": " << sink << ",\n"
        << "  \"samples_ms\": [";
    for (std::size_t index = 0; index < samples.size(); ++index)
        std::cout << (index ? ", " : "") << samples[index];
    std::cout << "]\n}\n";
    return 0;
}
''',
    ]

    # atomic_signal_fence needs <atomic>; keeping this include next to the
    # generated source list makes the dependency explicit.
    lines.insert(2, "#include <atomic>")
    output.write_text("\n".join(lines), encoding="utf-8")
    return {
        "schema": "one-stage-sparse-ea-vs-bch-build-v1",
        "status": "EXACT_CIRCUITS_VERIFIED_BY_GENERATORS",
        "sparse": {
            "rows_sha256": sparse_row_hash,
            "circuit_sha256": sparse_circuit_digest,
            "optimizer_seed": sparse_seed,
            "forward_xors": sparse_circuit.xor_count,
            "transpose_xors_before_accumulator": sparse_transpose_xors,
            "transpose_xors_with_accumulator": sparse_transpose_xors + 511,
            "schedule_seed": sparse_schedule_seed,
            "schedule_peak_live": sparse_peak,
            "schedule_live_area": sparse_area,
        },
        "bch": {
            "rows_sha256": bch.rows_sha256(bch_rows),
            "optimizer_seed": bch_seed,
            "forward_xors": bch_circuit.xor_count,
            "transpose_xors_per_constituent": paar.transpose_xor_count(
                bch_circuit
            ),
            "transpose_xors_per_512_to_256_capacity": 2
            * paar.transpose_xor_count(bch_circuit),
            "schedule_seed": bch_schedule_seed,
            "schedule_peak_live": bch_peak,
            "schedule_live_area": bch_area,
        },
        "benchmark": {
            "input_blocks": 1 << 21,
            "output_blocks": 1 << 20,
            "packing": "two independent constituent calls per AVX2 vector",
            "generated_source": str(output),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=Path,
        default=HERE / "one_stage_sparse_ea_vs_bch_bench.cpp",
    )
    parser.add_argument(
        "--receipt",
        type=Path,
        default=HERE / "one_stage_sparse_ea_vs_bch_build.json",
    )
    parser.add_argument("--schedule-trials", type=int, default=16)
    args = parser.parse_args()
    payload = generate(args.source, args.schedule_trials)
    args.receipt.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
