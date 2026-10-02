"""Emit isolated, exact-map compact-GL32/BCH scheduling experiments.

Usage: python outer_schedule_codegen.py BchCircuit.h --mode 10 > Outer10.cpp

This translation unit exports all three bchPackedCoeff entry points plus
packedCoeffTileMode.  Original and Aligned retain the established tile-mode-1
implementation; only Compact selects the new schedule.  The translation unit
replaces PackedCoeff1.o, not PackedMixer.o.  Modes 10--12 preserve the
canonical sampled GL32 coefficient layout and all BCH coordinates.  They
change only BCH register blocking, calls, and intermediate storage.

10: established 2-output/8-plane paired schedule, force-inlined tile calls.
11: 4-output/4-plane paired schedule, packed-output buffer and wide stores.
12: 4-output/4-plane single-input schedule, the same packed-output buffer.

The latter schedules halve dense-input vector loads, but double matrix
broadcasts and introduce an 8-KiB store/read pass.  They are hypotheses for
measurement, not asserted improvements.  Mode 12 lowers coefficient register
pressure relative to mode 11 at the cost of separate XOR accumulations.
"""

import argparse
import contextlib
import io
from pathlib import Path
import re
import subprocess
import sys

import packed_bch_tune_codegen as tune
import packed_coeff_codegen as coeff
import packed_mixer_codegen as mixer


MODES = {
    10: "force-inlined established 2-output/8-plane paired tiles",
    11: "4-output/4-plane paired tiles, packed-output buffer",
    12: "4-output/4-plane single-input tiles, packed-output buffer",
}


def verify_schedule(header):
    """Assert complete, duplicate-free GFNI and packed-output coverage."""
    tune.verify()
    for paired in (False, True):
        coverage = [[[0] * 16 for _ in range(8)] for _ in range(16)]
        stores = [[0] * 8 for _ in range(16)]
        for output in range(0, 16, 4):
            for plane in range(0, 8, 4):
                for p in range(4):
                    for j in range(4):
                        stores[output+p][plane+j] += 1
                step = 2 if paired else 1
                for inp in range(0, 16, step):
                    for p in range(4):
                        for j in range(4):
                            for i in range(inp, inp+step):
                                coverage[output+p][plane+j][i] += 1
        assert all(n == 1 for row in stores for n in row)
        assert all(n == 1 for row in coverage for plane in row for n in plane)
    # Independently check the blocked tile's sparse/dense/parity formula on
    # all 256 coordinate basis vectors, against the source generator rows.
    words = [int(x, 16) for x in re.findall(r"0x([0-9a-f]+)ULL", header.read_text())]
    assert len(words) == 512
    rows = [sum(words[4*i+j] << (64*j) for j in range(4)) for i in range(128)]
    matrices = [[sum(((rows[8*o+j] >> (128+8*i)) & 255) << (8*(7-j))
                     for j in range(8)) for i in range(16)] for o in range(16)]
    parity = [sum(((rows[8*o+j] >> 127) & 1) << j for j in range(8))
              for o in range(16)]

    def gfni(value, matrix):
        return sum(((value & (matrix >> (8*(7-j))) & 255).bit_count() & 1) << j
                   for j in range(8))

    for coordinate in range(256):
        octets = [((1 << coordinate) >> (8*i)) & 255 for i in range(32)]
        for output in range(16):
            value = octets[output] & (127 if output == 15 else 255)
            value ^= parity[output] if octets[15] & 128 else 0
            for inp in range(16):
                value ^= gfni(octets[16+inp], matrices[output][inp])
            expected = sum(((rows[8*output+j] >> coordinate) & 1) << j for j in range(8))
            assert value == expected


def retained_source(header):
    result = subprocess.run(
        [sys.executable, str(Path(mixer.__file__).resolve()),
         str(header)], check=True, capture_output=True, text=True)
    if result.stderr:
        print(result.stderr, file=sys.stderr, end="")
    return result.stdout


def emit_inline_tile():
    # Reuse the exact selected arithmetic/output generator.  The only change
    # to this definition is SPIN_NOINLINE -> SPIN_FORCEINLINE.
    captured = io.StringIO()
    with contextlib.redirect_stdout(captured):
        tune.emit_tile("packedScheduleInline", 2, True, False)
    text = captured.getvalue()
    assert text.count("SPIN_NOINLINE") == 1
    print(text.replace("SPIN_NOINLINE", "SPIN_FORCEINLINE"), end="")


def emit_blocked_tile(paired):
    print("template<unsigned output,unsigned plane> static SPIN_NOINLINE void "
          "packedScheduleBlocked(const __m512i* __restrict src,"
          "__m512i* __restrict packed) {")
    print("static_assert(output+4<=16 && plane+4<=8);")
    print("static_assert(output%4==0 && plane%4==0);")
    # Parity coordinate 127 remains separate.  Systematic coordinate 127
    # must be cleared; its generator column is not the 128th identity column.
    for p in range(4):
        print(f"constexpr auto parity{p}=[] {{std::uint64_t m=0;"
              f"for(unsigned j=0;j<8;++j)if(BchRows[8*(output+{p})+j][1]>>63)"
              "m|=std::uint64_t(0x80)<<(8*(7-j));return m;}();")
        for j in range(4):
            index = 4*p+j
            print(f"auto y{index}=_mm512_load_si512(src+8*(output+{p})+plane+{j});")
            print(f"if constexpr(output+{p}==15)y{index}="
                  f"_mm512_and_si512(y{index},_mm512_set1_epi8(0x7f));")
            print(f"y{index}=_mm512_xor_si512(y{index},"
                  f"_mm512_gf2p8affine_epi64_epi8("
                  f"_mm512_load_si512(src+120+plane+{j}),"
                  f"_mm512_set1_epi64(parity{p}),0));")
    step = 2 if paired else 1
    print(f"for(unsigned input=0;input<16;input+={step}) {{")
    for p in range(4):
        print(f"const auto m{p}=_mm512_set1_epi64(matrices[output+{p}][input]);")
        if paired:
            print(f"const auto n{p}=_mm512_set1_epi64(matrices[output+{p}][input+1]);")
    for j in range(4):
        print("{")
        print(f"const auto x=_mm512_load_si512(src+128+8*input+plane+{j});")
        if paired:
            print(f"const auto z=_mm512_load_si512(src+136+8*input+plane+{j});")
        for p in range(4):
            index = 4*p+j
            term = f"_mm512_gf2p8affine_epi64_epi8(x,m{p},0)"
            if paired:
                term2 = f"_mm512_gf2p8affine_epi64_epi8(z,n{p},0)"
                print(f"y{index}=_mm512_ternarylogic_epi64(y{index},{term},{term2},0x96);")
            else:
                print(f"y{index}=_mm512_xor_si512(y{index},{term});")
        print("}")
    print("}")
    for p in range(4):
        for j in range(4):
            print(f"_mm512_store_si512(packed+8*(output+{p})+plane+{j},y{4*p+j});")
    print("}")


def emit_writeback():
    print("static SPIN_FORCEINLINE void packedScheduleWriteback("
          "const __m512i* __restrict packed,block* __restrict out) {")
    print("for(unsigned output=0;output<16;output+=2) {")
    for p in range(2):
        for j in range(8):
            print(f"auto y{8*p+j}=_mm512_load_si512(packed+8*(output+{p})+{j});")
    tune.emit_output(2)
    print("}\n}")


def emit(header, mode):
    source = retained_source(header)
    print('#include "Spin.h"\n#include "generated/BchCircuit.h"')
    print("namespace spin::detail::kernel { namespace outer_schedule_detail {")
    print(coeff.definition(source, "alignas(64) static constexpr std::uint64_t matrices", True))
    print(coeff.definition(source, "static SPIN_FORCEINLINE void orthoBlend("))
    # These two unchanged variants let the existing basis/check harness keep
    # checking all coefficient representations.  Their selected BCH schedule
    # is exactly packed_coeff_codegen.py --tile-mode 1.
    tune.emit_tile("packedCoeffTile", 2, True, False)
    if mode == 10:
        emit_inline_tile()
    else:
        emit_blocked_tile(paired=mode == 11)
        emit_writeback()
    print("}")
    for storage in ("Original", "Aligned"):
        print(f"SPIN_NOINLINE void bchPackedCoeff{storage}(const block* __restrict a,"
              "block* __restrict out,const std::uint64_t* __restrict coeff) {")
        print("using namespace outer_schedule_detail;")
        coeff.emit_preparation(storage)
        for output in range(0, 16, 2):
            print(f"packedCoeffTile<{output}>(src,out);")
        print("}")
    print("SPIN_NOINLINE void bchPackedCoeffCompact(const block* __restrict a,"
          "block* __restrict out,const std::uint64_t* __restrict coeff) {")
    print("using namespace outer_schedule_detail;")
    coeff.emit_preparation("Compact")
    if mode == 10:
        for output in range(0, 16, 2):
            print(f"packedScheduleInline<{output}>(src,out);")
    else:
        print("alignas(64) __m512i packed[128];")
        for output in range(0, 16, 4):
            for plane in range(0, 8, 4):
                print(f"packedScheduleBlocked<{output},{plane}>(src,packed);")
        print("packedScheduleWriteback(packed,out);")
    print("}")
    print(f"unsigned packedCoeffTileMode() {{return {mode};}}")
    print("}")
    print(f"Outer schedule {mode}: {MODES[mode]}; exact GFNI coverage, parity "
          "semantics, packed layout, and wide output transpose checked", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("header", type=Path)
    parser.add_argument("--mode", type=int, choices=MODES, default=10)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    verify_schedule(args.header)
    if args.verify_only:
        retained_source(args.header)
        print("All generator checks passed for modes 10--12", file=sys.stderr)
    else:
        emit(args.header, args.mode)


if __name__ == "__main__":
    main()
