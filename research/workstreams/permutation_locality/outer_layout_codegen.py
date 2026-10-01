"""Emit exact-map GL32 preparation variants with BCH tile mode 1 or 5.

Usage: python outer_layout_codegen.py BchCircuit.h --mode 0..3 [--tile-mode 1|5] > Candidate.cpp

All variants export bchPackedCoeffCompact with the existing compact coefficient
layout. They do not change setup, the BCH map, or output layout. Tile mode 1
retains the original bit packing; tile mode 5 retains its byte/GFNI packing.
Only one variant is emitted per translation unit, so the existing encoder can
link it in place of PackedCoeff1.o or PackedCoeff5.o without runtime selection.
The default tile mode remains 1 and produces the same source as before.
"""

import argparse
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import random
import subprocess
import sys

import packed_bch_tune_codegen as tune
import packed_coeff_codegen as coeff_codegen
from packed_coeff_codegen import definition


VARIANTS = {
    0: "retained compact preparation and selected BCH schedule",
    1: "four GL32 terms with paired XOR/ternary reduction",
    2: "factored even/odd GL32 rotations, one plane at a time",
    3: "factored even/odd GL32 rotations, four-plane scheduling",
}


def rotate(values, count):
    return [values[(lane + count) % 4] for lane in range(4)]


def affine(byte, matrix):
    """Intel GFNI byte map; the most significant matrix byte is output row 0."""
    return sum((((matrix >> (8 * (7 - bit))) & byte).bit_count() & 1) << bit
               for bit in range(8))


def apply_lanes(values, matrices):
    return [affine(value, matrix) for value, matrix in zip(values, matrices)]


def xor_lanes(*terms):
    result = [0] * 4
    for term in terms:
        result = [a ^ b for a, b in zip(result, term)]
    return result


def original(values, matrices):
    return xor_lanes(*(apply_lanes(rotate(values, d), matrices[d])
                       for d in range(4)))


def factored(values, matrices):
    half = rotate(values, 2)
    odd = xor_lanes(apply_lanes(values, rotate(matrices[1], -1)),
                    apply_lanes(half, rotate(matrices[3], -1)))
    return xor_lanes(apply_lanes(values, matrices[0]),
                     apply_lanes(half, matrices[2]), rotate(odd, 1))


def verify():
    """Check exact linear-map identities; this is not a performance test."""
    # Track the input-lane and coefficient-lane indices symbolically. This
    # verifies the factoring for arbitrary byte maps, not only sampled GL32s.
    for lane in range(4):
        old = {(d, lane, (lane + d) % 4) for d in range(4)}
        new = {(0, lane, lane), (2, lane, (lane + 2) % 4),
               (1, ((lane + 1) - 1) % 4, (lane + 1) % 4),
               (3, ((lane + 1) - 1) % 4, (lane + 3) % 4)}
        assert old == new
    # Compact coefficient expansion duplicates each 64-bit map into both
    # qwords of a 128-bit payload lane. The odd expansion absorbs R^-1.
    assert [0, 0, 1, 1, 2, 2, 3, 3] == [lane for lane in range(4) for _ in range(2)]
    assert [3, 3, 0, 0, 1, 1, 2, 2] == [lane for lane in rotate(list(range(4)), -1)
                                        for _ in range(2)]
    # Exhaust every single-entry coefficient matrix, using the input basis
    # on which that entry acts. Matrices need not be invertible for the
    # identity, so this also covers arbitrary GL32 submatrices.
    for d in range(4):
        for lane in range(4):
            for output_bit in range(8):
                for input_bit in range(8):
                    matrices = [[0] * 4 for _ in range(4)]
                    matrices[d][lane] = 1 << (8 * (7 - output_bit) + input_bit)
                    values = [0] * 4
                    values[(lane + d) % 4] = 1 << input_bit
                    expected = [0] * 4
                    expected[lane] = 1 << output_bit
                    assert original(values, matrices) == expected
                    assert factored(values, matrices) == expected
    # Full dense examples and every input basis check interactions and XOR
    # accumulation. The fixed seed is solely for repeatable generator tests.
    rng = random.Random(0x474C3332)
    for _ in range(32):
        matrices = [[rng.getrandbits(64) for _ in range(4)] for _ in range(4)]
        for bit in range(32):
            values = [0] * 4
            values[bit // 8] = 1 << (bit % 8)
            assert original(values, matrices) == factored(values, matrices)


def emit_factored_planes(planes):
    # Keep the physical input, its half-turn, and the odd-term sum live for
    # this small compile-time batch. There are no runtime lane/plane loops.
    for j in planes:
        print(f'const auto u{j}=_mm512_shuffle_i32x4(v{j},v{j},0x4e);')
    for j in planes:
        print(f'const auto odd{j}=_mm512_xor_si512('
              f'_mm512_gf2p8affine_epi64_epi8(v{j},m1,0),'
              f'_mm512_gf2p8affine_epi64_epi8(u{j},m3,0));')
    for j in planes:
        print(f'const auto mixed{j}=_mm512_ternarylogic_epi64('
              f'_mm512_gf2p8affine_epi64_epi8(v{j},m0,0),'
              f'_mm512_gf2p8affine_epi64_epi8(u{j},m2,0),'
              f'_mm512_shuffle_i32x4(odd{j},odd{j},0x39),0x96);')
        print(f'_mm512_store_si512(src+8*group+{j},mixed{j});')


def emit_preparation(mode, packed_layout=False):
    print('alignas(64) __m512i src[256];')
    print('const auto duplicate=_mm512_setr_epi64(0,0,1,1,2,2,3,3);')
    if mode >= 2:
        print('const auto duplicateOdd=_mm512_setr_epi64(3,3,0,0,1,1,2,2);')
    print('for(unsigned group=0;group<32;++group) {')
    for j in range(8):
        print(f'auto v{j}=_mm512_loadu_si512(a+4*(8*group+{j}));')
    pack = 'packedTunePack5' if packed_layout else 'orthoBlend'
    print(pack + '(' + ','.join(f'v{j}' for j in range(8)) + ');')
    for d in range(4):
        print(f'const auto c{d}=_mm256_load_si256('
              f'reinterpret_cast<const __m256i*>(coeff+16*group+4*{d}));')
        index = 'duplicateOdd' if mode >= 2 and d % 2 else 'duplicate'
        # Only initialized low qwords of the cast are selected.
        print(f'const auto m{d}=_mm512_permutexvar_epi64('
              f'{index},_mm512_castsi256_si512(c{d}));')
    if mode >= 2:
        batch = 4 if mode == 3 else 1
        for start in range(0, 8, batch):
            print('{')
            emit_factored_planes(range(start, start + batch))
            print('}')
    else:
        for j in range(8):
            print('{')
            for d, immediate in ((0, None), (1, '0x39'), (2, '0x4e'), (3, '0x93')):
                source = f'v{j}' if d == 0 else f'_mm512_shuffle_i32x4(v{j},v{j},{immediate})'
                print(f'const auto term{d}=_mm512_gf2p8affine_epi64_epi8({source},m{d},0);')
            if mode == 1:
                value = '_mm512_ternarylogic_epi64(term0,term1,_mm512_xor_si512(term2,term3),0x96)'
            else:
                value = '_mm512_xor_si512(_mm512_xor_si512(_mm512_xor_si512(term0,term1),term2),term3)'
            print(f'_mm512_store_si512(src+8*group+{j},{value});')
            print('}')
    print('}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('header', type=Path, nargs='?')
    parser.add_argument('--mode', type=int, choices=VARIANTS, default=2)
    parser.add_argument('--tile-mode', type=int, choices=(1, 5), default=1)
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    verify()
    if args.verify_only:
        if args.tile_mode == 5:
            # This proves that tile 5 applies one common payload permutation
            # in every row lane, with an exact inverse. The GL32 identities
            # above hold for arbitrary byte inputs after either packing.
            tune.verify()
            print('Tile-5 common payload permutation and inverse verified on '
                  'all 1024 packing basis bits.', file=sys.stderr)
        print('GL32 rotation/coefficient indices, 1024 matrix-entry bases, and '
              '1024 dense-map input bases verified.', file=sys.stderr)
        return
    if args.header is None:
        parser.error('header is required unless --verify-only is supplied')
    # Locate the retained generator through its import, allowing this probe to
    # live elsewhere with PYTHONPATH pointing to the retained source directory.
    # The retained generator also runs tune.verify(), including the tile-5
    # common payload permutation and inverse check. GL32 acts independently
    # on each payload bit, so this permutation commutes with its lane maps.
    # Preserve Original, Aligned, and the selected BCH helpers byte-for-byte;
    # replace only the Compact export required by the measured encoder.
    retained = subprocess.run([sys.executable,
                              str(Path(coeff_codegen.__file__).resolve()),
                              str(args.header), '--tile-mode', str(args.tile_mode)],
                              check=True, capture_output=True, text=True)
    if retained.stderr:
        print(retained.stderr, file=sys.stderr, end='')
    old_compact = definition(retained.stdout, 'SPIN_NOINLINE void bchPackedCoeffCompact(')
    replacement = StringIO()
    with redirect_stdout(replacement):
        print('SPIN_NOINLINE void bchPackedCoeffCompact(const block* __restrict a,'
              'block* __restrict out,const std::uint64_t* __restrict coeff) {')
        print('using namespace coeff_probe_detail;')
        emit_preparation(args.mode, packed_layout=args.tile_mode == 5)
        for output in range(0, 16, tune.VARIANTS[args.tile_mode][0]):
            print(f'packedCoeffTile<{output}>(src,out);')
        print('}')
    print(retained.stdout.replace(old_compact, replacement.getvalue().rstrip()))
    print('namespace spin::detail::kernel {')
    print(f'unsigned outerLayoutMode() {{return {args.mode};}}\n}}')
    print(f'Outer layout mode {args.mode}: {VARIANTS[args.mode]}; exact map, '
          f'compact coefficients, retained BCH mode {args.tile_mode}.', file=sys.stderr)


if __name__ == '__main__':
    main()
