"""Emit exact-layout packed GL32/BCH schedule probes, not a new encoder.

Usage: python packed_bch_tune_codegen.py path/to/generated/BchCircuit.h
The generated translation unit replaces PackedMixer.o. Compile it with
-DSPIN_PACKED_TUNE=N; the existing bchPackedFullGl32 entry point statically
selects N=0 (retained), 1 (wide stores), 2 (wide stores/unrolled pairs),
3 (wide stores/one output group), 4 (wide stores/three output groups), or
5 (byte/GFNI packing, paired two-output tiles, wide stores).
All other packed-mixer entry points remain available from the old generator.
"""
import argparse
from pathlib import Path
import subprocess
import sys


VARIANTS = {
    1: (2, True, False, 'two output groups, paired inputs, wide output stores'),
    2: (2, True, True, 'two output groups, fully unrolled pairs, wide output stores'),
    3: (1, True, False, 'one output group, paired inputs, wide output stores'),
    4: (3, False, False, 'three output groups, single inputs, wide output stores'),
    5: (2, True, False, 'byte/GFNI packing, two output groups, paired inputs, wide stores'),
}


def shuffle128(left, right, immediate):
    """Symbolic VSHUFI64X2 selection, in 128-bit lanes."""
    return [left[immediate & 3], left[(immediate >> 2) & 3],
            right[(immediate >> 4) & 3], right[(immediate >> 6) & 3]]


def unpack128(left, right, bits, high=False):
    """Exact lane-local SIMD unpack on one 128-bit lane."""
    count = 128 // bits
    start = count//2 if high else 0
    mask = (1 << bits)-1
    return sum((((value >> (bits*(start+j))) & mask) << (bits*(2*j+k)))
               for j in range(count//2) for k, value in enumerate((left, right)))


def byte_transpose128(values):
    t = [unpack128(values[i], values[i+1], 8, high)
         for i in (0, 2, 4, 6) for high in (False, True)]
    u = [unpack128(t[i], t[j], 16, high)
         for i, j in ((0, 2), (1, 3), (4, 6), (5, 7)) for high in (False, True)]
    return [unpack128(u[i], u[i+4], 32, high)
            for i in range(4) for high in (False, True)]


def bit_transpose128(value, descending):
    """GFNI(identity, value, 0): the payload supplies the matrix operand."""
    result = 0
    for half in range(2):
        matrix = (value >> (64*half)) & ((1 << 64)-1)
        for byte in range(8):
            basis = 1 << (7-byte if descending else byte)
            output = sum((((matrix >> (8*(7-row))) & basis).bit_count() & 1) << row
                         for row in range(8))
            result |= output << (64*half+8*byte)
    return result


def pack128(values):
    return [bit_transpose128(value, True) for value in byte_transpose128(values[::-1])]


def unpack_packed128(values):
    t = byte_transpose128([bit_transpose128(value, False) for value in values])
    return [unpack128(t[i], t[i+4], 8, high)
            for i in range(4) for high in (False, True)]


def verify():
    # Each input vector holds the same column from all four row lanes.
    # Verify that every wide store instead holds four columns of one row.
    columns = [[(row, col) for row in range(4)] for col in range(4)]
    a = shuffle128(columns[0], columns[1], 0x44)
    b = shuffle128(columns[0], columns[1], 0xee)
    c = shuffle128(columns[2], columns[3], 0x44)
    d = shuffle128(columns[2], columns[3], 0xee)
    stores = [shuffle128(a, c, 0x88), shuffle128(a, c, 0xdd),
              shuffle128(b, d, 0x88), shuffle128(b, d, 0xdd)]
    assert stores == [[(row, col) for col in range(4)] for row in range(4)]
    # All 1024 basis bits of one 128-bit row lane. Old ortho packing puts
    # payload bit b at register b%8, byte b//8, coordinate bit c. The new
    # layout moves only that payload position; coordinate bits are unchanged.
    # The same permutation in every row lane commutes with all shared BCH
    # matrices and with the four-row GL32 maps, which act per payload bit.
    from gfni_bch import ortho
    mask64 = (1 << 64)-1
    for col in range(8):
        for bit in range(128):
            source = [0]*8
            source[col] = 1 << bit
            packed = pack128(source)
            old = [0]*8
            for half in range(2):
                lane = ortho([(x >> (64*half)) & mask64 for x in source])
                for j in range(8):
                    old[j] |= lane[j] << (64*half)
            expected = []
            for j in range(8):
                value = 0
                for byte in range(16):
                    payload = 16*j+8*(byte//8)+7-byte%8
                    value |= ((old[payload%8] >> (8*(payload//8))) & 255) << (8*byte)
                expected.append(value)
            assert packed == expected
            assert unpack_packed128(packed) == source
    # Check every dense 8x8 submatrix action is scheduled exactly once for
    # every output group and each of the eight packed element bitplanes.
    for width, paired, _, _ in VARIANTS.values():
        coverage = [[[0] * 16 for _ in range(8)] for _ in range(16)]
        for output in range(0, 16, width):
            for p in range(min(width, 16-output)):
                for step in range(0, 16, 2 if paired else 1):
                    for inp in range(step, step+(2 if paired else 1)):
                        for plane in range(8):
                            coverage[output+p][plane][inp] += 1
        assert all(n == 1 for group in coverage for plane in group for n in plane)


def emit_byte_transpose_helpers():
    args = ','.join(f'__m512i& x{j}' for j in range(8))
    print(f'static SPIN_FORCEINLINE void packedTuneByteTranspose5({args}) {{')
    for pair in range(4):
        for high, suffix in ((0, 'lo'), (1, 'hi')):
            print(f'const auto t{2*pair+high}=_mm512_unpack{suffix}_epi8(x{2*pair},x{2*pair+1});')
    for pair, (left, right) in enumerate(((0, 2), (1, 3), (4, 6), (5, 7))):
        for high, suffix in ((0, 'lo'), (1, 'hi')):
            print(f'const auto u{2*pair+high}=_mm512_unpack{suffix}_epi16(t{left},t{right});')
    for pair in range(4):
        print(f'x{2*pair}=_mm512_unpacklo_epi32(u{pair},u{pair+4});')
        print(f'x{2*pair+1}=_mm512_unpackhi_epi32(u{pair},u{pair+4});')
    print('}')
    print(f'static SPIN_FORCEINLINE void packedTunePack5({args}) {{')
    print('packedTuneByteTranspose5('+','.join(f'x{j}' for j in range(7, -1, -1))+');')
    # Passing references in reverse order puts each computed output in the
    # corresponding reversed source variable. Rename in reverse on readback.
    print('const auto basis=_mm512_set1_epi64(0x0102040810204080ULL);')
    for j in range(8):
        print(f'const auto y{j}=_mm512_gf2p8affine_epi64_epi8(basis,x{7-j},0);')
    for j in range(8):
        print(f'x{j}=y{j};')
    print('}')
    print(f'static SPIN_FORCEINLINE void packedTuneUnpack5({args}) {{')
    print('const auto basis=_mm512_set1_epi64(0x8040201008040201ULL);')
    for j in range(8):
        print(f'x{j}=_mm512_gf2p8affine_epi64_epi8(basis,x{j},0);')
    print('packedTuneByteTranspose5('+','.join(f'x{j}' for j in range(8))+');')
    for j in range(4):
        print(f'const auto y{2*j}=_mm512_unpacklo_epi8(x{j},x{j+4});')
        print(f'const auto y{2*j+1}=_mm512_unpackhi_epi8(x{j},x{j+4});')
    for j in range(8):
        print(f'x{j}=y{j};')
    print('}')


def emit_output(parallel, packed_layout=False):
    for p in range(parallel):
        unpack = 'packedTuneUnpack5' if packed_layout else 'orthoBlend'
        print(unpack+'(' + ','.join(f'y{8*p+j}' for j in range(8)) + ');')
        for j in (0, 4):
            at = 8*p+j
            print('{')
            print(f'const auto a=_mm512_shuffle_i64x2(y{at},y{at+1},0x44);')
            print(f'const auto b=_mm512_shuffle_i64x2(y{at},y{at+1},0xee);')
            print(f'const auto c=_mm512_shuffle_i64x2(y{at+2},y{at+3},0x44);')
            print(f'const auto d=_mm512_shuffle_i64x2(y{at+2},y{at+3},0xee);')
            for row, (left, right, immediate) in enumerate((
                    ('a', 'c', '0x88'), ('a', 'c', '0xdd'),
                    ('b', 'd', '0x88'), ('b', 'd', '0xdd'))):
                print(f'_mm512_storeu_si512(out+{128*row}+8*(output+{p})+{j},'
                      f'_mm512_shuffle_i64x2({left},{right},{immediate}));')
            print('}')


def emit_input_step(parallel, paired, index):
    print('{')
    for p in range(parallel):
        print(f'const auto m{p}=_mm512_set1_epi64(matrices[output+{p}][{index}]);')
        if paired:
            print(f'const auto n{p}=_mm512_set1_epi64(matrices[output+{p}][{index}+1]);')
    for j in range(8):
        print(f'const auto x{j}=_mm512_load_si512(src+128+8*({index})+{j});')
        if paired:
            print(f'const auto z{j}=_mm512_load_si512(src+128+8*({index})+8+{j});')
        for p in range(parallel):
            first = f'_mm512_gf2p8affine_epi64_epi8(x{j},m{p},0)'
            if paired:
                second = f'_mm512_gf2p8affine_epi64_epi8(z{j},n{p},0)'
                print(f'y{8*p+j}=_mm512_ternarylogic_epi64(y{8*p+j},{first},{second},0x96);')
            else:
                print(f'y{8*p+j}=_mm512_xor_si512(y{8*p+j},{first});')
    print('}')


def emit_tile(name, parallel, paired, unroll, packed_layout=False):
    print(f'template<unsigned output> static SPIN_NOINLINE void {name}('
          'const __m512i* __restrict src,block* __restrict out) {')
    print(f'static_assert(output+{parallel}<=16);')
    # This sparse half is exactly the retained implementation. Coordinate
    # 127 is parity, not a systematic coordinate, and is masked accordingly.
    for p in range(parallel):
        for j in range(8):
            print(f'auto y{8*p+j}=_mm512_load_si512(src+8*(output+{p})+{j});')
            print(f'if constexpr(output+{p}==15)y{8*p+j}='
                  f'_mm512_and_si512(y{8*p+j},_mm512_set1_epi8(0x7f));')
        print(f'constexpr auto parity{p}=[] {{std::uint64_t m=0;for(unsigned j=0;j<8;++j)'
              f'if(BchRows[8*(output+{p})+j][1]>>63)m|=std::uint64_t(0x80)<<(8*(7-j));return m;}}();')
        for j in range(8):
            print(f'y{8*p+j}=_mm512_xor_si512(y{8*p+j},_mm512_gf2p8affine_epi64_epi8('
                  f'_mm512_load_si512(src+120+{j}),_mm512_set1_epi64(parity{p}),0));')
    if unroll:
        for index in range(0, 16, 2 if paired else 1):
            emit_input_step(parallel, paired, str(index))
    else:
        print(f'for(unsigned input=0;input<16;input+={2 if paired else 1}) {{')
        emit_input_step(parallel, paired, 'input')
        print('}')
    emit_output(parallel, packed_layout)
    print('}')


def emit_preparation(packed_layout=False):
    # Keep the sampled 32x32 matrices, their cyclic-lane coefficient layout,
    # and all pack/mix operations identical to bchPackedFullGl32Retained.
    print('alignas(64) __m512i src[256];')
    print('for(unsigned group=0;group<32;++group) {')
    for j in range(8):
        print(f'auto v{j}=_mm512_loadu_si512(a+4*(8*group+{j}));')
    pack = 'packedTunePack5' if packed_layout else 'orthoBlend'
    print(pack+'(' + ','.join(f'v{j}' for j in range(8)) + ');')
    for d in range(4):
        print(f'const auto m{d}=_mm512_loadu_si512(coeff+32*group+8*{d});')
    for j in range(8):
        print(f'const auto z{j}=v{j};v{j}=_mm512_gf2p8affine_epi64_epi8(z{j},m0,0);')
        for d, immediate in ((1, '0x39'), (2, '0x4e'), (3, '0x93')):
            print(f'v{j}=_mm512_xor_si512(v{j},_mm512_gf2p8affine_epi64_epi8('
                  f'_mm512_shuffle_i32x4(z{j},z{j},{immediate}),m{d},0));')
        print(f'_mm512_store_si512(src+8*group+{j},v{j});')
    print('}')


def emit_variants():
    print('#ifndef SPIN_PACKED_TUNE\n#define SPIN_PACKED_TUNE 0\n#endif')
    print('#if SPIN_PACKED_TUNE < 0 || SPIN_PACKED_TUNE > 5\n#error Invalid SPIN_PACKED_TUNE\n#endif')
    print('namespace spin::detail::kernel {')
    emit_byte_transpose_helpers()
    for mode, (parallel, paired, unroll, description) in VARIANTS.items():
        name = f'packedTuneTile{mode}'
        print(f'// Mode {mode}: {description}.')
        emit_tile(name, parallel, paired, unroll, packed_layout=mode == 5)
        if 16 % parallel:
            emit_tile(name+'Tail', 16 % parallel, paired, unroll)
        print(f'SPIN_NOINLINE void bchPackedTune{mode}(const block* __restrict a,'
              'block* __restrict out,const std::uint64_t* __restrict coeff) {')
        emit_preparation(packed_layout=mode == 5)
        for output in range(0, 16, parallel):
            tile_name = name if output+parallel <= 16 else name+'Tail'
            print(f'{tile_name}<{output}>(src,out);')
        print('}')
    print('SPIN_NOINLINE void bchPackedFullGl32(const block* __restrict a,'
          'block* __restrict out,const std::uint64_t* __restrict coeff) {')
    print('#if SPIN_PACKED_TUNE == 0\nbchPackedFullGl32Retained(a,out,coeff);')
    for mode in VARIANTS:
        print(f'#elif SPIN_PACKED_TUNE == {mode}\nbchPackedTune{mode}(a,out,coeff);')
    print('#endif\n}\n}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('header', type=Path)
    args = parser.parse_args()
    verify()
    retained = subprocess.run([sys.executable, str(Path(__file__).with_name('packed_mixer_codegen.py')),
        str(args.header)], check=True, capture_output=True, text=True)
    before = 'void bchPackedFullGl32('
    if retained.stdout.count(before) != 1:
        raise ValueError('expected exactly one retained full-GL32 entry point')
    print(retained.stdout.replace(before, 'void bchPackedFullGl32Retained('))
    if retained.stderr:
        print(retained.stderr, file=sys.stderr, end='')
    emit_variants()
    print('Packed BCH tune: exact output-lane transpose and dense schedules checked; '
          '1024 packed-layout basis bits and inverse checked; '
          'SPIN_PACKED_TUNE=0..5, setup layout unchanged', file=sys.stderr)


if __name__ == '__main__':
    main()
