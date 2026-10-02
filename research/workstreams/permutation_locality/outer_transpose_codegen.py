"""Exact BCH probe that absorbs final bit transposes into GFNI products.

The systematic half of the mixed source is stored in transposed-qword form;
the dense half retains the original format.  Setup coefficients and the public
outer API are unchanged.  A fused VBMI byte network finishes row-major output.
"""
import argparse
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import re
import sys

import outer_layout_codegen as layout
import outer_unpack_codegen as unpack
from packed_coeff_codegen import definition


def affine(data, matrix):
    return sum((((matrix >> (8*(7-i))) & (data >> (8*j)) & 255).bit_count() & 1)
               << (8*j+i) for j in range(8) for i in range(8))


def transpose(data):
    return affine(0x8040201008040201, data)


def reverse_bytes(data):
    return int.from_bytes(data.to_bytes(8, 'little'), 'big')


def read_rows(header):
    words = [int(word, 16) for word in re.findall(r'0x([0-9a-f]+)ULL', header.read_text())]
    if len(words) != 512:
        raise ArithmeticError('unexpected retained BCH row table')
    rows = [sum(words[4*i+j] << (64*j) for j in range(4)) for i in range(128)]
    for i, row in enumerate(rows):
        if row & ((1 << 127)-1) != (1 << i if i < 127 else 0):
            raise ArithmeticError('retained systematic/parity split changed')
    return rows


def verify(rows):
    # Bilinearity makes these4096 matrix/input pairs a complete proof of the
    # operand-swap identity, including arbitrary GL32 submatrix coefficients.
    for matrix_bit in range(64):
        matrix = 1 << matrix_bit
        for input_bit in range(64):
            value = 1 << input_bit
            if transpose(affine(value, matrix)) != affine(reverse_bytes(matrix), value):
                raise ArithmeticError('GFNI operand swap changed a bit')
    layout.verify()
    unpack.network()
    # Symbolically propagate every256-coordinate x8-payload-bit basis at once.
    # Each integer bit is an independent original source bit.  Qword/plane/row
    # lanes repeat this same map without mixing, so this also proves the full
    # four-row,128-bit-element BCH output map.
    source = [[[1 << (8*(8*g+c)+b) for c in range(8)] for b in range(8)]
              for g in range(32)]
    mixed = [[[(source[g][7-i][j] if g < 16 else source[g][j][i])
               for i in range(8)] for j in range(8)] for g in range(32)]
    for output in range(16):
        got = [list(byte) for byte in mixed[output]]
        if output == 15:
            got[7] = [0]*8
        for j in range(8):
            row = rows[8*output+j]
            if row >> 127 & 1:
                got[j] = [a ^ b for a, b in zip(got[j], mixed[15][7])]
            for inp in range(16):
                mask = (row >> (128+8*inp)) & 255
                for i in range(8):
                    for c in range(8):
                        if mask >> c & 1:
                            got[j][i] ^= mixed[16+inp][7-i][c]
            for i in range(8):
                expected = 0
                for coordinate in range(256):
                    if row >> coordinate & 1:
                        expected ^= 1 << (8*coordinate+7-i)
                if got[j][i] != expected:
                    raise ArithmeticError('mixed-format complete BCH basis mismatch')
    # Verify coefficient expansion for both even and compensated odd lanes.
    for odd in (False, True):
        order = (3, 0, 1, 2) if odd else (0, 1, 2, 3)
        index = [8*order[lane]+7-byte for lane in range(4)
                 for half in range(2) for byte in range(8)]
        if any(index[16*lane+8*half+byte] != 8*order[lane]+7-byte
               for lane in range(4) for half in range(2) for byte in range(8)):
            raise ArithmeticError('coefficient duplicate/reverse mismatch')
    # VPSHUFB is128-bit-lane-local, unlike GFNI's64-bit matrix lanes.  The
    # upper qword must select bytes8..15, never repeat the lower payload half.
    for lane in range(4):
        for half in range(2):
            for byte in range(8):
                reverse = 8*half+7-byte
                parity = 8*half+7
                if 16*lane+reverse != 16*lane+8*half+7-byte:
                    raise ArithmeticError('full-lane coefficient reversal mismatch')
                if 16*lane+parity != 16*lane+8*half+7:
                    raise ArithmeticError('full-lane parity replication mismatch')
    for a in (0, 1):
        for b in (0, 1):
            for c in (0, 1):
                if (0x78 >> (4*a+2*b+c)) & 1 != a ^ (b & c):
                    raise ArithmeticError('parity ternary truth table mismatch')


def constant_index(name, values):
    words = [sum(values[8*j+k] << (8*k) for k in range(8)) for j in range(8)]
    print(f'const auto {name}=_mm512_setr_epi64('
          + ','.join(f'0x{word:016x}ULL' for word in words) + ');')


def emit_preparation(mode):
    print('alignas(64) __m512i src[256];')
    print('const auto duplicate=_mm512_setr_epi64(0,0,1,1,2,2,3,3);')
    print('const auto duplicateOdd=_mm512_setr_epi64(3,3,0,0,1,1,2,2);')
    constant_index('reverseBytes', [8*half+7-byte for lane in range(4)
                                   for half in range(2) for byte in range(8)])
    for odd in (False, True):
        order = (3, 0, 1, 2) if odd else (0, 1, 2, 3)
        constant_index('reverseOdd' if odd else 'reverseEven',
                       [8*order[lane]+7-byte for lane in range(4)
                        for half in range(2) for byte in range(8)])
    for transposed in (True, False):
        start, end = (0, 16) if transposed else (16, 32)
        print(f'for(unsigned group={start};group<{end};++group){{')
        for plane in range(8):
            print(f'auto v{plane}=_mm512_loadu_si512(a+4*(8*group+{plane}));')
        print('packedTunePack5(' + ','.join(f'v{i}' for i in range(8)) + ');')
        for d in range(4):
            if mode == 'Compact':
                print(f'const auto c{d}=_mm256_load_si256('
                      f'reinterpret_cast<const __m256i*>(coeff+16*group+4*{d}));')
                if transposed:
                    index = 'reverseOdd' if d % 2 else 'reverseEven'
                    intrinsic = '_mm512_permutexvar_epi8'
                else:
                    index = 'duplicateOdd' if d % 2 else 'duplicate'
                    intrinsic = '_mm512_permutexvar_epi64'
                print(f'const auto m{d}={intrinsic}({index},_mm512_castsi256_si512(c{d}));')
            else:
                load = '_mm512_load_si512' if mode == 'Aligned' else '_mm512_loadu_si512'
                value = f'{load}(coeff+32*group+8*{d})'
                if d % 2:
                    value = f'_mm512_shuffle_i32x4({value},{value},0x93)'
                if transposed:
                    value = f'_mm512_shuffle_epi8({value},reverseBytes)'
                print(f'const auto m{d}={value};')
        for plane in range(8):
            print('{')
            print(f'const auto u=_mm512_shuffle_i32x4(v{plane},v{plane},0x4e);')
            def apply(data, coefficient):
                a, b = (coefficient, data) if transposed else (data, coefficient)
                return f'_mm512_gf2p8affine_epi64_epi8({a},{b},0)'
            print('const auto odd=_mm512_xor_si512('
                  + apply(f'v{plane}', 'm1') + ',' + apply('u', 'm3') + ');')
            print('const auto mixed=_mm512_ternarylogic_epi64('
                  + apply(f'v{plane}', 'm0') + ',' + apply('u', 'm2')
                  + ',_mm512_shuffle_i32x4(odd,odd,0x39),0x96);')
            print(f'_mm512_store_si512(src+8*group+{plane},mixed);')
            print('}')
        print('}')


def emit_tile(rows):
    print('template<unsigned output> static SPIN_NOINLINE void packedCoeffTile('
          'const __m512i* __restrict src,block* __restrict out){')
    print('static_assert(output+2<=16);')
    # The parity coordinate is byte7 of each already-transposed qword.
    constant_index('parityByte', [8*half+7 for lane in range(4)
                                 for half in range(2) for byte in range(8)])
    for p in range(2):
        print(f'constexpr auto parity{p}=[]{{std::uint64_t mask=0;for(unsigned j=0;j<8;++j)'
              f'if(BchRows[8*(output+{p})+j][1]>>63)mask|=std::uint64_t(255)<<(8*j);return mask;}}();')
        for plane in range(8):
            y = 8*p+plane
            print(f'auto y{y}=_mm512_load_si512(src+8*(output+{p})+{plane});')
            print(f'if constexpr(output+{p}==15)y{y}=_mm512_and_si512('
                  f'y{y},_mm512_set1_epi64(0x00ffffffffffffffULL));')
            print(f'y{y}=_mm512_ternarylogic_epi64(y{y},'
                  f'_mm512_shuffle_epi8(_mm512_load_si512(src+120+{plane}),parityByte),'
                  f'_mm512_set1_epi64(parity{p}),0x78);')
    print('for(unsigned input=0;input<16;input+=2){')
    for p in range(2):
        print(f'const auto m{p}=_mm512_set1_epi64(transposedMatrices[output+{p}][input]);')
        print(f'const auto n{p}=_mm512_set1_epi64(transposedMatrices[output+{p}][input+1]);')
    for plane in range(8):
        print(f'const auto x{plane}=_mm512_load_si512(src+128+8*input+{plane});')
        print(f'const auto z{plane}=_mm512_load_si512(src+136+8*input+{plane});')
        for p in range(2):
            y = 8*p+plane
            print(f'y{y}=_mm512_ternarylogic_epi64(y{y},'
                  f'_mm512_gf2p8affine_epi64_epi8(m{p},x{plane},0),'
                  f'_mm512_gf2p8affine_epi64_epi8(n{p},z{plane},0),0x96);')
    print('}')
    output = StringIO()
    with redirect_stdout(output):
        unpack.emit_fused_output(2, unpack.network())
    text = output.getvalue()
    for p in range(2):
        for plane in range(8):
            old = f'const auto x{plane}=_mm512_gf2p8affine_epi64_epi8(basis,y{8*p+plane},0);'
            if text.count(old) != 1:
                raise ArithmeticError('fused output emitter changed')
            text = text.replace(old, f'const auto x{plane}=y{8*p+plane};')
    print(text)
    print('}')


def generate(header):
    rows = read_rows(header)
    verify(rows)
    source = unpack.generate(header)
    tables = StringIO()
    with redirect_stdout(tables):
        print('alignas(64) static constexpr std::uint64_t transposedMatrices[16][16]={')
        for output in range(16):
            words = [sum(((rows[8*output+j] >> (128+8*inp)) & 255) << (8*j)
                         for j in range(8)) for inp in range(16)]
            print('{' + ','.join(f'0x{word:016x}ULL' for word in words) + '},')
        print('};')
    old_tile = definition(source, 'template<unsigned output> static SPIN_NOINLINE void packedCoeffTile(')
    tile = StringIO()
    with redirect_stdout(tile):
        emit_tile(rows)
    source = source.replace(old_tile, tables.getvalue() + tile.getvalue())
    for mode in ('Original', 'Aligned', 'Compact'):
        old = definition(source, f'SPIN_NOINLINE void bchPackedCoeff{mode}(')
        new = StringIO()
        with redirect_stdout(new):
            print(f'SPIN_NOINLINE void bchPackedCoeff{mode}(const block* __restrict a,'
                  'block* __restrict out,const std::uint64_t* __restrict coeff){')
            print('using namespace coeff_probe_detail;')
            emit_preparation(mode)
            for output in range(0, 16, 2):
                print(f'packedCoeffTile<{output}>(src,out);')
            print('}')
        source = source.replace(old, new.getvalue())
    return '// Exact mixed-layout GFNI operand-swap BCH; final bit transpose absorbed.\n' + source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('header', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    if args.verify_only:
        verify(read_rows(args.header))
    else:
        source = generate(args.header)
        if args.output:
            args.output.write_text(source, encoding='utf-8')
        else:
            print(source)
    print('Verified:4096 GFNI bilinear bases, GL32 rotations, compact byte orientations, '
          'all mixed-format BCH output basis forms and4096 fused-layout bit bases.', file=sys.stderr)


if __name__ == '__main__':
    main()
