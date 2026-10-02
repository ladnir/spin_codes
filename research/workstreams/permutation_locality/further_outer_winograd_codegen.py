"""Exact one-level Winograd factorization of the mixed-layout dense BCH map.

The retained GL32 preparation, systematic/parity forms, and fused output byte
network are unchanged.  Only the fixed dense BCH multiplication is factored.
One tile uses 1792 rather than 2048 GFNI instructions, at the cost of extra
XORs, coefficient broadcasts, and L1 traffic.  This is a screening candidate,
not a selected kernel.  --planes 4 uses 28 accumulators; --planes 2 uses 14.
"""
import argparse
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import sys

import outer_transpose_codegen as retained
import outer_unpack_codegen as unpack
from packed_coeff_codegen import definition


def xor(*values):
    result = list(values[0])
    for value in values[1:]:
        result = [a ^ b for a, b in zip(result, value)]
    return result


def product(matrix, value):
    """Symbolic GFNI with constant first operand and variable second operand."""
    result = []
    for byte in range(8):
        row = (matrix >> (8 * byte)) & 255
        for bit in range(8):
            form = 0
            for source in range(8):
                if row >> source & 1:
                    form ^= value[8 * (7 - bit) + source]
            result.append(form)
    return result


def matrices(rows):
    return [[sum(((rows[8 * out + j] >> (128 + 8 * inp)) & 255)
                 << (8 * j) for j in range(8)) for inp in range(16)]
            for out in range(16)]


def factor(a, out, inp):
    a11, a12 = a[out][inp], a[out][inp + 8]
    a21, a22 = a[out + 8][inp], a[out + 8][inp + 8]
    return (a11, a12, a11 ^ a12 ^ a21 ^ a22, a22,
            a21 ^ a22, a11 ^ a21 ^ a22, a11 ^ a21)


def verify(rows):
    """Propagate all 8192 input-bit forms through both complete dense maps."""
    a = matrices(rows)
    source = [[[1 << (64 * (8 * group + plane) + bit) for bit in range(64)]
               for plane in range(8)] for group in range(16)]
    zero = [0] * 64
    for out in range(8):
        for plane in range(4):
            p = [zero for _ in range(7)]
            for inp in range(8):
                b11, b12 = source[inp][plane], source[inp][plane + 4]
                b21, b22 = source[inp + 8][plane], source[inp + 8][plane + 4]
                t1 = xor(b12, b11)
                t2 = xor(b22, t1)
                t3 = xor(b22, b12)
                t4 = xor(t2, b21)
                operands = (b11, b21, b22, t4, t1, t2, t3)
                for i, (coefficient, operand) in enumerate(zip(factor(a, out, inp), operands)):
                    p[i] = xor(p[i], product(coefficient, operand))
            u2 = xor(p[0], p[5])
            u3 = xor(u2, p[6])
            actual = (xor(p[0], p[1]), xor(u2, p[4], p[2]),
                      xor(u3, p[3]), xor(u3, p[4]))
            for value, group, component in zip(actual, (out, out, out + 8, out + 8),
                                                (plane, plane + 4, plane, plane + 4)):
                expected = zero
                for inp in range(16):
                    expected = xor(expected, product(a[group][inp], source[inp][component]))
                if value != expected:
                    raise ArithmeticError('complete Winograd dense BCH basis mismatch')


def emit_tables(rows):
    a = matrices(rows)
    print('alignas(64) static constexpr std::uint64_t furtherWinogradCoeff[8][7][8]={')
    for output in range(8):
        print('{')
        for term in range(7):
            print('{' + ','.join(f'0x{factor(a, output, inp)[term]:016x}ULL'
                                for inp in range(8)) + '},')
        print('},')
    print('};')


def operand(term, inp, plane):
    if term == 0:
        return f'src+128+8*({inp})+{plane}'
    if term == 1:
        return f'src+192+8*({inp})+{plane}'
    if term == 2:
        return f'src+196+8*({inp})+{plane}'
    temporary = {3: 3, 4: 0, 5: 1, 6: 2}[term]
    return f'temporary+{32 * temporary}+4*({inp})+{plane}'


def emit_tile(planes):
    print('template<unsigned output> static SPIN_NOINLINE void furtherWinogradTile('
          'const __m512i* __restrict src,const __m512i* __restrict temporary,'
          'block* __restrict out){')
    print('static_assert(output<8);alignas(64) __m512i result[16];')
    for first in range(0, 4, planes):
        print('{')
        for term in range(7):
            for plane in range(first, first + planes):
                print(f'auto p{term}_{plane}=_mm512_setzero_si512();')
        print('for(unsigned input=0;input<8;input+=2){')
        for term in range(7):
            print('{')
            print(f'const auto m=_mm512_set1_epi64(furtherWinogradCoeff[output][{term}][input]);')
            print(f'const auto n=_mm512_set1_epi64(furtherWinogradCoeff[output][{term}][input+1]);')
            for plane in range(first, first + planes):
                left, right = operand(term, 'input', plane), operand(term, 'input+1', plane)
                print(f'p{term}_{plane}=_mm512_ternarylogic_epi64(p{term}_{plane},'
                      f'_mm512_gf2p8affine_epi64_epi8(m,_mm512_load_si512({left}),0),'
                      f'_mm512_gf2p8affine_epi64_epi8(n,_mm512_load_si512({right}),0),0x96);')
            print('}')
        print('}')
        for plane in range(first, first + planes):
            print(f'const auto u2_{plane}=_mm512_xor_si512(p0_{plane},p5_{plane});')
            print(f'const auto u3_{plane}=_mm512_xor_si512(u2_{plane},p6_{plane});')
            forms = (f'_mm512_xor_si512(p0_{plane},p1_{plane})',
                     f'_mm512_ternarylogic_epi64(u2_{plane},p4_{plane},p2_{plane},0x96)',
                     f'_mm512_xor_si512(u3_{plane},p3_{plane})',
                     f'_mm512_xor_si512(u3_{plane},p4_{plane})')
            for offset, form in zip((0, 4, 8, 12), forms):
                print(f'_mm512_store_si512(result+{offset + plane},{form});')
        print('}')
    retained.constant_index('parityByte', [8 * half + 7 for lane in range(4)
                                         for half in range(2) for byte in range(8)])
    for p in range(2):
        offset = 8 * p
        print(f'constexpr auto parity{p}=[]{{std::uint64_t mask=0;for(unsigned j=0;j<8;++j)'
              f'if(BchRows[8*(output+{offset})+j][1]>>63)mask|=std::uint64_t(255)<<(8*j);return mask;}}();')
        for plane in range(8):
            y = 8 * p + plane
            print(f'auto y{y}=_mm512_load_si512(src+8*(output+{offset})+{plane});')
            print(f'if constexpr(output+{offset}==15)y{y}=_mm512_and_si512(y{y},'
                  '_mm512_set1_epi64(0x00ffffffffffffffULL));')
            print(f'y{y}=_mm512_ternarylogic_epi64(y{y},'
                  f'_mm512_shuffle_epi8(_mm512_load_si512(src+120+{plane}),parityByte),'
                  f'_mm512_set1_epi64(parity{p}),0x78);')
            print(f'y{y}=_mm512_xor_si512(y{y},_mm512_load_si512(result+{y}));')
    output = StringIO()
    with redirect_stdout(output):
        unpack.emit_fused_output(2, unpack.network())
    text = output.getvalue().replace('8*(output+1)', '8*(output+8)')
    for p in range(2):
        for plane in range(8):
            old = f'const auto x{plane}=_mm512_gf2p8affine_epi64_epi8(basis,y{8*p+plane},0);'
            if text.count(old) != 1:
                raise ArithmeticError('fused output emitter changed')
            text = text.replace(old, f'const auto x{plane}=y{8*p+plane};')
    print(text)
    print('}')


def emit_temporaries():
    print('alignas(64) __m512i temporary[128];')
    print('for(unsigned input=0;input<8;++input){')
    for plane in range(4):
        print('{')
        for label, base in (('b11', 128), ('b12', 132), ('b21', 192), ('b22', 196)):
            print(f'const auto {label}=_mm512_load_si512(src+{base}+8*input+{plane});')
        print('const auto t1=_mm512_xor_si512(b12,b11);')
        print('const auto t2=_mm512_xor_si512(b22,t1);')
        print('const auto t3=_mm512_xor_si512(b22,b12);')
        print('const auto t4=_mm512_xor_si512(t2,b21);')
        for temporary in range(4):
            print(f'_mm512_store_si512(temporary+{32*temporary}+4*input+{plane},t{temporary+1});')
        print('}')
    print('}')


def generate(header, planes):
    rows = retained.read_rows(header)
    verify(rows)
    source = retained.generate(header)
    old_tile = definition(source, 'template<unsigned output> static SPIN_NOINLINE void packedCoeffTile(')
    replacement = StringIO()
    with redirect_stdout(replacement):
        emit_tables(rows)
        emit_tile(planes)
    source = source.replace(old_tile, replacement.getvalue())
    for mode in ('Original', 'Aligned', 'Compact'):
        old = definition(source, f'SPIN_NOINLINE void bchPackedCoeff{mode}(')
        new = StringIO()
        with redirect_stdout(new):
            print(f'SPIN_NOINLINE void bchPackedCoeff{mode}(const block* __restrict a,'
                  'block* __restrict out,const std::uint64_t* __restrict coeff){')
            print('using namespace coeff_probe_detail;')
            retained.emit_preparation(mode)
            emit_temporaries()
            for output in range(8):
                print(f'furtherWinogradTile<{output}>(src,temporary,out);')
            print('}')
        source = source.replace(old, new.getvalue())
    return source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('header', type=Path)
    parser.add_argument('--planes', type=int, choices=(2, 4), default=4)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    if args.verify_only:
        rows = retained.read_rows(args.header)
        retained.verify(rows)
        verify(rows)
    else:
        source = generate(args.header, args.planes)
        if args.output:
            args.output.write_text(source, encoding='utf-8')
        else:
            print(source)
    print('Verified all8192 dense BCH input-bit forms under exact Winograd; '
          'retained mixed-layout, systematic/parity, and output-network proofs pass.', file=sys.stderr)


if __name__ == '__main__':
    main()
