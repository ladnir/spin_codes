"""Emit a bounded raw-message-basis BCH CLMUL experiment.

Usage: python outer_clmul_codegen.py BchCircuit.h > OuterClmul.cpp
       python outer_clmul_codegen.py BchCircuit.h --verify-only

Compile with AVX512F/BW/VL/VBMI, GFNI and VPCLMULQDQ. The ordinary coefficient
exports retain outer_transpose_codegen's exact map. RawPolynomial/RawClmul
instead return the explicitly different Q/P message coordinates, with the
same rank-128 BCH code image. The existing raw_outer_probe may link unchanged.

Each native qword represents 64 BCH coordinates for ONE payload bit. Four
quarter arrays therefore carry eight independent payload streams per ZMM.
The GL32 matrices differ between eight-coordinate groups, so that operation
still precedes this representation change. Six 64x64 carryless products per
payload stream compute the required middle product; the small top limb uses
four GFNI operations per ZMM. Five exceptional P coordinates use packed GFNI.
No performance claim is made: all conversion instructions are retained.
"""
from __future__ import annotations

import argparse
from contextlib import redirect_stdout
from io import StringIO
import json
from pathlib import Path
import random
import re
import sys

import outer_composed_codegen as polynomial
import outer_layout_codegen as layout
import outer_transpose_codegen as exact
import outer_unpack_codegen as unpack
import packed_bch_tune_codegen as tune

MASK = (1 << 64) - 1
R = int(f'{polynomial.Q:0133b}'[::-1], 2)
R0, R1, R2 = (R >> (64*i) & MASK for i in range(3))


def mul(a, b):
    result = 0
    while b:
        bit = b & -b
        result ^= a << (bit.bit_length() - 1)
        b ^= bit
    return result


def byte_matrix(function):
    return polynomial.matrix_from_rows([
        sum(((function(1 << j) >> i) & 1) << j for j in range(8))
        for i in range(8)])


SMALL_LO = byte_matrix(lambda x: mul(x, R2) & 255)
SMALL_HI = byte_matrix(lambda x: mul(x, R2) >> 8)


def word_affine(value, matrix):
    return sum(layout.affine(value >> (8*j) & 255, matrix) << (8*j)
               for j in range(8))


def native_q(x):
    """Literal qword form of the emitted middle-product arithmetic."""
    x0, x1, x2, x3 = x
    c1 = mul(x0, R1) ^ mul(x1, R0)
    c2 = mul(x1, R1) ^ mul(x2, R0)
    c3 = mul(x2, R1) ^ mul(x3, R0)
    lo = (c1 >> 64) ^ (c2 & MASK)
    hi = (c2 >> 64) ^ (c3 & MASK)
    a, b = word_affine(x0, SMALL_LO), word_affine(x0, SMALL_HI)
    c, d = word_affine(x1, SMALL_LO), word_affine(x1, SMALL_HI)
    lo ^= a ^ ((b << 8) & MASK)
    hi ^= c ^ ((d << 8) & MASK) ^ (b >> 56)
    parity = MASK if x3 >> 63 else 0
    return (((lo >> 4) ^ ((hi << 60) & MASK) ^ parity),
            ((hi >> 4) ^ parity) & ((1 << 59) - 1))


def transpose_bytes(values):
    t = [unpack.unpack_bytes(values[i], values[i+1], 1, high)
         for i in (0, 2, 4, 6) for high in (False, True)]
    u = [unpack.unpack_bytes(t[i], t[j], 2, high)
         for i, j in ((0, 2), (1, 3), (4, 6), (5, 7))
         for high in (False, True)]
    return [unpack.unpack_bytes(u[i], u[i+4], 4, high)
            for i in range(4) for high in (False, True)]


def inverse_network():
    original = [list(range(64*i, 64*i+64)) for i in range(8)]
    forward = transpose_bytes(original)
    inverse = [0] * 512
    for i, value in enumerate(sum(forward, [])):
        inverse[value] = i
    targets = [inverse[64*i:64*i+64] for i in range(8)]
    nodes = {f'x{i}': original[i] for i in range(8)}
    gates = []

    def gate(name, left, right, wanted):
        source = nodes[left] + nodes[right]
        lookup = {value: i for i, value in enumerate(source)}
        index = tuple(lookup[value] for value in wanted)
        nodes[name] = [source[i] for i in index]
        if nodes[name] != wanted or len(wanted) != 64:
            raise ArithmeticError('invalid inverse transpose gate')
        gates.append((name, left, right, index))

    for pair in range(4):
        for half in range(2):
            wanted = [v for target in targets[4*half:4*half+4]
                      for v in target if v // 128 == pair]
            gate(f'a{2*pair+half}', f'x{2*pair}', f'x{2*pair+1}', wanted)
    for half in range(2):
        for pair in range(4):
            wanted = [v for target in targets[2*pair:2*pair+2]
                      for v in target if v // 256 == half]
            gate(f'b{4*half+pair}', f'a{4*half+pair//2}',
                 f'a{4*half+2+pair//2}', wanted)
    for target in range(8):
        gate(f'c{target}', f'b{target//2}', f'b{4+target//2}', targets[target])
    if len(gates) != 24:
        raise ArithmeticError('inverse network size changed')
    for i in range(8):
        if [sum(forward, [])[j] for j in nodes[f'c{i}']] != original[i]:
            raise ArithmeticError('inverse transpose changed a byte')
    return gates


def p_matrices():
    rows = polynomial.raw_rows()[123:]
    return [polynomial.matrix_from_rows(
        [(r >> (8*g)) & 255 for r in rows] + [0]*3) for g in range(32)]


def verify_pipeline():
    """Replay the complete post-GL32 layout with independently random bits.

    The exhaustive component checks below prove linear-map identities. These
    dense cases additionally catch quarter/plane indexing and in-place writes
    in their composition, with distinct bits in all four physical row lanes.
    """
    rng = random.Random(0x434c4d554c)
    gates, matrices = inverse_network(), p_matrices()
    for case in range(3):
        values = [[rng.getrandbits(128) for _ in range(4)] for _ in range(256)]
        src = [[0]*4 for _ in range(256)]
        for group in range(32):
            for row in range(4):
                packed = tune.pack128([values[8*group+j][row] for j in range(8)])
                for plane, value in enumerate(packed):
                    src[64*(group//8)+8*plane+group%8][row] = value
        extra = [[0]*4 for _ in range(8)]
        for group, matrix in enumerate(matrices):
            if matrix:
                for plane in range(8):
                    for row in range(4):
                        value = src[64*(group//8)+8*plane+group%8][row]
                        extra[plane][row] ^= sum(layout.affine(
                            value >> (8*j) & 255, matrix) << (8*j) for j in range(16))
        for base in range(0, 256, 8):
            for row in range(4):
                native = tune.byte_transpose128([src[base+i][row] for i in range(8)])
                for i in range(8):
                    src[base+i][row] = native[i]
        for index in range(64):
            for row in range(4):
                source = [src[64*q+index][row] for q in range(4)]
                lo, hi = 0, 0
                for half in range(2):
                    a, b = native_q([x >> (64*half) & MASK for x in source])
                    lo |= a << (64*half)
                    hi |= b << (64*half)
                src[index][row], src[64+index][row] = lo, hi
        for base in range(0, 128, 8):
            nodes = {f'x{i}': [src[base+i][row] >> (8*j) & 255
                               for row in range(4) for j in range(16)] for i in range(8)}
            for name, left, right, index in gates:
                joined = nodes[left] + nodes[right]
                nodes[name] = [joined[j] for j in index]
            for i in range(8):
                for row in range(4):
                    src[base+i][row] = sum(nodes[f'c{i}'][16*row+j] << (8*j)
                                           for j in range(16))
        got = [[0]*128 for _ in range(4)]
        for group in range(16):
            for row in range(4):
                packed = [src[64*(group//8)+8*plane+group%8][row] for plane in range(8)]
                if group == 15:
                    packed = [value ^ (extra[plane][row] << 3)
                              for plane, value in enumerate(packed)]
                got[row][8*group:8*group+8] = tune.unpack_packed128(packed)
        for row in range(4):
            for output, form in enumerate(polynomial.raw_rows()):
                expected = 0
                for column in range(256):
                    if form >> column & 1:
                        expected ^= values[column][row]
                if got[row][output] != expected:
                    raise ArithmeticError(f'complete pipeline failed case{case}, row{row}, output{output}')


def verify(header):
    words = [int(word, 16) for word in re.findall(
        r'0x([0-9a-fA-F]+)ULL', header.read_text(encoding='utf-8'))]
    if len(words) != 512:
        raise ArithmeticError('expected exactly128 BCH rows')
    retained = [sum(words[4*i+j] << (64*j) for j in range(4)) for i in range(128)]
    rows = polynomial.raw_rows()
    raw_basis = polynomial.reduced_basis(rows)
    if len(raw_basis) != 128 or raw_basis != polynomial.reduced_basis(retained):
        raise ArithmeticError('raw message basis changed the BCH code image')
    if R2 != 0x1e or polynomial.Q.bit_count() % 2 != 1:
        raise ArithmeticError('small-limb/parity specialization no longer applies')
    matrices = p_matrices()
    if [g for g, m in enumerate(matrices) if m] != list(range(17)) + [31]:
        raise ArithmeticError('P sparse support changed')
    # All256 coordinate bases authenticate arithmetic and the P/Q seam.
    for column in range(256):
        value = 1 << column
        lo, hi = native_q([value >> (64*j) & MASK for j in range(4)])
        extra = 0
        for g, matrix in enumerate(matrices):
            extra ^= layout.affine(value >> (8*g) & 255, matrix)
        got = lo | ((hi | (extra << 59)) << 64)
        expected = sum(((r >> column) & 1) << i for i, r in enumerate(rows))
        if got != expected:
            raise ArithmeticError(f'CLMUL raw transpose failed column{column}')
    # The native transpose must concatenate groups in increasing order,
    # independently in each payload bit and each physical row lane.
    forward = transpose_bytes([list(range(64*i, 64*i+64)) for i in range(8)])
    for pair in range(8):
        for lane in range(4):
            for half in range(2):
                for group in range(8):
                    want = 64*group + 16*lane + 2*pair + half
                    if forward[pair][16*lane + 8*half + group] != want:
                        raise ArithmeticError('native qword coordinate/payload mapping changed')
    inverse_network()
    tune.verify()
    layout.verify()
    unpack.network()
    verify_pipeline()


def emit_inverse(gates):
    args = ','.join(f'__m512i& x{i}' for i in range(8))
    print(f'static SPIN_FORCEINLINE void rawClmulInverse({args}) {{')
    indices = list(dict.fromkeys(gate[3] for gate in gates))
    for i, index in enumerate(indices):
        words = [sum(index[8*j+k] << (8*k) for k in range(8)) for j in range(8)]
        print(f'const auto index{i}=_mm512_setr_epi64('
              + ','.join(f'0x{word:016x}ULL' for word in words) + ');')
    for name, left, right, index in gates:
        print(f'const auto {name}=_mm512_permutex2var_epi8('
              f'{left},index{indices.index(index)},{right});')
    for i in range(8):
        print(f'x{i}=c{i};')
    print('}')


def emit_q():
    print('static SPIN_FORCEINLINE void rawClmulQ(__m512i* src) {')
    print(f'const auto r0=_mm512_set1_epi64(0x{R0:016x}ULL);')
    print(f'const auto r1=_mm512_set1_epi64(0x{R1:016x}ULL);')
    print(f'const auto ml=_mm512_set1_epi64(0x{SMALL_LO:016x}ULL);')
    print(f'const auto mh=_mm512_set1_epi64(0x{SMALL_HI:016x}ULL);')
    print('for(unsigned index=0;index<64;++index) {')
    for quarter in range(4):
        print(f'const auto x{quarter}=_mm512_load_si512(src+{64*quarter}+index);')
    print('__m512i lo,hi;')
    for stage in range(3):
        print('{')
        for parity, immediate in (('e', '0x00'), ('o', '0x01')):
            print(f'const auto {parity}=_mm512_xor_si512('
                  f'_mm512_clmulepi64_epi128(x{stage},r1,{immediate}),'
                  f'_mm512_clmulepi64_epi128(x{stage+1},r0,{immediate}));')
        if stage == 0:
            print('lo=_mm512_unpackhi_epi64(e,o);')
        elif stage == 1:
            print('lo=_mm512_xor_si512(lo,_mm512_unpacklo_epi64(e,o));')
            print('hi=_mm512_unpackhi_epi64(e,o);')
        else:
            print('hi=_mm512_xor_si512(hi,_mm512_unpacklo_epi64(e,o));')
        print('}')
    print('const auto a=_mm512_gf2p8affine_epi64_epi8(x0,ml,0);')
    print('const auto b=_mm512_gf2p8affine_epi64_epi8(x0,mh,0);')
    print('const auto c=_mm512_gf2p8affine_epi64_epi8(x1,ml,0);')
    print('const auto d=_mm512_gf2p8affine_epi64_epi8(x1,mh,0);')
    print('lo=_mm512_ternarylogic_epi64(lo,a,_mm512_slli_epi64(b,8),0x96);')
    print('hi=_mm512_ternarylogic_epi64(hi,c,_mm512_slli_epi64(d,8),0x96);')
    print('hi=_mm512_xor_si512(hi,_mm512_srli_epi64(b,56));')
    print('const auto parity=_mm512_srai_epi64(x3,63);')
    print('const auto q0=_mm512_ternarylogic_epi64('
          '_mm512_srli_epi64(lo,4),_mm512_slli_epi64(hi,60),parity,0x96);')
    print('const auto q1=_mm512_and_si512('
          '_mm512_xor_si512(_mm512_srli_epi64(hi,4),parity),'
          '_mm512_set1_epi64(0x07ffffffffffffffULL));')
    print('_mm512_store_si512(src+index,q0);')
    print('_mm512_store_si512(src+64+index,q1);')
    print('}\n}')


def emit_raw():
    print('namespace spin::detail::kernel {\nnamespace coeff_clmul_detail {')
    emit_inverse(inverse_network())
    emit_q()
    print('alignas(64) static constexpr std::uint64_t pMatrix[32]={')
    print(','.join(f'0x{value:016x}ULL' for value in p_matrices()) + '};')
    print('}\nSPIN_NOINLINE void bchPackedCoeffRawPolynomial('
          'const block* __restrict a,block* __restrict out,'
          'const std::uint64_t* __restrict coeff) {')
    print('using namespace coeff_probe_detail;using namespace coeff_clmul_detail;')
    preparation = StringIO()
    with redirect_stdout(preparation):
        layout.emit_preparation(2, packed_layout=True)
    text = preparation.getvalue()
    # Preparation is plane-major within each64-coordinate quarter, allowing
    # every later conversion to overwrite exactly the eight vectors it read.
    for plane in range(8):
        old = f'src+8*group+{plane},mixed{plane}'
        if text.count(old) != 1:
            raise ArithmeticError('retained preparation emitter changed')
        text = text.replace(old, f'src+64*(group/8)+{8*plane}+group%8,mixed{plane}')
    print(text)
    print('alignas(64) __m512i extra[8];')
    for plane in range(8):
        print(f'auto p{plane}=_mm512_setzero_si512();')
    print('for(unsigned step=0;step<18;++step) {')
    print('const unsigned group=step==17?31:step;')
    print('const auto pm=_mm512_set1_epi64(pMatrix[group]);')
    for plane in range(8):
        print(f'p{plane}=_mm512_xor_si512(p{plane},'
              f'_mm512_gf2p8affine_epi64_epi8(_mm512_load_si512('
              f'src+64*(group/8)+{8*plane}+group%8),pm,0));')
    print('}')
    for plane in range(8):
        print(f'_mm512_store_si512(extra+{plane},p{plane});')
    print('for(unsigned base=0;base<256;base+=8) {')
    for i in range(8):
        print(f'auto v{i}=_mm512_load_si512(src+base+{i});')
    print('packedTuneByteTranspose5(' + ','.join(f'v{i}' for i in range(8)) + ');')
    for i in range(8):
        print(f'_mm512_store_si512(src+base+{i},v{i});')
    print('}\nrawClmulQ(src);')
    print('for(unsigned base=0;base<128;base+=8) {')
    for i in range(8):
        print(f'auto v{i}=_mm512_load_si512(src+base+{i});')
    print('rawClmulInverse(' + ','.join(f'v{i}' for i in range(8)) + ');')
    for i in range(8):
        print(f'_mm512_store_si512(src+base+{i},v{i});')
    print('}')
    print('for(unsigned output=0;output<16;++output) {')
    for plane in range(8):
        print(f'auto y{plane}=_mm512_load_si512('
              f'src+64*(output/8)+{8*plane}+output%8);')
        print(f'if(output==15)y{plane}=_mm512_xor_si512(y{plane},'
              f'_mm512_slli_epi64(_mm512_load_si512(extra+{plane}),3));')
    unpack.emit_fused_output(1, unpack.network())
    print('}\n}')
    print('SPIN_NOINLINE void bchPackedCoeffRawClmul(const block* a,block* out,'
          'const std::uint64_t* coeff) {bchPackedCoeffRawPolynomial(a,out,coeff);}')
    polynomial.emit_raw_reference(polynomial.raw_rows())
    print('}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('header', type=Path)
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    verify(args.header)
    print(json.dumps(dict(
        code_image_rank=128, message_basis='raw Q/P, NOT retained outputs',
        input_coordinate_bases=256, transpose_byte_labels=512, dense_four_row_tiles=3,
        native_qword_payload_streams=8, scratch_bytes=16896,
        per_tile=dict(clmul=768, small_limb_gfni=256, exceptional_p_gfni=144,
                      native_input_unpacks=768, native_inverse_vbmi=384,
                      output_gfni=128, output_vbmi=384),
        unchanged_gl32_gfni=1024, measured_speedup=False), sort_keys=True), file=sys.stderr)
    if not args.verify_only:
        print(exact.generate(args.header))
        print('#if !defined(__VPCLMULQDQ__)\n#error Raw CLMUL requires -mvpclmulqdq\n#endif')
        emit_raw()


if __name__ == '__main__':
    main()
