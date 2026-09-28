"""Emit bounded-register GFNI tiles for the exact existing BCH transpose.

The 128 dense input coordinates are bit-transposed eight at a time. Each
GFNI matrix implements one 8-by-8 binary submatrix of the fixed generator.
This generator verifies the matrix packing and bit transpose exhaustively.
"""
from pathlib import Path
import re
import sys


def ortho(x):
    x = x.copy()
    for shift,mask in ((1,0x5555555555555555),(2,0x3333333333333333),(4,0x0f0f0f0f0f0f0f0f)):
        for a in range(8):
            if not a & shift:
                b = a+shift
                t = ((x[a] >> shift)^x[b]) & mask
                x[a] ^= t << shift
                x[b] ^= t
    return x


def main():
    words = [int(v,16) for v in re.findall(r'0x([0-9a-f]+)ULL',Path(sys.argv[1]).read_text())]
    assert len(words)==512
    rows = [sum(words[4*i+j] << (64*j) for j in range(4)) for i in range(128)]
    for r,row in enumerate(rows):
        assert row & ((1 << 127)-1) == (1 << r if r<127 else 0)
    for row in range(8):
        for bit in range(64):
            x = [0]*8
            x[row] = 1 << bit
            y = [0]*8
            y[bit%8] = 1 << ((bit//8)*8+row)
            assert ortho(x)==y and ortho(y)==x
    matrices = []
    for output in range(16):
        line=[]
        for inp in range(16):
            matrix = sum(((rows[8*output+j] >> (128+8*inp))&255) << (8*(7-j)) for j in range(8))
            for value in range(256):
                actual = sum((((value & (matrix >> (8*(7-j)))).bit_count()&1) << j) for j in range(8))
                expected = sum((((value & (rows[8*output+j] >> (128+8*inp)) & 255).bit_count()&1) << j) for j in range(8))
                assert actual==expected
            line.append(matrix)
        matrices.append(line)
    print('GFNI: checked 65536 submatrix actions and 512 transpose basis vectors',file=sys.stderr)
    print('#include "Spin.h"\n#include "generated/BchCircuit.h"\nnamespace spin::detail::kernel {')
    print('alignas(64) static constexpr std::uint64_t matrices[16][16]={')
    for line in matrices:
        print('{'+','.join(f'0x{x:016x}ULL' for x in line)+'},')
    print('};')
    for blend in (False,True):
        suffix = 'Blend' if blend else ''
        print(f'static SPIN_FORCEINLINE void ortho{suffix}(__m512i& x0,__m512i& x1,__m512i& x2,__m512i& x3,__m512i& x4,__m512i& x5,__m512i& x6,__m512i& x7) {{')
        for shift,mask in ((1,0x5555555555555555),(2,0x3333333333333333),(4,0x0f0f0f0f0f0f0f0f)):
            print(f'const auto mask{shift}=_mm512_set1_epi64(0x{mask:016x}ULL);')
            for a in range(8):
                if not a & shift:
                    b=a+shift
                    if blend:
                        print(f'{{const auto lo=_mm512_srli_epi64(x{a},{shift});const auto hi=_mm512_slli_epi64(x{b},{shift});')
                        print(f'x{a}=_mm512_ternarylogic_epi64(x{a},hi,mask{shift},0xe4);x{b}=_mm512_ternarylogic_epi64(lo,x{b},mask{shift},0xe4);}}')
                    else:
                        print(f'{{const auto z=_mm512_and_si512(_mm512_xor_si512(_mm512_srli_epi64(x{a},{shift}),x{b}),mask{shift});')
                        print(f'x{a}=_mm512_xor_si512(x{a},_mm512_slli_epi64(z,{shift}));x{b}=_mm512_xor_si512(x{b},z);}}')
        print('}')
    for parallel,name,specialized,paired,unroll in ((1,'1',False,False,False),(2,'2',False,False,False),
            (2,'Static',True,False,False),(2,'Ternary',True,True,False),(2,'Unroll',True,True,True),(2,'Blend',True,True,False),
            (2,'Mapped',True,True,False),(2,'Packed',True,True,False),(1,'Mapped1',True,True,False)):
        mapped = name in ('Mapped','Packed','Mapped1')
        ortho_name = 'orthoBlend' if name in ('Blend','Mapped','Packed','Mapped1') else 'ortho'
        if specialized:
            print('template<unsigned output>')
        argument = '' if specialized else ',unsigned output'
        if mapped:
            argument += ',const unsigned* __restrict permutation'
        print(f'static SPIN_NOINLINE void tile{name}(const __m512i* __restrict src,const block* __restrict a,block* __restrict out{argument}) {{')
        for j in range(8*parallel):
            print(f'auto y{j}=_mm512_setzero_si512();')
        if unroll:
            print('#if defined(__GNUC__)\n#pragma GCC unroll 8\n#endif')
        print(f'for(unsigned input=0;input<16;input+={2 if paired else 1}) {{')
        for p in range(parallel):
            print(f'const auto m{p}=_mm512_set1_epi64(matrices[output+{p}][input]);')
            if paired:
                print(f'const auto n{p}=_mm512_set1_epi64(matrices[output+{p}][input+1]);')
        for j in range(8):
            print(f'const auto x{j}=_mm512_load_si512(src+8*input+{j});')
            if paired:
                print(f'const auto z{j}=_mm512_load_si512(src+8*input+8+{j});')
            for p in range(parallel):
                term = f'_mm512_gf2p8affine_epi64_epi8(x{j},m{p},0)'
                if paired:
                    term2 = f'_mm512_gf2p8affine_epi64_epi8(z{j},n{p},0)'
                    print(f'y{8*p+j}=_mm512_ternarylogic_epi64(y{8*p+j},{term},{term2},0x96);')
                else:
                    print(f'y{8*p+j}=_mm512_xor_si512(y{8*p+j},{term});')
        print('}')
        for p in range(parallel):
            print(ortho_name+'('+','.join(f'y{8*p+j}' for j in range(8))+');')
            for j in range(8):
                print(f'{{const unsigned row=8*(output+{p})+{j};')
                # These branches depend only on the fixed matrix, not data.
                offset = 'permutation[row]' if mapped else '4*row'
                parity = 'permutation[127]' if mapped else '4*127'
                print(f'if(row<127)y{8*p+j}=_mm512_xor_si512(y{8*p+j},_mm512_loadu_si512(a+{offset}));')
                print(f'if(BchRows[row][1]>>63)y{8*p+j}=_mm512_xor_si512(y{8*p+j},_mm512_loadu_si512(a+{parity}));')
                if name!='Packed':
                    for lane in range(4):
                        print(f'out[{128*lane}+row]=block(_mm512_extracti32x4_epi32(y{8*p+j},{lane}));')
                print('}')
            if name=='Packed':
                for j in (0,4):
                    at=8*p+j
                    print('{')
                    print(f'const auto a=_mm512_shuffle_i64x2(y{at},y{at+1},0x44);')
                    print(f'const auto b=_mm512_shuffle_i64x2(y{at},y{at+1},0xee);')
                    print(f'const auto c=_mm512_shuffle_i64x2(y{at+2},y{at+3},0x44);')
                    print(f'const auto d=_mm512_shuffle_i64x2(y{at+2},y{at+3},0xee);')
                    for lane,(left,right,imm) in enumerate((('a','c','0x88'),('a','c','0xdd'),('b','d','0x88'),('b','d','0xdd'))):
                        print(f'_mm512_storeu_si512(out+{128*lane}+8*(output+{p})+{j},_mm512_shuffle_i64x2({left},{right},{imm}));')
                    print('}')
        print('}')
        argument = ',const unsigned* __restrict permutation' if mapped else ''
        print(f'SPIN_NOINLINE void bchTranspose4Gfni{name}(const block* __restrict a,block* __restrict x{argument}) {{')
        print('alignas(64) __m512i src[128];\nfor(unsigned group=0;group<16;++group) {')
        for j in range(8):
            offset = f'permutation[128+8*group+{j}]' if mapped else f'4*(128+8*group+{j})'
            print(f'auto v{j}=_mm512_loadu_si512(a+{offset});')
        print(ortho_name+'('+','.join(f'v{j}' for j in range(8))+');')
        for j in range(8):
            print(f'_mm512_store_si512(src+8*group+{j},v{j});')
        print('}')
        if specialized:
            for output in range(0,16,parallel):
                arguments = ',permutation' if mapped else ''
                print(f'tile{name}<{output}>(src,a,x{arguments});')
        else:
            print(f'for(unsigned output=0;output<16;output+={parallel})tile{name}(src,a,x,output);')
        print('}')
    # Prepare the sparse/systematic inputs too, so output tiles never revisit
    # the scattered source or its permutation. Only 8 KiB of extra copying.
    print('SPIN_NOINLINE void bchTranspose4GfniPrepared(const block* __restrict a,block* __restrict x,const unsigned* __restrict permutation) {')
    print('alignas(64) __m512i src[128];alignas(64) block systematic[512];')
    print('for(unsigned group=0;group<16;++group) {')
    for j in range(8):
        print(f'auto v{j}=_mm512_loadu_si512(a+permutation[128+8*group+{j}]);')
    print('orthoBlend('+','.join(f'v{j}' for j in range(8))+');')
    for j in range(8):
        print(f'_mm512_store_si512(src+8*group+{j},v{j});')
    print('}\nfor(unsigned row=0;row<128;row+=4) {')
    for j in range(4):
        print(f'const auto v{j}=_mm512_loadu_si512(a+permutation[row+{j}]);')
    for j in range(4):
        print(f'_mm512_store_si512(systematic+4*(row+{j}),v{j});')
    print('}')
    for output in range(0,16,2):
        print(f'tileBlend<{output}>(src,systematic,x);')
    print('}')
    # Each row may use its own coordinate permutation. Gather only within
    # the current 16 KiB four-row tile, then reuse the unchanged GFNI tiles.
    print('static SPIN_FORCEINLINE __m512i independentColumn(const block* a,const std::uint16_t* offsets,unsigned c) {')
    print('auto x=_mm512_castsi128_si512(_mm_load_si128(reinterpret_cast<const __m128i*>(a+4*c)));')
    for lane in range(1,4):
        print(f'x=_mm512_inserti32x4(x,_mm_load_si128(reinterpret_cast<const __m128i*>(a+offsets[4*c+{lane}])),{lane});')
    print('return x;\n}')
    print('SPIN_NOINLINE void bchTranspose4GfniIndependent(const block* __restrict a,block* __restrict x,const std::uint16_t* __restrict offsets) {')
    print('alignas(64) __m512i src[128];alignas(64) block systematic[512];')
    print('for(unsigned group=0;group<16;++group) {')
    for j in range(8):
        print(f'auto v{j}=independentColumn(a,offsets,128+8*group+{j});')
    print('orthoBlend('+','.join(f'v{j}' for j in range(8))+');')
    for j in range(8):
        print(f'_mm512_store_si512(src+8*group+{j},v{j});')
    print('}\nfor(unsigned row=0;row<128;row+=4) {')
    for j in range(4):
        print(f'const auto v{j}=independentColumn(a,offsets,row+{j});')
    for j in range(4):
        print(f'_mm512_store_si512(systematic+4*(row+{j}),v{j});')
    print('}')
    for output in range(0,16,2):
        print(f'tileBlend<{output}>(src,systematic,x);')
    print('}\n}')


if __name__ == '__main__':
    main()
