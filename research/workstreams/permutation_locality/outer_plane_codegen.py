"""Bound BCH register pressure by processing payload planes in smaller batches.

Exact-map experiment. The outer preparation and packed layout are unchanged.
The temporary is per outer tile, not an allocation or full-codeword pass.
"""
import argparse
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import subprocess
import sys

import packed_bch_tune_codegen as tune
import outer_layout_codegen
from packed_coeff_codegen import definition


def emit(planes):
    print('template<unsigned output,unsigned plane> static SPIN_NOINLINE void planeCompute('
          'const __m512i* __restrict src,__m512i* __restrict result) {')
    for p in range(2):
        for j in range(planes):
            y = p * planes + j
            print(f'auto y{y}=_mm512_load_si512(src+8*(output+{p})+plane+{j});')
            print(f'if constexpr(output+{p}==15)y{y}=_mm512_and_si512(y{y},_mm512_set1_epi8(0x7f));')
        print(f'constexpr auto parity{p}=[] {{std::uint64_t m=0;for(unsigned j=0;j<8;++j)'
              f'if(BchRows[8*(output+{p})+j][1]>>63)m|=std::uint64_t(0x80)<<(8*(7-j));return m;}}();')
        for j in range(planes):
            y = p * planes + j
            print(f'y{y}=_mm512_xor_si512(y{y},_mm512_gf2p8affine_epi64_epi8('
                  f'_mm512_load_si512(src+120+plane+{j}),_mm512_set1_epi64(parity{p}),0));')
    print('for(unsigned input=0;input<16;input+=2) {')
    for p in range(2):
        print(f'const auto m{p}=_mm512_set1_epi64(matrices[output+{p}][input]);')
        print(f'const auto n{p}=_mm512_set1_epi64(matrices[output+{p}][input+1]);')
    for j in range(planes):
        print(f'const auto x{j}=_mm512_load_si512(src+128+8*input+plane+{j});')
        print(f'const auto z{j}=_mm512_load_si512(src+136+8*input+plane+{j});')
        for p in range(2):
            y = p * planes + j
            print(f'y{y}=_mm512_ternarylogic_epi64(y{y},'
                  f'_mm512_gf2p8affine_epi64_epi8(x{j},m{p},0),'
                  f'_mm512_gf2p8affine_epi64_epi8(z{j},n{p},0),0x96);')
    print('}')
    for p in range(2):
        for j in range(planes):
            print(f'_mm512_store_si512(result+{8*p}+plane+{j},y{p*planes+j});')
    print('}')
    print('template<unsigned output> static SPIN_NOINLINE void packedCoeffTile('
          'const __m512i* __restrict src,block* __restrict out) {')
    print('alignas(64) __m512i result[16];')
    for plane in range(0, 8, planes):
        print(f'planeCompute<output,{plane}>(src,result);')
    for j in range(16):
        print(f'auto y{j}=_mm512_load_si512(result+{j});')
    tune.emit_output(2, packed_layout=True)
    print('}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('header', type=Path)
    parser.add_argument('--planes', type=int, choices=(1, 2, 4, 8), required=True)
    parser.add_argument('--unrolled', action='store_true',
                        help='retain all8 planes but fully unroll constant BCH input pairs')
    args = parser.parse_args()
    retained = subprocess.run([sys.executable, str(Path(outer_layout_codegen.__file__).resolve()),
                               str(args.header), '--mode', '2', '--tile-mode', '5'],
                              check=True, capture_output=True, text=True)
    old = definition(retained.stdout, 'template<unsigned output> static SPIN_NOINLINE void packedCoeffTile(')
    replacement = StringIO()
    with redirect_stdout(replacement):
        if args.unrolled:
            if args.planes != 8:
                parser.error('--unrolled requires --planes8')
            tune.emit_tile('packedCoeffTile', 2, True, True, packed_layout=True)
        else:
            emit(args.planes)
    print(retained.stdout.replace(old, replacement.getvalue()))
    print(retained.stderr, file=sys.stderr, end='')
    print(f'Exact payload partition: {args.planes} planes per helper; all8 planes covered once.', file=sys.stderr)


if __name__ == '__main__':
    main()
