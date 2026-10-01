"""Emit exact-map GL32 coefficient-layout probes with one shared BCH schedule.

Usage: python packed_coeff_codegen.py BchCircuit.h [--tile-mode 0..5]
Link the resulting object alongside the retained PackedMixer.o. Only the three
bchPackedCoeff* entry points are exported; no retained symbol is redefined.
"""
import argparse
from pathlib import Path
import subprocess
import sys

import packed_bch_tune_codegen as tune


def definition(source, marker, semicolon=False):
    """Copy one retained generated definition, rejecting ambiguous markers."""
    if source.count(marker) != 1:
        raise ValueError(f'expected one retained definition: {marker}')
    start = source.index(marker)
    opening = source.index('{', start)
    depth = 1
    end = opening + 1
    while depth:
        if source[end] == '{':
            depth += 1
        elif source[end] == '}':
            depth -= 1
        end += 1
    if semicolon:
        if source[end] != ';':
            raise ValueError('missing retained declaration terminator')
        end += 1
    return source[start:end]


def emit_preparation(mode, packed_layout=False):
    print('alignas(64) __m512i src[256];')
    if mode == 'Compact':
        print('const auto duplicate=_mm512_setr_epi64(0,0,1,1,2,2,3,3);')
    print('for(unsigned group=0;group<32;++group) {')
    for j in range(8):
        print(f'auto v{j}=_mm512_loadu_si512(a+4*(8*group+{j}));')
    pack = 'packedTunePack5' if packed_layout else 'orthoBlend'
    print(pack+'(' + ','.join(f'v{j}' for j in range(8)) + ');')
    for d in range(4):
        if mode == 'Compact':
            print(f'const auto c{d}=_mm256_load_si256(reinterpret_cast<const __m256i*>(coeff+16*group+4*{d}));')
            # Only the initialized low four qwords of the cast are selected.
            print(f'const auto m{d}=_mm512_permutexvar_epi64(duplicate,_mm512_castsi256_si512(c{d}));')
        else:
            load = '_mm512_load_si512' if mode == 'Aligned' else '_mm512_loadu_si512'
            print(f'const auto m{d}={load}(coeff+32*group+8*{d});')
    for j in range(8):
        print(f'const auto z{j}=v{j};v{j}=_mm512_gf2p8affine_epi64_epi8(z{j},m0,0);')
        for d, immediate in ((1, '0x39'), (2, '0x4e'), (3, '0x93')):
            print(f'v{j}=_mm512_xor_si512(v{j},_mm512_gf2p8affine_epi64_epi8('
                  f'_mm512_shuffle_i32x4(z{j},z{j},{immediate}),m{d},0));')
        print(f'_mm512_store_si512(src+8*group+{j},v{j});')
    print('}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('header', type=Path)
    parser.add_argument('--tile-mode', type=int, choices=range(6), default=0)
    args = parser.parse_args()
    tune.verify()
    retained = subprocess.run([sys.executable, str(Path(__file__).with_name('packed_mixer_codegen.py')),
                              str(args.header)], check=True, capture_output=True, text=True)
    if retained.stderr:
        print(retained.stderr, file=sys.stderr, end='')
    # Keep the retained matrix packing and bit transpose verbatim. A private
    # namespace lets this object coexist with the original/tuned kernel object.
    print('#include "Spin.h"\n#include "generated/BchCircuit.h"')
    print('namespace spin::detail::kernel { namespace coeff_probe_detail {')
    print(definition(retained.stdout, 'alignas(64) static constexpr std::uint64_t matrices', True))
    print(definition(retained.stdout, 'static SPIN_FORCEINLINE void orthoBlend('))
    if args.tile_mode == 5:
        tune.emit_byte_transpose_helpers()
    if args.tile_mode == 0:
        print(definition(retained.stdout, 'template<unsigned output> static SPIN_NOINLINE void packedFullTile(')
              .replace('packedFullTile(', 'packedCoeffTile('))
        parallel = 2
    else:
        parallel, paired, unroll, _ = tune.VARIANTS[args.tile_mode]
        tune.emit_tile('packedCoeffTile', parallel, paired, unroll, packed_layout=args.tile_mode == 5)
        if 16 % parallel:
            tune.emit_tile('packedCoeffTileTail', 16 % parallel, paired, unroll)
    print('}')
    for mode in ('Original', 'Aligned', 'Compact'):
        print(f'SPIN_NOINLINE void bchPackedCoeff{mode}(const block* __restrict a,'
              'block* __restrict out,const std::uint64_t* __restrict coeff) {')
        print('using namespace coeff_probe_detail;')
        emit_preparation(mode, packed_layout=args.tile_mode == 5)
        for output in range(0, 16, parallel):
            name = 'packedCoeffTile' if output + parallel <= 16 else 'packedCoeffTileTail'
            print(f'{name}<{output}>(src,out);')
        print('}')
    print(f'unsigned packedCoeffTileMode() {{return {args.tile_mode};}}')
    print('}')
    print(f'Coefficient probes: original/aligned/compact; common BCH tile mode {args.tile_mode}; '
          'compact expansion selects [0,0,1,1,2,2,3,3]', file=sys.stderr)


if __name__ == '__main__':
    main()
