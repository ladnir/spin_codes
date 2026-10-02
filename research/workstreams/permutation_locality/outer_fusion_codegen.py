"""Emit exact GL32/BCH probes that fuse the one-use systematic preparation.

Usage: python outer_fusion_codegen.py BchCircuit.h --mode {20,21,22,23,24,25}

Only the compact-coefficient entry point changes. Modes 20 and 21 use two
and one BCH output groups per tile, respectively. Mode 22 prepares both
systematic groups before a shared per-plane parity pass; mode 23 shares one
runtime-indexed two-output tile to bound instruction footprint. Mode 24 reduces
each plane with an explicit ternary before its fourth term; mode 25 factors
the four terms into even/odd rotations. The original and aligned
entry points, coefficient representation, bit packing, and BCH map are copied
unchanged from packed_coeff_codegen.py --tile-mode 1.
"""
import argparse
import contextlib
import io
from pathlib import Path
import subprocess
import sys

import packed_coeff_codegen as retained
import packed_bch_tune_codegen as tune


def emit_mix(group, names, declare, store_base=None, style='retained'):
    """The retained compact GL32 preparation, with register destinations."""
    assert len(names) == 8 and len(set(names)) == 8
    assert style in ('retained', 'ternary', 'factored')
    for j, name in enumerate(names):
        prefix = 'auto ' if declare else ''
        print(f'{prefix}{name}=_mm512_loadu_si512(a+4*(8*({group})+{j}));')
    print('orthoBlend(' + ','.join(names) + ');')
    print('const auto duplicate=_mm512_setr_epi64(0,0,1,1,2,2,3,3);')
    if style == 'factored':
        print('const auto duplicateOdd=_mm512_setr_epi64(3,3,0,0,1,1,2,2);')
    for d in range(4):
        print(f'const auto c{d}=_mm256_load_si256('
              f'reinterpret_cast<const __m256i*>(coeff+16*({group})+{4*d}));')
        duplicate = 'duplicateOdd' if style == 'factored' and d % 2 else 'duplicate'
        print(f'const auto m{d}=_mm512_permutexvar_epi64('
              f'{duplicate},_mm512_castsi256_si512(c{d}));')
    for j, name in enumerate(names):
        if style == 'retained':
            print(f'const auto z{j}={name};'
                  f'{name}=_mm512_gf2p8affine_epi64_epi8(z{j},m0,0);')
            for d, immediate in ((1, '0x39'), (2, '0x4e'), (3, '0x93')):
                print(f'{name}=_mm512_xor_si512({name},'
                      f'_mm512_gf2p8affine_epi64_epi8('
                      f'_mm512_shuffle_i32x4(z{j},z{j},{immediate}),m{d},0));')
        elif style == 'ternary':
            print('{')
            print(f'const auto z={name};')
            for d, immediate in ((0, None), (1, '0x39'), (2, '0x4e')):
                value = 'z' if immediate is None else f'_mm512_shuffle_i32x4(z,z,{immediate})'
                print(f'const auto term{d}=_mm512_gf2p8affine_epi64_epi8({value},m{d},0);')
            print(f'{name}=_mm512_ternarylogic_epi64(term0,term1,term2,0x96);')
            print(f'{name}=_mm512_xor_si512({name},_mm512_gf2p8affine_epi64_epi8('
                  '_mm512_shuffle_i32x4(z,z,0x93),m3,0));')
            print('}')
        else:
            print('{')
            print(f'const auto z={name};')
            print('const auto u=_mm512_shuffle_i32x4(z,z,0x4e);')
            print('const auto odd=_mm512_xor_si512('
                  '_mm512_gf2p8affine_epi64_epi8(z,m1,0),'
                  '_mm512_gf2p8affine_epi64_epi8(u,m3,0));')
            print(f'{name}=_mm512_ternarylogic_epi64('
                  '_mm512_gf2p8affine_epi64_epi8(z,m0,0),'
                  '_mm512_gf2p8affine_epi64_epi8(u,m2,0),'
                  '_mm512_shuffle_i32x4(odd,odd,0x39),0x96);')
            print('}')
        if store_base is not None:
            # Match the retained preparation's source schedule: release each
            # mixed plane immediately rather than keeping eight planes live.
            print(f'_mm512_store_si512({store_base}+{j},{name});')


def emit_tile(parallel):
    print('template<unsigned output> static SPIN_NOINLINE void fusedCoeffTile('
          'const block* __restrict a,block* __restrict out,'
          'const std::uint64_t* __restrict coeff,'
          'const __m512i* __restrict dense,const __m512i* __restrict paritySrc) {')
    print(f'static_assert(output+{parallel}<=16);')
    for p in range(parallel):
        names = [f'y{8*p+j}' for j in range(8)]
        print('__m512i ' + ','.join(names) + ';')
        print(f'if constexpr(output+{p}==15) {{')
        for j, name in enumerate(names):
            print(f'{name}=_mm512_and_si512(_mm512_load_si512(paritySrc+{j}),'
                  '_mm512_set1_epi8(0x7f));')
        print('} else {')
        emit_mix(f'output+{p}', names, False)
        print('}')
        print(f'constexpr auto parity{p}=[] {{std::uint64_t m=0;for(unsigned j=0;j<8;++j)'
              f'if(BchRows[8*(output+{p})+j][1]>>63)'
              'm|=std::uint64_t(0x80)<<(8*(7-j));return m;}();')
        for j, name in enumerate(names):
            print(f'{name}=_mm512_xor_si512({name},_mm512_gf2p8affine_epi64_epi8('
                  f'_mm512_load_si512(paritySrc+{j}),_mm512_set1_epi64(parity{p}),0));')
    # Reuse the exact retained paired-input schedule; the dense half has moved
    # from src[128..255] into its own 128-vector array, with no layout change.
    body = io.StringIO()
    with contextlib.redirect_stdout(body):
        tune.emit_input_step(parallel, True, 'input')
    dense_body = body.getvalue().replace('src+128+', 'dense+')
    assert 'src+' not in dense_body
    print('for(unsigned input=0;input<16;input+=2) {')
    print(dense_body, end='')
    print('}')
    tune.emit_output(parallel)
    print('}')


def emit_compact(parallel):
    print('namespace fusion_probe_detail {')
    print('using coeff_probe_detail::orthoBlend;')
    print('using coeff_probe_detail::matrices;')
    emit_tile(parallel)
    print('}')
    print('SPIN_NOINLINE void bchPackedCoeffCompact(const block* __restrict a,'
          'block* __restrict out,const std::uint64_t* __restrict coeff) {')
    print('using namespace fusion_probe_detail;')
    print('alignas(64) __m512i dense[128];')
    print('alignas(64) __m512i paritySrc[8];')
    print('{')
    names = [f'v{j}' for j in range(8)]
    emit_mix('15', names, True)
    for j, name in enumerate(names):
        print(f'_mm512_store_si512(paritySrc+{j},{name});')
    print('}')
    print('for(unsigned group=16;group<32;++group) {')
    emit_mix('group', names, True)
    for j, name in enumerate(names):
        print(f'_mm512_store_si512(dense+8*(group-16)+{j},{name});')
    print('}')
    for output in range(0, 16, parallel):
        print(f'fusedCoeffTile<{output}>(a,out,coeff,dense,paritySrc);')
    print('}')


def emit_shared_parity_tile(runtime_output, style='retained'):
    """Two-output refinement with one parity load per bitplane."""
    if runtime_output:
        print('#if defined(__GNUC__) && !defined(__clang__)\n'
              '#define SPIN_FUSION_NOCLONE __attribute__((noclone))\n'
              '#else\n#define SPIN_FUSION_NOCLONE\n#endif')
        print('static constexpr std::array<std::uint64_t,16> fusedParityMatrices=[] {')
        print('std::array<std::uint64_t,16> result{};')
        print('for(unsigned output=0;output<16;++output)for(unsigned j=0;j<8;++j)'
              'if(BchRows[8*output+j][1]>>63)'
              'result[output]|=std::uint64_t(0x80)<<(8*(7-j));')
        print('return result;}();')
    template = '' if runtime_output else 'template<unsigned output> '
    argument = ',unsigned output' if runtime_output else ''
    no_clone = ' SPIN_FUSION_NOCLONE' if runtime_output else ''
    print(template + f'static SPIN_NOINLINE{no_clone} void fusedSharedParityTile('
          'const block* __restrict a,block* __restrict out,'
          'const std::uint64_t* __restrict coeff,'
          f'const __m512i* __restrict cached{argument}) {{')
    if not runtime_output:
        print('static_assert(output+2<=16);')
    print('const auto dense=cached+8;')
    # Complete both one-use GL32 groups before opening the parity live ranges.
    for p in range(2):
        names = [f'y{8*p+j}' for j in range(8)]
        print('__m512i ' + ','.join(names) + ';')
        condition = f'if(output+{p}==15)' if runtime_output else f'if constexpr(output+{p}==15)'
        if runtime_output and p == 0:
            # This helper is private and is called only with even output
            # indices 0,2,...,14. The first group therefore cannot be 15.
            print('{')
            emit_mix('output', names, False, style=style)
            print('}')
        else:
            print(condition + ' {')
            for j, name in enumerate(names):
                print(f'{name}=_mm512_and_si512(_mm512_load_si512(cached+{j}),'
                      '_mm512_set1_epi8(0x7f));')
            print('} else {')
            emit_mix(f'output+{p}', names, False, style=style)
            print('}')
    for p in range(2):
        if runtime_output:
            print(f'const auto parity{p}=_mm512_set1_epi64(fusedParityMatrices[output+{p}]);')
        else:
            print(f'constexpr auto parityWord{p}=[] {{std::uint64_t m=0;for(unsigned j=0;j<8;++j)'
                  f'if(BchRows[8*(output+{p})+j][1]>>63)'
                  'm|=std::uint64_t(0x80)<<(8*(7-j));return m;}();')
            print(f'const auto parity{p}=_mm512_set1_epi64(parityWord{p});')
    for j in range(8):
        print('{')
        print(f'const auto z=_mm512_load_si512(cached+{j});')
        for p in range(2):
            print(f'y{8*p+j}=_mm512_xor_si512(y{8*p+j},'
                  f'_mm512_gf2p8affine_epi64_epi8(z,parity{p},0));')
        print('}')
    body = io.StringIO()
    with contextlib.redirect_stdout(body):
        tune.emit_input_step(2, True, 'input')
    dense_body = body.getvalue().replace('src+128+', 'dense+')
    assert 'src+' not in dense_body
    print('for(unsigned input=0;input<16;input+=2) {')
    print(dense_body, end='')
    print('}')
    tune.emit_output(2)
    print('}')
    if runtime_output:
        print('#undef SPIN_FUSION_NOCLONE')


def emit_contiguous_compact(runtime_output, style='retained'):
    print('namespace fusion_probe_detail {')
    print('using coeff_probe_detail::orthoBlend;')
    print('using coeff_probe_detail::matrices;')
    emit_shared_parity_tile(runtime_output, style)
    print('}')
    print('SPIN_NOINLINE void bchPackedCoeffCompact(const block* __restrict a,'
          'block* __restrict out,const std::uint64_t* __restrict coeff) {')
    print('using namespace fusion_probe_detail;')
    print('alignas(64) __m512i cached[136];')
    print('for(unsigned group=15;group<32;++group) {')
    names = [f'v{j}' for j in range(8)]
    emit_mix('group', names, True, store_base='cached+8*(group-15)', style=style)
    print('}')
    for output in range(0, 16, 2):
        if runtime_output:
            print(f'fusedSharedParityTile(a,out,coeff,cached,{output});')
        else:
            print(f'fusedSharedParityTile<{output}>(a,out,coeff,cached);')
    print('}')


def verify_schedule(parallel):
    """Check preparation and contribution multiplicities independently."""
    mixes = [0] * 32
    mixes[15] += 1
    for group in range(16, 32):
        mixes[group] += 1
    dense = [[0] * 16 for _ in range(16)]
    parity = [0] * 16
    systematic = [0] * 16
    for output in range(0, 16, parallel):
        for p in range(parallel):
            group = output + p
            if group != 15:
                mixes[group] += 1
            systematic[group] += 1
            parity[group] += 1
            for step in range(0, 16, 2):
                for inp in (step, step + 1):
                    dense[group][inp] += 1
    assert mixes == [1] * 32
    assert systematic == [1] * 16 and parity == [1] * 16
    assert all(row == [1] * 16 for row in dense)
    # The final systematic byte masks coordinate 127; the unmasked cached
    # group 15 is used independently for every output's parity contribution.
    assert (0x7f & 0x80) == 0 and (0x7f | 0x80) == 0xff


def verify_factoring():
    # R rotates the four 128-bit lanes. Moving R through a lane-dependent
    # byte map requires applying R^-1 to that map's coefficient lanes.
    # Compare (map term, coefficient lane, source lane) for arbitrary maps.
    for lane in range(4):
        original = {(d, lane, (lane+d) % 4) for d in range(4)}
        factored = {(0, lane, lane), (2, lane, (lane+2) % 4),
                    (1, ((lane+1)-1) % 4, (lane+1) % 4),
                    (3, ((lane+1)-1) % 4, (lane+3) % 4)}
        assert original == factored
    assert [3,3,0,0,1,1,2,2] == [((lane-1) % 4) for lane in range(4) for _ in range(2)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('header', type=Path)
    parser.add_argument('--mode', type=int, choices=(20, 21, 22, 23, 24, 25), required=True)
    args = parser.parse_args()
    parallel = 1 if args.mode == 21 else 2
    verify_schedule(parallel)
    if args.mode == 25:
        verify_factoring()
    tune.verify()
    reference = subprocess.run(
        [sys.executable, str(Path(retained.__file__).resolve()), str(args.header),
         '--tile-mode', '1'], check=True, capture_output=True, text=True)
    if reference.stderr:
        print(reference.stderr, file=sys.stderr, end='')
    source = reference.stdout
    if args.mode >= 23:
        source = '#include <array>\n' + source
    old = retained.definition(source, 'SPIN_NOINLINE void bchPackedCoeffCompact(')
    replacement = io.StringIO()
    with contextlib.redirect_stdout(replacement):
        if args.mode in (20, 21):
            emit_compact(parallel)
        else:
            style = {24: 'ternary', 25: 'factored'}.get(args.mode, 'retained')
            emit_contiguous_compact(args.mode >= 23, style)
    source = source.replace(old, replacement.getvalue().rstrip(), 1)
    marker = 'unsigned packedCoeffTileMode() {return 1;}'
    assert source.count(marker) == 1
    source = source.replace(marker, f'unsigned packedCoeffTileMode() {{return {args.mode};}}')
    # These two public reference entry points must stay byte-for-byte intact.
    for mode in ('Original', 'Aligned'):
        marker = f'SPIN_NOINLINE void bchPackedCoeff{mode}('
        assert retained.definition(source, marker) == retained.definition(reference.stdout, marker)
    print(source, end='')
    print(f'Fusion mode {args.mode}: {parallel} output group(s), 15 one-use '
          'systematic groups prepared in registers, 16 dense groups plus one '
          'parity group cached; all 32 GL32 groups evaluated exactly once', file=sys.stderr)


if __name__ == '__main__':
    main()
