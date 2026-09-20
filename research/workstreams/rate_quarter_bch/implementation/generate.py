"""Generate and symbolically verify the fixed quarter-rate paired AVX2 outer."""
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path[:0] = [str(HERE.parent), str(ROOT/'workstreams/bare_bch_rm2sub')]
import smaller_outer
import generate as shared


def main():
    construction = smaller_outer.construction()
    raw = [int(v, 16) for v in construction['generator_rows_hex']]
    rows = shared.reduce_rows(raw)
    assert len(rows) == 32
    # Check both inclusions, not just rank, before changing the message basis.
    def contains(word):
        for row in rows:
            if word & (row & -row):
                word ^= row
        return word == 0
    assert all(contains(word) for word in raw)
    assert contains((1 << 128)-1)
    shared.paar.DIMENSION = 128
    circuit = shared.paar.synthesize(0, 2, rows)
    forms = [1 << i for i in range(128)]
    for a, b in circuit.gates:
        forms.append(forms[a] ^ forms[b])
    rebuilt = []
    for signals in circuit.output_signals:
        word = 0
        for signal in signals:
            word ^= forms[signal]
        rebuilt.append(word)
    assert rebuilt == rows
    generated = HERE/'generated'
    generated.mkdir(exist_ok=True)
    header = ['#pragma once', '#include <cryptoTools/Common/Defines.h>',
              '#include <cstdint>', 'namespace bare_spin {',
              'inline constexpr std::uint64_t QuarterRows[32][2] = {']
    header += ['{'+','.join(f'0x{(r>>(64*j)) & ((1<<64)-1):016x}ULL' for j in range(2))+'},'
               for r in rows]
    header += ['};', 'void quarterTranspose2(const osuCrypto::block*, const osuCrypto::block*, osuCrypto::block*, osuCrypto::block*);', '}']
    (generated/'QuarterCircuit.h').write_text('\n'.join(header)+'\n', encoding='utf-8', newline='\n')
    code = ['#include "QuarterCircuit.h"', '#include <immintrin.h>',
            'namespace bare_spin {', 'using osuCrypto::block;',
            'static inline __m256i vx(__m256i a,__m256i b) {return _mm256_xor_si256(a,b);}',
            'void quarterTranspose2(const block* a,const block* b,block* x,block* y) {']
    code += [f'const auto v{i}=_mm256_set_m128i(b[{i}].mData,a[{i}].mData);' for i in range(128)]
    code += [f'const auto v{i+128}=vx(v{a},v{b});' for i, (a,b) in enumerate(circuit.gates)]
    for i, signals in enumerate(circuit.output_signals):
        code += [f'const auto o{i}={shared.expression(signals)};',
                 f'x[{i}]=block(_mm256_castsi256_si128(o{i}));',
                 f'y[{i}]=block(_mm256_extracti128_si256(o{i},1));']
    code += ['}', '}']
    (generated/'QuarterCircuit.cpp').write_text('\n'.join(code)+'\n', encoding='utf-8', newline='\n')
    bare = ROOT/'workstreams/bare_bch_rm2sub'
    inner = json.loads((bare/'generated/MANIFEST.json').read_text())['t128_s19']
    files = [Path(__file__), HERE.parent/'smaller_outer.py', HERE.parent/'SMALLER_OUTER_AUDIT.json',
             bare/'generate.py', ROOT/'scripts/probe_bch_forward_xor_circuit.py',
             bare/'Inner.h', bare/'generated/SelectedMaps.h', bare/'Spin.h', bare/'Spin.cpp',
             bare/'CMakeLists.txt', bare/'benchmark.cpp', bare/'correctness.cpp', HERE/'correctness.cpp']
    manifest = dict(outer=construction, basis='RREF of certified generator, increasing coordinate pivots',
                    generator_rows_hex=[f'{r:032x}' for r in rows],
                    transpose_xors=circuit.xor_count, dense_transpose_xors=sum(r.bit_count()-1 for r in rows),
                    symbolic_circuit_verified=True, inner=inner,
                    source_sha256={str(p.relative_to(ROOT)): shared.digest(p) for p in files},
                    generated_sha256={p.name: shared.digest(p) for p in sorted(generated.glob('*')) if p.suffix in ('.h','.cpp')})
    (generated/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({k:manifest[k] for k in ('transpose_xors','dense_transpose_xors','symbolic_circuit_verified')}))


if __name__ == '__main__':
    main()
