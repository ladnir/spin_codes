"""Generate selected RM2Sub maps and a verified Q <= C <= P BCH XOR circuit."""
import hashlib
import json
import random
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BRIDGE = ROOT / 'workstreams/bch_rm2sub_bridge'
sys.path[:0] = [str(ROOT / 'scripts'), str(ROOT / 'bch_spectrum_work/bch_spectrum_codex_bundle/code')]
import bch_quotient as bch
import probe_bch_forward_xor_circuit as paar


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reduce_rows(rows):
    rows = list(rows)
    pivot = 0
    for j in range(256):
        found = next((i for i in range(pivot, len(rows)) if (rows[i] >> j) & 1), None)
        if found is None:
            continue
        rows[pivot], rows[found] = rows[found], rows[pivot]
        for i in range(len(rows)):
            if i != pivot and (rows[i] >> j) & 1:
                rows[i] ^= rows[pivot]
        pivot += 1
        if pivot == len(rows):
            return rows
    raise ValueError('Dependent rows')


def modulus(s):
    def rem(x, p):
        while x.bit_length() >= p.bit_length():
            x ^= p << (x.bit_length() - p.bit_length())
        return x
    def gcd(a, b):
        while b:
            a, b = b, rem(a, b)
        return a
    for low in range(3, 1 << s, 2):
        p = (1 << s) | low
        x, valid = 2, True
        for i in range(1, s + 1):
            x = rem(sum(((x >> j) & 1) << (2*j) for j in range(s)), p)
            if i <= s//2 and gcd(x ^ 2, p) != 1:
                valid = False
                break
        if valid and x == 2:
            return low
    raise ValueError('No irreducible polynomial')


def expression(signals, prefix='v'):
    level = [f'{prefix}{i}' for i in signals]
    if not level:
        return '_mm_setzero_si128()'
    while len(level) > 1:
        level = [f'vx({level[i]}, {level[i+1]})' if i+1 < len(level) else level[i]
                 for i in range(0, len(level), 2)]
    return level[0]


def grouping(columns, s):
    # Minimize nonzero nibble lookups in the fully unrolled A emission.
    # Cache each group's cost; swapping bits changes only two small groups.
    cache = {}
    def score(order):
        total = 0
        for start in range(0, s, 4):
            mask = sum(1 << b for b in order[start:start+4])
            if mask not in cache:
                cache[mask] = sum(bool(c & mask) for c in columns)
            total += cache[mask]
        return total
    best = list(range(s))
    initial = score(best)
    value = initial
    rng = random.Random(0)
    for restart in range(8):
        order = list(range(s))
        if restart:
            rng.shuffle(order)
        cost = score(order)
        for _ in range(2000):
            a, b = rng.sample(range(s), 2)
            order[a], order[b] = order[b], order[a]
            candidate = score(order)
            if candidate <= cost:
                cost = candidate
            else:
                order[a], order[b] = order[b], order[a]
            if cost < value:
                best, value = order.copy(), cost
    return best, initial, value


def main():
    generated = HERE / 'generated'
    generated.mkdir(exist_ok=True)
    records = {}
    code = ['#pragma once', '#include <array>', '#include <cstdint>',
            '#include <immintrin.h>', 'namespace bare_spin {']
    for t, s in ((64, 16), (64, 20), (128, 19), (256, 14)):
        name = f't{t}_s{s}'
        directory = BRIDGE / ('generated/larger_state_inputs_v1' if name in ('t64_s20', 't128_s19') else 'inputs')
        path = directory / f'{name}_selection.json'
        manifest = json.loads((directory / 'manifest.json').read_text())
        assert digest(path) == manifest['snapshot_sha256'][path.name]
        data = json.loads(path.read_text())
        m = data['parameters']['variables']
        selected = data['selected']
        columns = [int(c, 16) for c in selected['B_columns_hex']]
        masks = [int(q, 16) for q in selected['quadratic_masks_hex']]
        mono = [(1 << i) | (1 << j) for i in range(m) for j in range(i+1, m)]
        rebuilt = []
        for x in range(t):
            q = sum(int((x & mask) == mask) << i for i, mask in enumerate(mono))
            rebuilt.append(1 | (x << 1) | sum(((q & mask).bit_count() & 1) << (m+1+i)
                                                for i, mask in enumerate(masks)))
        assert columns == rebuilt
        rows = [sum(((c >> j) & 1) << i for i, c in enumerate(columns)) for j in range(s)]
        assert rows == [int(r, 16) for r in selected['A_generator_words_hex']]
        assert all(not ((a & b).bit_count() & 1) for a in rows for b in rows)
        # Symbolically audit the pruned zeta dependency rule used by Inner.h.
        z = [1 << i for i in range(t)]
        distance = t//2
        while distance:
            for base in range(0, t, 2*distance):
                if (base//(2*distance)).bit_count() <= 2:
                    for j in range(distance):
                        z[base+j] ^= z[base+j+distance]
            distance //= 2
        for monomial in [0, *[1 << i for i in range(m)], *mono]:
            assert z[monomial] == sum(1 << x for x in range(t) if x & monomial == monomial)
        paar.DIMENSION = len(mono)
        circuit = paar.synthesize(0, 2, masks)
        low = modulus(s)
        order, before, after = grouping(columns, s)
        grouped_columns = [sum(((c >> b) & 1) << j for j, b in enumerate(order)) for c in columns]
        assert sorted(order) == list(range(s))
        assert columns == [sum(((c >> j) & 1) << b for j, b in enumerate(order)) for c in grouped_columns]
        code += [f'struct Map{t}S{s} {{', f'static constexpr unsigned T={t}, S={s}, M={m};',
                 f'static constexpr std::uint32_t modulusLow=0x{low:x};',
                 f'static constexpr const char* name="{name}";',
                 'static constexpr std::array<std::uint32_t,T> columns{' + ','.join(hex(c) for c in columns) + '};',
                 'static constexpr std::array<unsigned,S> groupOrder{' + ','.join(map(str,order)) + '};',
                 'static constexpr std::array<std::uint32_t,T> groupedColumns{' + ','.join(hex(c) for c in grouped_columns) + '};',
                 'static inline __m128i vx(__m128i a,__m128i b) { return _mm_xor_si128(a,b); }',
                 'static inline void finish(const __m128i* z,__m128i* out) {', 'out[0]=z[0];']
        code += [f'out[{i+1}]=z[{1<<i}];' for i in range(m)]
        code += [f'const auto v{i}=z[{p}];' for i, p in enumerate(mono)]
        code += [f'const auto v{i+len(mono)}=vx(v{a},v{b});' for i, (a,b) in enumerate(circuit.gates)]
        code += [f'out[{i+m+1}]={expression(signals)};' for i, signals in enumerate(circuit.output_signals)]
        code += ['}', '};']
        records[name] = dict(selection_sha256=digest(path), modulus_low=hex(low),
                             quadratic_xors=circuit.xor_count, columns=columns,
                             A_group_order=order, A_table_lookups_before=before, A_table_lookups_after=after)
    code += ['}']
    (generated / 'SelectedMaps.h').write_text('\n'.join(code)+'\n', encoding='utf-8')

    p, q = bch.generator_polynomial(37), bch.generator_polynomial(39)
    assert p.bit_length() == 125 and q.bit_length() == 133
    assert bch.binary_poly_divmod(q, p)[1] == 0
    def extend(word):
        return word | ((word.bit_count() & 1) << 255)
    raw = [extend(q << i) for i in range(123)] + [extend(p << i) for i in range(5)]
    rows = reduce_rows(raw)
    assert len(rows) == 128
    assert all(bch.binary_poly_divmod(w & ((1 << 255)-1), p)[1] == 0 for w in rows)
    # RREF membership establishes Q containment and preservation of all-one word.
    def contains(w):
        for r in rows:
            pivot = (r & -r).bit_length()-1
            if (w >> pivot) & 1:
                w ^= r
        return w == 0
    assert all(contains(extend(q << i)) for i in range(123))
    assert contains((1 << 256)-1)
    paar.DIMENSION = 256
    circuit = paar.synthesize(0, 2, rows)
    code = ['#pragma once', '#include <immintrin.h>', '#include <cstdint>', 'namespace bare_spin {',
            'inline constexpr std::uint64_t BchRows[128][4] = {']
    code += ['{' + ','.join(f'0x{(r >> (64*j)) & ((1<<64)-1):016x}ULL' for j in range(4)) + '},' for r in rows]
    code += ['};', 'void bchTranspose2(const osuCrypto::block*, const osuCrypto::block*, osuCrypto::block*, osuCrypto::block*);', '}']
    (generated / 'BchCircuit.h').write_text('\n'.join(code)+'\n', encoding='utf-8')
    code = ['#include "Spin.h"', '#include "generated/BchCircuit.h"', 'namespace bare_spin {',
            'static inline __m256i vx(__m256i a,__m256i b) { return _mm256_xor_si256(a,b); }',
            'void bchTranspose2(const block* a,const block* b,block* x,block* y) {']
    code += [f'const auto v{i}=_mm256_set_m128i(b[{i}].mData,a[{i}].mData);' for i in range(256)]
    code += [f'const auto v{i+256}=vx(v{a},v{b});' for i,(a,b) in enumerate(circuit.gates)]
    for i, signals in enumerate(circuit.output_signals):
        code += [f'const auto o{i}={expression(signals)};', f'x[{i}]=block(_mm256_castsi256_si128(o{i}));',
                 f'y[{i}]=block(_mm256_extracti128_si256(o{i},1));']
    code += ['}', '}']
    (generated / 'BchCircuit.cpp').write_text('\n'.join(code)+'\n', encoding='utf-8')
    records['bch'] = dict(description='Parity-extended Q plus p*x^j, j=0..4; systematic basis',
        p_generator=hex(p), q_generator=hex(q), dimension=128, Q_containment=True, P_containment=True,
        generator_rows_sha256=hashlib.sha256(b''.join(r.to_bytes(32,'little') for r in rows)).hexdigest(),
        transpose_xors=circuit.xor_count)
    records['source_sha256'] = {str(p.relative_to(ROOT)): digest(p) for p in
        [Path(__file__), ROOT/'scripts/probe_bch_forward_xor_circuit.py',
         ROOT/'bch_spectrum_work/bch_spectrum_codex_bundle/code/bch_quotient.py',
         ROOT/'bch_spectrum_work/bch_spectrum_codex_bundle/code/affine_wambach.py']}
    records['generated_sha256'] = {p.name: digest(p) for p in sorted(generated.glob('*.h'))}
    records['generated_sha256']['BchCircuit.cpp'] = digest(generated/'BchCircuit.cpp')
    (generated / 'MANIFEST.json').write_text(json.dumps(records, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: {a:b for a,b in v.items() if a != 'columns'} for k,v in records.items()
                      if k not in ('generated_sha256', 'source_sha256')}, indent=2))


if __name__ == '__main__':
    main()
