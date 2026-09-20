"""Fixed certified K16 map and symbolically verified, allocation-free XOR circuits.

No search receipts, NumPy, or proof runtime are needed by the encoder build.
Columns use the forward convention y=x+Aq, q'=Fq+Bx.
"""
from functools import reduce
from operator import xor
from pathlib import Path
import sys

A = [820,821,822,823,816,2353,2866,307,828,2813,2110,511,696,377,442,635,
     804,1893,1638,551,2144,1569,1314,2915,1516,2157,2990,1583,3880,2217,2410,3819,
     788,3925,3798,663,3408,2833,2194,3795,1692,797,94,1503,2392,1753,1946,2075,
     3140,1093,1222,3271,2368,2881,2498,3011,3852,3789,3214,3407,2952,73,10,3019]
B = [3539,1337,2103,985,3260,3139,192,1268,3152,3087,3140,714,1251,1667,4091,1304,
     2598,1412,2608,696,2536,1015,2673,3951,156,153,3694,3958,774,831,3879,2212,
     1657,3189,2095,2363,2120,1133,93,1479,2535,2734,3115,3728,926,4090,587,4010,
     1998,893,1921,799,396,1937,1426,3892,1192,1449,745,446,593,3064,3785,868]


def circuit(targets, inputs):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
    import probe_bch_forward_xor_circuit as paar
    paar.DIMENSION = len(inputs)
    c = paar.synthesize(0, 2, targets)
    forms = [1 << j for j in range(len(inputs))]
    lines = [f'const auto v{j}={term};' for j, term in enumerate(inputs)]
    for j, (a, b) in enumerate(c.gates, len(inputs)):
        forms.append(forms[a] ^ forms[b])
        lines.append(f'const auto v{j}=vx(v{a},v{b});')
    if [reduce(xor, (forms[j] for j in out), 0) for out in c.output_signals] != targets:
        raise RuntimeError('incorrect XOR circuit')
    def expression(signals):
        level = [f'v{j}' for j in signals]
        if not level:
            return '_mm_setzero_si128()'
        while len(level) > 1:
            level = [f'vx({level[j]},{level[j+1]})' if j+1 < len(level) else level[j]
                     for j in range(0, len(level), 2)]
        return level[0]
    return lines, [expression(x) for x in c.output_signals]


def header():
    coefficients = A.copy()
    for bit in range(6):
        for x in range(64):
            if x >> bit & 1:
                coefficients[x] ^= coefficients[x ^ (1 << bit)]
    if any(c and x.bit_count() > 2 for x, c in enumerate(coefficients)):
        raise RuntimeError('expansion is not degree two')
    for x in range(64):
        if reduce(xor, (c for m, c in enumerate(coefficients) if m & x == m), 0) != A[x]:
            raise RuntimeError('incorrect expansion transform')
    monomials = [x for x in range(64) if x.bit_count() <= 2]
    targets = [sum(((coefficients[x] >> j) & 1) << i for i, x in enumerate(monomials)) for j in range(12)]
    finish, outs = circuit(targets, [f'z[{x}]' for x in monomials])
    emission, emits = circuit(B, [f'state[{j}]' for j in range(12)])
    feedback, feedback_out = circuit([sum(((c >> j) & 1) << p for p, c in enumerate(B))
                                     for j in range(12)], [f'x[{p}]' for p in range(64)])
    text = ['#pragma once', '// Generated from the fixed certified K16 map; forward A/B convention.',
            'namespace bare_spin {', 'struct Map64S12 {',
            'static OC_FORCEINLINE __m128i vx(__m128i a,__m128i b) {return _mm_xor_si128(a,b);}',
            'static constexpr unsigned T=64,S=12;',
            '// Unused by IMT; only needed to parse the legacy field-template fallback.',
            'static constexpr u32 modulusLow=0;',
            'static constexpr const char* name="imt_t64_s12_subspace_v1";']
    for name, values in [('columns', A), ('feedbackColumns', B), ('groupedColumns', A),
                         ('groupOrder', list(range(12)))]:
        text.append(f'static constexpr std::array<u32,{len(values)}> {name}{{'+','.join(map(str, values))+'};')
    text += ['static OC_FORCEINLINE void finish(const __m128i* z,__m128i* out) {'] + finish
    text += [f'out[{j}]={v};' for j, v in enumerate(outs)] + ['}']
    text += ['static OC_FORCEINLINE void feedback(const __m128i* x,__m128i* out) {'] + feedback
    text += [f'out[{j}]={v};' for j, v in enumerate(feedback_out)] + ['}']
    text += ['template<class Emit> static OC_FORCEINLINE void emitShared(const block* in,__m128i* raw,',
             'const __m128i* state,std::size_t base,Emit& emit) {'] + emission
    for p in reversed(range(64)):
        text += [f'raw[{p}]=in[{p}].mData;', f'emit(base+{p},block(vx(raw[{p}],{emits[p]})));']
    text += ['}', '};', '// Adapter for the historical transpose-only interface.',
             'struct Map64S12Transpose : Map64S12 {',
             'static constexpr auto columns=Map64S12::feedbackColumns;',
             'static constexpr auto feedbackColumns=Map64S12::columns;',
             '};', '}']
    return '\n'.join(text)+'\n'


if __name__ == '__main__':
    Path(sys.argv[1]).write_text(header(), newline='\n')
