"""Generate an authenticated, fixed-schedule packet t128/s16 inner circuit.

Only the new T128Map.h is generated; the pinned calibration map is read-only.
The two 64-output halves share one old state and one physical state update.
"""
from hashlib import sha256
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / 'research/workstreams/permutation_locality'))
import inner_packet_codegen as cse

SOURCE = ROOT / 'research/workstreams/rate_quarter_bch/inner_calibration/maps/t128_s20_nested.json'
MONOMIALS = tuple(m for m in range(128) if m.bit_count() <= 2)
HIGH = tuple(m for m in range(16) if m.bit_count() <= 2)


def anf(row):
    values = [(row >> p) & 1 for p in range(128)]
    for bit in (1, 2, 4, 8, 16, 32, 64):
        for p in range(128):
            if p & bit:
                values[p] ^= values[p ^ bit]
    assert not any(values[m] for m in range(128) if m.bit_count() > 2)
    return sum(values[m] << j for j, m in enumerate(MONOMIALS))


def xor(terms):
    assert terms
    value = terms[0]
    for term in terms[1:]:
        value = f'_mm512_xor_si512({value},{term})'
    return value


def half(upper):
    code = ['{']
    for h in HIGH:
        m = 4 * h
        base = f'c{m}'
        if upper and h.bit_count() <= 1:
            base = f'vx({base},c{64+m})'
        if h.bit_count() == 2:
            code.append(f'const auto d{h}=_mm512_broadcast_i32x4({base});')
        else:
            c1, c2 = f'c{m+1}', f'c{m+2}'
            if upper and h == 0:
                c1, c2 = f'vx({c1},c65)', f'vx({c2},c66)'
            code.append(f'const auto l{h}0={base};')
            code.append(f'const auto l{h}1=vx(l{h}0,{c1});')
            code.append(f'const auto l{h}2=vx(l{h}0,{c2});')
            last = f'vx(l{h}1,{c2})'
            if h == 0:
                last = f'vx({last},c3)'
            code.append(f'const auto l{h}3={last};')
            code.append(f'const auto d{h}=join(l{h}0,l{h}1,l{h}2,l{h}3);')
    for group in (3, 2, 1, 0):
        code += [f'__m512i g{group}0,g{group}1,g{group}2,g{group}3;', '{']
        for j in range(3):
            terms = [f'd{j}'] + [f'd{j+4*b}' for b in (1, 2) if group & b]
            if group == 3 and j == 0:
                terms.append('d12')
            code.append(f'const auto b{j}={xor(terms)};')
        code += ['const auto e1=_mm512_xor_si512(b0,b1);',
                 'const auto e2=_mm512_xor_si512(b0,b2);',
                 'const auto e3=_mm512_xor_si512(_mm512_xor_si512(e1,b2),d3);']
        for low in (3, 2, 1, 0):
            h = 16 * upper + 4 * group + low
            expansion = f'e{low}' if low else 'b0'
            code += [f'const auto y{low}=_mm512_xor_si512({expansion},_mm512_loadu_si512(raw+{4*h}));',
                     f'emit(packetBase+{h},y{low});']
        code += [f'g{group}1=_mm512_xor_si512(y1,y3);',
                 f'g{group}2=_mm512_xor_si512(y2,y3);',
                 f'g{group}0=_mm512_ternarylogic_epi64(y0,y2,g{group}1,0x96);',
                 f'g{group}3=y3;', '}']
        if group == 2:
            code += [f'const auto cd{j}=_mm512_xor_si512(g2{j},g3{j});' for j in range(4)]
    values = {}
    for j in range(3):
        values[j] = f'_mm512_ternarylogic_epi64(g0{j},g1{j},cd{j},0x96)'
        values[4+j] = f'_mm512_xor_si512(g1{j},g3{j})'
        values[8+j] = f'cd{j}'
    values[3] = '_mm512_ternarylogic_epi64(g03,g13,cd3,0x96)'
    values[12] = 'g30'
    for h in HIGH:
        if upper:
            code.append(f'_mm512_store_si512(upperHigh+{h},{values[h]});')
            if h.bit_count() <= 1:
                code.append(f'storeLowMoments<{4*h},{1-h.bit_count()}>({values[h]},moments+64);')
        else:
            code.append(f'storeLowMoments<{4*h},{2-h.bit_count()}>(_mm512_xor_si512({values[h]},_mm512_load_si512(upperHigh+{h})),moments);')
    code.append('}')
    return code


def main():
    raw = SOURCE.read_bytes()
    rows = [int(row, 16) for row in json.loads(raw)['generator_rows_hex'][:16]]
    assert len(rows) == 16
    assert not any((a & b).bit_count() & 1 for a in rows for b in rows)
    anfs = [anf(row) for row in rows]
    coefficients = {m: sum(((row >> j) & 1) << i for i, row in enumerate(anfs))
                    for j, m in enumerate(MONOMIALS)}
    columns = [sum(((row >> p) & 1) << j for j, row in enumerate(rows)) for p in range(128)]
    # Authenticate every emitted packet expression against all 16 state bits.
    for upper in (0, 1):
        d = {}
        for h in HIGH:
            terms = [coefficients.get(4*h+b, 0) ^ (coefficients.get(64+4*h+b, 0) if upper else 0)
                     for b in range(4)]
            d[h] = (terms[0], terms[0] ^ terms[1], terms[0] ^ terms[2],
                    terms[0] ^ terms[1] ^ terms[2] ^ terms[3])
        for packet in range(16):
            for lane in range(4):
                got = 0
                for h in HIGH:
                    if packet & h == h:
                        got ^= d[h][lane]
                assert got == columns[64*upper+4*packet+lane]
    # Authenticate the split-half moment combination on all 128 raw inputs.
    monomial_rows = []
    for m in MONOMIALS:
        h, low = (m & 63) >> 2, m & 3
        halves = (1,) if m & 64 else (0, 1)
        got = sum(1 << (64*upper+4*p+lane) for upper in halves for p in range(16)
                  for lane in range(4) if p & h == h and lane & low == low)
        want = sum(1 << p for p in range(128) if p & m == m)
        assert got == want
        monomial_rows.append(got)
    assert all(cse.compose(row, monomial_rows) == literal for row, literal in zip(anfs, rows))
    code = [f'// Generated by generate_t128.py; prefix16 of pinned t128_s20_nested.json.',
            f'// Input SHA256: {sha256(raw).hexdigest()}', '#pragma once',
            '#include "../../../src/packet/PacketLargeInner.h"',
            'namespace spin::research::k16codesign::t128 {',
            'using namespace detail::packet::large;',
            'inline constexpr std::uint16_t columns[128]={'+','.join(map(hex, columns))+'};',
            'template<class Emit> static SPIN_FORCEINLINE void step(const __m128i* state,',
            'const Block* raw,std::size_t packetBase,Emit& emit,__m128i* moments) {',
            'alignas(64) __m512i upperHigh[16];']
    generated, cost = cse.emit_cse([coefficients[m] for m in MONOMIALS],
                                  [f'state[{j}]' for j in range(16)], 'cs', [f'c{m}' for m in MONOMIALS])
    code.append(f'// ANF coefficient circuit: {cost} XMM XORs.')
    code += generated + half(1) + half(0) + ['}']
    code += ['static SPIN_FORCEINLINE void finish128(const __m128i* z,__m128i* out) {']
    generated, finish_cost = cse.emit_cse(anfs, [f'z[{m}]' for m in MONOMIALS],
                                         'fs', [f'f{j}' for j in range(16)])
    code += [f'// Literal A^T feedback: {finish_cost} XMM XORs.'] + generated
    code += [f'out[{j}]=f{j};' for j in range(16)] + ['}', '}']
    (HERE / 'T128Map.h').write_text('\n'.join(code)+'\n', encoding='utf-8', newline='\n')
    print(f'Validated 128 expansion columns, 29 feedback moments, CA=0; coefficient XORs={cost}, feedback XORs={finish_cost}.')


if __name__ == '__main__':
    main()
