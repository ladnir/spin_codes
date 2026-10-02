"""Generate exact diagnostic-map circuits; outputs are build artifacts only."""
import argparse
from functools import reduce
from operator import xor
from pathlib import Path
import sys
import json
import model

ROOT = model.base.ROOT
sys.path.insert(0,str(ROOT/'scripts'))
import probe_bch_forward_xor_circuit as paar


def circuit(targets,inputs):
    paar.DIMENSION = len(inputs)
    c = paar.synthesize(0,2,targets)
    forms = [1 << j for j in range(len(inputs))]
    lines = [f'const auto v{j}={term};' for j,term in enumerate(inputs)]
    for j,(a,b) in enumerate(c.gates,len(inputs)):
        forms.append(forms[a]^forms[b])
        lines.append(f'const auto v{j}=vx(v{a},v{b});')
    assert [reduce(xor,(forms[j] for j in out),0) for out in c.output_signals] == targets
    def expression(signals):
        level = [f'v{j}' for j in signals]
        if not level:
            return '_mm_setzero_si128()'
        while len(level)>1:
            level = [f'vx({level[j]},{level[j+1]})' if j+1<len(level) else level[j]
                     for j in range(0,len(level),2)]
        return level[0]
    return lines,[expression(x) for x in c.output_signals],c.xor_count


def run(t,s,output,record=None):
    record = model.Engine(t,s).record if record is None else record
    assert (record['t'],record['s'])==(t,s)
    a,b = record['expansion_columns'],record['feedback_columns']
    # Boolean Möbius inversion expresses each expansion coordinate as a
    # degree-at-most-two polynomial. The existing pruned zeta produces its sums.
    coefficients = a.copy()
    for bit in range(t.bit_length()-1):
        for x in range(t):
            if x >> bit & 1:
                coefficients[x] ^= coefficients[x ^ (1 << bit)]
    assert all(not c or x.bit_count()<=2 for x,c in enumerate(coefficients))
    for x in range(t):
        assert reduce(xor,(c for m,c in enumerate(coefficients) if m & x == m),0) == a[x]
    monomials = [x for x in range(t) if x.bit_count()<=2]
    targets = [sum(((coefficients[x]>>j)&1)<<i for i,x in enumerate(monomials)) for j in range(s)]
    finish,outs,finish_cost = circuit(targets,[f'z[{x}]' for x in monomials])
    emission,emits,emit_cost = circuit(b,[f'state[{j}]' for j in range(s)])
    text = ['#pragma once','#include "Inner.h"','namespace bare_spin {',
            'struct AsymmetricMap : Map128S19 {',f'static constexpr unsigned T={t},S={s};']
    for name,values in [('columns',b),('feedbackColumns',a),('groupedColumns',b),('groupOrder',list(range(s)))]:
        text.append(f'static constexpr std::array<u32,{len(values)}> {name}{{'+','.join(map(str,values))+'};')
    text += ['static inline void finish(const __m128i* z,__m128i* out) {']+finish
    text += [f'out[{j}]={v};' for j,v in enumerate(outs)]+['}']
    text += ['template<class Emit> static OC_FORCEINLINE void emitShared(const block* in,__m128i* raw,',
             'const __m128i* state,std::size_t base,Emit& emit) {']+emission
    for p in reversed(range(t)):
        text += [f'raw[{p}]=in[{p}].mData;',f'emit(base+{p},block(vx(raw[{p}],{emits[p]})));']
    text += ['}','};','}']
    output.mkdir(parents=True,exist_ok=True)
    (output/f'Map{t}S{s}.h').write_text('\n'.join(text)+'\n')
    (output/f'Map{t}S{s}.json').write_text(json.dumps(dict(inner=record,
        finish_xors=finish_cost,emission_xors=emit_cost,
        unshared_emission_xors=sum(c.bit_count()-1 for c in b)),indent=2)+'\n')
    print(t,s,'finish XORs',finish_cost,'shared emission XORs',emit_cost,flush=True)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--t',type=int,required=True);p.add_argument('--s',type=int,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--maps',type=Path)
    a=p.parse_args()
    record=next(r['inner'] for r in json.loads(a.maps.read_text())['rows'] if r['s']==a.s and r['found']) if a.maps else None
    run(a.t,a.s,a.output,record)
