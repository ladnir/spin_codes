"""Compile the exactly audited balanced-image map into isolated mixer variants.

The state coordinates are unchanged. Nibble grouping only permutes lookup-table
inputs. The feedback circuit is checked symbolically against every map column;
the dense C++ oracle uses raw columns and raw transvection samples instead.
"""
import hashlib
import json
from pathlib import Path
import sys

import generate_mixer as mixer

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
BARE=ROOT/'workstreams/bare_bch_rm2sub'
sys.path.insert(0,str(BARE))
import generate as shared


def header(columns):
    t,s=128,19
    assert len(columns)==t and all(0<c<1<<s for c in columns)
    coefficients=columns.copy()
    for bit in range(7):
        for x in range(t):
            if x&(1<<bit):coefficients[x]^=coefficients[x^(1<<bit)]
    assert all(not c or x.bit_count()<=2 for x,c in enumerate(coefficients))
    monomials=[x for x in range(t) if x.bit_count()<=2]
    targets=[sum(((coefficients[x]>>j)&1)<<i for i,x in enumerate(monomials)) for j in range(s)]
    # Audit the actual pruned zeta schedule, not just the unpruned transform.
    z=[1<<x for x in range(t)]
    distance=t//2
    while distance:
        for base in range(0,t,2*distance):
            if (base//(2*distance)).bit_count()<=2:
                for j in range(distance):z[base+j]^=z[base+j+distance]
        distance//=2
    shared.paar.DIMENSION=len(monomials)
    circuit=shared.paar.synthesize(0,2,targets)
    symbolic=[z[x] for x in monomials]
    for a,b in circuit.gates:symbolic.append(symbolic[a]^symbolic[b])
    for j,signals in enumerate(circuit.output_signals):
        rebuilt=0
        for i in signals:rebuilt^=symbolic[i]
        assert rebuilt==sum(((c>>j)&1)<<x for x,c in enumerate(columns))
    order,before,after=shared.grouping(columns,s)
    grouped=[sum(((c>>b)&1)<<j for j,b in enumerate(order)) for c in columns]
    assert sorted(order)==list(range(s))
    assert columns==[sum(((c>>j)&1)<<b for j,b in enumerate(order)) for c in grouped]
    code=['#pragma once','#include "Inner.h"','namespace bare_spin {',
          'struct BalancedMap : Map128S19 {',
          'static constexpr std::array<std::uint32_t,T> columns{'+','.join(map(hex,columns))+'};',
          'static constexpr std::array<std::uint32_t,T> groupedColumns{'+','.join(map(hex,grouped))+'};',
          'static constexpr std::array<unsigned,S> groupOrder{'+','.join(map(str,order))+'};',
          'static inline void finish(const __m128i* z,__m128i* out) {']
    code += [f'const auto v{i}=z[{x}];' for i,x in enumerate(monomials)]
    code += [f'const auto v{i+len(monomials)}=vx(v{a},v{b});' for i,(a,b) in enumerate(circuit.gates)]
    code += [f'out[{i}]={shared.expression(signals)};' for i,signals in enumerate(circuit.output_signals)]
    code += ['}','};','}']
    return '\n'.join(code)+'\n',dict(feedback_xors=circuit.xor_count,
          emission_lookups_before=before,emission_lookups=after,group_order=order,
          exact_feedback_circuit_checked=True,exact_lookup_grouping_checked=True)


def main():
    map_path=HERE/'NO_CONSTANT_MAP.json'
    record=json.loads(map_path.read_text());assert (record['t'],record['s'])==(128,19)
    code,stats=header(record['columns'])
    old_manifest=HERE/'MIXER_IMPLEMENTATION.json'
    old=json.loads(old_manifest.read_text())
    for name,digest in old['source_sha256'].items():assert shared.digest(ROOT/name)==digest,name
    records=[]
    sources=[Path(__file__),map_path,old_manifest,Path(shared.__file__),Path(shared.paar.__file__)]
    for style in ('sparse','masked'):
        name='balanced_'+style;directory=HERE/'generated'/name;directory.mkdir(parents=True,exist_ok=True)
        source_path=HERE/'generated'/('mixer_'+style)/'CandidateSpin.cpp'
        expected=next(v['source_sha256'] for v in old['candidates'] if v['name']=='mixer_'+style)
        assert shared.digest(source_path)==expected
        sources.append(source_path)
        source=mixer.once(source_path.read_text(),'#include "Inner.h"','#include "Inner.h"\n#include "BalancedMap.h"')
        source=mixer.once(source,'"t128_s19_transvection_r1"','"t128_s19_balanced_transvection_r1"')
        source=mixer.once(source,'case Configuration::T128S19: setupInner<Map128S19>(coefficientSeed); break;',
            'case Configuration::T128S19:\n'
            '            if(mOuter==Outer::Bch128x32) setupInner<BalancedMap>(coefficientSeed);\n'
            '            else setupInner<Map128S19>(coefficientSeed);\n'
            '            break;')
        for layout in ('true','false'):
            source=mixer.once(source,f'run<Map128S19,{layout},true>',f'run<BalancedMap,{layout},true>')
        source=mixer.once(source,'oracle<Map128S19,true>','oracle<BalancedMap,true>')
        (directory/'BalancedMap.h').write_text(code,encoding='utf-8',newline='\n')
        (directory/'CandidateSpin.cpp').write_text(source,encoding='utf-8',newline='\n')
        records.append(dict(name=name,**stats,source_sha256={p.name:shared.digest(p) for p in
                       (directory/'BalancedMap.h',directory/'CandidateSpin.cpp')}))
    payload=dict(status='SYMBOLICALLY_CHECKED_IMPLEMENTATION_PENDING_RUNTIME_TESTS',
                 map='NO_CONSTANT_MAP.json',rounds=1,candidates=records,
                 source_sha256={p.relative_to(ROOT).as_posix():shared.digest(p) for p in sources})
    (HERE/'BALANCED_IMPLEMENTATION.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(records,indent=2))


if __name__=='__main__':main()
