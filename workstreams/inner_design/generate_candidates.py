"""Generate isolated exact-conjugacy variants without changing certified sources."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import search_basis as search

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
BARE=ROOT/'workstreams/bare_bch_rm2sub'
sys.path.insert(0,str(BARE))
import generate as shared


def replace_once(text,old,new):
    assert text.count(old)==1,old
    return text.replace(old,new)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('names',nargs='*')
    parser.add_argument('--search',default='BASIS_SEARCH.json')
    args=parser.parse_args()
    data=json.loads((HERE/args.search).read_text())
    names=args.names or ['grouped','feedback_rref','emission_rref','joint_w1_r5','joint_w2_r5','joint_w0_r5']
    records={r['name']:r for r in data['candidates']}
    built=[]
    for name in names:
        r=records[name];directory=HERE/'generated'/name;directory.mkdir(parents=True,exist_ok=True)
        a=[int(x,16) for x in r['emission_rows_hex']]
        columns=[sum(((row>>p)&1)<<j for j,row in enumerate(a)) for p in range(128)]
        targets=[int(x,16) for x in r['feedback_monomials_hex']]
        shared.paar.DIMENSION=len(data['monomials'])
        circuit=shared.paar.synthesize(0,2,targets)
        symbolic=[1<<j for j in range(len(data['monomials']))]
        for u,v in circuit.gates:symbolic.append(symbolic[u]^symbolic[v])
        for target,signals in zip(targets,circuit.output_signals):
            rebuilt=0
            for s in signals:rebuilt^=symbolic[s]
            assert rebuilt==target
        code=['#pragma once','#include "Inner.h"','namespace bare_spin {',
              'struct CandidateMap : Map128S19 {',
              'static constexpr std::array<std::uint32_t,T> columns{'+','.join(map(hex,columns))+'};',
              'static constexpr auto groupedColumns=columns;',
              'static constexpr std::array<unsigned,S> groupOrder{'+','.join(map(str,range(19)))+'};',
              'static inline void finish(const __m128i* z,__m128i* out) {']
        code += [f'const auto v{i}=z[{m}];' for i,m in enumerate(data['monomials'])]
        code += [f'const auto v{i+len(data["monomials"])}=vx(v{a},v{b});' for i,(a,b) in enumerate(circuit.gates)]
        code += [f'out[{i}]={shared.expression(signals)};' for i,signals in enumerate(circuit.output_signals)]
        code += ['}', 'static void conjugate(u32* rows) {',
                 'constexpr u32 V[19]={'+','.join(map(hex,r['V']))+'};',
                 'constexpr u32 inverse[19]={'+','.join(map(hex,r['V_inverse']))+'};',
                 'u32 left[19]{},result[19]{};',
                 'for(unsigned i=0;i<19;++i) for(unsigned j=0;j<19;++j) if((V[i]>>j)&1) left[i]^=rows[j];',
                 'for(unsigned i=0;i<19;++i) for(unsigned j=0;j<19;++j) if((left[i]>>j)&1) result[i]^=inverse[j];',
                 'std::memcpy(rows,result,sizeof(result));','}','};','}']
        (directory/'CandidateMap.h').write_text('\n'.join(code)+'\n',encoding='utf-8',newline='\n')
        source=(BARE/'Spin.cpp').read_text()
        source=replace_once(source,'#include "Inner.h"','#include "Inner.h"\n#include "CandidateMap.h"')
        source=replace_once(source,'#include "../rate_quarter_bch/implementation/generated/QuarterCircuit.h"','#include "QuarterCircuit.h"')
        for layout in ('true','false'):
            source=replace_once(source,f'run<Map128S19,{layout},true>',f'run<CandidateMap,{layout},true>')
        old='for(unsigned j=0;j<Map::S;++j) mFieldRows[e*Map::S+j]=fieldMultiply<Map>(1U<<j,a);'
        source=replace_once(source,old,old+'\n        if constexpr(Map::T==128 && Map::S==19) {\n            if(mOuter==Outer::Bch128x32) CandidateMap::conjugate(mFieldRows.data()+e*Map::S);\n        }')
        (directory/'CandidateSpin.cpp').write_text(source,encoding='utf-8',newline='\n')
        built.append(dict(name=name,emission_lookups=r['emission_lookups'],feedback_xors=circuit.xor_count,
                          source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in directory.glob('*.h')},
                          cpp_sha256=hashlib.sha256((directory/'CandidateSpin.cpp').read_bytes()).hexdigest()))
    (HERE/(Path(args.search).stem+'_CANDIDATES.json')).write_text(json.dumps(dict(baseline_spin_sha256=shared.digest(BARE/'Spin.cpp'),
        search_sha256=shared.digest(HERE/args.search),candidates=built),indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(built,indent=2))


if __name__=='__main__':main()
