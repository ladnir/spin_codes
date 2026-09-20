"""Exact S19 expansion circuit; only change the forward evaluation schedule."""
from pathlib import Path
from functools import reduce
from operator import xor
import json
import re
import sys

root=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(root/'scripts'))
import probe_bch_forward_xor_circuit as paar

p=Path(sys.argv[1])
source=(p/'generated/SelectedMaps.h').read_text().split('struct Map128S19 {',1)[1]
match=re.search(r'columns\{([^}]+)\}',source)
assert match
targets=[int(v,0) for v in match[1].split(',')]
assert len(targets)==128 and all(0<v<1<<19 for v in targets)
paar.DIMENSION=19
circuit=paar.synthesize(0,2,targets)
forms=[1<<i for i in range(19)]
for a,b in circuit.gates: forms.append(forms[a]^forms[b])
assert [reduce(xor,(forms[i] for i in ids),0) for ids in circuit.output_signals]==targets
code=['#pragma once','#include "Spin.h"','namespace bare_spin {',
      'template<class Gather> OC_FORCEINLINE void forwardExpansionS19(Gather& gather,block* out,__m128i* values,const __m128i* state,std::size_t base) {']
emitted=set()
def visit(i):
    if i in emitted: return
    if i<19: code.append(f'const auto e{i}=state[{i}];')
    else:
        a,b=circuit.gates[i-19];visit(a);visit(b)
        code.append(f'const auto e{i}=_mm_xor_si128(e{a},e{b});')
    emitted.add(i)
for j,ids in enumerate(circuit.output_signals):
    for i in ids: visit(i)
    terms=[f'e{i}' for i in ids]
    while len(terms)>1:
        terms=[f'_mm_xor_si128({terms[i]},{terms[i+1]})' if i+1<len(terms) else terms[i] for i in range(0,len(terms),2)]
    assert terms
    code.extend([f'const auto raw{j}=gather(base+{j}).mData;',f'values[{j}]=raw{j};',
                 f'out[base+{j}]=block(_mm_xor_si128(raw{j},{terms[0]}));'])
code+=['}','}']
(p/'ForwardExpansion.h').write_text('\n'.join(code)+'\n',newline='\n')
path=p/'Inner.h';text=path.read_text()
text='#include "ForwardExpansion.h"\n'+text
start=text.index('template<class Map,class Gather> OC_FORCEINLINE void innerForward(')
end=text.index('template<u32 Mask> OC_FORCEINLINE __m128i imtSparseSum',start)
body=text[start:end]
anchor='''            if constexpr(SPIN_GROUPED_A && isImtMap<Map>) {'''
assert body.count(anchor)==1
body=body.replace(anchor,'''            if constexpr(std::is_same_v<Map,Map128S19>) {
                forwardExpansionS19(gather,out,values,state,base);
            } else {
            if constexpr(SPIN_GROUPED_A && isImtMap<Map>) {''')
anchor='            forwardPoints<Map>(gather,out,values,table,base,std::make_index_sequence<Map::T>{});'
assert body.count(anchor)==1
body=body.replace(anchor,anchor+'\n            }')
path.write_text(text[:start]+body+text[end:],newline='\n')
record={'map':'Map128S19.columns','expansion_xors':circuit.xor_count,'exact_outputs_verified':128,
        'raw_input_xors':128,'seed':0,'minimum_overlap':2}
(p/'FORWARD_EXPANSION.json').write_text(json.dumps(record,indent=2)+'\n',newline='\n')
print(record)
