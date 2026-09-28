"""Emit an exact BCH circuit with a different output schedule (build tool).

All XOR gates and output expressions are retained. Check every output as a
formal 256-bit linear form against BchRows before emitting C++ to stdout.
"""
from collections import Counter
from functools import reduce
from operator import xor
from pathlib import Path
import re
import sys


def main():
    source,header,mode = Path(sys.argv[1]),Path(sys.argv[2]),sys.argv[3]
    text = source.read_text()
    definitions = {int(i):line for line,i in re.findall(r'(const auto v(\d+)=.*;)',text)}
    gates = {int(i):(int(a),int(b)) for i,a,b in re.findall(r'const auto v(\d+)=vx\(v(\d+),v(\d+)\);',text)}
    expressions = {int(i):expression for i,expression in re.findall(r'const auto o(\d+)=(.*);',text)}
    outputs = {i:list(map(int,re.findall(r'v(\d+)',expression))) for i,expression in expressions.items()}
    assert set(outputs) == set(range(128))
    forms = {i:1 << i for i in range(256)}
    for i,(a,b) in sorted(gates.items()):
        forms[i] = forms[a]^forms[b]
    words = [int(v,16) for v in re.findall(r'0x([0-9a-f]+)ULL',header.read_text())]
    assert len(words) == 512
    targets = [sum(words[4*i+j] << (64*j) for j in range(4)) for i in range(128)]
    assert [reduce(xor,(forms[v] for v in outputs[i]),0) for i in range(128)] == targets
    use = Counter(v for pair in gates.values() for v in pair)
    use.update(v for output in outputs.values() for v in output)
    computed,live = set(),set()

    def simulate(output,computed,live,use,emit=None):
        peak = len(live)
        def consume(v):
            use[v] -= 1
            if use[v] == 0:
                live.remove(v)
        def visit(v):
            nonlocal peak
            if v in computed:
                return
            if v in gates:
                for dependency in gates[v]:
                    visit(dependency)
                for dependency in gates[v]:
                    consume(dependency)
            computed.add(v)
            live.add(v)
            peak = max(peak,len(live))
            if emit is not None:
                emit.append(definitions[v])
        for v in outputs[output]:
            visit(v)
        for v in outputs[output]:
            consume(v)
        if emit is not None:
            emit.append(f'const auto o{output}={expressions[output]};')
            for lane in range(4):
                emit.append(f'x[{128*lane+output}]=block(_mm512_extracti32x4_epi32(o{output},{lane}));')
        return peak,len(live)

    order=[]
    remaining=set(range(128))
    body=[]
    peak=0
    while remaining:
        if mode in ('Greedy','Release'):
            def cost(i):
                p,l=simulate(i,computed.copy(),live.copy(),use.copy())
                return (p,l,i) if mode=='Greedy' else (l,p,i)
            selected=min(remaining,key=cost)
        elif mode=='Reverse':
            selected=max(remaining)
        elif mode=='Restrict':
            selected=min(remaining)
        else:
            raise ValueError(mode)
        p,_=simulate(selected,computed,live,use,body)
        peak=max(peak,p)
        order.append(selected)
        remaining.remove(selected)
    assert not live and all(value==0 for value in use.values())
    assert set(computed)==set(definitions)
    print(f'BCH {mode}: verified 128 forms; symbolic peak live signals={peak}',file=sys.stderr)
    print('#include "Spin.h"\n#include "generated/BchCircuit.h"\nnamespace spin::detail::kernel {')
    print('static inline __m512i vx(__m512i a,__m512i b) {return _mm512_xor_si512(a,b);}')
    restrict=' __restrict' if mode=='Restrict' else ''
    print(f'SPIN_NOINLINE void bchTranspose4{mode}(const block*{restrict} a,block*{restrict} x) {{')
    print('\n'.join(body))
    print('}\n}')


if __name__ == '__main__':
    main()
