"""Build-only BCH resynthesis with less sharing; emits C++ to stdout.

Reuse the existing Paar synthesizer, verify every formal output, and retain
four-row SIMD packing. No production sources or generated files are edited.
"""
from functools import reduce
from operator import xor
from pathlib import Path
import re
import sys


def expression(signals):
    level = [f'v{i}' for i in signals]
    assert level
    while len(level)>1:
        level = [f'vx({level[i]},{level[i+1]})' if i+1<len(level) else level[i]
                 for i in range(0,len(level),2)]
    return level[0]


def main():
    header,helper,overlap = Path(sys.argv[1]),Path(sys.argv[2]),int(sys.argv[3])
    sys.path.insert(0,str(helper.parent))
    import probe_bch_forward_xor_circuit as paar
    paar.DIMENSION = 256
    words = [int(v,16) for v in re.findall(r'0x([0-9a-f]+)ULL',header.read_text())]
    assert len(words)==512
    targets = [sum(words[4*i+j] << (64*j) for j in range(4)) for i in range(128)]
    circuit = paar.synthesize(0,overlap,targets)
    forms = [1 << i for i in range(256)]
    for a,b in circuit.gates:
        forms.append(forms[a]^forms[b])
    assert [reduce(xor,(forms[v] for v in row),0) for row in circuit.output_signals] == targets
    print(f'BCH Share{overlap}: verified 128 forms; XORs={circuit.xor_count}; shared gates={len(circuit.gates)}',file=sys.stderr)
    print('#include "Spin.h"\n#include "generated/BchCircuit.h"\nnamespace spin::detail::kernel {')
    print('static inline __m512i vx(__m512i a,__m512i b) {return _mm512_xor_si512(a,b);}')
    print(f'SPIN_NOINLINE void bchTranspose4Share{overlap}(const block* __restrict a,block* __restrict x) {{')
    computed = set()
    def visit(v):
        if v in computed:
            return
        if v < 256:
            print(f'const auto v{v}=_mm512_loadu_si512(a+4*{v});')
        else:
            a,b = circuit.gates[v-256]
            visit(a)
            visit(b)
            print(f'const auto v{v}=vx(v{a},v{b});')
        computed.add(v)
    for i,row in enumerate(circuit.output_signals):
        for v in row:
            visit(v)
        print(f'const auto o{i}={expression(row)};')
        for lane in range(4):
            print(f'x[{128*lane+i}]=block(_mm512_extracti32x4_epi32(o{i},{lane}));')
    print('}\n}')


if __name__ == '__main__':
    main()
