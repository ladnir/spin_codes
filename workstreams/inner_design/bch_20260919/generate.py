"""Exact-map BCH scheduling experiments; never edits the pinned implementation."""
from pathlib import Path
import re
import sys
import json
import hashlib
from functools import reduce
from operator import xor

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
import probe_bch_forward_xor_circuit as paar


def expression(signals):
    level = [f'v{i}' for i in signals]
    assert level
    while len(level) > 1:
        level = [f'vx({level[i]},{level[i+1]})' if i+1 < len(level) else level[i]
                 for i in range(0, len(level), 2)]
    return level[0]


def verify(gates, outputs, targets):
    forms = [1 << i for i in range(256)]
    for a, b in gates:
        forms.append(forms[a] ^ forms[b])
    assert [reduce(xor, (forms[i] for i in out), 0) for out in outputs] == targets


def emit(name, gates, outputs, targets, chunk, synth=False, width=256):
    code = ['#include "Spin.h"', '#include "generated/BchCircuit.h"', 'namespace bare_spin {']
    vec = {128:'__m128i',256:'__m256i',512:'__m512i'}[width]
    op = {128:'_mm_xor_si128',256:'_mm256_xor_si256',512:'_mm512_xor_si512'}[width]
    code += [f'static inline {vec} vx({vec} a,{vec} b) {{return {op}(a,b);}}']
    count = 0
    for start in range(0, 128, chunk):
        localgates, localouts = gates, outputs[start:start+chunk]
        if synth:
            paar.DIMENSION = 256
            circuit = paar.synthesize(0, 2, targets[start:start+chunk])
            localgates, localouts = circuit.gates, circuit.output_signals
        verify(localgates, localouts, targets[start:start+chunk])
        namefn = f'part{start}'
        args = 'const block* a,const block* b,block* x,block* y' if width == 256 else 'const block* a,block* x'
        code += [f'__attribute__((noinline)) static void {namefn}({args}) {{']
        computed = set()

        def visit(i):
            nonlocal count
            if i in computed:
                return
            if i < 256:
                load = f'_mm256_set_m128i(b[{i}].mData,a[{i}].mData)' if width == 256 else f'a[{i}].mData'
                if width == 512:
                    load = f'_mm512_castsi128_si512(a[{i}].mData)'
                    for lane in (1,2,3):
                        load = f'_mm512_inserti32x4({load},a[{256*lane+i}].mData,{lane})'
                code.append(f'const auto v{i}={load};')
            else:
                a, b = localgates[i-256]
                visit(a); visit(b)
                code.append(f'const auto v{i}=vx(v{a},v{b});')
                count += 1
            computed.add(i)

        for j, signals in enumerate(localouts, start):
            for i in signals:
                visit(i)
            count += len(signals)-1
            code.append(f'const auto o{j}={expression(signals)};')
            if width == 256:
                code += [f'x[{j}]=block(_mm256_castsi256_si128(o{j}));',
                         f'y[{j}]=block(_mm256_extracti128_si256(o{j},1));']
            elif width == 512:
                for lane in range(4):
                    code += [f'x[{128*lane+j}]=block(_mm512_extracti32x4_epi32(o{j},{lane}));']
            else:
                code += [f'x[{j}]=block(o{j});']
        code += ['}']
    code += ['void bchTranspose4(const block* a,block* x) {' if width == 512 else
             'void bchTranspose2(const block* a,const block* b,block* x,block* y) {']
    for start in range(0, 128, chunk):
        if width == 256:
            code += [f'part{start}(a,b,x,y);']
        elif width == 512:
            code += [f'part{start}(a,x);']
        else:
            code += [f'part{start}(a,x);part{start}(b,y);']
    code += ['}', '}']
    path = HERE / 'generated' / f'{name}.cpp'
    path.write_text('\n'.join(code)+'\n', newline='\n')
    return dict(xors=count, width=width, chunk=chunk, resynthesized=synth,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    src = ROOT / 'workstreams/bare_bch_rm2sub/generated/BchCircuit.cpp'
    header = ROOT / 'workstreams/bare_bch_rm2sub/generated/BchCircuit.h'
    text = src.read_text()
    gates = [(int(a), int(b)) for _, a, b in re.findall(r'const auto v(\d+)=vx\(v(\d+),v(\d+)\);', text)]
    outputs = [[int(i) for i in re.findall(r'v(\d+)', ex)]
               for _, ex in re.findall(r'const auto o(\d+)=(.*);', text)]
    words = [int(i, 16) for i in re.findall(r'0x([0-9a-f]+)ULL', header.read_text())]
    targets = [sum(words[4*i+j] << (64*j) for j in range(4)) for i in range(128)]
    verify(gates, outputs, targets)
    (HERE / 'generated').mkdir(exist_ok=True)
    records = {'original': dict(xors=len(gates)+sum(len(s)-1 for s in outputs),
                               sha256=hashlib.sha256(src.read_bytes()).hexdigest())}
    for name, chunk, synth, width in [('dfs',128,False,256),('chunk32',32,False,256),
            ('synth32',32,True,256),('synth16',16,True,256),('single',128,False,128),('four',128,False,512)]:
        records[name] = emit(name, gates, outputs, targets, chunk, synth, width)
        print(name, records[name], flush=True)
    for overlap in (3,4,8):
        paar.DIMENSION = 256
        circuit = paar.synthesize(0, overlap, targets)
        name = f'share{overlap}'
        records[name] = emit(name, circuit.gates, circuit.output_signals, targets, 128)
        print(name, records[name], flush=True)
        if overlap in (3,4):
            qname = f'packedshare{overlap}'
            record = emit(qname, circuit.gates, circuit.output_signals, targets, 128, width=512)
            path = HERE / 'generated' / f'{qname}.cpp'
            code, count = re.subn(r'const auto v(\d+)=_mm512_inserti32x4\(.*;',
                lambda m: f'const auto v{m[1]}=_mm512_loadu_si512(a+4*{m[1]});', path.read_text())
            assert count == 256
            path.write_text(code, newline='\n')
            records[qname] = dict(record, sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                                 internal_layout='four interleaved rows', minimum_overlap=overlap)
    spin = (ROOT / 'workstreams/inner_design/asymmetric/bch256/weight5/implementation/Weight5Spin.cpp').read_text()
    spin = spin.replace('namespace bare_spin {','namespace bare_spin {\nvoid bchTranspose4(const block*,block*);',1)
    spin = spin.replace('    mK=std::size_t{1}<<exponent;',
        '    if(tileRows && tileRows<4) throw std::invalid_argument("four-row BCH requires tileRows>=4");\n'
        '    mK=std::size_t{1}<<exponent;')
    old = 'for(std::size_t j=0;j<tileSize;j+=512)\n                bchTranspose2(tile+j,tile+j+256,out+base/2+j/2,out+base/2+j/2+128);'
    assert spin.count(old) == 1
    spin = spin.replace(old, 'for(std::size_t j=0;j<tileSize;j+=1024)\n                bchTranspose4(tile+j,out+base/2+j/2);')
    (HERE / 'generated/SpinFour.cpp').write_text(spin, newline='\n')
    # Exact internal workspace relabeling: four rows become interleaved lanes.
    # The external route, recurrence, and output row order remain unchanged.
    quad = (HERE / 'generated/four.cpp').read_text()
    quad, count = re.subn(r'const auto v(\d+)=_mm512_inserti32x4\(.*;',
        lambda m: f'const auto v{m[1]}=_mm512_loadu_si512(a+4*{m[1]});', quad)
    assert count == 256
    (HERE / 'generated/packedfour.cpp').write_text(quad, newline='\n')
    mapping = '((local&~1023U)|((local&255U)<<2)|((local>>8)&3U))'
    old = 'mSlots32[inner]=slot; mOffsets32[slot]=local;\n        pack(mSlots24.data()+3*inner,slot); pack(mOffsets24.data()+3*slot,local);'
    assert spin.count(old) == 1
    spin = spin.replace(old, f'mSlots32[inner]=slot; const u32 packedLocal={mapping}; mOffsets32[slot]=packedLocal;\n'
                            '        pack(mSlots24.data()+3*inner,slot); pack(mOffsets24.data()+3*slot,packedLocal);')
    old = 'mOffsets32[slot]!=outer'
    assert spin.count(old) == 1
    spin = spin.replace(old, '((mOffsets32[slot]&~1023U)|((mOffsets32[slot]&3U)<<8)|((mOffsets32[slot]>>2)&255U))!=outer')
    (HERE / 'generated/SpinPackedFour.cpp').write_text(spin, newline='\n')
    for local in range(1 << 20):
        packed = (local&~1023)|((local&255)<<2)|((local>>8)&3)
        assert (packed&~1023)|((packed&3)<<8)|((packed>>2)&255) == local
    records['packedfour'] = dict(records['four'], internal_layout='four interleaved rows',
        sha256=hashlib.sha256(quad.encode()).hexdigest())
    (HERE / 'generated/manifest.json').write_text(json.dumps(records, indent=2)+'\n', newline='\n')


if __name__ == '__main__': main()
