"""Four-row BCH output-local synthesis; prove every emitted output exactly."""
from pathlib import Path
import importlib.util
import re
import sys

root=Path(__file__).resolve().parents[2]
out=Path(sys.argv[1]);mode=sys.argv[2]
assert mode in ('local8','local16','local32','single','reload')
if mode=='reload':
    # Reload immutable input leaves on demand instead of keeping 256 vectors live.
    # The GNU unaligned/may-alias SIMD type also serves the loadu intrinsic.
    # Volatile prevents cross-output load commoning; it does not add CPU fences.
    p=out/'generated/BchAvx512.cpp';text=p.read_text()
    pattern=r'^const auto v(\d+)=_mm512_loadu_si512\(a\+4\*\d+\);\n'
    ids={int(x) for x in re.findall(pattern,text,re.M)}
    assert ids==set(range(256)),ids
    text=re.sub(pattern,'',text,flags=re.M)
    text=re.sub(r'\bv(\d+)\b',lambda m:f'loadInput(a+4*{m[1]})' if int(m[1])<256 else m[0],text)
    anchor='namespace bare_spin {\n';assert text.count(anchor)==1
    text=text.replace(anchor,anchor+'static inline __m512i loadInput(const block* p) {return *reinterpret_cast<const volatile __m512i_u*>(p);}\n')
    p.write_text(text,newline='\n')
    print('Four-row reload schedule: unchanged XOR DAG; 256 immutable input leaves reloaded on demand')
    sys.exit(0)
spec=importlib.util.spec_from_file_location('bch_emit',root/'workstreams/inner_design/bch_20260919/generate.py')
emit=importlib.util.module_from_spec(spec);spec.loader.exec_module(emit)
emit.HERE=out
words=[int(x,16) for x in re.findall(r'0x([0-9a-f]+)ULL',(out/'generated/BchCircuit.h').read_text())]
targets=[sum(words[4*i+j]<<(64*j) for j in range(4)) for i in range(128)]
emit.paar.DIMENSION=256
circuit=emit.paar.synthesize(0,4,targets)
chunk=1 if mode=='single' else int(mode[5:])
stats=emit.emit('BchAvx512',circuit.gates,circuit.output_signals,targets,chunk,
                synth=mode!='single',width=512)
p=out/'generated/BchAvx512.cpp'
text,count=re.subn(r'const auto v(\d+)=_mm512_inserti32x4\(.*;',
    lambda m:f'const auto v{m[1]}=_mm512_loadu_si512(a+4*{m[1]});',p.read_text())
assert count>=256,count
p.write_text(text,newline='\n')
print('Four-row local BCH',mode,'exact output matrix verified;',stats['xors'],'XORs;',count,'input loads')
