"""Reschedule the existing BCH forward DAG; verify all 256 linear forms exactly."""
from pathlib import Path
import json
import re
import sys
from functools import reduce
from operator import xor

p=Path(sys.argv[1]);mode=sys.argv[2]
assert mode in ('dfs','local16','local32')
path=p/'generated/BchForward512.cpp'
original=path.read_text()
prefix,body=original.split('void bchForward4(const block* a,block* x) {',1)
declarations=re.findall(r'const auto ([ogz]\d+)=(.*);',body)
definitions=dict(declarations)
assert len(definitions)==len(declarations)
inputs={f'o{i}' for i in range(128)}
assert {name for name in definitions if name.startswith('o')}==inputs
assert {name for name in definitions if name.startswith('z')}=={f'z{j}' for j in range(256)}
for i in range(128):
    assert definitions[f'o{i}']==f'join(a[{i}],a[128+{i}],a[256+{i}],a[384+{i}])'
for name,expr in definitions.items():
    if name not in inputs:
        assert re.fullmatch(r'(?:vx|[ogz]\d+|[(),])+',expr), (name,expr)
dependencies={name:re.findall(r'\b[ogz]\d+\b',expr) if name not in inputs else []
              for name,expr in definitions.items()}
forms={f'o{i}':1<<i for i in range(128)}
def form(name):
    if name not in forms:
        forms[name]=reduce(xor,(form(dep) for dep in dependencies[name]),0)
    return forms[name]

# The independent, committed generator matrix is row-major 128 x 256 bits.
words=[int(x,16) for x in re.findall(r'0x([0-9a-fA-F]+)ULL',(p/'generated/BchCircuit.h').read_text())]
assert len(words)==128*4
targets=[sum(((words[4*i+j//64]>>(j%64))&1)<<i for i in range(128)) for j in range(256)]
assert [form(f'z{j}') for j in range(256)]==targets

chunk=256 if mode=='dfs' else int(mode[5:])
code=[prefix];count=0
for start in range(0,256,chunk):
    fn='void bchForward4' if mode=='dfs' else f'__attribute__((noinline)) static void forwardPart{start}'
    code.append(f'{fn}(const block* a,block* x) {{')
    emitted=set()
    def visit(name):
        global count
        if name in emitted: return
        for dep in dependencies[name]: visit(dep)
        expr=definitions[name]
        code.append(f'const auto {name}={expr};')
        count+=expr.count('vx(')
        emitted.add(name)
    for j in range(start,start+chunk):
        visit(f'z{j}')
        code.append(f'_mm512_storeu_si512(x+4*{j},z{j});')
    code.append('}')
if mode!='dfs':
    code.append('void bchForward4(const block* a,block* x) {')
    code.extend(f'forwardPart{j}(a,x);' for j in range(0,256,chunk))
    code.append('}')
code.append('}')
path.write_text('\n'.join(code)+'\n',newline='\n')
record={'schedule':mode,'xors':count,'original_xors':sum(e.count('vx(') for e in definitions.values()),
        'outputs_verified_against_generator_matrix':256}
(p/'FORWARD_SCHEDULE.json').write_text(json.dumps(record,indent=2)+'\n',newline='\n')
print(record)
