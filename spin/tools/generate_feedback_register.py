from pathlib import Path
from collections import Counter
from itertools import combinations
import re,json
p=Path(__file__).resolve().parents[1]/'src/kernels/wide';root=p
source=(p/'generated/WideCircuit.h').read_text()
base=(Path(__file__).with_name('generate_bch_register.py')).read_text()
suffix=base[base.index('def canon(v):'):].replace("'generated/BchRegisterSchedule.h'","'feedbackRegister.h'")
def emit(name,defs,outs,basis):
 global root
 env=dict(root=root,defs=defs,outs=outs,basis=basis)
 exec('from collections import Counter,defaultdict,deque\nimport json\n'+suffix,env)
 h=(root/'feedbackRegister.h').read_text().replace('bchForwardRegister',name).replace('BCH linear map','feedback linear map')
 (root/(name+'.h')).write_text(h)
 return env['results'][0]
defs={};outs=[];basis={}
body=source.split('static SPIN_FORCEINLINE void feedback(const Vec* z,Vec* out) {',1)[1].split('\n}',1)[0]
for line in body.splitlines():
 m=re.fullmatch(r'const auto (f\d+)=(.*);',line)
 if m:
  name,expr=m.groups()
  if expr.startswith('z['):defs[name]=int(re.search(r'\d+',expr)[0]);basis[name]=1<<defs[name]
  else:
   defs[name]=re.findall(r'f\d+',expr);basis[name]=0
   for t in defs[name]:basis[name]^=basis[t]
  continue
 m=re.fullmatch(r'out\[(\d+)\]=(.*);',line)
 if m:
  idx=int(m[1]);name='y'+str(idx);defs[name]=re.findall(r'f\d+',m[2]);basis[name]=0
  for t in defs[name]:basis[name]^=basis[t]
  outs.append((idx,name));continue
 assert not line.strip(),line
columns=[int(v,16) for v in re.search(r'feedbackColumns\{([^}]+)',source)[1].split(',')]
for i,name in outs:assert basis[name]==sum(1<<j for j,c in enumerate(columns) if c>>i&1)
stats={'feedback':emit('feedbackRegister',defs,outs,basis)}

print(stats)
