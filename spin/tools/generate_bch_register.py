from pathlib import Path
from collections import Counter,defaultdict,deque
import re,json
root=Path(__file__).resolve().parents[1]/'src/kernels/wide'
s=(root/'generated/WideCircuit.h').read_text()
body=s.split('static void bchForward(const Vec* a,Vec* x) {',1)[1].split('\n}',1)[0]
defs={};outs=[];basis={}
for line in body.splitlines():
 m=re.fullmatch(r'const auto (\w+)=(.*);',line)
 if m:
  name,expr=m.groups()
  if 'Ops::load' in expr:defs[name]=int(re.search(r'a\+(\d+)',expr)[1]);basis[name]=1<<defs[name]
  else:
   terms=re.findall(r'\b[ogz]\d+\b',expr);parity=Counter(terms)
   defs[name]=[t for t,c in parity.items() if c%2];v=0
   for t in defs[name]:v^=basis[t]
   basis[name]=v
  continue
 m=re.fullmatch(r'Ops::store\(x\+(\d+),(\w+)\);',line)
 if m:outs.append((int(m[1]),m[2]));continue
 assert not line.strip(),line
assert len(outs)==256
def canon(v):
 while not isinstance(defs[v],int) and len(defs[v])==1:v=defs[v][0]
 return v
outs=[(i,canon(k)) for i,k in outs];outset={k for _,k in outs}
defs={k:([canon(v) for v in d] if not isinstance(d,int) else d) for k,d in defs.items() if canon(k)==k}
uses=Counter(v for d in defs.values() if isinstance(d,list) for v in d);uses.update(k for _,k in outs)
def expanded(k):
 if isinstance(defs[k],int) or uses[k]!=1 or k in outset:return [k]
 return [x for v in defs[k] for x in expanded(v)]
expr={}
for k,d in defs.items():
 if isinstance(d,int):expr[k]=d
 elif uses[k]!=1 or k in outset:
  c=Counter(x for v in d for x in expanded(v));expr[k]=[x for x,n in c.items() if n%2]

def schedule(strategy):
 emitted=set();order=[]
 def emit(k):
  if k in emitted or isinstance(expr[k],int):return
  for v in expr[k]:emit(v)
  order.append(('node',k));emitted.add(k)
 def missing(k,seen):
  if k in emitted or k in seen or isinstance(expr[k],int):return
  seen.add(k)
  for v in expr[k]:missing(v,seen)
 remaining=list(outs)
 if strategy=='original':
  for k in expr:emit(k)
  order.extend(('store',(i,k)) for i,k in outs)
 else:
  while remaining:
   costs=[]
   for i,k in remaining:
    cone=set();missing(k,cone);costs.append((len(cone),i,k))
   _,i,k=min(costs);emit(k);order.append(('store',(i,k)));remaining.remove((i,k))
 return order

def generate(strategy):
 order=schedule(strategy);future=defaultdict(deque)
 for pos,(op,v) in enumerate(order):
  for a in (expr[v] if op=='node' else [v[1]]):future[a].append(pos)
 reg={};owner={};mem={k:('in',d*32) for k,d in expr.items() if isinstance(d,int)}
 slots=[];slot_count=0;max_slots=0;lines=[];stats=Counter();values={};memory={('in',d*32):basis[k] for k,d in expr.items() if isinstance(d,int)}
 def rr(r):return f'%%ymm{r}'
 def loc(v):
  a,b=mem[v];return f'{b}(%[{a}])'
 def asm(inst):lines.append(inst);stats[inst.split()[0]]+=1
 def free_mem(v):
  if v in mem and mem[v][0]=='tmp':slots.append(mem[v][1]//32)
  mem.pop(v,None)
 def alloc(pos,pins=()):
  nonlocal slot_count,max_slots
  unused=set(range(32))-set(owner)
  if unused:return min(unused)
  choices=[(future[v][0] if future[v] else 10**9,r,v) for r,v in owner.items() if r not in pins]
  _,r,v=max(choices)
  if future[v] and v not in mem:
   slot=slots.pop() if slots else slot_count
   if slot==slot_count:slot_count+=1;max_slots=max(max_slots,slot_count)
   mem[v]=('tmp',slot*32);asm(f'vmovdqa64 {rr(r)}, {loc(v)}');memory[mem[v]]=values[r];stats['spill_store']+=1
  del reg[v];del owner[r]
  return r
 def load(v,pos,pins=()):
  if v in reg:return reg[v]
  r=alloc(pos,pins);asm(f'vmovdqu64 {loc(v)}, {rr(r)}');values[r]=memory[mem[v]]
  reg[v]=r;owner[r]=v;return r
 def consume(v,pos,preserve=None):
  assert future[v].popleft()==pos
  if not future[v]:
   if v in reg:
    r=reg.pop(v)
    if r!=preserve:owner.pop(r)
   free_mem(v)
 for pos,(op,k) in enumerate(order):
  if op=='store':
   i,v=k;r=load(v,pos);asm(f'vmovdqu64 {rr(r)}, {i*32}(%[out])');assert values[r]==basis[v]
   memory[('out',i*32)]=values[r]
   if len(future[v])>1:
    free_mem(v);mem[v]=('out',i*32)
   consume(v,pos);continue
  terms=list(expr[k]);assert terms
  # Reuse a dying register where possible; all operands are unique after XOR cancellation.
  first=min(terms,key=lambda v:(not(v in reg and len(future[v])==1),v not in reg))
  terms.remove(first)
  if first in reg and len(future[first])==1:
   acc=reg[first];consume(first,pos,preserve=acc);owner[acc]=k;reg[k]=acc
  elif strategy=='dfs2' and first in reg and terms:
   acc=alloc(pos,pins=[reg[first]])
   second=terms.pop(0);operand=rr(reg[second]) if second in reg else loc(second)
   vsecond=values[reg[second]] if second in reg else memory[mem[second]]
   asm(f'vpxorq {operand}, {rr(reg[first])}, {rr(acc)}');values[acc]=values[reg[first]]^vsecond
   consume(first,pos);consume(second,pos);owner[acc]=k;reg[k]=acc
  else:
   acc=alloc(pos,pins=([reg[first]] if first in reg else []))
   if first in reg:asm(f'vmovdqa64 {rr(reg[first])}, {rr(acc)}');values[acc]=values[reg[first]]
   else:asm(f'vmovdqu64 {loc(first)}, {rr(acc)}');values[acc]=memory[mem[first]]
   consume(first,pos);owner[acc]=k;reg[k]=acc
  while terms:
   available=[v for v in terms if v in reg]
   if len(terms)>=2 and available:
    a=available[0];b=next(v for v in terms if v!=a)
    operand=rr(reg[b]) if b in reg else loc(b)
    vb=values[reg[b]] if b in reg else memory[mem[b]]
    asm(f'vpternlogq $150, {operand}, {rr(reg[a])}, {rr(acc)}');values[acc]^=values[reg[a]]^vb
    consume(a,pos);consume(b,pos);terms.remove(a);terms.remove(b)
   else:
    a=terms.pop(0);operand=rr(reg[a]) if a in reg else loc(a)
    va=values[reg[a]] if a in reg else memory[mem[a]]
    asm(f'vpxorq {operand}, {rr(acc)}, {rr(acc)}');values[acc]^=va;consume(a,pos)
  assert values[acc]==basis[k],k
 for i,k in outs:assert memory[('out',i*32)]==basis[k]
 text='// Generated fixed-register evaluation of the unchanged BCH linear map.\n#pragma once\n#include <immintrin.h>\n'
 text+=f'__attribute__((noinline,target("avx512f,avx512vl"))) static void bch_asm_{strategy}(const __m256i* __restrict a,__m256i* __restrict x) {{\n alignas(64) __m256i scratch[{max_slots or 1}];\n asm volatile(\n'
 text+=''.join('  "'+l+'\\n\\t"\n' for l in lines)
 text+='  : : [in] "r"(a), [out] "r"(x), [tmp] "r"(scratch)\n  : '+','.join('"ymm'+str(i)+'"' for i in range(32))+',"memory");\n}\n'
 text=text.replace('bch_asm_dfs2','bchForwardRegister')
 text=text.replace('#include <immintrin.h>\n','#include <immintrin.h>\nnamespace spin::detail::kernel {\n')+'}\n'
 (root/'generated/BchRegisterSchedule.h').write_text(text)
 return dict(strategy=strategy,instructions=len(lines),scratch_bytes=max_slots*32,counts=dict(stats),symbolic_basis_check=True)
results=[generate('dfs2')]
print(json.dumps(results,indent=2))
