"""Propose the first occupation witness using exact outward integer arithmetic.
Generated Lean checks are the proof; this Python computation is untrusted.
"""
from pathlib import Path
from fractions import Fraction as Q
from decimal import Decimal, localcontext
import json, math, argparse
args=argparse.ArgumentParser()
args.add_argument("index",type=int)
index=args.parse_args().index
assert 1<=index<=132
tag=f"W{index:03d}"
ROOT=Path(__file__).resolve().parents[1]
S=10**30
zero=(0,0)
def sc(n):return (n,n)
def oi(n):return sc(n*S)
def frac(n,d):assert d>0;return (n*S//d,-((-n*S)//d))
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def sub(a,b):return (a[0]-b[1],a[1]-b[0])
def mul(a,b):
    p=[x*y for x in a for y in b];return (min(p)//S,-((-max(p))//S))
def div(a,d):assert d>0;return (a[0]//d,-((-a[1])//d))
def mn(a,b):return (min(a[0],b[0]),min(a[1],b[1]))
def mx(a,b):return (max(a[0],b[0]),max(a[1],b[1]))
def isum(xs):
    a=zero
    for b in xs:a=add(a,b)
    return a
M=json.loads((ROOT/'scripts/sparse_data/model.json').read_text())
D=json.loads((ROOT.parent/'workstreams/inner_design/imt_asymptotic/d11/DENSE_REPLAY.json').read_text())
seen=set();unique=[]
for box in D['leaves']:
    if box['rational_witness']['family']!='occupation': continue
    key=json.dumps(box['rational_witness'],sort_keys=True)
    if key not in seen:seen.add(key);unique.append(box)
b=unique[index]
w=b['rational_witness'];q=Q(w['p'])*Q(w['y'])
with localcontext() as ctx:
    ctx.prec=70
    ell=Q(w['log_lam']); ev=(Decimal(ell.numerator)/Decimal(ell.denominator)).exp()
    zn=int((-ev).exp()*Decimal(10**18))*10**12
z=sc(zn);v=[sc(round(Q(x)*10**20)*10**10) for x in w['vector']]
powers=[oi(1)]
for _ in range(208):powers.append(mul(powers[-1],z))
levels=M['levels'];counts=M['counts'];den=q.denominator**128
avg=div(isum(mul(oi(c),x) for c,x in zip(counts,v[2:])),524287)
def choose(n,k):return math.comb(n,k) if 0<=k<=n else 0
def pat(j,p):return isum(mul(frac(c,choose(128,j)),powers[int(out)]) for out,c in p)
def moment(w,j):return isum(mul(frac(choose(w,h)*choose(128-w,j-h),choose(128,j)),powers[w+j-2*h]) for h in range(129) if h<=j and choose(w,h)*choose(128-w,j-h))
cols=[zero]*7
for j in range(129):
    total=choose(128,j);ell=frac(total-M['kernel'][j],total)
    ms=[moment(w,j) for w in levels];a=zero
    for m in ms:a=mx(m,a)
    cd=mn(mn(mn(a,ell),mul(frac(M['caps'][j],total),powers[min(abs(w-j) for w in levels)])),a)
    cs=[mn(m,mul(frac(min(total-M['kernel'][j],c*M['caps'][j]),c*total),powers[abs(w-j)])) for w,c,m in zip(levels,counts,ms)]
    if str(j) in M['low']:
        low=M['low'][str(j)];lowmax=zero
        for pp in low['patterns']:lowmax=mx(pat(j,pp),lowmax)
        cd=mn(cd,lowmax)
        cs=[mn(c,div(pat(j,[(int(o),n) for o,n in low['by_weight'][str(w)].items()]),nc)) for c,w,nc in zip(cs,levels,counts)]
    f=[add(mul(mul(sub(oi(1),ell),powers[j]),v[0]),mul(mul(ell,powers[j]),v[1])),
       add(mul(add(div(cd,2),div(mn(a,ell),1048574)),v[0]),div(mul(a,add(v[1],avg)),2))]
    for i,m in enumerate(ms):f.append(add(mul(add(div(cs[i],2),div(mn(m,ell),1048574)),v[0]),div(mul(m,add(avg,v[i+2] if total==M['kernel'][j] else v[1])),2)))
    pr=frac(total*q.numerator**j*(q.denominator-q.numerator)**(128-j),den)
    cols=[add(c,mul(pr,fj)) for c,fj in zip(cols,f)]
rn=max(-((-a[1]*10**18)//vv[0]) for a,vv in zip(cols,v))*10**12
radius=sc(rn)
assert all(a[1]<=mul(radius,vv)[0] for a,vv in zip(cols,v))
prefactor=-((-S)//min(x[0] for x in v))
print(tag,'z',zn/S,'radius',rn/S,'prefactor',prefactor)
fmt=lambda a:f'⟨{a[0]}, {a[1]}⟩'
lines=['import SpinCodes.Structured.DenseOccupationFixedColumnDefs','', f'namespace Spin.Structured.DenseOccupationFixed.{tag}','open Spin.Numeric','set_option maxRecDepth 1000000','set_option maxHeartbeats 0',f'def qn : Int := {q.numerator}',f'def qd : Int := {q.denominator}',f'def z : Fix := Fix.sc {zn}',f'def radius : Fix := Fix.sc {rn}',f'def v : FCoords := ⟨Fix.sc {v[0][0]}, Fix.sc {v[1][0]}, fun i => ([{", ".join("Fix.sc "+str(a[0]) for a in v[2:])}]).getD i Fix.zero⟩',f'def zPowers : List Fix := [{", ".join(map(fmt,powers))}]','theorem powers_checked : checkPowers z zPowers 208 = true := by decide','def col : FCoords := column qn qd (powers z zPowers 208) v', 'theorem column_Z_checked : col.Z.hi ≤ (Fix.mul radius v.Z).lo := by decide','theorem column_D_checked : col.D.hi ≤ (Fix.mul radius v.D).lo := by decide','theorem column_S_checked : ∀ i : Fin 5, (col.S i).hi ≤ (Fix.mul radius (v.S i)).lo := by\n  intro i\n  match i with\n  | ⟨0, _⟩ => change (col.S (0:Fin 5)).hi ≤ (Fix.mul radius (v.S (0:Fin 5))).lo; decide\n  | ⟨1, _⟩ => change (col.S (1:Fin 5)).hi ≤ (Fix.mul radius (v.S (1:Fin 5))).lo; decide\n  | ⟨2, _⟩ => change (col.S (2:Fin 5)).hi ≤ (Fix.mul radius (v.S (2:Fin 5))).lo; decide\n  | ⟨3, _⟩ => change (col.S (3:Fin 5)).hi ≤ (Fix.mul radius (v.S (3:Fin 5))).lo; decide\n  | ⟨4, _⟩ => change (col.S (4:Fin 5)).hi ≤ (Fix.mul radius (v.S (4:Fin 5))).lo; decide',f'end Spin.Structured.DenseOccupationFixed.{tag}']
(ROOT/f'SpinCodes/Structured/DenseOccupationFixed{tag}Data.lean').write_text('\n'.join(lines)+'\n',encoding='utf-8')
(ROOT/f'scripts/map_data/dense_occupation_fixed_{tag}_candidate.json').write_text(json.dumps({'status':'UNTRUSTED_CANDIDATE','index':index,'prefactor':prefactor,'source_box':b,'q':str(q),'z_scaled':zn,'radius_scaled':rn,'v_scaled':[a[0] for a in v],'columns':cols},indent=2)+'\n')


semantic=(ROOT/'SpinCodes/Structured/DenseOccupationFixedFirst.lean').read_text(encoding='utf-8-sig')
semantic=semantic.replace('DenseOccupationFixedFirstData',f'DenseOccupationFixed{tag}Data')
semantic=semantic.replace('DenseOccupationFixed.First',f'DenseOccupationFixed.{tag}')
semantic=semantic.replace('372',str(prefactor))
semantic=semantic.replace('first numerical occupation witness',f'numerical occupation witness {index}')
(ROOT/f'SpinCodes/Structured/DenseOccupationFixed{tag}.lean').write_text(semantic,encoding='utf-8')
