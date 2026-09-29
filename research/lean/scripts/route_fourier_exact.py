"""Untrusted proposals for exact kernel-checked Fourier integer witnesses."""
from pathlib import Path
from fractions import Fraction as Q
from decimal import Decimal,localcontext
import json,math,argparse,sys
sys.set_int_max_str_digits(0)
R=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('index',type=int);args=p.parse_args()
D=json.loads((R.parent/'workstreams/inner_design/imt_asymptotic/d11/DENSE_REPLAY.json').read_text())
seen=set();uni=[]
for b in D['leaves']:
 w=b['rational_witness']
 if w['family']!='fourier':continue
 key=json.dumps(w,sort_keys=True)
 if key not in seen:seen.add(key);uni.append(b)
b=uni[args.index];w=b['rational_witness'];q=Q(w['p'])*Q(w['y'])
with localcontext() as c:
 c.prec=110;ell=Q(w['log_lam']);z=Q(int((-(Decimal(ell.numerator)/Decimal(ell.denominator)).exp()).exp()*10**30),10**30)
qn,qd,zn,zd=q.numerator,q.denominator,z.numerator,z.denominator
A=(qd-qn)*zd+qn*zn;B=qn*zd+(qd-qn)*zn;U=abs((qd-qn)*zd-qn*zn);V=abs(qn*zd-(qd-qn)*zn)
assert A>0 and B>0
arrays=[[a**i for i in range(129)] for a in [A,B,U,V]]
pa,pb,pu,pv=arrays
M=json.loads((R/'scripts/sparse_data/model.json').read_text())
import re
source=(R/'SpinCodes/Structured/FiberNumericsData/Tables.lean').read_text()
spectrum=json.loads(re.search(r'def spectrum : List Nat := (\[[^\n]+\])',source).group(1))
vv=[round(Q(a)*10**20) for a in w['vector']]
weighted=vv[0]+sum(c*v for c,v in zip(M['counts'],vv[2:]))
def mon(d,t,h):return pa[128-d-(t-h)]*pb[d-h]*pu[t-h]*pv[h]
def cap(d):return 524287*sum(c*max(mon(d,t,max(0,d+t-128)),mon(d,t,min(d,t))) for t,c in enumerate(spectrum) if c)+524288*pa[128-d]*pb[d]
targets=[cap(d) for d in M['levels']]
zero=sum((k*vv[0]+(math.comb(128,j)-k)*vv[1])*(qn*zn)**j*((qd-qn)*zd)**(128-j) for j,k in enumerate(M['kernel']))
den=(qd*zd)**128;factor=2*524287*524288
values=[Q(zero,den*vv[0]),Q(max(targets)*weighted,factor*den*vv[1])]+[Q(a*weighted,factor*den*v) for a,v in zip(targets,vv[2:])]
exact=max(values);radius=Q(-((-exact.numerator*10**100)//exact.denominator),10**100)
assert 0<radius<1
prefactor=-((-vv[0])//min(vv))
tag=f'W{args.index:03d}'
lines=['import SpinCodes.Structured.DenseFourierExactColumnDefs','',f'namespace Spin.Structured.DenseFourierExact.{tag}','set_option maxHeartbeats 0','set_option maxRecDepth 1000000',f'def qn : Int := {qn}',f'def qd : Int := {qd}',f'def zn : Int := {zn}',f'def zd : Int := {zd}',f'def rn : Int := {radius.numerator}',f'def rd : Int := {radius.denominator}',f'def v : ICoords := ⟨{vv[0]}, {vv[1]}, fun i => ([{", ".join(map(str,vv[2:]))}]).getD i 0⟩','theorem check_Z : checkZ qn qd zn zd rn rd v=true := by decide','theorem check_D : checkD qn qd zn zd rn rd v=true := by decide','theorem check_S : ∀ i:Fin 5,checkS qn qd zn zd rn rd v i=true := by\n  intro i\n  match i with'+''.join(f'\n  | ⟨{i}, _⟩ => change checkS qn qd zn zd rn rd v ({i}:Fin 5)=true; decide' for i in range(5)),f'end Spin.Structured.DenseFourierExact.{tag}']
(R/f'SpinCodes/Structured/DenseFourierExact{tag}Data.lean').write_text('\n'.join(lines)+'\n',encoding='utf-8')
(R/f'scripts/map_data/dense_fourier_exact_{tag}_candidate.json').write_text(json.dumps({'status':'UNTRUSTED_CANDIDATE','source_box':b,'q':str(q),'z':str(z),'radius':str(radius),'vector':vv,'prefactor':prefactor,'unique_witnesses':len(uni)},indent=2)+'\n')
print(tag,'proposed','radius',float(radius),'z',float(z),'prefactor',prefactor,'unique',len(uni))
sem=[f'import SpinCodes.Structured.DenseFourierExact{tag}Data','import SpinCodes.Structured.DenseFourierExactCheck','', 'noncomputable section',f'namespace Spin.Structured.DenseFourierExact.{tag}', 'open Spin.Imt',f'def wmin : ℝ := {min(vv)}',f'def prefactor : ℝ := {prefactor}', 'theorem collatz : ((ConcreteFourier.matrix ((qn:ℝ)/qd) ((zn:ℝ)/zd)).applyCol v.real).le', '    (Coords.smul ((rn:ℝ)/rd) v.real) :=', '  checked_collatz (by decide) (by decide) (by decide) (by decide) (by decide)', '    v (by decide) check_Z check_D check_S', 'theorem witness_floor : 0 < wmin ∧ wmin ≤ v.real.Z ∧ wmin ≤ v.real.D ∧ ∀ i,wmin ≤ v.real.S i := by', '  refine ⟨by norm_num [wmin], by norm_num [wmin,v,ICoords.real], by norm_num [wmin,v,ICoords.real], ?_⟩', '  intro i', '  fin_cases i <;> norm_num [wmin,v,ICoords.real]', 'theorem witness_prefactor : v.real.Z/wmin ≤ prefactor := by norm_num [v,ICoords.real,wmin,prefactor]', '#print axioms collatz', '#print axioms witness_floor', '#print axioms witness_prefactor', f'end Spin.Structured.DenseFourierExact.{tag}']
(R/f'SpinCodes/Structured/DenseFourierExact{tag}.lean').write_text('\n'.join(sem)+'\n',encoding='utf-8')
