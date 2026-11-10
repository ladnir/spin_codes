"""Bounded q2 common-marker proof diagnostic; NOT a Lean certificate or benchmark."""
import argparse, hashlib, json, subprocess, sys
from fractions import Fraction as Q
from math import comb
from pathlib import Path
from time import monotonic
import numpy as np
from flint import arb, arb_mat, ctx
from turbo_coarse_physical import CoarsePhysical, aq, up, binary_target_floor
ROOT = Path(__file__).resolve().parents[3]
BASE = Path('C:/Users/peter/repo/permute_conv-github-bch/tmp')

def endpoint(x): return [int(v) for v in up(x).man_exp()]
def dyadic(x):
    if isinstance(x, dict): return arb(int(x['numerator_hex'],16)) * arb(2)**x['exponent']
    return arb(x[0]) * arb(2)**x[1]
def rounded(M): return arb_mat([[up(M[i,j]) for j in range(M.ncols())] for i in range(M.nrows())])
def poly_mul(A,B):
    n=A[0].nrows(); out=[]
    for k in range(3):
        M=arb_mat(n,n)
        for i in range(k+1): M += A[i]*B[k-i]
        out.append(rounded(M))
    return out
def regions(family,R):
    n=family[0].nrows()
    A=[family[0], 16*family[1],120*family[2]]
    V=[arb_mat([[int(i==j) for j in range(n)] for i in range(n)]),arb_mat(n,n),arb_mat(n,n)]
    e=R
    while e:
        if e&1: V=poly_mul(V,A)
        e//=2
        if e: A=poly_mul(A,A)
    return [rounded(V[j]/comb(16*R,j)) for j in range(3)]
def collatz(M):
    n=M.nrows(); mf=np.array([[float(M[i,j]) for j in range(n)] for i in range(n)])
    eig=max(np.linalg.eigvals(mf).real)
    best=None
    for gap in (1/256,1/64,1/16):
        f=np.linalg.solve((eig+gap)*np.eye(n)-mf,np.ones(n))
        if not np.all(f>0): continue
        f=f/f[0]
        v=[Q(format(float(x),'.17g')) for x in f]
        ratios=[up(sum((M[i,j]*aq(v[j]) for j in range(n)),arb(0))/aq(v[i])) for i in range(n)]
        rho=max(ratios)
        c=Q(v[0],min(v))
        bound=up(aq(c)*rho**64)
        rec=dict(rho=endpoint(rho),vector=list(map(str,v)),prefactor=str(c),moment_upper=endpoint(bound))
        if best is None or bound<best[0]: best=(bound,rec)
    if best is None: raise ArithmeticError('no positive Collatz witness')
    return best

def run(output,k20_scaled=False):
    if output.exists(): raise FileExistsError('fresh output required')
    ctx.prec=256; started=monotonic()
    cp=ROOT/'research/lean/.lake/turbo/coarse_packet_profiles.json'
    model=CoarsePhysical(json.loads(cp.read_text()))
    records={}
    for name,receipt,diag,L,start,bits in (
      ('k18','k18-rs8-gf32-s22-v1.json','k18_coarse_full_diagnostic.json',2048,0,Q('58.24')),
      ('k20','k20-two-ms-rs8-transfer-v1.json','k20_coarse_extended_diagnostic.json',8192,9,Q('51.82'))):
        if k20_scaled and name!='k20': continue
        ref=json.loads((BASE/receipt).read_text()); tail=json.loads((ROOT/'research/lean/.lake/turbo'/diag).read_text())
        counts=list(map(Q,ref['exact_expected_shells']))
        if len(counts)!=65 or sum(counts)!=(1<<128)-1: raise ArithmeticError('count normalization')
        q1=dyadic(ref['component_uppers']['q1'])
        rare=dyadic(ref['component_uppers']['rare']) if name=='k20' else arb(0)
        tkey='exact_tail_union' if 'exact_tail_union' in tail else 'exact_tail'
        if tkey not in tail: raise KeyError(str(tail.keys()))
        tv=dyadic(tail[tkey])
        records[name]=dict(groups=L,R=L//16,cutoff=(L*256)//10,min_support=start,counts=counts,
          q1=q1,rare=rare,tail=tv,target=aq(binary_target_floor(bits)),best=None)
    result=dict(schema='turbo-q2-common-marker-diagnostic-1',status='partial',
      scope='Untrusted compact-count numerical diagnostic; complete2145shell marker bound, not Lean authenticated.',
      candidate_sha256=hashlib.sha256(cp.read_bytes()).hexdigest(),trials=[],best={})
    def save():
        result['elapsed_seconds']=monotonic()-started
        output.write_text(json.dumps(result,indent=2)+'\n')
    save()
    tilt_grid=('1/5000','3/10000','1/2500','3/5000','1/1250') if k20_scaled else ('1/1250','3/2500','1/625','3/1250','2/625')
    for tilt in map(Q,tilt_grid):
        family=model.family(tilt)
        for name,r in records.items():
            R=regions(family,r['R'])
            factor=up((aq(tilt)*r['cutoff']).exp())
            for x in (8,15,24):
                marker=R[0]+(2*x)*R[1]+(x*x)*R[2]
                moment,witness=collatz(marker)
                outer=sum((r['counts'][v]/(comb(64,v)*x**v) for v in range(r['min_support'],65)),Q())
                q2=up(comb(r['groups'],2)*factor*moment*aq(outer)**2)
                whole=up(q2+r['q1']+r['rare']+r['tail'])
                rec=dict(tilt=str(tilt),marker=x,q2_bits=float(-q2.log()/arb(2).log()),
                  whole_bits=float(-whole.log()/arb(2).log()),target_met=bool(whole<r['target']))
                result['trials'].append(dict(geometry=name,**rec))
                if r['best'] is None or q2<r['best']:
                    r['best']=q2
                    result['best'][name]=dict(**rec,q2_upper=endpoint(q2),whole_upper=endpoint(whole),
                      outer_weight_sum=str(outer),covered_shells=65,covered_triangular_pairs=2145,
                      retained_min_support=r['min_support'],collatz=witness)
                save()
                print(name,str(tilt),x,rec,flush=True)
                if monotonic()-started>105: return
    result['status']='complete'; save()
    print(json.dumps(dict(status=result['status'],best={k:{z:v[z] for z in ('tilt','marker','q2_bits','whole_bits','target_met')} for k,v in result['best'].items()},elapsed_seconds=result['elapsed_seconds'])),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--worker',action='store_true');p.add_argument('--k20-scaled',action='store_true');a=p.parse_args()
    if a.worker: run(a.output,a.k20_scaled)
    else:
      try: raise SystemExit(subprocess.run([sys.executable,__file__,'--output',str(a.output),'--worker']+(['--k20-scaled'] if a.k20_scaled else []),timeout=115).returncode)
      except subprocess.TimeoutExpired: print('Bounded diagnostic timed out; partial output retained',file=sys.stderr);raise SystemExit(2)

