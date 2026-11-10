"""Bounded all65-shell q1 replay of the new compact physical family; NOT Lean proof."""
import argparse,json,subprocess,sys
from fractions import Fraction as Q
from math import comb
from pathlib import Path
from time import monotonic
from flint import arb,arb_mat,ctx
from turbo_coarse_physical import CoarsePhysical,aq,up,binary_target_floor
from diagnose_turbo_q2_common_marker import regions,endpoint,dyadic,rounded
ROOT=Path(__file__).resolve().parents[3]; D=ROOT/'research/lean/.lake/turbo'
BASE=Path('C:/Users/peter/repo/permute_conv-github-bch/tmp')
def moments(R0,R1):
    n=R0.nrows(); prev=[arb_mat([[1]+[0]*(n-1)])]
    for r in range(1,65):
        nxt=[]
        for v in range(r+1):
            row=arb_mat(1,n)
            if v<r:row+=prev[v]*R0
            if v:row+=prev[v-1]*R1
            nxt.append(rounded(row))
        prev=nxt
    return [up(sum((row[0,j] for j in range(n)),arb(0))/comb(64,v)) for v,row in enumerate(prev)]
def run(output):
    if output.exists():raise FileExistsError('fresh output required')
    ctx.prec=256;started=monotonic()
    model=CoarsePhysical(json.loads((D/'coarse_packet_profiles.json').read_text()))
    result=dict(schema='turbo-q1-coarse-complete65-diagnostic-1',status='partial',scope='New compact physical candidates, complete65-shell replay; not Lean authenticated.',trials=[],geometries={})
    def save():
        result['elapsed_seconds']=monotonic()-started;output.write_text(json.dumps(result,indent=2)+'\n')
    save()
    for name,receipt,sparse,diag,q2file,L,bits in (
      ('k18','k18-rs8-gf32-s22-v1.json','rs16-s22-k262144-whole-fresh-p256-q1.json','k18_coarse_full_diagnostic.json','q2_common_marker_diagnostic.json',2048,Q('58.24')),
      ('k20','k20-two-ms-rs8-transfer-v1.json','k20-small-outer-q1-p192-v1.json','k20_coarse_extended_diagnostic.json','q2_common_marker_k20_full_recombination.json',8192,Q('51.82'))):
        reference=json.loads((BASE/receipt).read_text());oldq1=json.loads((BASE/sparse).read_text())
        counts=list(map(Q,reference['exact_expected_shells']))
        tailrec=json.loads((D/diag).read_text());tail=dyadic(tailrec['exact_tail_union'])
        q2rec=json.loads((D/q2file).read_text())
        q2=dyadic(q2rec['best']['k18']['q2_upper'] if name=='k18' else q2rec['q2_upper'])
        rare=dyadic(reference['component_uppers']['rare']) if name=='k20' else arb(0)
        threshold=(L*256)//10;best=[arb(1) for _ in range(65)];choices=[None]*65
        for tilt in map(Q,oldq1['tilts']):
            R=regions(model.family(tilt),L//16)
            m=moments(R[0],R[1]);factor=up((aq(tilt)*threshold).exp())
            for v in range(65):
                bound=up(factor*m[v])
                if bound<best[v]:best[v]=bound;choices[v]=str(tilt)
            q1=up(L*sum((aq(counts[v])*best[v] for v in range(65)),arb(0)))
            whole=up(q1+q2+rare+tail)
            rec=dict(geometry=name,tilt=str(tilt),q1_bits=float(-q1.log()/arb(2).log()),whole_bits=float(-whole.log()/arb(2).log()),target_met=bool(whole<aq(binary_target_floor(bits))))
            result['trials'].append(rec)
            result['geometries'][name]=dict(**rec,q1_upper=endpoint(q1),whole_upper=endpoint(whole),covered_shells=65,
              support_probability_uppers=list(map(endpoint,best)),support_choices=choices,
              actual_retained_tilt_count=len(oldq1['tilts']),retained_tilts=oldq1['tilts'])
            print(json.dumps(rec),flush=True);save()
            if monotonic()-started>105:return
    result['status']='complete';save()
    print(json.dumps(dict(status='complete',geometries={k:{z:v[z] for z in ('q1_bits','whole_bits','target_met')} for k,v in result['geometries'].items()},elapsed_seconds=result['elapsed_seconds'])),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--worker',action='store_true');a=p.parse_args()
    if a.worker:run(a.output)
    else:
      try:raise SystemExit(subprocess.run([sys.executable,__file__,'--output',str(a.output),'--worker'],timeout=115).returncode)
      except subprocess.TimeoutExpired:print('Bounded q1 replay timed out; partial output retained',file=sys.stderr);raise SystemExit(2)
