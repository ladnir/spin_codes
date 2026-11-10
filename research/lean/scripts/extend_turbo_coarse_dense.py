"""Bounded scalar-only diagnostic extending existing dense coefficient witnesses."""
import json,sys
from pathlib import Path
from fractions import Fraction as Q
from math import comb
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'workstreams/k18_margin60'))
import replay_fractional as replay
from turbo_coarse_physical import aq,binary_target_floor
from flint import arb,ctx
ctx.prec=256
for label,L,bits,path in [('k18',2048,Q(1456,25),'k18_coarse_full_diagnostic.json'),('k20',8192,Q(2591,50),'k20_coarse_full_diagnostic.json')]:
 src=Path(__file__).resolve().parents[1]/'.lake/turbo'/path
 x=json.loads(src.read_text()); original=x['occupancy_uppers']; chosen=dict(original); changes={}
 beta=Q(1<<256,65535**8)
 qs=range(3,385) if L==2048 else range(248,616)
 logs={q:arb(comb(L,q)).log() for q in qs}
 for w in x['witnesses']:
  if not w['kind'].startswith('iid'):continue
  p=Q(w['marker_probability'])
  if p==1:continue
  alpha=Q(w.get('alpha','1')); q0=w['minimum']; lb0=arb(comb(L,q0)).log()
  initial=aq(replay.fraction(original[str(q0)])).log()
  slope=aq(alpha)*aq(beta).log()-64*(aq(p).log()-aq(1-p).log())
  for q in qs:
   ep=replay.endpoint((initial-63*(logs[q]-lb0)+(q-q0)*slope).exp())
   if replay.fraction(ep)<replay.fraction(chosen[str(q)]):chosen[str(q)]=ep;changes[str(q)]=dict(w)
 vals=list(chosen.values());emin=min(e for n,e in vals)
 summed=Q(sum(n<<(e-emin) for n,e in vals))*Q(2)**emin
 if L==2048:
  total=summed;tail=total-replay.fraction(chosen['1'])-replay.fraction(chosen['2'])
 else:
  tail=summed;total=tail+sum((replay.fraction(v) for v in x['retained_components'].values()),Q())
 result=dict(schema='turbo-coarse-extended-diagnostic-2',scope='Outward Python diagnostic only; existing dense witnesses extended by coefficient-ratio identity.',source_checkpoint=str(src),occupancy_uppers=chosen,changed_witnesses=changes,exact_tail_union=replay.encode_exact(tail),exact_union=replay.encode_exact(total),tail_margin_bits_display=replay.bits(tail),margin_bits_display=replay.bits(total),rounded_paper_target_met=total<binary_target_floor(bits))
 out=src.with_name(label+'_coarse_extended_diagnostic.json');out.write_text(json.dumps(result,indent=2)+'\n')
 print(label,'changed',len(changes),'range',([min(map(int,changes)),max(map(int,changes))] if changes else []),'tailbits',replay.bits(tail),'wholebits',replay.bits(total),'target',result['rounded_paper_target_met'],flush=True)
