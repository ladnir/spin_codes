"""Exact recombination of saved common-marker bounds; no new numeric search."""
from fractions import Fraction as Q
from pathlib import Path
import json
from flint import arb, ctx
from turbo_coarse_physical import binary_target_floor
ROOT=Path(__file__).resolve().parents[3]
D=ROOT/'research/lean/.lake/turbo'
BASE=Path('C:/Users/peter/repo/permute_conv-github-bch/tmp')
def dyadic(x):
    if isinstance(x,dict): n,e=int(x['numerator_hex'],16),x['exponent']
    else: n,e=x
    return Q(n*(1<<e)) if e>=0 else Q(n,1<<(-e))
def ceil_endpoint(x):
    n=(x.numerator*(1<<256)+x.denominator-1)//x.denominator
    return [n,-256]
def display_bits(x):
    v=arb(x.numerator)/x.denominator
    return float(-v.log()/arb(2).log())
ctx.prec=256
full=json.loads((D/'q2_common_marker_diagnostic.json').read_text())['best']['k18']
heavy=json.loads((D/'q2_common_marker_k20_scaled_diagnostic.json').read_text())['best']['k20']
a=json.loads((BASE/'k18-rs8-gf32-s22-v1.json').read_text())
b=json.loads((BASE/'k20-two-ms-rs8-transfer-v1.json').read_text())
if list(map(Q,a['exact_expected_shells']))!=list(map(Q,b['exact_expected_shells'])):raise ArithmeticError('current Cv mismatch')
if full['marker']!=heavy['marker']:raise ArithmeticError('marker mismatch')
ratio=Q(full['outer_weight_sum'])/Q(heavy['outer_weight_sum'])
q2=dyadic(heavy['q2_upper'])*ratio**2
other=dyadic(b['component_uppers']['q1'])+dyadic(b['component_uppers']['rare'])
tail=json.loads((D/'k20_coarse_extended_diagnostic.json').read_text())
other+=dyadic(tail['exact_tail_union'])
whole=q2+other
result=dict(schema='turbo-q2-common-marker-k20-full-recombination-1',
 scope='Exact rational recombination of untrusted saved diagnostic upper bounds; not a Lean certificate.',
 marker=heavy['marker'],tilt=heavy['tilt'],covered_shells=65,covered_triangular_pairs=2145,
 min_support=0,truncated_q2_le_unconditional_full=True,rare_paid_once=True,
 outer_weight_ratio=str(ratio),q2_upper=ceil_endpoint(q2),whole_upper=ceil_endpoint(whole),
 q2_margin_bits_display=display_bits(q2),whole_margin_bits_display=display_bits(whole),
 target_51_82_met=whole<binary_target_floor(Q('51.82')))
p=D/'q2_common_marker_k20_full_recombination.json'
if p.exists():raise FileExistsError('fresh output required')
p.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ('q2_margin_bits_display','whole_margin_bits_display','target_51_82_met','covered_shells','covered_triangular_pairs','rare_paid_once')}))
