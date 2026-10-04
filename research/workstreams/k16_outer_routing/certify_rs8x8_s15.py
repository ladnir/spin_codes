"""Fresh sparse proof plus exact authenticated tail transfer for RS8x8/s15.

Each128-to256bit group has8 GF16 RS[8,4] rows and independent transitive
32-bit symbol maps. Its64four-bit packets and all regional slots are shuffled
independently. The inner is exactly the completed paired-s15-drop10 map.
The old numerical tail is reused only after its receipt, all106 source pins,
geometry, inner map, and exact dyadic sum have been authenticated.
For q>=3 the new pointwise envelope differs by (65535/65537)^(4q).
This is a full bound composition, not an all-components-fresh replay.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys
from time import monotonic

HERE=Path(__file__).resolve().parent
PROOF=HERE.parent/'k16_codesign_100us'/'proof'
sys.path.insert(0,str(PROOF))
import screen_paired_restriction as restriction
import reproduce_disjoint_pairs as shared
from flint import arb,ctx
from rs_outer import expected_group_support_counts
from rs_uniform_envelope import UniformInputEnvelope

screen=shared.screen
OLD_RECEIPT=PROOF/'paired-s15-drop10-whole-p256.json'
OLD_SHA256='89bed0c4442349a520d54b1aa6e914eac326beb3530123f8273f78fc764b867d'
Q1_TILTS=('.00256','.004','.00512','.006','.01024')
Q2_TILTS=('.00512','.01024','.01536')


def fraction(endpoint):
    if not isinstance(endpoint,(tuple,list)) or len(endpoint)!=2:
        raise ValueError('dyadic endpoint required')
    m,e=endpoint
    if type(m) is not int or type(e) is not int or m<=0:
        raise ValueError('positive integer dyadic endpoint required')
    return Q(m)*Q(2)**e


def ceil_fraction(value,precision=256):
    """Least precision-bit positive dyadic endpoint above an exact rational."""
    if not isinstance(value,Q) or value<=0:raise ValueError('positive rational required')
    p,q=value.numerator,value.denominator
    k=p.bit_length()-q.bit_length()
    if (p<q<<k) if k>=0 else (p<<(-k)<q):k-=1
    exponent=k-precision+1
    num,den=(p,q<<exponent) if exponent>=0 else (p<<(-exponent),q)
    mantissa=(num+den-1)//den
    if mantissa==(1<<precision):mantissa>>=1;exponent+=1
    if not (1<<(precision-1))<=mantissa<(1<<precision):
        raise ArithmeticError('normalization failed')
    if fraction((mantissa,exponent))<value or fraction((mantissa-1,exponent))>=value:
        raise ArithmeticError('not the least upward endpoint')
    return [mantissa,exponent]


def authenticated_old():
    raw=OLD_RECEIPT.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=OLD_SHA256:raise ArithmeticError('old receipt changed')
    old=json.loads(raw)
    pins=old['source_sha256']
    if len(pins)!=106 or any(hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h for p,h in pins.items()):
        raise ArithmeticError('old mathematical source changed')
    required=dict(K=65536,N=131072,threshold=13107,groups=512,group_dimension=128,
        regions=64,packet_bits=4,physical_t=64,state_bits=15,physical_steps=2048,
        physical_steps_per_region=32,zero_initial_state=True,final_flush=False,
        all_occupancies_covered=True,whole_code_certificate=True,target_met=True,precision=256)
    if any(old.get(k)!=v for k,v in required.items()):raise ArithmeticError('old premises mismatch')
    if set(old['occupancy_uppers'])!={str(q) for q in range(1,513)}:
        raise ArithmeticError('incomplete old occupancies')
    total=sum(map(fraction,old['occupancy_uppers'].values()),Q())
    saved=Q(int(old['union_exact_mantissa_hex'],16))*Q(2)**old['union_exact_exponent']
    if total!=saved or fraction(old['union_upper'])<total or not fraction(old['union_upper'])<Q(1,2**40):
        raise ArithmeticError('invalid old exact union')
    if fraction(ceil_fraction(total))!=fraction(old['union_upper']):
        raise ArithmeticError('old final upward rounding mismatch')
    return old


def certify(output):
    output=Path(output).resolve()
    if output.exists():raise FileExistsError('fresh output required')
    ctx.prec=256;start=monotonic();old=authenticated_old()
    data,map_record=restriction.prepare((10,))
    map_record=json.loads(json.dumps(map_record))
    if map_record!=old['map_record']:raise ArithmeticError('fresh inner map changed')
    env=UniformInputEnvelope(8,4,4,8);env.verify_shell_domination()
    old_env=UniformInputEnvelope(16,8,4,4);ratio=env.beta/old_env.beta
    if ratio!=Q(65535,65537)**4 or Q(old['beta'])!=old_env.beta:
        raise ArithmeticError('exact pointwise envelope transfer mismatch')
    counts=expected_group_support_counts(8,4,4,8)
    if sum(counts)!=(1<<128)-1:raise ArithmeticError('outer count is incomplete')
    actual=shared.sources()
    new_pins={str(Path(__file__).resolve()):hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    for path,digest in actual.items():
        if path in old['source_sha256']:
            if old['source_sha256'][path]!=digest:raise ArithmeticError('loaded old source mismatch')
        else:new_pins[path]=digest
    values={};record=dict(schema='rs8x8-paired15-complete-composition-1',
        K=65536,N=131072,threshold=13107,minimum_distance_outside_failure=13108,
        target_margin_bits=40,groups=512,group_dimension=128,regions=64,packet_bits=4,
        physical_t=64,state_bits=15,physical_steps=2048,physical_steps_per_region=32,
        zero_initial_state=True,final_flush=False,state_continuity='retained across every step and region',
        outer='eight parallel GF16 RS[8,4] rows; systematic evaluations0..3 and parity4..7',
        outer_randomness='independent invertible32-bit symbol maps with uniform nonzero fixed-input law',
        coordinate_order='symbol-major; coordinate4*row+bit inside each32-bit symbol; eight contiguous4-bit packets',
        inner_randomness='independent uniform GL15 at every physical t64 step',
        routing_randomness='independent uniform per-group packet permutations and regional slot permutations',
        map_record=map_record,outer_envelope=env.metadata(),precision=256,
        q1_tilts=Q1_TILTS,q2_tilts=Q2_TILTS,q1_fresh=True,q2_fresh=True,
        fresh_replay=False,proof_method='fresh q1/q2 plus exact pointwise-envelope transfer of authenticated q>=3 endpoints',
        old_receipt=dict(path=str(OLD_RECEIPT.resolve()),sha256=OLD_SHA256,verified_source_count=106),
        new_source_sha256=new_pins,all_actual_math_sources_verified=True,
        tail_transfer=dict(beta_old=str(old_env.beta),beta_new=str(env.beta),ratio=str(ratio),
            exponent='q',occupancies=[3,512],same_physical_and_regional_operators=True),
        all_occupancies_covered=False,whole_code_certificate=False,
        scope='Ideal independent setup ensemble; not a deterministic-seed guarantee. Tail endpoints are authenticated and transferred, not recomputed.')

    def checkpoint():
        if ctx.prec!=256 or authenticated_old()!=old:
            raise ArithmeticError('precision or authenticated premises changed')
        if any(hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h for p,h in new_pins.items()):
            raise ArithmeticError('new mathematical source changed')
        record.update(occupancy_uppers={str(q):v for q,v in values.items()},
            occupancy_margin_bits={str(q):shared.margin(arb(v[0])*arb(2)**v[1]) for q,v in values.items()},
            elapsed_seconds=monotonic()-start)
        output.write_text(json.dumps(record,indent=2)+'\n')

    one=[arb(1)]*65
    for tilt in Q1_TILTS:
        local=screen.q1.kernel_t64.local_operators(data,Q(tilt),Q(1,2))
        regional=screen.q1.placement(local,epochs=16,windows=32,rounding=screen.q1.rounded,maximum_groups=1)
        factor=(screen.aq(Q(tilt))*13107).exp()
        moments=screen.q1.support_moments(regional[0],regional[1],64)
        one=[min(a,screen.up(factor*b)) for a,b in zip(one,moments)]
    values[1]=screen.q1.endpoint(screen.up(512*sum((screen.aq(c)*m for c,m in zip(counts,one)),arb(0))))
    record['q1_support_uppers']=[screen.q1.endpoint(x) for x in one]
    checkpoint();print('q1='+record['occupancy_margin_bits']['1'],flush=True)
    two=[[arb(1)]*(v+1) for v in range(65)]
    for tilt in Q2_TILTS:
        local=screen.q1.kernel_t64.local_operators(data,Q(tilt),Q(1,2))
        regional=screen.q1.placement(local,epochs=16,windows=32,rounding=screen.q1.rounded,maximum_groups=2)
        factor=(screen.aq(Q(tilt))*13107).exp()
        moments=screen.q2.pair_support_moments(regional,64)
        for v,column in enumerate(moments):
            for u,m in enumerate(column):two[v][u]=min(two[v][u],screen.up(factor*m))
    values[2]=screen.q1.endpoint(screen.q2.fold_shell_pairs(counts,two,512))
    checkpoint();print('q2='+record['occupancy_margin_bits']['2'],flush=True)
    for q in range(3,513):values[q]=ceil_fraction(fraction(old['occupancy_uppers'][str(q)])*ratio**q)
    if set(values)!=set(range(1,513)):raise ArithmeticError('incomplete new coverage')
    total=shared.exact_sum.dyadic_sum(values.values())
    endpoint=shared.exact_sum.rounded_endpoint(total,256)
    passed=fraction(endpoint)<Q(1,2**40)
    record.update(all_occupancies_covered=True,occupancy_covered=[1,512],target_met=passed,
        whole_code_certificate=passed,union_exact_mantissa_hex=hex(total[0]),union_exact_exponent=total[1],
        union_upper=endpoint,margin_bits=shared.margin(arb(endpoint[0])*arb(2)**endpoint[1]),
        arithmetic='Exact rational tail scaling and upward256-bit rounding; exact dyadic total and final upward rounding; exact target comparison')
    checkpoint();print('Complete union='+record['margin_bits'],flush=True)
    if not passed:raise ArithmeticError('new complete union did not meet target')
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    certify(parser.parse_args().output)
