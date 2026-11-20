"""Exact shared-row parity mixing for the actual BCH tail pair experiment.

The result concerns the 256 column parities, not all B-syndrome checks.
It is uniform over both occupied-row supports of size Q=2620.
"""
import math
from fractions import Fraction as F
from pathlib import Path
import bridge as base
import screen_exponential_modes as cap_source


def kraw_values(n,w):
    values=[1,n-2*w]
    for j in range(1,n):
        value,remainder=divmod((n-2*w)*values[j]-(n-j+1)*values[j-1],j+1)
        assert remainder==0;values.append(value)
    return values


def run():
    caps,sources=cap_source.latest_caps()
    tail_path=base.HERE/'generated/joint_tail_lower_80/lower.json';tail=base.read(tail_path)
    for name,digest in tail['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    for name,digest in tail['outer_sha256'].items():assert base.sha(base.BCH/name)==digest
    assert tail['rational_primal_dual_checks_passed'];a0=tail['lower']
    close_path=base.HERE/'generated/bch_tail_overlap_outward.json';close=base.read(close_path)
    for name,digest in close['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    close_upper=base.decode(close['close_pair_probability_upper']);assert close_upper<F(1,128)
    table={w:kraw_values(256,w) for w in range(38,161,2)}
    single=max(F(sum(caps[w]*abs(table[w][j]) for w in range(38,81,2)),a0*math.comb(256,j)) for j in range(1,256))
    far=max(F(abs(table[w][j]),math.comb(256,j)) for w in range(82,161,2) for j in range(1,256))
    assert single<F(2,5) and far==F(23,64)
    second=close_upper+(1-close_upper)*far;assert second<F(3,8)
    q=2620
    epsilon=((1<<256)-2)*F(2,5)**q+F(((1<<256)-2)**2,4)*F(3,8)**q
    assert epsilon<F(1,1<<3100)
    result=dict(status='EXACT_SHARED_ROW_PAIR_PARITY_MIXING',occupation=q,
        single_character_absolute_upper=base.encode(single),nontrivial_character_second_moment_upper=base.encode(second),
        parity_pair_uniform_relative_error_upper=base.encode(epsilon),relative_error_below_2_to_minus_3100=True,
        parity_pair_reference_probability=base.encode(F(1,1<<510)),uniform_over_occupied_supports=True,
        full_second_moment_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__),tail_path,close_path]+sources})
    path=base.HERE/'generated/tail_parity_mixing.json'
    if path.exists():assert result==base.read(path)
    else:base.write_new(path,result)
    print('Exact parity mixing: single gap',float(single),'pair gap',float(second),
        'relative error <2^-3100; diagnostic error bits',math.log2(epsilon.denominator)-math.log2(epsilon.numerator),flush=True)


if __name__=='__main__':run()
