"""Exact shell caps from the already certified strength-29 outer array.

For a uniform C word, polynomials of its weight of degree <=29 have
the Bin(256,1/2) expectation. Squaring a degree-14 reproducing kernel
therefore bounds each shell without solving another linear program.
"""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path
import bridge as base
import occupation_two as q2


def kraw(n,j,w):
    return sum((-1)**i*math.comb(w,i)*math.comb(n-w,j-i)
               for i in range(max(0,j-(n-w)),min(w,j)+1))


def shell_bound(n,k,degree,w):
    kernel=sum((F(kraw(n,j,w)**2,math.comb(n,j)) for j in range(degree+1)),F(0))
    return F(1<<k)/kernel


def deterministic_caps():
    original=q2.deterministic_caps()
    return {w:min(cap,math.floor(shell_bound(256,128,14,w))) if cap else 0 for w,cap in original.items()}


def build(verify=False):
    from certify_bch_shift_rank_split import build as rank_build
    proof=rank_build();source=base.BCH/'generated/bch256_shift_rank_q30_refined.json'
    assert proof==base.read(source)
    assert proof['Qdual_extended_distance_lower']>=30
    caps=deterministic_caps();old=q2.deterministic_caps()
    path=base.HERE/'generated/christoffel_oa29_caps.json'
    local,outer=q2.source_paths();local.append(Path(__file__));outer.append(source)
    outer += [base.BCH/name for name in proof['source_sha256']]
    payload=dict(status='EXACT_OA29_CHRISTOFFEL_SHELL_CAPS',outer_length=256,outer_dimension=128,
        polynomial_degree=14,required_moment_degree=28,certified_OA_strength=29,
        premise='The fixed Q is contained in C, so C dual is contained in Q dual. This reuses the fixed-code containment audit.',
        rows=[dict(weight=w,cap=caps[w],old_cap=old[w],christoffel_upper=base.encode(shell_bound(256,128,14,w))) for w in base.WEIGHTS],
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in local},
        outer_sha256={str(p.relative_to(base.BCH)):base.sha(p) for p in outer})
    if verify:assert payload==base.read(path)
    else:base.write_new(path,payload)
    print('Exact OA29 caps replayed' if verify else 'Exact OA29 caps written',flush=True)
    print([(w,round(math.log2(old[w]),3),round(math.log2(caps[w]),3)) for w in (50,54,58,62,66,70,78,90,110,128)],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();build(a.verify)
