"""Directed quiet-class coverage for the dimension-four SUM-RANK candidate.

The unchanged FFT64 construction's receipts are separate. The sparse
unweighted and dense inverse-output receipts transfer by their map-agnostic
hypotheses; the newly replayed middle intervals use exact sum-rank moments.
"""

from fractions import Fraction as F
from functools import lru_cache
from math import comb

from flint import arb,ctx

from hd_randomized_outer_joint_input import interval_proposal,directed_interval
from hd_randomized_outer_quiet_exact import _arb_fraction
from hd_sumrank_transfer import profiles,fractional_packets


@lru_cache(maxsize=20)
def enclosures(z,alpha):
    exponent=_arb_fraction(alpha)
    quiet,free=[arb(0) for _ in range(65)],[arb(0) for _ in range(65)]
    for profile,j,multiplicity,k,t,v in profiles(F(z)):
        quiet[j]+=multiplicity*_arb_fraction(k)**exponent
        free[j]+=multiplicity*_arb_fraction(t)**exponent
    return (tuple(value/comb(64,j) for j,value in enumerate(quiet)),
            tuple(value/comb(64,j) for j,value in enumerate(free)))


def receipt():
    # Every integer in each interval is covered by convexity. The center
    # merely selects a positive marker; no floating sign decides acceptance.
    selected=((960,F(1,2),F(1,8)),(1280,F(1,4),F(1,16)),
              (1560,F(1,6),F(1,16)),(2050,F(1,12),F(1,16)),
              (2740,F(1,32),F(1,16)),
              (1800,F(1,8),F(1,16)),(2400,F(1,16),F(1,16)),
              (3100,F(1,64),F(1,16)),(4200,F(1,256),F(1,8)),
              (5500,F(1,1024),F(1,8)),(6700,F(1,4096),F(1,8)))
    with ctx.workprec(192):
        receipts=[]
        covered=set(range(1,954))|set(range(7548,8193))
        for center,z,a in selected:
            moments=fractional_packets(z,a,quiet_only=True)
            proposal=interval_proposal(center,z,a,F(1),threshold=-100,moments=moments)
            if proposal is None:
                continue
            r=directed_interval(proposal['first'],proposal['last'],z,a,F(1),F(proposal['marker']),
                                enclosures=enclosures(z,a))
            if not r['below_minus_40']:
                raise ArithmeticError('Middle interval did not certify')
            receipts.append(r)
            covered.update(range(proposal['first'],proposal['last']+1))
        # These previous class bounds use only C's MDS rank and invertibility.
        inherited=arb(2)**F(-1421972,10000)+arb(2)**F(-3266848,10000)
        union=inherited+sum((arb(2)**arb(r['union_log2_upper']) for r in receipts),arb(0))
        log_union=union.log()/arb(2).log()
        missing=[r for r in range(1,8193) if r not in covered]
        return dict(construction='sumrank-s4',intervals=receipts,uncovered_quiet_rows=missing,
                    quiet_union_log2_upper=str(log_union),all_quiet_rows_certified=not missing and bool(log_union < -40),
                    directed=True,whole_code_certificate=False)


if __name__=='__main__':
    import json
    print(json.dumps(receipt()),flush=True)
