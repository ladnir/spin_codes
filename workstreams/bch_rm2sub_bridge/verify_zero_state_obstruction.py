"""Exact simplified obstruction ledger, after the 512-bit lower replay."""
import math
from fractions import Fraction as F
import bridge as base
from replay_zero_state_lower import replay


def verify():
    replay()
    conditional=base.read(base.HERE/'generated/zero_state_q2620_conditional_lower.json')
    tail=base.read(base.HERE/'generated/joint_tail_lower_80/lower.json')
    q=2620;assert conditional['occupation']==q
    positions=math.comb(8192,q)
    assert base.decode(conditional['factor_lower'])>=F(positions,1<<245760)
    lower=F(positions*tail['lower']**q,1<<245760)
    assert lower>F(1<<19600)
    print('Simplified exact lower: binom(8192,2620) * tail_lower^2620 / 2^245760 >2^19600.')
    print('Diagnostic log2 lower:',math.log2(lower.numerator)-math.log2(lower.denominator))
    print('This rules out a small unconditional first-moment bound, NOT a small failure probability by itself.')


if __name__=='__main__':verify()
