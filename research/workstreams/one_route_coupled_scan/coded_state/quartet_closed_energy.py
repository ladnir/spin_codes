"""Exact telescoping score for the seven-type closed-set transfer.

Integer identities and bounded FLOAT location tests. Construction fixed;
no new code-distance claim. Serial Peach90s, <=128MiB, no generated files.
"""
from fractions import Fraction as F
import itertools
import math
from quartet_closed_quotient import kernel,coefficient_table
from quartet_closed_joint import bound


WIDTH=512
UPDATES=4096
THRESHOLD=81074176


def energy(population,regional_inactive=2100):
    """Numerator of the output-state score, in units1/2560."""
    z,v0,v1,blocked,l0,l1,l2=population
    return 707*(l0+l1+l2)-1905*l0+3072*z-192*regional_inactive


def tests():
    checked=0
    for old in itertools.combinations_with_replacement(range(7),4):
        old_selected=sum(t>=4 for t in old)
        for (new,u,v,a,d),ways in kernel(old):
            for paid in range(v+1):
                original=u+paid
                left=707*a-1905*d-3072*paid
                right=707*a-1905*new.count(4)+3072*u-3072*original+1905*(a-old_selected)
                assert left==right;checked+=1
    assert F(THRESHOLD,2560)==F(713983,5)-111127
    # The free terminal contribution is retained, not discarded.
    assert 1905*WIDTH<THRESHOLD
    # Every member of this low-score25-state set has negative score.
    center=(104,176,72,32,16,16,96)
    scores=[energy((104,176-a-c,72,32,16+a,16+c,96))
            for a,c in itertools.product(range(-2,3),repeat=2)]
    assert max(scores)==-19886
    print(dict(exact_score_identities=checked,global_threshold_numerator=THRESHOLD,
               maximum_terminal_credit=1905*WIDTH,low_box_maximum_energy=max(scores)),flush=True)


def diagnostic():
    table=coefficient_table()
    profiles=((104,120,64,32,16,16,160),(104,88,48,48,24,24,176),
              (104,56,32,64,32,32,192),(104,40,24,56,40,40,208),
              (104,20,12,56,48,48,224),(104,16,8,24,16,16,328))
    for pop in profiles:
        assert sum(pop)==WIDTH and pop[4]==pop[5]
        assert energy(pop)*UPDATES>THRESHOLD
        values=[bound(table,pop,pop,i) for i in (131,132)]
        regional=12*values[0][0]+4*values[1][0]+math.log2(math.comb(16,4))-math.log2(math.comb(8192,2100))
        print(dict(population=pop,score_per_step=str(F(energy(pop),2560)),
                   conditional_region_log2=float(regional),maximum_fit_residual=max(v[2] for v in values),
                   full_certificate=False),flush=True)


if __name__=='__main__':
    tests();diagnostic()
