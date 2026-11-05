"""Directed burst test for closed paid sets; not a code-distance bound.

The canonical support covers every guaranteed crossing-burst input, so
choosing paid positions only there gives a subset of actual closed sets.
"""
from fractions import Fraction as F
import math
from quartet_quota_audit import choose_log,linear_log


def audit():
    lam=F(707,2560);tau=F(381,512);B=F(165556);rowset=choose_log(8192,2100)
    for half in range(6,13):
        v=512*half;A=(2*half-1)*512*128;D=512*128;score=lam*A-tau*D
        closed_terms=[(1,choose_log(256,44))];open_terms=list(closed_terms)
        for c,multiplicity in ((435,44),(436,212)):
            shared=[(multiplicity,choose_log(8192-v,2100-c)),(-multiplicity,choose_log(8192,2100))]
            closed_terms.extend(shared+[(multiplicity,choose_log(v,c))])
            open_terms.extend(shared+[(multiplicity,choose_log(6092+c,c))])
        closed=linear_log(closed_terms);opened=linear_log(open_terms)
        reduction=linear_log(((1,opened),(-1,closed)))
        current=linear_log(((1,closed),(1,rowset)));current=tuple(x+4*(score-B) for x in current)
        threshold=None
        if score>B:
            threshold=(-40-rowset[0]-closed[0])/(score-B)
        print(dict(full_updates_per_half=half,guaranteed_activity=A,score=str(score),
                   directed_log2_closed_event_lower=str(closed[0]),closed_event_log2=float(closed[0]),
                   unused_payment_log2_reduction_lower=float(reduction[0]),
                   directed_b16_expression_lower=str(current[0]),b16_expression_lower=float(current[0]),
                   necessary_log2_b_upper=None if threshold is None else float(threshold),
                   positive_is_moment_obstruction=current[0]>0,not_failure_lower=True),flush=True)


if __name__=='__main__':audit()
