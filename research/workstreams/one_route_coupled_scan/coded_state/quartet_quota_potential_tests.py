"""Exact toy checks of a regional boundary potential and initial row."""
from collections import defaultdict
from fractions import Fraction as F
import math
import numpy as np
import quartet_fresh_regional as regional
from quartet_quota_potential import evaluate


def main():
    P=[[F(1,8),F(1,16),F(1,8),F(1,8),F(1,16),F(1,8)],
       [F(1,16),F(1,8),F(1,8),F(1,16),F(1,8),F(1,8)],
       [F(1,4),F(1,4),F(0),F(0),F(1,8),F(1,8)]]
    targets=np.array([0,1,0,1,0,1]);labels=np.array([408,408,409,409,410,410])
    coefficients=[]
    for source in range(3):
        current={(0,source):F(1)}
        for _ in range(2):
            after=defaultdict(F)
            for (total,i),value in current.items():
                for e in range(6):after[total+int(labels[e])-408,int(targets[e])]+=value*P[i][e]
            current=dict(after)
        coefficients.append(current)
    qs=np.arange(5);quota=[F(1),F(2),F(3),F(2),F(1)];checks=0;lower_checks=0
    pair=np.array([[float(P[i][e]*P[int(targets[e])][f]) for e in range(6) for f in range(6)] for i in range(3)])
    pair_labels=np.array([labels[e]+labels[f] for e in range(6) for f in range(6)])
    pair_targets=np.tile(targets,6);paired_checks=0
    for w in ((F(1),F(1)),(F(1),F(2)),(F(4),F(1))):
        rows=regional.coefficient_rows(np.array(P,dtype=float),labels,steps=2,
               terminal=np.array([float(w[0]),float(w[1]),1.]),all_totals=True,target_map=targets)
        pair_rows=regional.coefficient_rows(pair,pair_labels,steps=1,center=816,
               terminal=np.array([float(w[0]),float(w[1]),1.]),all_totals=True,target_map=pair_targets)
        assert rows.keys()==pair_rows.keys()
        for q in rows:assert np.array_equal(rows[q],pair_rows[q]);paired_checks+=1
        logrows=np.log2(np.array([rows[q] for q in qs]))+np.log2(np.array(quota,dtype=float))[:,None]
        logw=np.log2(np.array(w,dtype=float))
        for y in (F(1,2),F(1),F(2)):
            R=[[sum(value*quota[q]*y**q for (q,j),value in c.items() if j==destination)
                for destination in range(2)] for c in coefficients]
            output,rho,_,_=evaluate(logrows,qs,logw,math.log2(y),0.)
            exact=[sum(R[i][j]*w[j] for j in range(2)) for i in range(3)]
            assert all(abs(2**output[i]-float(exact[i]))<=1e-12*max(1,float(exact[i])) for i in range(3))
            ratio=max(exact[i]/w[i] for i in range(2))
            assert abs(rho-math.log2(ratio))<1e-12
            terminal=[F(1),F(3)];v=terminal
            for _ in range(2):v=[sum(R[i][j]*v[j] for j in range(2)) for i in range(2)]
            actual=sum(R[2][j]*v[j] for j in range(2))
            bound=exact[2]*ratio**2*max(terminal[i]/w[i] for i in range(2))
            assert actual<=bound
            numeric=2**(output[-1]+2*rho+max(math.log2(terminal[i])-logw[i] for i in range(2)))
            assert abs(numeric-float(bound))<=1e-12*max(1,float(bound));checks+=1
            C={q:[[sum(value*quota[q]*y**q for (b,j),value in c.items() if b==q and j==destination)
                   for destination in range(2)] for c in coefficients] for q in (1,2)}
            initial=sum(C[2][2][j]*w[j] for j in range(2))
            floors={q:min(sum(C[q][i][j]*w[j] for j in range(2))/w[i] for i in range(2)) for q in (1,2)}
            lower=2*initial*floors[1]*floors[2]*min(terminal[i]/w[i] for i in range(2))
            actual=F(0)
            for a,b in ((1,2),(2,1)):
                v=[sum(C[b][i][j]*terminal[j] for j in range(2)) for i in range(2)]
                v=[sum(C[a][i][j]*v[j] for j in range(2)) for i in range(2)]
                actual+=sum(C[2][2][j]*v[j] for j in range(2))
            assert lower<=actual;lower_checks+=1
    # The all-covered-count marker charges each path once, not once for
    # every possible count. Compare direct coefficient summation and the
    # regional generating-function product at b=4,beta=3/2,y=1/8.
    current={(0,2):F(1)}
    for _ in range(3):
        after=defaultdict(F)
        for (total,i),weight in current.items():
            for (q,j),value in coefficients[i].items():after[total+q,j]+=weight*value*quota[q]
        current=dict(after)
    y=F(1,8);terminal=[F(1),F(3)]
    direct=sum(value*y**q*terminal[j] for (q,j),value in current.items())
    R=[[sum(value*quota[q]*y**q for (q,j),value in c.items() if j==destination)
        for destination in range(2)] for c in coefficients]
    v=terminal
    for _ in range(2):v=[sum(R[i][j]*v[j] for j in range(2)) for i in range(2)]
    assert direct==sum(R[2][j]*v[j] for j in range(2))
    print(dict(exact_three_region_boundary_bounds=checks,numeric_initial_and_terminal_replays=checks,
               exact_balanced_coefficient_lower_bounds=lower_checks,exact_all_counts_marker_identity=True,
               exact_blocked_counter_rows=paired_checks,non_eigenvectors=True,files_written=0),flush=True)


if __name__=='__main__':main()
