"""Actual-setup mixed-minor exception bound; optional FLOAT diagnostic.

No rejection rule or construction changes. Old frozen helper files are only
imported. The optional process-local callback substitutes the stronger local
count; its finite screens exclude the reported exception overhead and are
not certificates. Serial, file-free, <=60 seconds, >=20 GiB free.
"""
import argparse
from fractions import Fraction as F
import itertools
import math
import os
from pathlib import Path
import shutil
import threading

import quartet_two_update as q
import quartet_joint_quota as joint

OLD_CHOICES=q.choices
TAIL_UPPER=F(-53488800536375,1099511627776)
UNION_OVERHEAD=F(28760690499381,274877906944)


def mixed_choices(g,before,after,target,strict):
    if not strict:return OLD_CHOICES(g,before,after,target,strict)
    active=before.bit_count();overlap=before&after
    if not active:return OLD_CHOICES(g,before,after,target,strict)
    counts=[0]*(g+1)
    for zeros in range(1<<g):
        if zeros&~overlap:continue
        # Input-zero rows outside the next support use A. Input-zero rows
        # at selected output zeros use A+delta I. Keep their union.
        if (target&(~after|zeros)).bit_count()>=active:continue
        counts[zeros.bit_count()]+=1
    return tuple(counts)


def exception_bound():
    q.THETA=F(7,40);q.ROOTS.clear()
    count=degree=0
    for k in range(1,5):
        subsets=[set(s) for s in itertools.combinations(range(4),k)]
        for rows,columns in itertools.product(subsets,repeat=2):
            common=len(rows&columns)
            # Empty and full diagonal shifts are already accepted. With
            # no intersection they coincide and leave no additional case.
            extra=max(0,2**common-2)
            count+=extra;degree+=k*extra
    assert (count,degree)==(74,224)
    field=1<<32;positions=1<<19
    probability=F(degree,field-280)
    tail=math.comb(positions,7)*probability**7
    tail_log=q.log_upper(1,tail)
    ell=F(13592,664659)
    overhead=sum(math.comb(positions,j)*q.root_upper((1+ell)**(4*j)) for j in range(7))
    overhead_log=q.log_upper(1,overhead)
    print(dict(additional_minors=count,total_degree=degree,
               maximum_exceptional_quartets=6,tail_log2_upper=str(tail_log),
               tail_float=float(tail_log),union_log2_overhead=str(overhead_log),
               overhead_float=float(overhead_log),directed=True,whole_code=False),flush=True)
    assert tail_log==TAIL_UPPER and overhead_log==UNION_OVERHEAD
    assert tail_log < -48 and overhead_log < 105


def witness_check():
    import inner_mechanisms as field
    alpha=field.verify_base();delta=field.embedded_delta()
    _,multiply,_,_=field.field_ops(alpha)
    key=(0xc423f3e3,0xb01ba9e8,0x5aa87c13,0x915276c9)
    matrix=field.matrix_for(key,multiply)
    assert all(value for _,_,value in field.minors(matrix))
    checks=0
    for k in range(1,5):
        for rows in itertools.combinations(range(4),k):
            for columns in itertools.combinations(range(4),k):
                common=sorted(set(rows)&set(columns))
                for mask in range(1,(1<<len(common))-1):
                    shifted={common[i] for i in range(len(common)) if mask>>i&1}
                    minor=[[matrix[i][j]^(delta if i==j and i in shifted else 0)
                            for j in columns] for i in rows]
                    assert field.determinant(minor)!=0;checks+=1
    assert checks==74
    print(dict(actual_field_additional_minor_checks=checks,frozen_key_unchanged=True),flush=True)


def tests():
    comparisons=removed=0
    for before,after,target in itertools.product(range(16),repeat=3):
        old=OLD_CHOICES(4,before,after,target,True)
        new=mixed_choices(4,before,after,target,True)
        assert all(x<=y for x,y in zip(new,old))
        assert new[0]==old[0]
        removed+=sum(old)-sum(new);comparisons+=1
    # For arbitrary ell>0, identical feasibility and constant coefficient
    # one imply P_old <= (1+ell)^4 P_mixed. Check exact representative values.
    ratios=0
    for ell in (F(1,64),F(1,4),F(1),F(3)):
        for a,after,target in itertools.product(range(5),range(16),range(16)):
            before=(1<<a)-1
            old=OLD_CHOICES(4,before,after,target,True)
            new=mixed_choices(4,before,after,target,True)
            left=sum(c*ell**j for j,c in enumerate(old))
            right=sum(c*ell**j for j,c in enumerate(new))
            assert left<=(1+ell)**4*right;ratios+=1
    print(dict(local_comparisons=comparisons,removed_selections=removed,
               exact_polynomial_ratio_checks=ratios),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(60,lambda:os._exit(124));timer.daemon=True;timer.start()
    try:
        parser=argparse.ArgumentParser();parser.add_argument('--screen',action='store_true')
        parser.add_argument('--tau',type=float,default=.1)
        parser.add_argument('--x',type=float,default=97905/432928)
        args=parser.parse_args();tests();witness_check();exception_bound()
        if args.screen:
            q.choices=mixed_choices
            joint.run(args.tau,112627/167924,args.x)
    finally:
        q.choices=OLD_CHOICES;timer.cancel()


if __name__=='__main__':main()
