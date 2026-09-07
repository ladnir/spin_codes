"""Bounded SAT counterexample search, with exact independent membership checks.

UNSAT/UNKNOWN is not exported as a distance certificate. A SAT word must pass
direct binary orthogonality, weight, and field-syndrome checks.
"""
import argparse
from pathlib import Path
import time
import z3
import bridge as base
import pdual_rank_refinement as rank
from affine_wambach import gf_pow


def problem():
    data=rank.code_data(37)
    g=int(data['generator_hex'],16)
    checks=[g<<i for i in range(data['dimension'])]
    values=[gf_pow(2,7*i) for i in range(255)]
    return data,checks,values


def verify_word(word):
    data,checks,values=problem()
    assert 0<word<1<<255 and word.bit_count()==30
    assert all((word&row).bit_count()%2==0 for row in checks)
    syndrome=0
    for i,value in enumerate(values):
        if word>>i&1:syndrome ^= value
    assert syndrome==1
    return dict(weight=30,F7=syndrome,punctured_Pdual_membership_checked=True,
                extended_word_hex=hex(word),Pdual_distance32_refuted=True)


def xor_tree(values):
    if not values:return z3.BoolVal(False)
    while len(values)>1:
        values=[z3.Xor(values[i],values[i+1]) if i+1<len(values) else values[i]
                for i in range(0,len(values),2)]
    return values[0]


def run(args):
    assert not args.output.exists()
    data,checks,values=problem()
    bits=[z3.Bool(f'c_{i}') for i in range(255)]
    solver=z3.Solver();solver.set(timeout=args.seconds*1000,random_seed=args.seed)
    for row in checks:
        solver.add(z3.Not(xor_tree([bits[i] for i in range(255) if row>>i&1])))
    for j in range(8):
        solver.add(xor_tree([bits[i] for i in range(255) if values[i]>>j&1])==bool(1>>j&1))
    solver.add(z3.PbEq([(v,1) for v in bits],30))
    begun=time.monotonic();status=solver.check();elapsed=time.monotonic()-begun
    result=dict(status=str(status),elapsed_seconds=elapsed,z3_version=z3.get_version_string(),
        timeout_seconds=args.seconds,seed=args.seed,proof_of_absence=False,
        normalization='F7=1; cyclic coordinate scaling is transitive on its nonzero values because gcd(7,255)=1.',
        source_sha256={p.relative_to(base.ROOT).as_posix():base.sha(p) for p in
            (Path(__file__),Path(rank.__file__),base.BCH/'code/certify_bch_shift_rank.py',
             base.BCH/'code/affine_wambach.py',base.BCH/'code/bch_quotient.py')})
    if status==z3.sat:
        model=solver.model();word=sum(1<<i for i,b in enumerate(bits) if z3.is_true(model.eval(b)))
        result['word_hex']=hex(word);result['exact_witness_audit']=verify_word(word)
    else:
        result['reason_unknown']=solver.reason_unknown() if status==z3.unknown else None
    base.write_new(args.output,result)
    print({k:v for k,v in result.items() if k!='source_sha256'},flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seconds',type=int,default=45)
    p.add_argument('--seed',type=int,default=20260907)
    p.add_argument('--output',type=Path,required=True)
    run(p.parse_args())
