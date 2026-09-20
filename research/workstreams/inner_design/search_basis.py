"""Search reverse-state bases, preserving the exact selected transposed inner.

For z=V r, emission is A V^-1, feedback is V A^T, transition is V M^T V^-1.
The cost proxy jointly scores emission-table loads and feedback XOR density.
It is a screening score, never a performance claim.
"""
import argparse
import hashlib
import json
import random
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import probe_bch_forward_xor_circuit as paar
S=19


def identity(): return [1<<j for j in range(S)]


def apply(rows, word):
    return sum(((r & word).bit_count() & 1)<<j for j,r in enumerate(rows))


def compose(a,b):
    def row(mask):
        value=0
        while mask:
            bit=mask & -mask; mask^=bit
            value^=b[bit.bit_length()-1]
        return value
    return [row(r) for r in a]


def inverse(rows):
    a=list(rows);b=identity()
    for j in range(S):
        pivot=next(i for i in range(j,S) if a[i]>>j&1)
        a[j],a[pivot]=a[pivot],a[j];b[j],b[pivot]=b[pivot],b[j]
        for i in range(S):
            if i!=j and a[i]>>j&1: a[i]^=a[j];b[i]^=b[j]
    assert a==identity()
    return b


def transpose(rows):
    return [sum(((r>>j)&1)<<i for i,r in enumerate(rows)) for j in range(S)]


def load():
    path=ROOT/'workstreams/bare_bch_rm2sub/generated/MANIFEST.json'
    record=json.loads(path.read_text())['t128_s19']
    columns=record['columns']
    emission=[sum(((c>>j)&1)<<x for x,c in enumerate(columns)) for j in range(S)]
    coefficients=columns.copy()
    for b in range(7):
        for x in range(128):
            if x>>b&1: coefficients[x]^=coefficients[x^(1<<b)]
    monomials=[x for x in range(128) if x.bit_count()<=2]
    feedback=[sum(((coefficients[x]>>j)&1)<<i for i,x in enumerate(monomials)) for j in range(S)]
    assert all(c==0 for x,c in enumerate(coefficients) if x.bit_count()>2)
    return record,emission,feedback,monomials


def emission_cost(rows):
    return sum(group_cost(rows,g) for g in range(5))


def group_cost(rows,g):
    union=0
    for r in rows[4*g:4*g+4]: union|=r
    return union.bit_count()


def combine(rows, coefficients):
    result=[]
    for row in rows:
        value=0
        for j in range(S):
            if row>>j&1:value^=coefficients[j]
        result.append(value)
    return result


def anneal(start,original_a,original_b,weight,seed,steps):
    rng=random.Random(seed)
    v=start.copy()
    a=combine(transpose(inverse(v)),original_a)
    b=combine(v,original_b)
    cost=emission_cost(a)+weight*sum(r.bit_count()-1 for r in b)
    best=(cost,v.copy())
    for step in range(steps):
        i,j=rng.sample(range(S),2)
        # z_i <- z_i + z_j changes feedback row i and emission row j.
        old_group=group_cost(a,j//4)
        next_a=a[j]^a[i];next_b=b[i]^b[j]
        delta_b=weight*(next_b.bit_count()-b[i].bit_count())
        old_a=a[j];a[j]=next_a
        delta=group_cost(a,j//4)-old_group+delta_b
        temperature=1.25*(1-step/steps)**2
        accept=delta<=0 or (temperature>0 and rng.random()<2.718281828459045**(-delta/temperature))
        if accept:
            b[i]=next_b;v[i]^=v[j];cost+=delta
            if cost<best[0]:best=(cost,v.copy())
        else:a[j]=old_a
    return best[1]


def reduced_basis(rows):
    a=rows.copy();v=identity();pivot=0
    for column in range(max(r.bit_length() for r in rows)):
        found=next((i for i in range(pivot,S) if a[i]>>column&1),None)
        if found is None:continue
        a[pivot],a[found]=a[found],a[pivot];v[pivot],v[found]=v[found],v[pivot]
        for i in range(S):
            if i!=pivot and a[i]>>column&1:a[i]^=a[pivot];v[i]^=v[pivot]
        pivot+=1
        if pivot==S:break
    assert pivot==S
    return v


def audit(v,original_a,original_b):
    inv=inverse(v)
    assert compose(v,inv)==identity() and compose(inv,v)==identity()
    a=combine(transpose(inv),original_a);b=combine(v,original_b)
    assert combine(transpose(v),a)==original_a
    assert combine(inv,b)==original_b
    # BA=0 is preserved, including nonsymmetric coordinate transformations.
    b_words=combine(v,original_a)
    assert all((r&s).bit_count()%2==0 for r in b_words for s in a)
    return a,b


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--steps',type=int,default=30000)
    parser.add_argument('--restarts',type=int,default=6)
    parser.add_argument('--seed-offset',type=int,default=0)
    parser.add_argument('--mixed-starts',action='store_true')
    parser.add_argument('--output',default='BASIS_SEARCH.json');args=parser.parse_args()
    record,a,b,monomials=load()
    group=[1<<j for j in record['A_group_order']]
    bases={'identity':identity(),'grouped':group,'feedback_rref':reduced_basis(b),
           'emission_rref':transpose(inverse(reduced_basis(a)))}
    for weight in (0,0.25,0.5,1,2):
        for seed in range(args.seed_offset,args.seed_offset+args.restarts):
            starts=[group,bases['feedback_rref'],bases['emission_rref']]
            start=starts[seed%3] if args.mixed_starts else group
            bases[f'joint_w{weight}_r{seed}']=anneal(start,a,b,weight,seed,args.steps)
    records=[];seen=set();paar.DIMENSION=len(monomials)
    for name,v in bases.items():
        if tuple(v) in seen:continue
        seen.add(tuple(v));aa,bb=audit(v,a,b)
        circuit=paar.synthesize(0,2,bb)
        records.append(dict(name=name,V=v,V_inverse=inverse(v),
                            emission_lookups=emission_cost(aa),feedback_xors=circuit.xor_count,
                            feedback_dense_xors=sum(r.bit_count()-1 for r in bb),
                            joint_proxy=emission_cost(aa)+circuit.xor_count,
                            emission_rows_hex=[hex(r) for r in aa],feedback_monomials_hex=[hex(r) for r in bb]))
    records.sort(key=lambda r:r['joint_proxy'])
    output=dict(parameters=dict(t=128,s=19,steps=args.steps,restarts=args.restarts,
                                seed_offset=args.seed_offset,mixed_starts=args.mixed_starts),
                convention='reverse state z=V r; emission A V^-1; feedback V A^T; transition V M^T V^-1',
                monomials=monomials,selection_sha256=record['selection_sha256'],
                producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),candidates=records)
    (HERE/args.output).write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps([{k:r[k] for k in ('name','emission_lookups','feedback_xors','joint_proxy')} for r in records],indent=2))


if __name__=='__main__':main()
