"""Exact Boolean-degree description of the residual F7 cosets in P dual."""
import argparse
from pathlib import Path
import bridge as base
import pdual_low_weight as low
from affine_wambach import gf_pow


def anf(word,coords):
    a=[0]*256
    for i,x in enumerate(coords):a[x]=(word>>i)&1
    for j in range(8):
        for x in range(256):
            if x>>j&1:a[x]^=a[x^(1<<j)]
    return sum(v<<x for x,v in enumerate(a))


def binary_rank(rows):
    pivots={}
    for v in rows:
        while v:
            j=v.bit_length()-1
            if j not in pivots:pivots[j]=v;break
            v^=pivots[j]
    return len(pivots)


def build():
    checks,rows,values=low.data();coords=[gf_pow(2,i) for i in range(255)]+[0]
    tops=[];pivots={};records=[]
    for word in rows:
        poly=anf(word,coords)
        assert all(not(poly>>m&1) for m in range(256) if m.bit_count()>5)
        top=sum(1<<m for m in range(256) if m.bit_count()==5 and poly>>m&1)
        f=low.syndrome(word,values);tops.append(top)
        records.append(dict(word_hex=hex(word),F7=f,quintic_mask_hex=hex(top)))
        while f:
            j=f.bit_length()-1
            if j not in pivots:pivots[j]=(f,top);break
            f^=pivots[j][0];top^=pivots[j][1]
        if not f:assert top==0
    assert len(pivots)==binary_rank(tops)==8
    def top_of(f):
        out=0
        while f:
            j=f.bit_length()-1;v,top=pivots[j];f^=v;out^=top
        return out
    masks=[top_of(f) for f in range(256)]
    assert len(set(masks))==256 and masks[0]==0
    # Independently check the normalized field-scaling permutation on all basis
    # functions: it acts transitively on nonzero F7, hence on nonzero top parts.
    assert {gf_pow(a,7) for a in range(1,256)}==set(range(1,256))
    return dict(status='EXACT_QUINTIC_STRUCTURE_NOT_A_NEW_DISTANCE_BOUND',
        length=256,Pdual_dimension=125,maximum_Boolean_degree=5,quintic_image_dimension=8,
        Pdual_intersection_RM4_equals_F7_kernel=True,kernel_dimension=117,
        normalized_F7_one_quintic_mask_hex=hex(masks[1]),
        normalized_F7_one_quintic_monomial_masks=[m for m in range(256) if masks[1]>>m&1],
        basis_checks=records,
        meaning='F7=1 fixes the entire degree-five part; allowed lower-degree perturbations form the specific 117-dimensional kernel, not all RM(4,8).',
        source_sha256={p.relative_to(base.ROOT).as_posix():base.sha(p) for p in
            (Path(__file__),Path(low.__file__),base.HERE/'pdual_rank_refinement.py',
             base.BCH/'code/certify_bch_shift_rank.py',base.BCH/'code/affine_wambach.py',base.BCH/'code/bch_quotient.py')})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--verify',action='store_true');a=p.parse_args();result=build()
    if a.verify:assert result==base.read(a.output)
    else:base.write_new(a.output,result)
    print({k:v for k,v in result.items() if k not in ('source_sha256','basis_checks')},flush=True)
