"""512-bit replay with exactly verified coefficient-hull pruning.

The assignment count remains the original band count, not the hull size.
This verifier writes no files and preserves every producer source hash.
"""
import argparse
import math
import numpy as np
from flint import arb,ctx
import bridge as base
import activation_bridge as q1
import tightened_occupancy as tight
import polynomial_regions as poly
import scaled_adaptive as scaled
import positive_line_hull as hull
from general_batch_certificate import density_cost


def replay(screen_name,certificate_name):
    screen_path=base.HERE/'generated'/screen_name;receipt=base.read(base.HERE/'generated'/certificate_name)
    assert receipt['status']=='OUTWARD_LISTED_CONSTANT_SPLIT_CELLS' and receipt['configuration']==tight.NAME
    assert str(screen_path.relative_to(base.HERE)) in receipt['local_sha256']
    for name,digest in receipt['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    screen=base.read(screen_path);bands=screen['bands'];assert sorted(w for band in bands for w in band)==list(base.WEIGHTS[:-1])
    caps={int(w):v for w,v in screen['used_caps'].items()}
    assert len(screen['rows'])==len(receipt['rows'])
    ctx.prec=512;t,s,spectrum=base.load_map(tight.NAME);cached={}
    from audit_bch_q1_full_arb import rational
    for row,stored in zip(screen['rows'],receipt['rows']):
        d=row['ordinary_rows'];h=row['all_one_rows'];tilt=row['witness_tenth']
        assert (d,h)==(stored['ordinary_rows'],stored['all_one_rows']) and 0<=d and 0<=h and 1<=d+h<=8192
        if tilt not in cached:
            lam=(arb(tilt)/10).exp()
            cached[tilt]=(poly.regions(t,s,spectrum,tight.kernel_spectrum(),(-lam).exp(),8192),(209716*lam).exp())
        region,correction=cached[tilt]
        if d:
            ps=[base.decode(v) for v in row['p']];assert len(ps)==len(bands) and all(0<p<1 for p in ps)
            gamma=[density_cost(band,p,caps) for band,p in zip(bands,ps)]
            roots=[((arb(v.numerator)/v.denominator).log()/256).exp().upper() for v in gamma]
            probabilities=[arb(p.numerator)/p.denominator for p in ps]
            left=np.nextafter(np.array([float((r*(1-p)).upper()) for r,p in zip(roots,probabilities)]),np.inf)
            right=np.nextafter(np.array([float((r*p).upper()) for r,p in zip(roots,probabilities)]),np.inf)
            keep=hull.indices(left,right)
            print('Exact coefficient hull retains',len(keep),'of',len(bands),'bands; assignment count unchanged',flush=True)
            mantissas,exponents=scaled.initial(region[h:h+d+1])
            for depth,matrix,exponent in scaled.matrices(mantissas,exponents,left[keep],right[keep]):
                if depth%2048==0:print('Pruned replay depth',depth,flush=True)
            assert depth==d;value=tuple(arb(float(v)) for v in matrix.flat)
        else:value=region[h];exponent=0
        for _ in range(8):value=q1.positive_mul(value,value)
        bound=(math.comb(8192,d)*math.comb(8192-d,h)*len(bands)**d
               *rational((sum(value[:3],arb(0))*correction*arb(2)**(256*exponent)).upper()))
        assert 0<bound<=base.decode(stored['upper'])
        print('512-bit exact-hull replay passed cell',d,h,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--screen',required=True);p.add_argument('--certificate',required=True)
    a=p.parse_args();replay(a.screen,a.certificate)
