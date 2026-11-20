"""Reuse a proved Fourier numerator for an entire integer type box.

Only exact factorial/probability ratios change between types. No endpoint
interpolation or new Fourier sampling is used.
"""
import math
from fractions import Fraction as F
from pathlib import Path
from flint import arb,arb_poly,ctx
import bridge as base
import tightened_occupancy as tight


def factorial_ratio(n,m):
    return F(math.prod(range(m+1,n+1))) if n>=m else F(1,math.prod(range(n+1,m+1)))


def run():
    source=base.HERE/'generated/pair_type_central_tracked_outward.json';saved=base.read(source)
    assert saved['status']=='OUTWARD_ONE_SHARED_REGION_PAIR_TYPE' and saved['configuration']==tight.NAME
    for name,digest in saved['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    original=saved['region_type'];ps=[base.decode(v) for v in saved['tilt_probabilities']]
    bound=base.decode(saved['pair_region_upper']);assert sum(ps)==1
    ctx.prec=512;kernel=tight.kernel_spectrum();poly=arb_poly([arb(kernel.get(j,0)) for j in range(129)])**64
    from audit_bch_q1_full_arb import rational
    singles={j:rational((poly[j]/math.comb(8192,j)).lower()) for j in range(814,823,2)}
    rows=[]
    for x in range(814,823,2):
        for y in range(814,823,2):
            for overlap in range(80,85):
                m=[8192-x-y+overlap,y-overlap,x-overlap,overlap]
                change=math.prod(factorial_ratio(n,old)*p**(old-n) for n,old,p in zip(m,original,ps))
                ratio=bound*change/(singles[x]*singles[y])
                rows.append(dict(type=m,ratio_upper=base.encode(ratio)))
    maximum=max(base.decode(row['ratio_upper']) for row in rows)
    result=dict(status='OUTWARD_SHARED_REGION_PAIR_TYPE_BOX',configuration=tight.NAME,
        first_weights=list(range(814,823,2)),second_weights=list(range(814,823,2)),overlaps=list(range(80,85)),
        rows=rows,maximum_ratio_upper=base.encode(maximum),all_below_eleven_tenths=maximum<F(11,10),
        full_second_moment_certified=False,local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__),source,Path(tight.__file__)]})
    path=base.HERE/'generated/pair_type_central_box_outward.json'
    if path.exists():assert result==base.read(path)
    else:base.write_new(path,result)
    print('Exact type-box extension:',len(rows),'types; max ratio',float(maximum),'all <1.1:',maximum<F(11,10),flush=True)


if __name__=='__main__':run()
