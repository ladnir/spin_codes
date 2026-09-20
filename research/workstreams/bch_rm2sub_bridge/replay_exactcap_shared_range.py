"""Read-only replay of the frozen first improved-cap shared receipt.

The producer's embedded verifier compares a nonzero-weight snapshot against
a dictionary also containing weight zero. This separate verifier compares
the same nonzero domains explicitly; it changes neither source nor receipt.
"""
import math
import numpy as np
from flint import arb,ctx
import bridge as base
import christoffel_caps as cap_base
import activation_bridge as q1
import tightened_occupancy as tight
import polynomial_regions as poly
import scaled_adaptive as scaled
import positive_line_hull as hull
from general_batch_certificate import density_cost


def replay():
    saved=base.read(base.HERE/'generated/exactcap_shared_q1025_q1280_outward.json')
    assert saved['status']=='OUTWARD_EXACTCAP_SHARED_RANGE'
    assert saved['configuration']=='t128_s15' and saved['occupancy_range']==[1025,1280]
    for name,digest in saved['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    caps=cap_base.deterministic_caps()
    for name in saved['cap_receipts']:
        receipt=base.read(base.HERE/name)
        for p,digest in receipt['local_sha256'].items():assert base.sha(base.HERE/p)==digest
        for p,digest in receipt['outer_sha256'].items():assert base.sha(base.BCH/p)==digest
        assert receipt['rational_primal_dual_checks_passed']
        w=receipt['weight'];caps[w]=caps[256-w]=min(caps[w],receipt['cap'])
    assert set(map(int,saved['used_caps']))==set(base.WEIGHTS)
    assert all(caps[w]==saved['used_caps'][str(w)] for w in base.WEIGHTS)
    screen=base.read(base.HERE/'generated/independent_p_anchors_screen.json')
    for name,digest in screen['source_sha256'].items():assert base.sha(base.HERE/name)==digest
    bands=screen['bands'];assert sorted(w for band in bands for w in band)==list(base.WEIGHTS)
    witness=next(r for r in screen['rows'] if r['occupation']==1024)
    t,s,spectrum=base.load_map(tight.NAME);ctx.prec=512
    assert saved['parameters']==dict(message_bits=1<<20,output_bits=1<<21,outer_rows=8192,outer_length=256,
        step_bits=t,state_bits=s,distance_cutoff=209716)
    lam=(arb(witness['witness_tenth'])/10).exp();correction=(209716*lam).exp()
    region=poly.regions(t,s,spectrum,tight.kernel_spectrum(),(-lam).exp(),1280)
    ps=[base.decode(v) for v in witness['p']];assert len(ps)==len(bands) and all(0<p<1 for p in ps)
    costs=[density_cost(band,p,caps) for band,p in zip(bands,ps)]
    roots=[((arb(c.numerator)/c.denominator).log()/256).exp().upper() for c in costs]
    probabilities=[arb(p.numerator)/p.denominator for p in ps]
    left=np.nextafter(np.array([float((r*(1-p)).upper()) for r,p in zip(roots,probabilities)]),np.inf)
    right=np.nextafter(np.array([float((r*p).upper()) for r,p in zip(roots,probabilities)]),np.inf)
    keep=hull.indices(left,right);mantissas,exponents=scaled.initial(region)
    from audit_bch_q1_full_arb import rational
    from fractions import Fraction as F
    total=F(0);assert len(saved['rows'])==256
    for q,matrix,exponent in scaled.matrices(mantissas,exponents,left[keep],right[keep]):
        if q<1025:continue
        value=tuple(arb(float(v)) for v in matrix.flat)
        for _ in range(8):value=q1.positive_mul(value,value)
        bound=math.comb(8192,q)*len(bands)**q*rational((sum(value[:3],arb(0))*correction*arb(2)**(256*exponent)).upper())
        stored=saved['rows'][q-1025];assert stored['occupation']==q
        assert 0<bound<=base.decode(stored['upper']);total+=bound
    assert total<=base.decode(saved['range_upper'])
    print('512-bit replay passed all 256 rows, including retained nonpassing rows.')


if __name__=='__main__':replay()
