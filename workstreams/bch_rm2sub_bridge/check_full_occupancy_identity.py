"""Independent all-row Bernoulli identity without region coefficient extraction.

At Q=L with a common p, each region is iid Bernoulli(p), so its moment
matrix is the 64th power of the binomially averaged epoch matrix.
This is a pure-group diagnostic, not a proof over all group mixtures.
"""
import argparse
import math
from pathlib import Path
from flint import arb,ctx
import bridge as base
import tightened_occupancy as tight
import activation_bridge as q1
import christoffel_caps as outer
import screen_exponential_modes as latest
from general_batch_certificate import density_cost


def run(screen_name,tag):
    screen=base.read(base.HERE/'generated'/screen_name);ctx.prec=512
    t,s,spectrum=base.load_map(tight.NAME);kernel=tight.kernel_spectrum()
    old=outer.deterministic_caps();new,sources=latest.latest_caps();rows=[]
    for witness in screen['rows']:
        assert witness['occupation']==8192
        lam=(arb(witness['witness_tenth'])/10).exp()
        epoch=tight.epoch_matrices(t,s,spectrum,kernel,(-lam).exp(),arb,128)
        for band,encoded in zip(screen['bands'],witness['p']):
            p=base.decode(encoded);prob=arb(p.numerator)/p.denominator
            weights=[math.comb(t,j)*prob**j*(1-prob)**(t-j) for j in range(t+1)]
            averaged=tuple(sum((mass*row[k] for mass,row in zip(weights,epoch)),arb(0)) for k in range(9))
            value=averaged
            for _ in range(14):value=q1.positive_mul(value,value)
            moment=sum(value[:3],arb(0))
            margins=[]
            for caps in (old,new):
                gamma=density_cost(band,p,caps)
                logarithm=moment.log()+8192*(arb(gamma.numerator)/gamma.denominator).log()+209716*lam
                margins.append(-float(logarithm/arb(2).log()))
            row=dict(band=band,p=encoded,old_cap_margin_bits=margins[0],new_cap_margin_bits=margins[1],
                     epoch_identity_log_moment_diagnostic=float(moment.log()))
            rows.append(row)
            print('Direct iid band',band[0],band[-1],'margins old,new',margins,flush=True)
    base.write_new(base.HERE/'generated'/f'full_iid_identity_{tag}.json',dict(status='PURE_GROUP_IDENTITY_DIAGNOSTIC',rows=rows,
        source_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in
          [Path(__file__),Path(tight.__file__),Path(latest.__file__),Path(outer.__file__),base.HERE/'generated'/screen_name]+sources}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--screen',required=True);p.add_argument('--tag',required=True)
    a=p.parse_args();run(a.screen,a.tag)
