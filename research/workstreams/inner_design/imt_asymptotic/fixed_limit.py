"""Outward inequalities for the proposed conservative IMT continuum limit.

Conditional local result only: the finite-region-to-continuum reduction is
not proved by this script, and Q=1,2 need separate treatment.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys

from flint import ctx
import screen
import certify_imt_dense as cert
sys.path.insert(0,str(screen.FROZEN))
import certify_golay_ba_rm2sub_weight_coupled_fixed as old


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    ctx.prec = 256
    majorant = screen.FROZEN/'golay_ba3_concave_majorant.json'
    segments = screen.load_segments(majorant)
    sigma,tau,v = F(127,250),F(133,125),F(3,1600)
    reset,p0 = F(1,(1<<19)-1),F(1<<18,(1<<19)-1)
    lines = [(F(1),sigma*v),(sigma,reset/v+sigma)]
    d = (1-sigma)/sigma
    critical = F(4,3)/tau-1/d
    assert 0<critical<1
    num = cert.number
    spacing = max(cert.upper(-num(tau*y)+num(F(4,3))*(1+num(d*y)).log()) for y in (F(0),critical,F(1)))
    support = num(2).log()/num(F(39,4))
    delta_term = num(F(11,100)*tau/p0)
    rows = []
    for segment in segments:
        u = old.propose_fugacity((segment.lower+segment.upper)/2)
        norm = max(a+b*u for a,b in lines)
        for x in (segment.lower,segment.upper):
            entropy = -num(x)*num(x).log()-(1-num(x))*(1-num(x)).log()
            row = num(segment.slope*x+segment.intercept)-entropy-num(x)*num(u).log()+num(norm).log()+num(spacing)+delta_term
            value = cert.upper(row+support)
            assert value<0
            rows.append(dict(segment=segment.index,x=str(x),u=str(u),exponent_upper=str(value)))
    result = dict(status='OUTWARD_CONDITIONAL_IMT_CONTINUUM_Q_GE_3_INEQUALITIES',
                  precision_bits=256,delta='11/100',block_constant='39/4',
                  impulse_upper=[['0','1'],[str(reset),'1']],p0=str(p0),
                  sigma=str(sigma),tau=str(tau),v=str(v),checks=rows,
                  maximum_upper=max(float(F(r['exponent_upper'])) for r in rows),
                  full_asymptotic_theorem=False,
                  missing=['Uniform finite-region-to-continuum comparison for IMT.','Separate Q=1,2 proof.'],
                  source_sha256={q.relative_to(screen.ROOT).as_posix():screen.sha(q) for q in [Path(__file__),Path(cert.__file__),Path(old.__file__),majorant]})
    with args.output.open('x',encoding='utf-8') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('checks','source_sha256')}))


if __name__=='__main__':
    main()
