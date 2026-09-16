"""Outward Q=1,2 inequalities for the conservative IMT continuum envelope."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path

from flint import arb,ctx
import screen
import certify_imt_dense as dense


def evaluate(precision=256):
    ctx.prec=precision
    num=dense.number
    r=num(F(1,(1<<19)-1))
    p0=num(F(1<<18,(1<<19)-1))
    ln2=num(2).log()
    rows=[]
    for q in (1,2):
        gamma=q*ln2
        theta=gamma/p0
        a=num(F(1,1<<q))
        f=(1-a)/gamma
        gap=2*(gamma-1+a)/gamma**2
        edge=2*(1-(1+gamma)*a)/gamma**2
        k0=[[arb(1),arb(0)],[arb(0),a]]
        k1=[[arb(0),f],[r*f,a]]
        k2=[[r*gap,edge],[r*edge,r*edge+a]]
        matrix=[[k0[i][j]+q*k1[i][j]+(k2[i][j] if q==2 else 0) for j in range(2)] for i in range(2)]
        x,y=matrix[0];z,w=matrix[1]
        radius=(x+w+((x-w)**2+4*y*z).sqrt())/2
        exponent=radius.log()+num(F(11,100))*theta-q*ln2/2+q*num(F(1281,100000))+q*ln2/num(F(39,4))
        upper=dense.upper(exponent)
        assert upper<0
        rows.append(dict(q=q,gamma=f'{q}*ln(2)',exponent_per_outer_bit_upper=str(upper),
                         exponent_per_active_row_bit_upper=str(upper/q)))
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    rows=evaluate(256)
    replay=evaluate(512)
    assert all(F(b['exponent_per_outer_bit_upper'])<=F(a['exponent_per_outer_bit_upper']) for a,b in zip(rows,replay))
    result=dict(status='OUTWARD_IMT_Q1_Q2_CONTINUUM_INEQUALITIES_REPLAYED',
                precision_bits=256,replay_precision_bits=512,delta='11/100',block_constant='39/4',
                outer_likelihood_excess='1281/100000',checks=rows,
                source_sha256={q.relative_to(screen.ROOT).as_posix():screen.sha(q) for q in [Path(__file__),Path(dense.__file__)]},
                full_asymptotic_theorem=False,
                scope='Local continuum inequalities, to be combined with the IMT homogenization lemma and imported outer event.')
    with args.output.open('x',encoding='utf-8') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({**result,'source_sha256':{},'decimal_exponents':[float(F(r['exponent_per_active_row_bit_upper'])) for r in rows]},indent=2))


if __name__=='__main__':main()
