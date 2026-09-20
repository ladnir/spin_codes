"""Outward positive-occupancy IMT inequalities, NOT a full asymptotic theorem.

The IMT-specific module name avoids legacy certificate-module collisions.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
import math
from pathlib import Path
import sys

from flint import arb, ctx
import numpy as np
import screen


def number(x):
    x = F(x)
    return arb(x.numerator)/x.denominator


def upper(x,bits=128):
    m,e = map(int,x.upper().man_exp())
    exact = F(m)*F(2)**e
    scaled = exact*(1<<bits)
    return F(-(-scaled.numerator//scaled.denominator),1<<bits)


def kl(x,y):
    if x==0:
        return -(1-number(y)).log()
    if x==1:
        return -number(y).log()
    a,b = number(x),number(y)
    return a*(a/b).log()+(1-a)*((1-a)/(1-b)).log()


def check_geometry(leaves,segments):
    assert {b['segment'] for b in leaves}==set(range(len(segments)))
    for segment in segments:
        group = [b for b in leaves if b['segment']==segment.index]
        area = F(0)
        for i,b in enumerate(group):
            a,c = map(F,b['alpha']);x,y = map(F,b['row_density'])
            assert F(1,10000)<=a<c<=1
            assert segment.lower<=x<y<=segment.upper
            area += (c-a)*(y-x)
            for other in group[:i]:
                d,e = map(F,other['alpha']);z,w = map(F,other['row_density'])
                assert max(a,d)>=min(c,e) or max(x,z)>=min(y,w),'overlapping boxes'
        assert area==F(9999,10000)*(segment.upper-segment.lower),'coverage gap'


class Checker:
    def __init__(self):
        self.model = screen.Model()
        self.engine = self.model.engine

    def propose(self,w):
        p,y = F.from_float(w['p']),F.from_float(w['y'])
        ell = F.from_float(math.log(w['lam']))
        family = w['family']
        if family=='statefree':
            vector = [F(1)]
        else:
            beta,lam = float(p*y),math.exp(float(ell))
            matrix = self.model.occupation(beta,lam) if family=='occupation' else self.model.fourier(beta,lam)
            values,vectors = np.linalg.eig(matrix)
            i = int(np.argmax(values.real))
            assert abs(values[i].imag)<1e-12
            v = np.abs(vectors[:,i].real)
            v /= v[0]
            assert np.all(v>0)
            vector = [F(format(float(a),'.17g')) for a in v]
        return dict(p=str(p),y=str(y),log_lam=str(ell),family=family,vector=list(map(str,vector)))

    @lru_cache(maxsize=2000)
    def rate(self,p,y,ell,family,vector):
        p,y,ell = F(p),F(y),F(ell)
        theta = number(p*y)
        lam = number(ell).exp()
        z = (-lam).exp()
        if family=='statefree':
            g0,g1 = 1-theta+theta*z,theta+(1-theta)*z
            radius = max(upper(g0**(128-w)*g1**w) for w in [0,*self.engine.levels])
        else:
            n = self.engine.n
            if family=='occupation':
                epochs = self.engine.epoch(ell)
                probabilities = [math.comb(128,j)*theta**j*(1-theta)**(128-j) for j in range(129)]
                matrix = [sum((probabilities[j]*epochs[j][k] for j in range(129)),arb(0)) for k in range(n*n)]
            else:
                assert family=='fourier'
                matrix = self.engine.bernoulli(theta,lam)
            vec = list(map(number,vector))
            assert len(vec)==n and all(v>0 for v in vec)
            radius = max(upper(sum((matrix[i*n+j]*vec[j] for j in range(n)),arb(0))/vec[i]) for i in range(n))
        assert radius>0
        return number(radius).log()/128,lam

    def box(self,box,segments,delta):
        w = box['rational_witness']
        rate,lam = self.rate(w['p'],w['y'],w['log_lam'],w['family'],tuple(w['vector']))
        p,y = F(w['p']),F(w['y'])
        assert 0<p<1 and 0<y<1
        segment = segments[box['segment']]
        values = []
        for a in map(F,box['alpha']):
            for x in map(F,box['row_density']):
                val = number(a*(segment.slope*x+segment.intercept))+kl(a,p)+number(a)*kl(x,y)+rate+number(delta)*lam
                values.append(upper(val))
        return max(values)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cover',type=Path,default=screen.HERE/'COVER_D1099.json')
    p.add_argument('--verify',type=Path)
    p.add_argument('--output',type=Path,required=True)
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    ctx.prec = 512 if args.verify else 256
    majorant = screen.FROZEN/'golay_ba3_concave_majorant.json'
    segments = screen.load_segments(majorant)
    saved = json.loads(args.verify.read_text()) if args.verify else None
    if saved:
        for name,digest in saved['source_sha256'].items():
            assert screen.sha(screen.ROOT/name)==digest,name
        leaves = saved['leaves']
        delta = F(saved['delta'])
    else:
        cover = json.loads(args.cover.read_text())
        assert cover['unresolved']==0
        leaves = cover['leaves']
        delta = F(str(cover['delta']))
    check_geometry(leaves,segments)
    checker = Checker()
    if saved:
        assert saved['inner']==checker.engine.identity()['inner']
    checked = []
    for i,b in enumerate(leaves):
        box = {k:b[k] for k in ('segment','alpha','row_density')}
        box['rational_witness'] = b['rational_witness'] if saved else checker.propose(b['witness'])
        value = checker.box(box,segments,delta)
        assert value<0,(i,float(value))
        if saved:
            assert value<=F(b['exponent_upper']),(i,'replay exceeds retained upper')
        checked.append(dict(**box,exponent_upper=str(value)))
        if (i+1)%100==0:
            print(ctx.prec,'checked',i+1,'/',len(leaves),flush=True)
    paths = [Path(__file__),Path(screen.__file__),majorant,screen.HERE.parent/'NO_CONSTANT_MAP.json']
    paths += [Path(module.__file__).resolve() for module in list(sys.modules.values())
              if getattr(module,'__file__',None) and Path(module.__file__).resolve().is_relative_to(screen.ROOT)
              and Path(module.__file__).suffix=='.py']
    paths += [(args.verify if saved else args.cover).resolve()]
    result = dict(status='OUTWARD_IMT_POSITIVE_OCCUPANCY_INEQUALITIES' if not saved else 'HIGHER_PRECISION_IMT_POSITIVE_OCCUPANCY_REPLAY_PASSED',
                  precision_bits=ctx.prec,delta=str(delta),alpha_range=['1/10000','1'],row_density_range=['13/125','112/125'],
                  inner=checker.engine.identity()['inner'],exact_geometry_checked=True,boxes=len(checked),
                  maximum_exponent_upper=max(float(F(b['exponent_upper'])) for b in checked),leaves=checked,
                  source_sha256={q.relative_to(screen.ROOT).as_posix():screen.sha(q) for q in set(paths)},
                  full_asymptotic_theorem=False,
                  scope='Positive occupancy only, conditional on the existing Golay--BA-3 concave spectrum majorant and routing reduction.')
    with args.output.open('x',encoding='utf-8') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('leaves','source_sha256','inner')}))


if __name__=='__main__':
    main()
