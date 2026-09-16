"""Bounded binary64 box coverage for IMT. No outward certificate is claimed."""
import argparse
from collections import Counter
from fractions import Fraction as F
import json
import math
from pathlib import Path
import time

import numpy as np
import screen


def geometry(box,segments):
    segment = segments[box['segment']]
    vertices = [(float(F(a)),float(F(x))) for a in box['alpha'] for x in box['row_density']]
    outer = [a*(float(segment.slope)*x+float(segment.intercept)) for a,x in vertices]
    probabilities = np.array([[a*x,a*(1-x),1-a] for a,x in vertices])
    entropy = np.sum(np.where(probabilities>0,probabilities*np.log(np.maximum(probabilities,1e-300)),0),axis=1)
    return vertices,probabilities,np.array(outer)+entropy


class Witnesses:
    def __init__(self,model,delta):
        self.model,self.delta = model,delta
        self.rows,self.logs,self.rates = [],[],[]
        self.keys = set()

    def add(self,w):
        p,y,lam,family = w['p'],w['y'],w['lam'],w['family']
        key = p,y,lam,family
        if key in self.keys:
            return
        self.keys.add(key)
        self.rows.append(dict(p=p,y=y,lam=lam,family=family))
        self.logs.append(np.log([p*y,p*(1-y),1-p]))
        self.rates.append(self.model.rate(p*y,lam,family)+self.delta*lam)

    def best(self,box,segments):
        vertices,probs,base = geometry(box,segments)
        values = base[:,None]-probs @ np.asarray(self.logs).T+np.asarray(self.rates)
        worst = values.max(axis=0)
        i = int(np.argmin(worst))
        return float(worst[i]),self.rows[i],list(map(float,values[:,i]))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--grid',type=Path,default=screen.HERE/'GRID_D11_STATEFREE.json')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--nodes',type=int,default=1200)
    p.add_argument('--seconds',type=float,default=180)
    p.add_argument('--delta',type=float,default=.11)
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    model = screen.Model()
    segments = screen.load_segments(screen.FROZEN/'golay_ba3_concave_majorant.json')
    old = json.loads((screen.FROZEN/'golay_ba3_rm2sub_joint_interval_d11.json').read_text())
    grid = json.loads(args.grid.read_text())
    bank = Witnesses(model,args.delta)
    for box in old['boxes']:
        p,y,z = (float(F(box[k])) for k in ('candidate_probability','value_probability','z'))
        for family in ('occupation','fourier','statefree'):
            bank.add(dict(p=p,y=y,lam=-math.log(z),family=family))
    for w in grid['rows']:
        bank.add(w)
    pending = [dict(segment=b['segment'],alpha=b['alpha'],row_density=b['row_density'],root=i,path='') for i,b in enumerate(old['boxes'])]
    accepted,unresolved = [],[]
    start = time.monotonic()
    splits = optimized = 0
    while pending:
        box = pending.pop()
        value,w,vertices = bank.best(box,segments)
        if value < -1e-10:
            accepted.append(dict(**box,exponent=value,witness=w,vertices=vertices))
            continue
        if optimized>=args.nodes or time.monotonic()-start>args.seconds or len(box['path'])>=12:
            unresolved.append(dict(**box,exponent=value,witness=w,vertices=vertices))
            continue
        a0,a1 = map(F,box['alpha'])
        x0,x1 = map(F,box['row_density'])
        a,x = float((a0+a1)/2),float((x0+x1)/2)
        segment = segments[box['segment']]
        fit = screen.optimize(model,a,x,float(segment.slope)*x+float(segment.intercept),args.delta,extra=[(w['p'],w['y'],w['lam'])])
        bank.add(fit)
        optimized += 1
        value,w,vertices = bank.best(box,segments)
        if value < -1e-10:
            accepted.append(dict(**box,exponent=value,witness=w,vertices=vertices))
        else:
            # A nonnegative point result is unresolved, not an impossibility.
            corners,_,_ = geometry(box,segments)
            ca,cx = corners[int(np.argmax(vertices))]
            point = screen.optimize(model,ca,cx,float(segment.slope)*cx+float(segment.intercept),args.delta,
                                    extra=[(w['p'],w['y'],w['lam'])])
            bank.add(point)
            if point['exponent']>=0:
                unresolved.append(dict(**box,exponent=value,witness=w,vertices=vertices,point_diagnostic=point))
                continue
            # Compare two exact bisections using the current witness bank.
            options = []
            for axis,lo,hi in (('alpha',a0,a1),('row_density',x0,x1)):
                mid = (lo+hi)/2
                children = [dict(box,**{axis:[str(lo),str(mid)]},path=box['path']+'0'),
                            dict(box,**{axis:[str(mid),str(hi)]},path=box['path']+'1')]
                score = max(bank.best(child,segments)[0] for child in children)
                options.append((score,axis,children))
            _,_,children = min(options,key=lambda v:v[0])
            pending.extend(children)
            splits += 1
        if optimized%20==0:
            print('optimized',optimized,'accepted',len(accepted),'pending',len(pending),'unresolved',len(unresolved),flush=True)
    # Each original leaf is refined by a binary prefix tree; verify coverage
    # independently using exact areas and pairwise interiors within each root.
    leaves = accepted+unresolved
    for root,original in enumerate(old['boxes']):
        group = [b for b in leaves if b['root']==root]
        a0,a1 = map(F,original['alpha']);x0,x1 = map(F,original['row_density'])
        area = F(0)
        for i,b in enumerate(group):
            bl,bh = map(F,b['alpha']);xl,xh = map(F,b['row_density'])
            assert a0<=bl<bh<=a1 and x0<=xl<xh<=x1
            assert b['segment']==original['segment']
            area += (bh-bl)*(xh-xl)
            for c in group[:i]:
                cl,ch = map(F,c['alpha']);yl,yh = map(F,c['row_density'])
                assert max(bl,cl)>=min(bh,ch) or max(xl,yl)>=min(xh,yh)
        assert area==(a1-a0)*(x1-x0)
    result = dict(status='BINARY64_IMT_DENSE_BOX_SCREEN_NOT_CERTIFICATE',delta=args.delta,
                  accepted=len(accepted),unresolved=len(unresolved),splits=splits,optimized=optimized,
                  exact_refinement_geometry_checked=True,original_boxes=len(old['boxes']),
                  transfer_counts=dict(Counter(b['witness']['family'] for b in accepted)),
                  worst=max(leaves,key=lambda b:b['exponent']),leaves=leaves,
                  source_sha256={p.relative_to(screen.ROOT).as_posix():screen.sha(p) for p in [Path(__file__),Path(screen.__file__),args.grid,
                    screen.FROZEN/'golay_ba3_rm2sub_joint_interval_d11.json',screen.FROZEN/'golay_ba3_concave_majorant.json']},
                  limitations=['Bounds use binary64 spectral radii, not outward Collatz witnesses.',
                    'Only positive occupancy alpha>=1e-4; not a full asymptotic theorem.'])
    with args.output.open('x',encoding='utf-8') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('leaves','source_sha256')}),flush=True)


if __name__=='__main__':
    main()
