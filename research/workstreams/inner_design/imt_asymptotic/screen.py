"""IMT asymptotic witness discovery. Binary64 diagnostics, NOT a certificate."""
import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
FROZEN = ROOT/'workstreams/paper_architecture/certificates/single_sampled_ba_rm2sub'
sys.path.insert(0, str(HERE.parent/'asymmetric/bch256/weight5'))
import candidate
sys.path.insert(0, str(HERE.parent/'asymmetric'))
import screen_dense
sys.path.insert(0, str(FROZEN))
from certify_golay_ba_rm2sub_joint_interval import load_segments


def kl(x, y):
    if x == 0:
        return -math.log1p(-y)
    if x == 1:
        return -math.log(y)
    return x*math.log(x/y)+(1-x)*math.log((1-x)/(1-y))


class Model:
    """Vectorized version of the existing fixed-weight envelope, unchanged bounds."""
    def __init__(self, name='weight5_seed0'):
        self.engine = candidate.Engine(20, name)
        e = self.engine
        self.levels = np.array(e.levels)
        self.counts = np.array([e.spectrum[v] for v in e.levels], dtype=float)
        self.m = e.m
        self.j = np.arange(129)
        totals = [math.comb(128, j) for j in range(129)]
        self.logcombs = np.log(np.array(totals, dtype=float))
        self.zero = np.array([k/t for k,t in zip(e.kernel, totals)])
        self.live = np.array([(t-k)/t for k,t in zip(e.kernel, totals)])
        self.cap = np.array([int(c['cap'])/t for c,t in zip(e.caps, totals)])
        self.mindiff = np.min(abs(self.levels[:, None]-self.j), axis=0)
        self.diff = abs(self.levels[:, None]-self.j)
        self.shellcap = np.minimum(self.live[None, :]/self.counts[:, None], self.cap)
        self.hyper = np.zeros((len(e.levels), 129, 129))
        for i,w in enumerate(e.levels):
            for j in range(129):
                for h in range(max(0,j-128+w), min(w,j)+1):
                    self.hyper[i,j,w+j-2*h] = math.comb(w,h)*math.comb(128-w,j-h)/totals[j]
        self.low_patterns = {}
        self.low_shells = {}
        for j,data in e.low.items():
            rows = np.zeros((len(data['patterns']),129))
            for i,pattern in enumerate(data['patterns']):
                for w,count in pattern:
                    rows[i,w] = count/totals[j]
            self.low_patterns[j] = rows
            rows = np.zeros((len(e.levels),129))
            for i,v in enumerate(e.levels):
                for w,count in data['by_weight'].get(v,{}).items():
                    rows[i,w] = count/(int(self.counts[i])*totals[j])
            self.low_shells[j] = rows
        self.bweights = np.array(list(e.b_spectrum))
        self.bcounts = np.array(list(e.b_spectrum.values()), dtype=float)

    def occupation(self, beta, lam):
        powers = np.exp(-lam*self.j)
        moments = self.hyper @ powers
        arbitrary = np.max(moments, axis=0)
        cd = np.minimum(np.minimum(arbitrary,self.live), self.cap*np.exp(-lam*self.mindiff))
        cs = np.minimum(moments,self.shellcap*np.exp(-lam*self.diff))
        for j,rows in self.low_patterns.items():
            cd[j] = min(cd[j],np.max(rows @ powers))
            cs[:,j] = np.minimum(cs[:,j],self.low_shells[j] @ powers)
        probs = np.exp(self.logcombs+self.j*math.log(beta)+(128-self.j)*math.log1p(-beta))
        n = self.engine.n
        out = np.zeros((n,n))
        out[0,0] = probs @ (self.zero*powers)
        out[0,1] = probs @ (self.live*powers)
        out[1,0] = probs @ (cd/2+np.minimum(arbitrary,self.live)/(2*self.m))
        out[1,1] = (probs @ arbitrary)/2
        out[1,2:] = (probs @ arbitrary)*self.counts/(2*self.m)
        out[2:,0] = (cs/2+np.minimum(moments,self.live)/(2*self.m)) @ probs
        out[2:,1] = (moments @ (probs*(self.live>0)))/2
        out[2:,2:] = (moments @ probs)[:,None]*self.counts[None,:]/(2*self.m)
        out[np.arange(2,n),np.arange(2,n)] += (moments @ (probs*(self.live==0)))/2
        return out

    def fourier(self, beta, lam):
        z = math.exp(-lam)
        g0,g1 = 1-beta+beta*z, beta+(1-beta)*z
        r0,r1 = abs(1-2*beta*z/g0),abs(1-2*(1-beta)*z/g1)
        w = self.bweights[None,:]
        v = self.levels[:,None]
        h = np.minimum(v,w) if r1>=r0 else np.maximum(0,v+w-128)
        caps = (1+(np.power(r0,w-h)*np.power(r1,h)) @ self.bcounts)/(self.m+1)
        entries = (g0**(128-self.levels)*g1**self.levels)*(caps/2+1/(2*self.m))
        out = np.zeros((self.engine.n,self.engine.n))
        weighted = np.exp(self.logcombs+self.j*(math.log(beta)-lam)+(128-self.j)*math.log1p(-beta))
        out[0,0] = weighted @ self.zero
        out[0,1] = weighted @ self.live
        out[1:,0] = np.r_[entries.max(),entries]
        out[1:,2:] = out[1:,0,None]*self.counts
        return out

    @lru_cache(maxsize=20000)
    def rate(self,beta,lam,family):
        if family=='statefree':
            z = math.exp(-lam)
            g0,g1 = 1-beta+beta*z,beta+(1-beta)*z
            return max(128*math.log(g0),max((128-w)*math.log(g0)+w*math.log(g1) for w in self.levels))/128
        matrix = self.occupation(beta,lam) if family=='occupation' else self.fourier(beta,lam)
        radius = float(np.max(np.abs(np.linalg.eigvals(matrix))))
        return math.log(radius)/128

    def evaluate(self,alpha,x,outer,p,y,lam,delta,family):
        return alpha*outer+kl(alpha,p)+alpha*kl(x,y)+self.rate(p*y,lam,family)+delta*lam


def optimize(model,alpha,x,outer,delta,extra=()):
    # At beta=1/2 the state-independent bound is exact. Minimize the chain
    # KL cost over p,y subject to p*y=1/2 analytically.
    q = alpha*x
    p = min(1-1e-12,.5+.5*alpha*(1-x)/(1-q))
    y,lam = .5/p,math.log((1-delta)/delta)
    best = dict(alpha=alpha,x=x,outer=outer,delta=delta,
                exponent=model.evaluate(alpha,x,outer,p,y,lam,delta,'statefree'),
                p=p,y=y,lam=lam,family='statefree',optimizer_success=True)
    starts = [(alpha,x,max(1e-5,alpha)),(min(.999,2*alpha),.4,max(1e-5,alpha)),(.8,.5,1.8)]
    starts.extend(extra)
    for family in ('occupation','fourier'):
        for p,y,lam in starts:
            p,y = np.clip([p,y],1e-10,1-1e-10)
            start = np.array([math.log(p/(1-p)),math.log(y/(1-y)),math.log(lam)])
            def decode(v):
                return float(np.clip(expit(v[0]),1e-12,1-1e-12)),float(np.clip(expit(v[1]),1e-12,1-1e-12)),math.exp(float(np.clip(v[2],-25,2.4)))
            def objective(v):
                return model.evaluate(alpha,x,outer,*decode(v),delta,family)
            fit = minimize(objective,start,method='Nelder-Mead',options={'maxiter':350,'xatol':1e-8,'fatol':1e-12})
            if float(fit.fun)<best['exponent']:
                p,y,lam = decode(fit.x)
                best = dict(alpha=alpha,x=x,outer=outer,delta=delta,exponent=float(fit.fun),p=p,y=y,lam=lam,family=family,optimizer_success=bool(fit.success))
    return best


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('replay','grid','points'),default='replay')
    parser.add_argument('--delta',type=float,default=.11)
    parser.add_argument('--alphas',nargs='+',type=float,default=[.0001,.001,.01,.1,.5,1.])
    parser.add_argument('--xs',nargs='+',type=float,default=[.104,.16,.25,.4,.5,.6,.75,.84,.896])
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    model = Model()
    majorant = FROZEN/'golay_ba3_concave_majorant.json'
    segments = load_segments(majorant)
    def outer(x):
        return min(float(s.slope)*x+float(s.intercept) for s in segments)
    rows = []
    if args.mode=='replay':
        old = json.loads((FROZEN/'golay_ba3_rm2sub_joint_interval_d11.json').read_text())
        for idx,box in enumerate(old['boxes']):
            p,y,z = (float(Fraction(box[k])) for k in ('candidate_probability','value_probability','z'))
            segment = segments[box['segment']]
            vertices = [(float(Fraction(a)),float(Fraction(x))) for a in box['alpha'] for x in box['row_density']]
            witnesses = []
            for family in ('occupation','fourier'):
                values = [model.evaluate(a,x,float(segment.slope)*x+float(segment.intercept),p,y,-math.log(z),args.delta,family) for a,x in vertices]
                witnesses.append((max(values),family,values))
            value,family,values = min(witnesses)
            rows.append(dict(index=idx,alpha=box['alpha'],x=box['row_density'],exponent=value,vertices=values,family=family,p=p,y=y,lam=-math.log(z)))
            if (idx+1)%100==0:
                print('replay',idx+1,'nonnegative',sum(r['exponent']>=0 for r in rows),flush=True)
    else:
        for alpha in args.alphas:
            for x in args.xs:
                row = optimize(model,alpha,x,outer(x),args.delta)
                rows.append(row)
                print(json.dumps(row),flush=True)
    source_paths = [Path(__file__),majorant,FROZEN/'golay_ba3_rm2sub_joint_interval_d11.json',
                    HERE.parent/'NO_CONSTANT_MAP.json',Path(candidate.__file__),Path(screen_dense.__file__),
                    Path(candidate.model.independent.g.__file__)]
    result = dict(status='BINARY64_IMT_ASYMPTOTIC_SCREEN_NOT_CERTIFICATE',mode=args.mode,delta=args.delta,
                  inner=model.engine.identity()['inner'],rows=rows,nonnegative=sum(r['exponent']>=0 for r in rows),
                  worst=max(rows,key=lambda r:r['exponent']),
                  source_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in source_paths},
                  limitations=['No outward arithmetic or IMT asymptotic theorem.',
                    'Replay covers existing boxes only in binary64; optimized grid covers points only.',
                    'Reuses the certified concave outer majorant, not a guessed ensemble exponent.'])
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x',encoding='utf-8') as stream:
        json.dump(result,stream,indent=2)
        stream.write('\n')
    print(json.dumps({k:result[k] for k in ('status','nonnegative','worst')}),flush=True)


if __name__=='__main__':
    main()
