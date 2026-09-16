"""Retain feedback-syndrome density when a Bernoulli epoch leaves zero."""
import argparse
import math
from pathlib import Path

import numpy as np
import parameter_no_constant as maps

base = maps.full


def activate(matrix,spectrum,b_spectrum,t,theta,lam):
    log_g0 = math.log1p(theta*math.expm1(-lam))
    tilted = theta*math.exp(-lam-log_g0)
    rho = abs(1-2*tilted)
    log_rho = math.log(rho) if rho else -math.inf
    log_cap = t*log_g0+float(np.logaddexp.reduce(
        [0.]+[math.log(count)+w*log_rho for w,count in b_spectrum.items()]))-math.log(sum(spectrum.values())+1)
    out = matrix.copy()
    out[0,1] = -math.inf
    for k,w in enumerate(sorted(spectrum)):
        out[0,k+2] = log_cap+math.log(spectrum[w])
    return out


class Dense(maps.dense.Dense):
    def moment(self,theta,tilt):
        old = super().moment(theta,tilt)
        lam = math.exp(tilt)
        direct = base.direct.bernoulli(self.record['spectrum'],self.bs,self.kernel,float(theta),lam)
        probabilities = np.array([math.log(math.comb(self.t,j))+j*math.log(theta)
                                  +(self.t-j)*math.log1p(-theta) for j in range(self.t+1)])
        occupation = np.logaddexp.reduce(self.epoch_cache[tilt]+probabilities[:,None,None],axis=0)
        values = [base.g.terminal(activate(matrix,self.record['spectrum'],self.bs,self.t,theta,lam),
                                 self.block*self.length//self.t) for matrix in (direct,occupation)]
        return min(old,*values)


def run(output,minimum,nodes):
    assert not output.exists()
    with maps.use():
        checker = Dense(64,20,128,20,sorted(set(np.arange(-7.,1.51,.5))|{.75}))
        result = checker.search(minimum,nodes,60)
        sources = base.grid.ladder.candidate.sources()
        _,path = base.grid.outer(128)
        sources[path.relative_to(base.grid.ROOT).as_posix()] = base.grid.maps.sha(path)
        base.grid.ladder.model.base.write_new(output,dict(status='BINARY64_IMT_ACTIVATED_DENSE_COVER',
            geometry=[128,64,20,20],inner=checker.record,dense=result,
            margin_bits=-result['log_union_upper']/math.log(2),full_distance_proved=False,source_sha256=sources))
        print('activated dense margin',-result['log_union_upper']/math.log(2),flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--minimum',type=int,default=65)
    p.add_argument('--nodes',type=int,default=1023)
    a = p.parse_args()
    run(a.output.resolve(),a.minimum,a.nodes)
