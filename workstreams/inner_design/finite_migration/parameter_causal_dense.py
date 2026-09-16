"""Add the causal XOR moment bound to the IMT dense diagnostics."""
import argparse
import math
from pathlib import Path

import numpy as np
import parameter_mixed_extended as prior


def causal_moment(theta,lam,length):
    if not 0 <= theta <= 1 or lam < 0 or length < 0:
        raise ValueError('Invalid Bernoulli input or moment geometry')
    # Condition on every preceding input epoch and on the sampled setup.
    # Each fresh input bit is XORed with a now-fixed bit. Its exponential
    # moment is at most 1-p+p*exp(-lam), p=min(theta,1-theta).
    return length*math.log1p(min(theta,1-theta)*math.expm1(-lam))


class Dense(prior.Dense):
    def moment(self,theta,tilt):
        generic = super().moment(theta,tilt)
        causal = causal_moment(theta,math.exp(tilt),self.block*self.length)
        return min(generic,causal)


def run(output,b,t,s,exponent,minimum,nodes):
    assert not output.exists()
    dense = Dense(t,s,b,exponent,np.arange(-7.,1.51,0.5))
    result = dense.search(minimum,nodes,60)
    base = prior.prior.base
    sources = base.grid.ladder.candidate.sources()
    _,path = base.grid.outer(b)
    sources[path.relative_to(base.grid.ROOT).as_posix()] = base.grid.maps.sha(path)
    base.grid.ladder.model.base.write_new(output,dict(status='BINARY64_IMT_CAUSAL_DENSE_COVER',
        geometry=[b,t,s,exponent],inner=dense.record,dense=result,
        margin_bits=-result['log_union_upper']/math.log(2),full_distance_proved=False,source_sha256=sources))
    print('causal dense margin',-result['log_union_upper']/math.log(2),flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--minimum',type=int,default=65)
    p.add_argument('--nodes',type=int,default=1023)
    a = p.parse_args()
    run(a.output.resolve(),128,64,20,20,a.minimum,a.nodes)
