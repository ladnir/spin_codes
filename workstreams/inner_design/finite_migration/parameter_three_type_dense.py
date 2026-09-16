"""Coarsen ordinary shells to reduce the dimension of the IMT dense cover."""
import argparse
import math
from pathlib import Path

import numpy as np
import parameter_causal_dense as prior


class Dense(prior.Dense):
    def __init__(self,*args):
        super().__init__(*args)
        counts,_ = prior.prior.prior.base.grid.outer(self.block)
        ordinary = sorted(w for w in counts if w != self.block)
        self.bands = [ordinary,[self.block]]
        self.ps = np.array([0.,.5,1.])
        cost = max(math.log(counts[w])-math.log(math.comb(self.block,w))
                   +self.block*math.log(2) for w in ordinary)
        self.log_gammas = np.array([0.,cost,math.log(counts[self.block])])
        self.probability_banks = [(0.,self.ps,self.log_gammas)]


def run(output,b,t,s,exponent,minimum,nodes):
    assert not output.exists()
    dense = Dense(t,s,b,exponent,sorted(set(np.arange(-7.,1.51,0.5))|{.75}))
    result = dense.search(minimum,nodes,60)
    base = prior.prior.prior.base
    sources = base.grid.ladder.candidate.sources()
    _,path = base.grid.outer(b)
    sources[path.relative_to(base.grid.ROOT).as_posix()] = base.grid.maps.sha(path)
    base.grid.ladder.model.base.write_new(output,dict(status='BINARY64_IMT_THREE_TYPE_DENSE_COVER',
        geometry=[b,t,s,exponent],inner=dense.record,dense=result,
        margin_bits=-result['log_union_upper']/math.log(2),full_distance_proved=False,source_sha256=sources))
    print('three-type dense margin',-result['log_union_upper']/math.log(2),flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--minimum',type=int,default=65)
    p.add_argument('--nodes',type=int,default=1023)
    a = p.parse_args()
    run(a.output.resolve(),128,64,20,20,a.minimum,a.nodes)
