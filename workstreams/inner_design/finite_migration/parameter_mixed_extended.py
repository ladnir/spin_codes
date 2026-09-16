"""Extend the IMT dense pilot's tilt and Bernoulli-reference search banks."""
import argparse
import math
from pathlib import Path

import numpy as np
import parameter_mixed_dense as prior


class Dense(prior.Dense):
    def __init__(self,*args):
        super().__init__(*args)
        counts,_ = prior.base.grid.outer(self.block)
        probabilities = np.array([0.]+[1. if band == [self.block] else .5 for band in self.bands])
        costs = np.array([0.]+[math.log(counts[self.block]) if band == [self.block] else
            max(math.log(counts[w])-math.log(math.comb(self.block,w))+self.block*math.log(2)
                for w in band) for band in self.bands])
        self.probability_banks.append((0.,probabilities,costs))


def run(output,b,t,s,exponent,minimum,nodes):
    assert not output.exists()
    dense = Dense(t,s,b,exponent,np.arange(-7.,1.51,0.5))
    result = dense.search(minimum,nodes,60)
    base = prior.base
    sources = base.grid.ladder.candidate.sources()
    _,path = base.grid.outer(b)
    sources[path.relative_to(base.grid.ROOT).as_posix()] = base.grid.maps.sha(path)
    base.grid.ladder.model.base.write_new(output,dict(status='BINARY64_IMT_EXTENDED_DENSE_COVER',
        geometry=[b,t,s,exponent],inner=dense.record,dense=result,
        margin_bits=-result['log_union_upper']/math.log(2),full_distance_proved=False,source_sha256=sources))
    print('extended dense margin',-result['log_union_upper']/math.log(2),flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--minimum',type=int,default=65)
    p.add_argument('--nodes',type=int,default=1023)
    a = p.parse_args()
    run(a.output.resolve(),128,64,20,20,a.minimum,a.nodes)
