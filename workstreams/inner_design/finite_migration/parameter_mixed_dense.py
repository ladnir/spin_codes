"""Complete dense type boxes with both independent-map IMT moment bounds."""
import argparse
import math
from pathlib import Path

import numpy as np
import parameter_full as base


class Dense(base.Dense):
    def __init__(self,*args):
        super().__init__(*args)
        self.epoch_cache = {}
        _,_,_,self.caps,self.low = base.prepare(self.t,self.record['s'])

    def moment(self,theta,tilt):
        direct = super().moment(theta,tilt)
        if tilt not in self.epoch_cache:
            self.epoch_cache[tilt] = base.g.epochs(self.record['spectrum'],self.kernel,
                self.record['feedback_columns'],self.caps,self.low,math.exp(tilt),maximum=self.t)
        epoch = self.epoch_cache[tilt]
        probabilities = np.array([math.log(math.comb(self.t,j))+j*math.log(theta)
                                  +(self.t-j)*math.log1p(-theta) for j in range(self.t+1)])
        matrix = np.logaddexp.reduce(epoch+probabilities[:,None,None],axis=0)
        occupation = base.g.terminal(matrix,self.block*self.length//self.t)
        return min(direct,occupation)


def run(output,b,t,s,exponent,minimum,nodes):
    assert not output.exists()
    dense = Dense(t,s,b,exponent,np.arange(-7.,0.01,0.5))
    result = dense.search(minimum,nodes,60)
    sources = base.grid.ladder.candidate.sources()
    _,path = base.grid.outer(b)
    sources[path.relative_to(base.grid.ROOT).as_posix()] = base.grid.maps.sha(path)
    base.grid.ladder.model.base.write_new(output,dict(status='BINARY64_IMT_MIXED_DENSE_COVER',
        geometry=[b,t,s,exponent],inner=dense.record,dense=result,
        margin_bits=-result['log_union_upper']/math.log(2),full_distance_proved=False,source_sha256=sources))
    print('mixed dense margin',-result['log_union_upper']/math.log(2),flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--minimum',type=int,default=65)
    p.add_argument('--nodes',type=int,default=1023)
    a = p.parse_args()
    run(a.output.resolve(),128,64,20,20,a.minimum,a.nodes)
