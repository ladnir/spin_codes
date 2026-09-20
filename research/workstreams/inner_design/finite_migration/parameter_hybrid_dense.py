"""Combine joint Fourier tuning with the complete fixed-weight IMT bound."""
import argparse
import math
from pathlib import Path

import numpy as np
import parameter_joint_dense as joint

maps,base = joint.maps,joint.base


class Dense(joint.Dense):
    def __init__(self,*args):
        super().__init__(*args)
        self.mixed = maps.dense.Dense(*args)
        assert self.mixed.record == self.record
        np.testing.assert_array_equal(self.mixed.ps,self.ps)
        np.testing.assert_array_equal(self.mixed.log_gammas,self.log_gammas)

    def evaluate(self,lower,upper):
        mixed = self.mixed.evaluate(lower,upper)
        if mixed is None:
            return None
        mixed['witness']['family'] = 'mixed'
        if mixed['own_log_bound'] < -80*math.log(2):
            return mixed
        tuned = super().evaluate(lower,upper)
        tuned['witness']['family'] = 'joint_fourier'
        return min((mixed,tuned),key=lambda row:row['own_log_bound'])

    def replay_moment(self,theta,witness):
        assert witness['family'] in ('mixed','joint_fourier')
        checker = self.mixed if witness['family'] == 'mixed' else self
        return checker.moment(theta,witness['tilt'])


def run(output,b,t,s,exponent,minimum,nodes):
    assert not output.exists()
    with maps.use():
        checker = Dense(t,s,b,exponent,sorted(set(np.arange(-7.,1.51,.5))|{.75}))
        result = checker.search(minimum,nodes,60)
        sources = base.grid.ladder.candidate.sources()
        _,path = base.grid.outer(b)
        sources[path.relative_to(base.grid.ROOT).as_posix()] = base.grid.maps.sha(path)
        base.grid.ladder.model.base.write_new(output,dict(status='BINARY64_IMT_HYBRID_DENSE_COVER',
            geometry=[b,t,s,exponent],inner=checker.record,dense=result,
            margin_bits=-result['log_union_upper']/math.log(2),full_distance_proved=False,source_sha256=sources))
        print('hybrid dense margin',-result['log_union_upper']/math.log(2),flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--b',type=int,choices=(64,128),default=128)
    p.add_argument('--t',type=int,choices=(64,128,256),default=64)
    p.add_argument('--s',type=int,default=20)
    p.add_argument('--m',type=int,default=20)
    p.add_argument('--minimum',type=int,default=65)
    p.add_argument('--nodes',type=int,default=1023)
    a = p.parse_args()
    run(a.output.resolve(),a.b,a.t,a.s,a.m,a.minimum,a.nodes)
