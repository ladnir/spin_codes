"""Jointly tune both IMT moment families only where the hybrid cover is weak."""
import argparse
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
import parameter_hybrid_safe as safe
import parameter_activated_dense as activated

maps,base = activated.maps,activated.base


class Refined(activated.Dense):
    def moment(self,theta,tilt):
        if not 0 < theta < 1:
            return math.inf
        if len(self.epoch_cache) > 4:
            self.epoch_cache.clear()
        return super().moment(theta,tilt)


class Dense(safe.Dense):
    def __init__(self,*args):
        super().__init__(*args)
        self.refined = Refined(*args)

    def evaluate(self,lower,upper):
        row = super().evaluate(lower,upper)
        if row is None or row['own_log_bound'] < -80*math.log(2):
            return row
        w = row['witness']
        initial = np.array(w['proposal'])
        anchor = int(np.argmax(initial))
        free = [j for j in range(len(initial)) if j != anchor]
        penalty = base.typed.lattice_log_count(lower,upper)
        def objective(coordinates,detail=False):
            tilt = float(coordinates[-1])
            if not -12 <= tilt <= 2 or max(abs(coordinates[:-1])) > 40:
                return 1e6
            logits = np.zeros(len(initial)); logits[free] = coordinates[:-1]
            proposal = np.exp(logits-logsumexp(logits))
            moment = self.refined.moment(float(proposal@self.ps),tilt)
            values = base.typed.point_logs(row['corners'],self.length,self.block,self.log_gammas,
                proposal,moment,self.cutoff,math.exp(tilt))
            value = float(max(values))+penalty
            return (value,dict(tilt=tilt,probability_bank=0,proposal=proposal.tolist(),family='joint_mixed')) if detail else value/(self.block*self.length)
        start = np.r_[np.log(initial/initial[anchor])[free],w['tilt']]
        result = minimize(objective,start,method='Nelder-Mead',
            options=dict(maxiter=300,xatol=1e-6,fatol=1e-10))
        value,witness = objective(result.x,True)
        if value < row['own_log_bound']:
            row.update(own_log_bound=value,best_log_bound=value,witness=witness)
        return row

    def replay_moment(self,theta,witness):
        if witness['family'] == 'joint_mixed':
            return self.refined.moment(theta,witness['tilt'])
        return super().replay_moment(theta,witness)


def run(output,minimum,nodes):
    assert not output.exists()
    with maps.use():
        checker = Dense(64,20,128,20,sorted(set(np.arange(-7.,1.51,.5))|{.75}))
        result = checker.search(minimum,nodes,60)
        sources = base.grid.ladder.candidate.sources()
        _,path = base.grid.outer(128)
        sources[path.relative_to(base.grid.ROOT).as_posix()] = base.grid.maps.sha(path)
        base.grid.ladder.model.base.write_new(output,dict(status='BINARY64_IMT_JOINT_MIXED_DENSE_COVER',
            geometry=[128,64,20,20],inner=checker.record,dense=result,
            margin_bits=-result['log_union_upper']/math.log(2),full_distance_proved=False,source_sha256=sources))
        print('joint mixed dense margin',-result['log_union_upper']/math.log(2),flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--minimum',type=int,default=65)
    p.add_argument('--nodes',type=int,default=1023)
    a = p.parse_args()
    run(a.output.resolve(),a.minimum,a.nodes)
