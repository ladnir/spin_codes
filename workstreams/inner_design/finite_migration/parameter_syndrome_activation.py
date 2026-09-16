"""Numerical zero-state activation using every feedback syndrome.

Fourier inversion computes P(CX=q) for iid tilted input. A shell coefficient
is its maximum point mass times the number of states in that A-weight shell.
No independence between the maps is assumed. Binary64 diagnostics only.
"""
import argparse
import math
from pathlib import Path
import numpy as np
import parameter_activated_dense as prior


def image_weights(columns, s):
    assert len(columns) <= 64
    words = prior.base.grid.maps.generator_words(columns,s)
    values = np.zeros(1 << s,dtype=np.uint64)
    for j,word in enumerate(words):
        width = 1 << j
        values[width:2*width] = values[:width] ^ np.uint64(word)
    return np.bitwise_count(values)


def syndrome_probabilities(dual_weights, probability):
    n = len(dual_weights)
    assert n and n & (n-1) == 0 and 0 <= probability <= 1
    values = np.power(1-2*probability,dual_weights.astype(np.int64))
    width = 1
    while width < n:
        blocks = values.reshape(-1,2*width)
        left = blocks[:,:width].copy()
        right = blocks[:,width:].copy()
        blocks[:,:width] = left+right
        blocks[:,width:] = left-right
        width *= 2
    return values/n


class Dense(prior.Dense):
    def __init__(self,*args):
        super().__init__(*args)
        s = self.record['s']
        self.a_weights = image_weights(self.record['expansion_columns'],s)
        self.b_weights = image_weights(self.record['feedback_columns'],s)
        self.shell_indices = [np.flatnonzero(self.a_weights == w) for w in sorted(self.record['spectrum'])]

    def moment(self,theta,tilt):
        old = super().moment(theta,tilt)
        lam = math.exp(tilt)
        log_g0 = math.log1p(theta*math.expm1(-lam))
        probability = theta*math.exp(-lam-log_g0)
        distribution = syndrome_probabilities(self.b_weights,probability)
        # Guard cancellation in the numerical transform. This deliberately
        # loose floor is not a substitute for an outward certificate replay.
        guard = 32*(self.t+self.record['s'])*np.finfo(float).eps
        caps = [max(0.,float(distribution[indices].max()))+guard for indices in self.shell_indices]
        direct = prior.base.direct.bernoulli(self.record['spectrum'],self.bs,self.kernel,theta,lam)
        probabilities = np.array([math.log(math.comb(self.t,j))+j*math.log(theta)
                                 +(self.t-j)*math.log1p(-theta) for j in range(self.t+1)])
        occupation = np.logaddexp.reduce(self.epoch_cache[tilt]+probabilities[:,None,None],axis=0)
        values = []
        for matrix in (direct,occupation):
            matrix[0,1] = -math.inf
            for k,indices in enumerate(self.shell_indices):
                matrix[0,k+2] = self.t*log_g0+math.log(caps[k]*len(indices))
            values.append(prior.base.g.terminal(matrix,self.block*self.length//self.t))
        return min(old,*values)


def run(seed,output):
    assert not output.exists()
    model = prior.base.grid.ladder.model
    saved = model.base.read(seed); model.authenticate(saved)
    assert saved['status'] == 'BINARY64_IMT_DENSE_FACE_MULTISTART'
    with prior.maps.use():
        checker = Dense(64,20,128,20,[])
        point = np.array(saved['point'])
        rows = []
        for selected in sorted(saved['refined'],key=lambda r:r['margin_bits'],reverse=True)[:2]:
            proposal = np.array(selected['proposal']); tilt = selected['tilt']
            moment = checker.moment(float(proposal@checker.ps),tilt)
            bound = float(prior.base.typed.point_logs(point[None],checker.length,128,
                checker.log_gammas,proposal,moment,checker.cutoff,math.exp(tilt))[0])
            row = dict(previous_margin_bits=selected['margin_bits'],margin_bits=-bound/math.log(2),
                       proposal=proposal.tolist(),tilt=tilt,moment_log=moment)
            rows.append(row); print('exact syndrome activation',row,flush=True)
        sources = prior.base.grid.ladder.candidate.sources()
        _,outer_path = prior.base.grid.outer(128)
        for path in (seed,outer_path):
            sources[path.relative_to(prior.base.grid.ROOT).as_posix()] = model.base.sha(path)
        model.base.write_new(output,dict(status='BINARY64_IMT_SYNDROME_ACTIVATION_POINT',
            geometry=saved['geometry'],inner=checker.record,point=point.tolist(),rows=rows,
            full_distance_proved=False,source_sha256=sources))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a = p.parse_args()
    run(a.seed.resolve(),a.output.resolve())
