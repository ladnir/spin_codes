"""Keep a pointwise density envelope after lazy IMT updates.

For a source uniform on an A-weight shell of size a_v, each target determines
the source q=y+CX. Nonzero syndrome fibers give weighted target mass at most
min(nu_j,a_v*h_j)/a_v * z**abs(v-j). Zero syndrome preserves the source shell
and costs beta_j*z**abs(v-j). This bounds a complete measure representation;
we compare complete moments, never entrywise minima of representations.
"""
import argparse
import math
from pathlib import Path
import numpy as np
import parameter_syndrome_activation as previous

prior = previous.prior


def diffuse(epochs,spectrum,kernel,caps,t,lam):
    result = epochs.copy()
    levels = sorted(spectrum)
    m = sum(spectrum.values())
    for j in range(t+1):
        count = math.comb(t,j)
        beta = kernel[j]/count
        nu = (count-kernel[j])/count
        if not nu:
            continue  # The existing lazy branch already preserves its shell.
        h = int(caps[j]['cap'])/count
        for i,v in enumerate(levels):
            a = spectrum[v]
            # The old S_i -> D entry is half its marginal output moment.
            marginal = epochs[j,i+2,1]+math.log(2)
            point = min(marginal,math.log(min(nu,a*h)/a)-abs(v-j)*lam) if h else -math.inf
            result[j,i+2,1] = -math.inf
            for k,w in enumerate(levels):
                fresh = marginal+math.log(spectrum[w]/(2*m))
                lazy = point+math.log(spectrum[w]/2)
                result[j,i+2,k+2] = np.logaddexp(fresh,lazy)
            if beta:
                preserve = math.log(beta/2)-abs(v-j)*lam
                result[j,i+2,i+2] = np.logaddexp(result[j,i+2,i+2],preserve)
    return result


class Dense(previous.Dense):
    def moment(self,theta,tilt):
        old = super().moment(theta,tilt)
        lam = math.exp(tilt)
        epoch = diffuse(self.epoch_cache[tilt],self.record['spectrum'],self.kernel,self.caps,self.t,lam)
        probabilities = np.array([math.log(math.comb(self.t,j))+j*math.log(theta)
                                 +(self.t-j)*math.log1p(-theta) for j in range(self.t+1)])
        matrix = np.logaddexp.reduce(epoch+probabilities[:,None,None],axis=0)
        log_g0 = math.log1p(theta*math.expm1(-lam))
        distribution = previous.syndrome_probabilities(self.b_weights,theta*math.exp(-lam-log_g0))
        guard = 32*(self.t+self.record['s'])*np.finfo(float).eps
        matrix[0,1] = -math.inf
        for k,indices in enumerate(self.shell_indices):
            cap = max(0.,float(distribution[indices].max()))+guard
            matrix[0,k+2] = self.t*log_g0+math.log(cap*len(indices))
        return min(old,prior.base.g.terminal(matrix,self.block*self.length//self.t))


def run(seed,output):
    assert not output.exists()
    model = prior.base.grid.ladder.model
    saved = model.base.read(seed); model.authenticate(saved)
    assert saved['status'] == 'BINARY64_IMT_DENSE_FACE_MULTISTART'
    with prior.maps.use():
        checker = Dense(64,20,128,20,[])
        point = np.array(saved['point'])
        selected = max(saved['refined'],key=lambda r:r['margin_bits'])
        proposal = np.array(selected['proposal']);tilt = selected['tilt']
        moment = checker.moment(float(proposal@checker.ps),tilt)
        bound = float(prior.base.typed.point_logs(point[None],checker.length,128,checker.log_gammas,
            proposal,moment,checker.cutoff,math.exp(tilt))[0])
        row = dict(previous_margin_bits=selected['margin_bits'],margin_bits=-bound/math.log(2),
                   proposal=proposal.tolist(),tilt=tilt,moment_log=moment)
        print('diffuse lazy',row,flush=True)
        sources = prior.base.grid.ladder.candidate.sources()
        _,outer_path = prior.base.grid.outer(128)
        for path in (seed,outer_path):
            sources[path.relative_to(prior.base.grid.ROOT).as_posix()] = model.base.sha(path)
        model.base.write_new(output,dict(status='BINARY64_IMT_DIFFUSE_LAZY_POINT',
            geometry=saved['geometry'],inner=checker.record,point=point.tolist(),result=row,
            full_distance_proved=False,source_sha256=sources))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a = p.parse_args()
    run(a.seed.resolve(),a.output.resolve())
