"""Check old dense bottleneck rectangles under changed transvection laws.

Only the coupled Bernoulli bound is used. No old-law bound is retained as a
fallback, and these isolated rectangles never imply complete coverage.
"""
import argparse
from pathlib import Path

from flint import arb, ctx
import coupled_input_dense as coupled
import budget_dense
import mixing_rounds as study

model = study.model
HERE = Path(__file__).resolve().parent


class Engine(study.ladder.Engine):
    def __init__(self, rounds):
        super().__init__(16)
        self.rounds = str(rounds)
        self.eps = study.epsilon(rounds)

    def identity(self):
        result = super().identity()
        result['inner']['transvection_rounds'] = None if self.rounds == 'refresh' else int(self.rounds)
        result['inner']['mixing_law'] = self.rounds
        result['setup'] = 'independent row/region permutations; mixer rounds='+self.rounds+'; zero start; output before update; persistent state; no flush'
        return result

    def bernoulli(self, theta, lam):
        assert theta > 0 and theta < 1 and lam > 0
        z = (-lam).exp()
        g0, g1 = 1-theta+theta*z, theta+(1-theta)*z
        rho0 = min(arb(1),model.up(abs(1-2*theta*z/g0)))
        rho1 = min(arb(1),model.up(abs(1-2*(1-theta)*z/g1)))
        eps, eta = model.number(self.eps), model.number(1-self.eps)
        entries = []
        for v in self.levels:
            cap = arb(1)
            for w,count in self.b_spectrum.items():
                ends = (max(0,v+w-128),min(v,w))
                cap += count*max(model.up(rho0**(w-h)*rho1**h) for h in ends)
            moment = g0**(128-v)*g1**v
            entries.append(model.up(moment*(eps*cap/(self.m+1)+eta/self.m)))
        n = self.n
        matrix = [arb(0)]*(n*n)
        matrix[0] = sum((self.kernel[j]*(theta*z)**j*(1-theta)**(128-j) for j in range(129)),arb(0))
        matrix[1] = sum(((study.math.comb(128,j)-self.kernel[j])*(theta*z)**j*(1-theta)**(128-j)
                         for j in range(129)),arb(0))
        for i,entry in [(1,max(entries))]+[(i+2,v) for i,v in enumerate(entries)]:
            matrix[i*n] = entry
            for k,w in enumerate(self.levels):
                matrix[i*n+k+2] = entry*self.spectrum[w]
        return tuple(map(model.up,matrix))


def run(output, verify=False):
    ctx.prec = 512 if verify else 256
    if not verify and output.exists():
        raise FileExistsError('Use a fresh output path')
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
    source = HERE/'DENSE_M16_coupled_v1.json'
    dense = model.base.read(source)
    keys = sorted(dense['leaves'],key=lambda k:dense['leaves'][k]['power'],reverse=True)[:3]
    results = []
    for rounds in ('1','2','3','refresh'):
        checker = coupled.Checker(16)
        checker.engine = Engine(rounds)
        for key in keys:
            node = dense['leaves'][key]
            # Do not call checker.bound: its fallback includes old-law transfers.
            value = checker.coupled_bound(*budget_dense.geometry.geometry(node),node['witness'])
            assert value is not None
            upper = model.pack(model.up(value.exp()))
            result = dict(rounds=rounds,key=key,instance=checker.engine.identity(),upper=upper)
            if saved:
                old = saved['rectangles'][len(results)]
                assert (old['rounds'],old['key'],old['instance']) == (rounds,key,checker.engine.identity())
                assert model.unpack(upper) <= model.unpack(old['upper'])
                result = old
            results.append(result)
            print(ctx.prec,rounds,key,study.audit.bits(model.exact(model.unpack(result['upper']))),flush=True)
    if saved:
        model.base.write_new(output.with_name(output.stem+'_replay.json'),dict(
            status='IMT_MIXING_DENSE_PROBE_512_BIT_REPLAY_PASSED',producer_sha256=model.base.sha(output),
            rectangles_checked=len(results),full_distance_proved=False))
    else:
        model.base.write_new(output,dict(status='OUTWARD_IMT_MIXING_DENSE_RECTANGLE_PROBE',
            rectangles=results,full_distance_proved=False,source_sha256={**study.ladder.candidate.sources(),
                source.relative_to(model.ROOT).as_posix():model.base.sha(source)}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--verify',action='store_true')
    args = p.parse_args()
    run(args.output.resolve(),args.verify)
