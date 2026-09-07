"""Outward replay of a complete fixed-type tail, optionally clipped at its start.

Binary64 witnesses specify exact rational probabilities, normalized exactly.
Neither floating-point counting costs nor reported moments are trusted.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import math
from pathlib import Path

from flint import arb, ctx

import activation_density_arb as transfer
import bch256_dense_types as search
from audit_bch_q1_full_arb import rational

core, base = search.core, search.base


def aa(value):
    value = F(value)
    return arb(value.numerator)/value.denominator


def probabilities(witness):
    proposal = [F(p) for p in witness['proposal']]
    ps = tuple(F(p) for p in witness['probabilities'])
    core.require(len(proposal) == len(ps) and all(p > 0 for p in proposal) and
                 ps[0] == 0 and all(0 <= p <= 1 for p in ps), 'Invalid witness probabilities')
    total = sum(proposal)
    return tuple(p/total for p in proposal), ps


class Checker:
    def __init__(self, spec, bands):
        self.spec, self.bands = spec, bands
        self.t, self.s, self.ac, self.kernel = core.inputs.load(spec['configuration'])
        self.caps = core.inputs.caps_module.caps()
        core.require(sorted(w for band in bands for w in band) == sorted(self.caps), 'Bands do not partition BCH caps')

    @lru_cache(maxsize=128)
    def epoch(self, log_tilt):
        lam = aa(log_tilt).exp()
        return transfer.epoch(self.t, self.s, self.ac, self.kernel, (-lam).exp())

    @lru_cache(maxsize=256)
    def moment(self, log_tilt, theta):
        lam = aa(log_tilt).exp()
        r = aa(theta)
        epoch = self.epoch(log_tilt)
        weights = [math.comb(self.t, j)*r**j*(1-r)**(self.t-j) for j in range(self.t+1)]
        matrix = tuple(sum((weights[j]*row[k] for j, row in enumerate(epoch)), arb(0)) for k in range(16))
        result = transfer.power(matrix, 256*self.spec['rows']//self.t)
        return (sum(result[:4], arb(0)).log()+self.spec['cutoff']*lam).upper()

    @lru_cache(maxsize=256)
    def costs(self, ps):
        core.require(len(ps) == len(self.bands)+1, 'Wrong number of bands')
        result = [arb(0)]
        for band, p in zip(self.bands, ps[1:]):
            if band == [256] and p == 1:
                value = F(self.caps[256])
            else:
                core.require(0 < p < 1, 'Nonconstant band needs an interior probability')
                value = max(F(self.caps[w], math.comb(256, w))/(p**w*(1-p)**(256-w)) for w in band)
            result.append(aa(value).log().upper())
        return tuple(result)

    def bound(self, box):
        rows = self.spec['rows']
        vertices = search.initial.engine.typed.vertices(box['lower'], box['upper'], rows)
        count = search.direct.lattice_count(box['lower'], box['upper'], rows)
        core.require(count > 0 and len(vertices) > 0, 'Empty type box')
        witness = box['witness']
        pi, ps = probabilities(witness)
        theta = sum(p*r for p, r in zip(pi, ps))
        core.require(0 < theta < 1, 'Degenerate moment density')
        common = self.moment(F(witness['log_surprisal']), theta)
        costs = self.costs(ps)
        terms = [g-256*aa(p).log() for g, p in zip(costs, pi)]
        values = []
        for vertex in vertices:
            multinomial = math.factorial(rows)//math.prod(math.factorial(int(n)) for n in vertex)
            value = sum((int(n)*g for n, g in zip(vertex, terms)), arb(0))-255*arb(multinomial).log()
            values.append(value.upper())
        return (common+max(values)+arb(count).log()).upper()


def clipped_cover(saved, minimum):
    rows = saved['instance']['rows']
    core.require(saved['minimum'] <= minimum <= rows, 'Requested range not in source cover')
    boxes = []
    for box in saved['leaves']:
        lo, hi = list(box['lower']), list(box['upper'])
        hi[0] = min(hi[0], rows-minimum)
        if any(a > b for a, b in zip(lo, hi)) or not search.direct.lattice_count(lo, hi, rows):
            continue
        boxes.append(dict(lower=lo, upper=hi, witness=box['witness']))
    search.direct.check_cover(boxes, rows, minimum)
    return boxes


def run(source, output, minimum=None, verify=False):
    receipt = base.read(output) if verify else None
    if receipt:
        core.authenticate(receipt, base.ROOT)
        source = base.ROOT/receipt['source_cover']
        minimum = receipt['minimum']
        core.require(base.sha(source) == receipt['source_cover_sha256'], 'Source cover changed')
    else:
        core.require(source is not None and not output.exists(), 'Provide a source and fresh output')
    saved = base.read(source)
    core.authenticate(saved, base.ROOT)
    spec = saved['instance']
    core.require(spec == core.instance(spec['configuration'], spec['message_exponent']), 'Wrong instance')
    minimum = saved['minimum'] if minimum is None else minimum
    boxes = clipped_cover(saved, minimum)
    ctx.prec = 512 if verify else 256
    checker = Checker(spec, saved['bands'])
    powers = []
    for index, box in enumerate(boxes):
        log_bound = checker.bound(box)
        bits = rational((log_bound/arb(2).log()).upper())
        power = max(-200, -(-bits.numerator//bits.denominator))
        if receipt:
            core.require(index < len(receipt['upper_powers']) and power <= receipt['upper_powers'][index], 'Box replay failed')
            power = receipt['upper_powers'][index]
        powers.append(power)
        if index % 25 == 0:
            print('box', index+1, '/', len(boxes), 'retained exponent', power, flush=True)
    total = sum((F(2)**p for p in powers), F(0))
    status = 'DENSE_TYPE_RANGE_CERTIFIED_50_BITS' if total <= F(1, 1 << 50) else 'DENSE_TYPE_RANGE_WEAK_UPPER_BOUND'
    if receipt:
        core.require(len(powers) == len(receipt['upper_powers']) and total == base.decode(receipt['union_upper']) and
                     status == receipt['status'], 'Union replay mismatch')
        base.write_new(output.with_name(output.stem+'_replay.json'), dict(status='DENSE_TYPE_512_BIT_REPLAY_PASSED',
            producer_sha256=base.sha(output), source_sha256=search.sources()))
    else:
        base.write_new(output, dict(status=status, instance=spec, minimum=minimum,
            source_cover=source.resolve().relative_to(base.ROOT).as_posix(), source_cover_sha256=base.sha(source),
            upper_powers=powers, union_upper=base.encode(total), source_sha256=search.sources(), full_distance_proved=False))
    print(status, 'margin', math.log2(total.denominator)-math.log2(total.numerator), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--minimum', type=int)
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    run(a.source, a.output, a.minimum, a.verify)
