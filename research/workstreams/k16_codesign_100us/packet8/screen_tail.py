"""All-occupancy packet8 proposals from fresh outward local operators.

The regional calculation is logarithmic floating point, not a certificate.
It proposes tilts for later outward replay. No positive entry is discarded
by an absolute underflow threshold.
"""
from __future__ import annotations

import argparse
from math import comb, log
from time import monotonic

import numpy as np
from scipy.special import gammaln, logsumexp
from flint import arb, arb_mat, ctx

import screen_packet8 as base
import packet_regional_log as regional_log
import kernel_birth_density as density


def walsh(values):
    result = np.asarray(values, dtype=np.int64).copy()
    length = 1
    while length < len(result):
        blocks = result.reshape(-1, 2*length)
        left, right = blocks[:, :length].copy(), blocks[:, length:].copy()
        blocks[:, :length], blocks[:, length:] = left+right, left-right
        length *= 2
    return result


def prepare_full(data):
    b, W = data['packet_bits'], data['windows']
    images = base.images_from_rows(data['rows'])
    profiles = [(W, *([0]*b)), *data['profiles']]
    lookup = {profile: i for i, profile in enumerate(profiles)}
    indices, weights = [], []
    mask = (1 << b)-1
    for image in images:
        profile = [0]*(b+1)
        for p in range(0, data['width'], b):
            profile[((image >> p) & mask).bit_count()] += 1
        indices.append(lookup[tuple(profile)])
        weights.append(image.bit_count())
    indices, weights = np.array(indices), np.array(weights)
    census = []
    for level in (0, *data['levels']):
        transformed = walsh(weights == level)
        grouped = np.zeros(len(profiles), dtype=np.int64)
        np.add.at(grouped, indices, transformed)
        census.append(grouped)
    assert sum(census).tolist() == [len(images), *([0]*(len(profiles)-1))]
    return dict(data, character_profiles=profiles, census=np.array(census),
                records=profiles, multiplicities=np.bincount(indices, minlength=len(profiles)),
                birth_character_indices=indices, birth_state_weights=weights)


def polynomial(profile, inactive, active):
    values = [arb(1)]
    for w, count in enumerate(profile):
        for _ in range(count):
            following = [arb(0)]*(len(values)+1)
            for j, value in enumerate(values):
                following[j] += inactive[w]*value
                following[j+1] += active[w]*value
            values = following
    W = sum(profile)
    return [value/comb(W, j) for j, value in enumerate(values)]


def full_local(data, tilt):
    b, W, S = data['packet_bits'], data['windows'], 1 << data['bits']
    z = (-base.aq(base.Q(tilt))).exp()
    fourier = [((1+z)**(b-r)*(1-z)**r-1)/((1 << b)-1) for r in range(b+1)]
    character_rows = [polynomial(profile, [arb(1)]*(b+1), fourier)
                      for profile in data['character_profiles']]
    character = arb_mat(character_rows)
    births = arb_mat(data['census'].tolist())*character/S
    inactive = [z**w for w in range(b+1)]
    active = [sum(((comb(b, k)-int(k == w))*z**k for k in range(b+1)), arb(0))/((1 << b)-1)
              for w in range(b+1)]
    moments = [polynomial(profile, inactive, active) for profile in data['profiles']]
    weights = [sum(w*c for w, c in enumerate(profile)) for profile in data['profiles']]
    levels = data['levels']
    n = 2+len(levels)
    result = []
    for occupied in range(W+1):
        rows = [[arb(0)]*n for _ in range(n)]
        rows[0][0] = max(arb(0), base.up(births[0, occupied]))
        for i in range(len(levels)):
            rows[0][2+i] = max(arb(0), base.up(births[1+i, occupied]))
        mean = base.up(sum((count*row[occupied] for count, row in zip(data['profiles'].values(), moments)), arb(0))/(S-1))
        bounds = [mean, *[max(base.up(row[occupied]) for row, weight in zip(moments, weights) if weight == level)
                         for level in levels]]
        for i, bound in enumerate(bounds, 1):
            # The return requires Cx!=0. Dropping that condition is a valid upper bound.
            rows[i][0] = base.up(bound/(S-1)) if occupied else arb(0)
            rows[i][1] = bound
        # Injectivity gives exact structural zeros at occupancies0 and1.
        if occupied == 0:
            rows[0] = [arb(1), *([arb(0)]*(n-1))]
        if occupied == 1:
            rows[0][0] = arb(0)
        result.append(arb_mat([[base.up(value) for value in row] for row in rows]))
    # Each cap gives a complete representation of the same birth measure.
    # Select whole rows, never entrywise minima across different caps.
    laws, denominators = density.upper_laws(data, character_rows)
    arrays = np.array([[[float(m[i, j]) for j in range(n)] for i in range(n)] for m in result])
    potential = density.continuation(sum(comb(W, j)*.5**W*arrays[j] for j in range(W+1)))
    state_weights = data['birth_state_weights'][1:]
    for j, matrix in enumerate(result):
        baseline = [matrix[0, k] for k in range(n)]
        options = [baseline]
        for cap in density.thresholds(laws[1:, j]):
            row = baseline[:]
            row[1] = base.up(arb((S-1)*int(cap))/denominators[j])
            residual = np.maximum(0, laws[1:, j]-int(cap))
            for k, level in enumerate(levels, 2):
                row[k] = min(baseline[k], base.up(arb(int(residual[state_weights == level].sum()))/denominators[j]))
            options.append(row)
        selected = min(options, key=lambda row: sum(float(x)*v for x, v in zip(row, potential)))
        for k, value in enumerate(selected):
            matrix[0, k] = value
    return result


def mixture(regional, q):
    js = np.arange(q+1)
    weights = gammaln(q+1)-gammaln(js+1)-gammaln(q-js+1)+js*log(255)-q*log(256)
    return logsumexp(regional[:q+1]+weights[:, None, None], axis=0)


def run(tilts, q_min=2, q_max=512):
    ctx.prec = 256
    data, _ = base.construction()
    data = prepare_full(data)
    best = np.full(q_max+1, -np.inf)
    choices = {}
    # Unchanged pointwise MDS outer majorant, with the8-bit packet geometry.
    beta_log = 256*log(2)-8*log(65535)
    start = monotonic()
    for tilt in tilts:
        local = full_local(data, tilt)
        regional = regional_log.log_placement(regional_log.upper_logs(local), q_max, epochs=64, windows=8)
        for q in range(q_min, q_max+1):
            moment = regional_log.log_power_matrix(mixture(regional, q), 32)
            margin = -(log(comb(512, q))+q*beta_log+float(tilt)*13107+moment)/log(2)
            if margin > best[q]:
                best[q], choices[q] = margin, tilt
        worst = sorted(range(q_min, q_max+1), key=lambda q: best[q])[:10]
        print(f'PROPOSAL tilt={tilt} elapsed={monotonic()-start:.2f}s worst={[(q,round(best[q],6),choices[q]) for q in worst]}', flush=True)
    print(f'PROPOSAL union margin={-logsumexp(-best[q_min:]*log(2))/log(2):.10f}; whole_code_certificate=false', flush=True)
    return best, choices


def outward_point(occupancy, tilt):
    """Fresh rigorous component endpoint, still not a whole-code certificate."""
    if not 1 <= occupancy <= 512:
        raise ValueError('occupancy from1 through512 required')
    ctx.prec = 256
    data, _ = base.construction()
    local = full_local(prepare_full(data), tilt)
    regional = base.retained.placement(local, epochs=64, windows=8,
                   maximum_groups=occupancy, rounding=base.retained.rounded)
    n = local[0].nrows()
    mixed = sum((regional[j]*base.aq(base.Q(comb(occupancy, j)*255**j, 256**occupancy))
                 for j in range(occupancy+1)), arb_mat(n, n))
    power = mixed**32
    moment = sum((power[0, j] for j in range(n)), arb(0))
    beta = base.aq(base.Q(2**256, 65535**8))
    upper = base.up(comb(512, occupancy)*beta**occupancy*
                    (base.aq(base.Q(tilt))*13107).exp()*moment)
    margin = -upper.log()/arb(2).log()
    print(f'OUTWARD q={occupancy} tilt={tilt} margin={margin}; whole_code_certificate=false', flush=True)
    return upper


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tilts', nargs='+', default=['.00512', '.01024', '.0256', '.0512', '.1024', '.2048', '.4096', '.8192', '1.6'])
    parser.add_argument('--q-min', type=int, default=2)
    parser.add_argument('--q-max', type=int, default=512)
    parser.add_argument('--outward-q', type=int)
    args = parser.parse_args()
    if args.outward_q is None:
        run(args.tilts, args.q_min, args.q_max)
    else:
        for tilt in args.tilts:
            outward_point(args.outward_q, tilt)


if __name__ == '__main__':
    main()
