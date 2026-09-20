"""Outward pair-kernel tiles for actual overlaps0..39, including zero cells."""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from flint import arb, arb_poly, ctx
import bridge as base
import tightened_occupancy as tight
import pair_fourier_error as tracked
import outward_radix2_fft as fft
import pair_positive_tilt as tilt
from screen_pair_type_bound import load_pairs
from certify_pair_type_fourier_tracked import roots


def run(tile, verify=False):
    assert 0 <= tile < 6
    output = base.HERE/'generated'/f'low_overlap_tile{tile:02d}.json'
    old = base.read(output) if verify else None
    if old:
        for name, digest in old['local_sha256'].items():
            assert base.sha(base.HERE/name) == digest
    else:
        assert not output.exists()
    ctx.prec = 512 if verify else 256
    ranges = [list(range(740, 819, 2)), list(range(820, 901, 2))]
    i, j = [(0, 0), (0, 1), (1, 1)][tile//2]
    xs, ys = ranges[i], ranges[j]
    ks = list(range(0, 20)) if tile % 2 == 0 else list(range(20, 40))
    if old:
        ps = list(map(base.decode, old['tilt_probabilities']))
    else:
        cx, cy, ck = 780+80*i, 780+80*j, 8 if tile % 2 == 0 else 28
        proposed, _, _ = tilt.PairTilt().fit([8192-cx-cy+ck, cy-ck, cx-ck, ck])
        denominator = 1 << 44
        integers = [round(v*denominator) for v in proposed[1:]]
        integers = [denominator-sum(integers)]+integers
        ps = [F(v, denominator) for v in integers]
    assert min(ps) > 0 and sum(ps) == 1
    p = np.array([float(v) for v in ps])
    assert all(F.from_float(float(v)) == exact for v, exact in zip(p, ps))
    powers, counts = load_pairs()
    grid, total = (256, 256, 64), 256*256*64
    table = roots(256)
    values, errors = np.empty(total, dtype=np.complex128), np.empty(total)
    for start in range(0, total, 2048):
        indices = np.arange(start, min(start+2048, total))
        a, b, c = indices//16384, indices//64 % 256, indices % 64
        inputs = np.vstack([np.full(len(indices), p[0], dtype=np.complex128),
                            p[1]*table[b], p[2]*table[a], p[3]*table[(a+b+4*c) % 256]])
        z, error = tracked.evaluate(inputs, powers, counts)
        z *= 2**30
        error = fft.up(error*2**30)
        for _ in range(6):
            z, error = fft.square(z, error)
        values[start:start+len(indices)], errors[start:start+len(indices)] = z, error
        if (start+len(indices)) % (total//4) == 0:
            print('Low-overlap tile', tile, 'grid', (start+len(indices))*100//total, 'percent', flush=True)
    aliases, radii = fft.transform(values.reshape(grid), errors.reshape(grid))
    del values, errors
    assert np.all(np.isfinite(aliases)) and np.all(np.isfinite(radii))
    assert np.all(np.abs(aliases.imag) <= radii)
    spectrum = tight.kernel_spectrum()
    kernel = arb_poly([spectrum.get(j, 0) for j in range(129)])**64
    singles = {j:(kernel[j]/math.comb(8192, j)).log() for j in set(xs+ys)}
    needed = {8192}
    for k in ks:
        needed.add(k)
        needed.update(x-k for x in xs)
        needed.update(y-k for y in ys)
        needed.update(8192-x-y+k for x in xs for y in ys)
    factorials = {n:arb(math.factorial(n)).log() for n in needed}
    logp = [(arb(v.numerator)/v.denominator).log() for v in ps]
    constant = -1920*arb(2).log()-factorials[8192]
    quadratic = F(3, 25000)
    from audit_bch_q1_full_arb import rational
    rows, maximum_relative_error = [], 0.
    for k in ks:
        maximum, worst = None, None
        for x in xs:
            for y in ys:
                m = [8192-x-y+k, y-k, x-k, k]
                index = (x % 256, y % 256, k % 64)
                value, radius = float(aliases[index].real), float(radii[index])
                atom_upper = float(fft.up(value+radius))
                assert atom_upper > 0
                delta = (F(k)-F(x*y, 8192))**2*quadratic
                adjusted = (arb(atom_upper).log()+constant-singles[x]-singles[y]+
                            sum(factorials[n]-n*lp for n, lp in zip(m, logp))-
                            arb(delta.numerator)/delta.denominator).upper()
                if maximum is None or adjusted > maximum:
                    maximum, worst = adjusted, m
                maximum_relative_error = max(maximum_relative_error, radius/value if value > 0 else math.inf)
        exact = rational(maximum)
        rounded = F(math.ceil(exact*(1 << 28)), 1 << 28)
        if old:
            assert exact <= base.decode(old['rows'][k-ks[0]]['intercept_upper'])
        rows.append(dict(overlap=k, intercept_upper=base.encode(rounded), diagnostic_worst_type=worst))
    intercept = max(base.decode(r['intercept_upper']) for r in rows)
    if old:
        print('512-bit low-overlap tile replay passed:', tile, flush=True)
        return
    local = [Path(__file__), Path(fft.__file__), Path(tracked.__file__), Path(tilt.__file__),
             base.HERE/'screen_pair_type_bound.py', base.HERE/'certify_pair_type_fourier_tracked.py',
             base.HERE/'generated/t128_s15_pair_spectrum.json', Path(tight.__file__),
             base.HERE/'inputs/t128_s15_b_kernel_spectrum.json']
    base.write_new(output, dict(status='OUTWARD_ACTUAL_PAIR_LOW_OVERLAP_TILE', configuration=tight.NAME,
        tile=tile, grid=list(grid), first_weights=xs, second_weights=ys, overlaps=ks,
        tilt_probabilities=[base.encode(v) for v in ps], type_count=len(xs)*len(ys)*len(ks),
        quadratic_coefficient=base.encode(quadratic), intercept_upper=base.encode(intercept),
        common_quadratic_envelope_passed=intercept<=0, rows=rows,
        diagnostic_maximum_relative_fft_error=maximum_relative_error, full_second_moment_certified=False,
        local_sha256={str(path.relative_to(base.HERE)):base.sha(path) for path in local}))
    print('Low-overlap tile', tile, 'intercept <=', float(intercept), 'common envelope passed:', intercept<=0,
          'relative FFT error', maximum_relative_error, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--tile', type=int, required=True)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    run(args.tile, args.verify)
