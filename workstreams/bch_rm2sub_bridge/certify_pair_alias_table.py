"""One outward Fourier transform certifies 73600 actual shared-region types.

Only compact per-overlap maxima are retained. The large arrays are temporary.
Input evaluation, complex squaring, twiddles, butterflies, and normalization
all have explicit absolute error enclosures. Aliases remain upper bounds.
"""
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
from screen_pair_type_bound import load_pairs
from certify_pair_type_fourier_tracked import roots


def run(verify=False):
    output = base.HERE/'generated/pair_alias_central_outward.json'
    old = base.read(output) if verify else None
    if old:
        for name, digest in old['local_sha256'].items():
            assert base.sha(base.HERE/name) == digest
    else:
        assert not output.exists()
    ctx.prec = 512 if verify else 256
    source = base.HERE/'generated/pair_type_central_tracked_outward.json'
    saved = base.read(source)
    for name, digest in saved['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    ps = [base.decode(v) for v in saved['tilt_probabilities']]
    p = np.array([float(v) for v in ps])
    assert sum(ps) == 1 and all(F.from_float(float(v)) == exact for v, exact in zip(p, ps))
    powers, counts = load_pairs()
    grid = (256, 256, 128)
    total = math.prod(grid)
    values = np.empty(total, dtype=np.complex128)
    errors = np.empty(total)
    tables = [roots(n) for n in grid]
    for start in range(0, total, 2048):
        indices = np.arange(start, min(start+2048, total))
        i = indices//(grid[1]*grid[2])
        j = indices//grid[2] % grid[1]
        k = indices % grid[2]
        inputs = np.vstack([np.full(len(indices), p[0], dtype=np.complex128),
                            p[1]*tables[0][i], p[2]*tables[1][j], p[3]*tables[2][k]])
        z, error = tracked.evaluate(inputs, powers, counts)
        z *= 2**30
        error = fft.up(error*2**30)
        for _ in range(6):
            z, error = fft.square(z, error)
        values[start:start+len(indices)] = z
        errors[start:start+len(indices)] = error
        if (start+len(indices)) % (total//4) == 0:
            print('Outward grid', (start+len(indices))*100//total, 'percent complete', flush=True)
    aliases, radii = fft.transform(values.reshape(grid), errors.reshape(grid))
    del values, errors
    assert np.all(np.isfinite(aliases)) and np.all(np.isfinite(radii))
    assert np.all(np.abs(aliases.imag) <= radii)
    weights = list(range(780, 859, 2))
    spectrum = tight.kernel_spectrum()
    kernel = arb_poly([spectrum.get(j, 0) for j in range(129)])**64
    singles = {j:(kernel[j]/math.comb(8192, j)).log() for j in weights}
    needed = {8192}
    for overlap in range(60, 106):
        needed.add(overlap)
        needed.update(y-overlap for y in weights)
        needed.update(8192-x-y+overlap for x in weights for y in weights)
    factorials = {n:arb(math.factorial(n)).log() for n in needed}
    logp = [(arb(v.numerator)/v.denominator).log() for v in ps]
    constant = -1920*arb(2).log()-factorials[8192]
    from audit_bch_q1_full_arb import rational
    rows = []
    max_relative_error = 0.
    for overlap in range(60, 106):
        maximum, worst = None, None
        for x in weights:
            for y in weights:
                m = [8192-x-y+overlap, y-overlap, x-overlap, overlap]
                index = (m[1] % grid[0], m[2] % grid[1], m[3] % grid[2])
                value, radius = aliases[index].real, radii[index]
                assert value > 0
                atom_upper = float(fft.up(value+radius))
                logratio = (arb(atom_upper).log()+constant-singles[x]-singles[y]+
                            sum(factorials[n]-n*lp for n, lp in zip(m, logp))).upper()
                if maximum is None or logratio > maximum:
                    maximum, worst = logratio, m
                max_relative_error = max(max_relative_error, float(radius/value))
        bound = rational(maximum.exp().upper())
        # Dyadic outward rounding leaves ample replay headroom for tiny
        # differences between root enclosures at256 and512 bits.
        rounded = F(math.ceil(bound*(1 << 30)), 1 << 30)
        assert rounded < F(41, 40)
        row = dict(overlap=overlap, maximum_ratio_upper=base.encode(rounded),
                   diagnostic_maximum_type=worst)
        if old:
            assert bound <= base.decode(old['rows'][overlap-60]['maximum_ratio_upper'])
        rows.append(row)
    maximum = max(base.decode(row['maximum_ratio_upper']) for row in rows)
    if old:
        print('512-bit root/final replay passed all73600 shared-region types', flush=True)
        return
    local = [Path(__file__), Path(fft.__file__), Path(tracked.__file__), source,
             base.HERE/'screen_pair_type_bound.py', base.HERE/'certify_pair_type_fourier_tracked.py',
             base.HERE/'generated/t128_s15_pair_spectrum.json', Path(tight.__file__),
             base.HERE/'inputs/t128_s15_b_kernel_spectrum.json']
    base.write_new(output, dict(status='OUTWARD_SHARED_REGION_PAIR_ALIAS_TABLE', configuration=tight.NAME,
        grid=list(grid), first_weights=weights, second_weights=weights, overlaps=list(range(60, 106)),
        type_count=len(weights)**2*46, rows=rows, maximum_ratio_upper=base.encode(maximum),
        all_ratios_below_41_over_40=True, diagnostic_max_relative_fft_error=max_relative_error,
        full_second_moment_certified=False,
        local_sha256={str(path.relative_to(base.HERE)):base.sha(path) for path in local}))
    print('Outward alias table:', len(weights)**2*46, 'types; maximum ratio <=', float(maximum),
          'maximum relative FFT error', max_relative_error, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    run(parser.parse_args().verify)
