"""Screen many actual pair types with one Fourier coefficient-alias table.

The floating FFT is diagnostic only. Nonnegative coefficient aliases give
upper bounds mathematically, but this screen does not enclose roundoff.
"""
import argparse
import math
from pathlib import Path
import numpy as np
from scipy.special import gammaln
from flint import arb, arb_poly, ctx
import bridge as base
import tightened_occupancy as tight
from screen_pair_type_bound import load_pairs
from certify_pair_type_fourier_tracked import evaluate


def run(grid, tag):
    output = base.HERE/'generated'/f'pair_alias_{tag}_screen.json'
    assert not output.exists()
    saved_path = base.HERE/'generated/pair_type_central_tracked_outward.json'
    saved = base.read(saved_path)
    for name, digest in saved['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    p = np.array([float(base.decode(v)) for v in saved['tilt_probabilities']])
    powers, counts = load_pairs()
    total = math.prod(grid)
    values = np.empty(total, dtype=np.complex128)
    roots = [np.exp(2j*np.pi*np.arange(n)/n) for n in grid]
    for start in range(0, total, 2048):
        indices = np.arange(start, min(start+2048, total))
        i = indices//(grid[1]*grid[2])
        j = indices//grid[2] % grid[1]
        k = indices % grid[2]
        inputs = np.vstack([np.full(len(indices), p[0], dtype=np.complex128),
                            p[1]*roots[0][i], p[2]*roots[1][j], p[3]*roots[2][k]])
        z = evaluate(inputs, powers, counts)*2**30
        for _ in range(6):
            z = z*z
        values[start:start+len(indices)] = z
    aliases = np.fft.fftn(values.reshape(grid))/total
    ctx.prec = 192
    kernel = tight.kernel_spectrum()
    poly = arb_poly([kernel.get(j, 0) for j in range(129)])**64
    weights = range(780, 859, 2)
    singles = {j:float((poly[j]/math.comb(8192, j)).log()) for j in weights}
    x, y = np.meshgrid(np.array(list(weights)), np.array(list(weights)), indexing='ij')
    pair_single = np.array([[singles[int(a)]+singles[int(b)] for b in weights] for a in weights])
    rows = []
    for overlap in range(60, 106):
        m = [8192-x-y+overlap, y-overlap, x-overlap, np.full(x.shape, overlap)]
        atom = aliases[m[1] % grid[0], m[2] % grid[1], m[3] % grid[2]]
        assert np.all(atom.real > 0)
        logden = gammaln(8193)-sum(gammaln(v+1) for v in m)+sum(v*math.log(prob) for v, prob in zip(m, p))
        ratio_bits = (np.log(atom.real)-1920*math.log(2)-logden-pair_single)/math.log(2)
        index = np.unravel_index(int(np.argmax(ratio_bits)), x.shape)
        rows.append(dict(overlap=overlap, maximum_ratio_bits=float(ratio_bits[index]),
                         maximum_type=[int(v[index]) for v in m],
                         max_relative_imaginary=float(np.max(np.abs(atom.imag)/atom.real))))
    maximum = max(rows, key=lambda row:row['maximum_ratio_bits'])
    result = dict(status='FLOATING_FOURIER_ALIAS_TABLE_SCREEN_ONLY', grid=grid,
                  first_weights=list(weights), second_weights=list(weights), overlaps=list(range(60, 106)),
                  type_count=len(weights)**2*46, rows=rows, maximum=maximum,
                  local_sha256={str(p.relative_to(base.HERE)):base.sha(p)
                                for p in [Path(__file__), saved_path, base.HERE/'screen_pair_type_bound.py',
                                          base.HERE/'certify_pair_type_fourier_tracked.py', Path(tight.__file__)]})
    base.write_new(output, result)
    print('Alias table screen:', result['type_count'], 'types; maximum ratio',
          2**maximum['maximum_ratio_bits'], 'at', maximum['maximum_type'],
          'max relative imaginary', max(row['max_relative_imaginary'] for row in rows), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--grid', type=int, nargs=3, default=[128, 128, 64])
    parser.add_argument('--tag', required=True)
    args = parser.parse_args()
    run(args.grid, args.tag)
