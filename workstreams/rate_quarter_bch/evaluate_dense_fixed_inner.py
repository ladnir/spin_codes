"""Bound all remaining occupations with the unchanged t128_s19 inner.

Uses the existing typed-box inequality with separate zero, ordinary, and
all-one row counts. Numerical search is binary64, not an outward certificate.
"""
import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
import evaluate_fixed_inner as fixed
import typed_dense_boxes as typed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exponent', type=int, default=20)
    parser.add_argument('--distance', type=Fraction, default=Fraction(19, 100))
    parser.add_argument('--minimum', type=int, default=65)
    parser.add_argument('--nodes', type=int, default=127)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    a, kernel, sources = fixed.load_inner()
    spectrum = fixed.outer.read_spectrum(fixed.HERE / 'BCH256_64.wd')
    fixed.outer.audit_spectrum(spectrum, 64, 62)
    counts = {w: n for w, n in enumerate(spectrum) if w and n}
    geo = fixed.geometry(args.exponent, args.distance)
    # The ordinary band is dominated by one optimized Bernoulli measure.
    # The all-one word is an exact Bernoulli(1) atom, never put in that band.
    bands = [sorted(w for w in counts if w != 256), [256]]
    evaluator = typed.TypedDense(counts, 256, 128, 19, a, kernel, geo['outer_rows'],
        np.arange(-12., 1.01, .5),
        probability_scales=(1.,), bands=bands)
    evaluator.cutoff = geo['bad_weight']  # base driver defaults to distance 1/10
    result = evaluator.search(args.minimum, args.nodes, target_bits=60)
    sources += [Path(__file__), Path(fixed.__file__), Path(typed.__file__),
        Path(typed.general.__file__), Path(typed.dense.__file__),
        Path(fixed.q1.__file__), Path(fixed.composition.__file__), Path(fixed.outer.__file__),
        Path(fixed.maps.__file__), fixed.HERE/'BCH256_64.wd']
    payload = dict(status='BINARY64_ALL_DENSE_OCCUPATIONS_DIAGNOSTIC', geometry=geo,
        inner='existing optimized t128_s19', inner_reoptimized=False,
        margin_bits=-float(result['log_union_upper'])/math.log(2),
        log_surprisals=list(evaluator.tilts), result=result,
        source_sha256={p.relative_to(fixed.ROOT).as_posix(): fixed.sha(p) for p in sources})
    args.output.write_text(json.dumps(payload, indent=2, default=lambda x: x.tolist())+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(dict(geometry=geo, margin_bits=payload['margin_bits'], nodes=result['nodes_evaluated'])), flush=True)


if __name__ == '__main__':
    main()
