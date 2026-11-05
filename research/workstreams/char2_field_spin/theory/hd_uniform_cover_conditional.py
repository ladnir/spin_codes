"""Diagnostic only: suppose every H[:,A] has a uniform information cover.

The cover hypothesis is NOT proved. This checks its numerical usefulness
before investing in a universal algebraic checker. Cached regional atoms
are FLOAT proposals, not directed bounds. No output is a certificate.

For |A|=a the hypothetical Holder bound is product_i h_a(t_i), where
h_a(t)=(1+255*t**(16/a))**(a/16). Keep the two boundary states. Compute
the nonconstant part of each column by positive flag propagation, avoiding
subtraction of nearly equal floating-point matrices. Enforce at least
nine nonzero columns per active group. All labelled A have the same bound;
their binomial multiplicity is counted once outside the word product.
"""
from fractions import Fraction
from hashlib import sha256
import json
from math import comb, expm1, log, log1p, log2
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from hd_bit_syndrome_gate import ENDPOINTS
from hd_post_refresh_quiet_rows import log_matmul, log_matrix_power, moment_atoms

DIRECTORY = Path(__file__).resolve().parent
INPUT = 'hd_mds_exact_empty_mixed789_9927_pinned_float.json'
INPUT_SHA = 'c89ef4973de29007fe293817d2e3b9f582246a45d19ca19c7fcce2554e7d33b8'


def column_parts(measure, arity, groups):
    if not 1 <= arity <= 16 or groups < 1:
        raise ValueError('Arity1..16 and positive group count required')
    masses, baseline = moment_atoms(measure)
    atoms = (-np.inf, measure['log2_atom00'], measure['log2_atom01'],
             measure['log2_atom10'], 0.)
    c = arity / 16
    excess = np.full((2, 2), -np.inf)
    for (i, j), mass, atom in zip(ENDPOINTS, masses, atoms):
        if not np.isfinite(atom):
            continue
        # Natural log of h_a(atom)^groups. No cancellation at a small atom.
        exponent = groups * c * log1p(255 * 2.**(atom/c))
        log_excess = (exponent + log1p(-np.exp(-exponent)))/log(2) if exponent > 1 else log2(expm1(exponent))
        excess[i, j] = np.logaddexp2(excess[i, j], mass+log_excess)
    total = np.logaddexp2(baseline, excess)
    constant = np.full((2, 2), -np.inf)
    constant[0, 0] = constant[1, 1] = 0.
    nonconstant = np.full((2, 2), -np.inf)
    for _ in range(16):
        nonconstant = np.logaddexp2(log_matmul(nonconstant, total),
                                   log_matmul(constant, excess))
        constant = log_matmul(constant, baseline)
    return constant, nonconstant


def norm(parts, arity, groups, log2_marker):
    constant, nonconstant = parts
    marked = np.logaddexp2(constant, nonconstant+groups*log2_marker)
    word = log_matrix_power(marked, 16)
    value = np.logaddexp2.reduce(word[0])
    return (value + groups*log2(comb(16, arity))
            -8*arity*groups*log2(255) -9*groups*log2_marker)


def run():
    raw = (DIRECTORY/INPUT).read_bytes()
    if sha256(raw).hexdigest() != INPUT_SHA:
        raise ValueError('Frozen input changed')
    saved = json.loads(raw)
    groups = saved['groups']
    prefix = saved['common_prefix_log2']
    rows = []
    for arity in range(1, 17):
        parts = column_parts(saved['state_measure'], arity, groups)
        objective = lambda y: norm(parts, arity, groups, y)
        optimized = minimize_scalar(objective, bounds=(0., 16.), method='bounded',
                                    options={'xatol': 1e-9})
        options = [(objective(y), y) for y in (0., 16., float(optimized.x))]
        best, marker = min(options)
        rows.append(dict(arity=arity, log2_marker=marker,
                         unmarked_log2=prefix+objective(0.), conditional_log2=prefix+best))
    mixed = prefix+groups*np.logaddexp2.reduce(
        [(row['conditional_log2']-prefix)/groups for row in rows])
    return dict(hypothesis='Each H[:,A] admits the uniform |A|/16 information cover',
        hypothesis_proved=False, directed=False, whole_code_certificate=False,
        groups=groups, source_class=saved['source_class'], mark=saved['mark'],
        conditional_selected_class_only=True, arities=rows, conditional_mixed_log2=mixed,
        input_receipt=INPUT, input_sha256=INPUT_SHA,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        regional_provider_runs=0, native_encoder_changed=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite saved evidence')
    result = run()
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(result, indent=2), flush=True)
