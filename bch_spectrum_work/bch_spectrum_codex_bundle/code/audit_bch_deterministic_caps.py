"""Recheck exact dual witnesses without overwriting any existing certificate."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

sys.dont_write_bytecode = True
from export_scaled_rational_lp import objective_spec, variable_scales
from verify_scaled_rational_solution import normalized_rows, parse_assignments
from export_constant_weight_delsarte import attach_qsopt_certificate
from run_higher_endpoint_preflight import sha, write_new

ROOT = Path(__file__).resolve().parents[1]
GEN = ROOT / 'generated'


def main():
    output = GEN / 'bch256_deterministic_caps_reaudit.json'
    assert not output.exists()
    model_path = GEN / 'coupled_lp_exact.json'
    model = json.loads(model_path.read_text())
    scales = variable_scales(model)
    rows = normalized_rows(model, scales)
    cap_path = GEN / 'exact_bch_lp_caps.json'
    stored_caps = {r['weight']: r['exact_Aw_C_cap']
                   for r in json.loads(cap_path.read_text())['shells']}
    paths = [model_path, cap_path, Path(__file__)]
    results = {}
    for weight in range(38, 52, 2):
        objective = 'h_38' if weight == 38 else f'c_{weight}'
        base = GEN / f'max_{objective}_scaled_rational'
        solution = base.with_suffix('.sol')
        certificate_path = base.with_suffix('.certificate.json')
        receipt = json.loads(certificate_path.read_text())
        coefficients, multiplier = objective_spec(objective, model['metadata']['variables'], scales)
        assert receipt['objective_scaled_coefficients'] == coefficients
        assert receipt['objective_physical_multiplier'] == multiplier
        prices = parse_assignments(solution.read_text(), 'PI:', 'SLACK:')
        dual = [prices.get(f'c{i}', Fraction(0)) for i in range(1, len(rows) + 1)]
        assert set(prices) <= {f'c{i}' for i in range(1, len(rows) + 1)}
        for y, row in zip(dual, rows):
            assert row['sense'] == 'eq' or (y >= 0 if row['sense'] == 'le' else y <= 0)
        for name in model['metadata']['variables']:
            lhs = sum((y * row['coeffs'].get(name, 0) for y, row in zip(dual, rows)), Fraction(0))
            assert lhs >= coefficients.get(name, 0)
        physical_upper = multiplier * sum((y * row['rhs'] for y, row in zip(dual, rows)), Fraction(0))
        assert physical_upper == Fraction(int(receipt['physical_optimum']['numerator']),
                                          int(receipt['physical_optimum']['denominator']))
        cap = physical_upper.numerator // physical_upper.denominator
        if weight == 38:
            cap = 31 * (cap - cap % 128)
        assert cap == stored_caps[weight]
        assert sha(solution) == receipt['sha256']['solution']
        assert sha(base.with_suffix('.lp')) == receipt['sha256']['lp']
        results[str(weight)] = dict(cap=cap, dual_rows=len(rows),
                                    dual_columns=len(model['metadata']['variables']))
        paths.extend([solution, certificate_path, base.with_suffix('.lp')])
    for weight in (52, 54):
        path = GEN / f'johnson_n256_w{weight}_d38.json'
        stored = json.loads(path.read_text())
        assert (stored['n'], stored['weight'], stored['minimum_distance']) == (256, weight, 38)
        solution = path.with_suffix('.sol')
        actual = attach_qsopt_certificate(stored.copy(), 256, weight, 19, solution)
        assert actual == stored
        results[str(weight)] = dict(cap=actual['integer_shell_upper'],
                                    exact_primal_and_dual_rechecked=True)
        paths.extend([path, solution])
    paths.extend(ROOT / 'code' / name for name in (
        'export_lp.py', 'export_scaled_rational_lp.py', 'verify_scaled_rational_solution.py',
        'export_constant_weight_delsarte.py', 'run_higher_endpoint_preflight.py'))
    result = dict(classification='Exact rational dual-witness recheck and cap-to-application mapping',
                  checks=results,
                  limitation='Rechecks LP arithmetic against the retained BCH-sandwich model; does not rederive every algebraic or orbit constraint in that model.',
                  source_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths})
    write_new(output, result)
    print(json.dumps(result['checks'], indent=2))


if __name__ == '__main__':
    main()
