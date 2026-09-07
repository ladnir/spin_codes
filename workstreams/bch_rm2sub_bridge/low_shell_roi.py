"""Exact low-shell sensitivity and screening at the frozen K28 endpoint.

Hypothetical shell caps and moment zeros are NOT new code facts. All optimization
here is a three-item fractional knapsack, solved over exact rational numbers.
"""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
import sys
from flint import fmpq
import bridge as base

sys.path.insert(0, str(base.BCH/'code'))
from verify_scaled_rational_solution import normalized_rows, parse_assignments
from audit_bch_closure_envelope import kraw_table


def fraction(x):
    return F(int(x.numerator), int(x.denominator))


def fq(x):
    x = F(x)
    return fmpq(x.numerator, x.denominator)


def margin(x):
    return math.log2(x.denominator)-math.log2(x.numerator)


def knapsack(cost, value, caps, budget):
    """Maximize sum value[w]*x[w], 0<=x<=caps, sum cost[w]*x[w]<=budget."""
    assert set(cost) == set(value) == set(caps)
    assert budget >= 0 and all(cost[w] > 0 and value[w] >= 0 and caps[w] >= 0 for w in cost)
    left = budget
    allocation = {}
    for w in sorted(cost, key=lambda w: value[w]/cost[w], reverse=True):
        allocation[w] = min(caps[w], left/cost[w])
        left -= cost[w]*allocation[w]
    score = sum((value[w]*allocation[w] for w in cost), F(0))
    # Exact dual witness: lambda for the budget, nonnegative prices for caps.
    slopes = [F(0)]+[value[w]/cost[w] for w in cost]
    dual = min(lam*budget+sum((caps[w]*max(F(0), value[w]-lam*cost[w]) for w in cost), F(0)) for lam in slopes)
    assert score == dual and left >= 0
    return score, allocation


def replay_primal(folder):
    model = base.read(folder/'model.json')
    scales = base.read(folder/'scales.json')
    variables = model['metadata']['variables']
    solution = (folder/'h_38.sol').read_text()
    assert 'status = OPTIMAL' in solution
    raw = parse_assignments(solution, 'VARS:', 'REDUCED COST:')
    assert set(raw) <= set(variables)
    x = {v: fq(raw.get(v, F(0))) for v in variables}
    assert all(v >= 0 for v in x.values())
    rows = normalized_rows(model, scales)
    for row in rows:
        lhs = sum((fq(a)*x[v] for v, a in row['coeffs'].items()), fmpq(0))
        rhs = fq(row['rhs'])
        assert {'eq': lhs == rhs, 'le': lhs <= rhs, 'ge': lhs >= rhs}[row['sense']], row['name']
    # Replay the old weighted dual, not merely its stored success flag.
    objective = base.read(folder/'objective.json')
    prices = parse_assignments(solution, 'PI:', 'SLACK:')
    assert set(prices) <= {f'c{i}' for i in range(1, len(rows)+1)}
    combined = {v: fmpq(0) for v in variables}
    bound = fmpq(0)
    for i, row in enumerate(rows, 1):
        price = fq(prices.get(f'c{i}', F(0)))
        assert row['sense'] == 'eq' or (price >= 0 if row['sense'] == 'le' else price <= 0)
        bound += price*fq(row['rhs'])
        for v, a in row['coeffs'].items():
            combined[v] += price*fq(a)
    target = {v: fq(base.decode(a)) for v, a in objective['scaled_objective'].items()}
    assert all(combined[v] >= target.get(v, fmpq(0)) for v in variables)
    assert bound == sum((a*x[v] for v, a in target.items()), fmpq(0))
    audit = base.read(folder/'audit.json')
    assert fraction(bound)*base.decode(objective['physical_objective_normalization']) == base.decode(audit['paired_shells_upper'])
    return model, {v: fraction(a)*scales[v] for v, a in x.items()}, len(rows)


def build():
    folder = base.BCH/'generated/shift_rank_oa29_joint'
    model, x, checked = replay_primal(folder)
    q1path = base.HERE/'generated/frontier_k28_q1_v1.json'
    fullpath = base.HERE/'generated/curve_k28_full_retained_v1.json'
    q1, full = base.read(q1path), base.read(fullpath)
    assert q1['message_exponent'] == full['message_exponent'] == 28
    co = {int(w): base.decode(v) for w, v in q1['coefficient_upper'].items()}
    upper, factor, rest = base.bch_bound(co)
    higher = base.decode(full['higher_occupancy_upper'])
    assert upper+higher == base.decode(full['failure_upper'])
    objective = base.read(folder/'objective.json')
    audit = base.read(folder/'audit.json')
    cost = {w: base.decode(objective['dyadic_upper_pair_coefficients'][str(w)]) for w in (38, 40, 42)}
    value = {w: co[w]+co[256-w] for w in cost}
    budget = base.decode(audit['paired_shells_upper'])
    primal_counts = {w: x[f'q_{w}']+31*x[f'h_{w}'] for w in cost}
    lower = sum((primal_counts[w]*value[w] for w in cost), F(0))
    caps, cap_paths = {}, []
    for w in cost:
        cap_path = (base.BCH/'generated/shift_rank_oa29_hull_probe/audit.json' if w == 38
                    else base.HERE/f'generated/joint_shell_exact_w{w}/cap.json')
        record = base.read(cap_path)
        caps[w] = F(record['cap'] if w != 38 else record['A38_cap_without_lattice_rounding'])
        cap_paths.append(cap_path)
    rows = []
    for changed in ((38,), (38, 40), (38, 40, 42)):
        for bits in (0, 1, 2, 3, 4, 6, 8, 12, 20):
            trial = {w: F(caps[w]//(1 << bits)) if w in changed else caps[w] for w in caps}
            score, allocation = knapsack(cost, value, trial, budget)
            total = min(upper, score+rest)+higher
            rows.append(dict(tightened_weights=list(changed),cap_reduction_bits=bits,
                hypothetical_caps={str(w): int(v) for w, v in trial.items()},
                conditional_full_upper=base.encode(total),conditional_margin_bits=margin(total),
                allocation={str(w): base.encode(v) for w, v in allocation.items()},
                improvement_bits=margin(total)-margin(upper+higher)))
    kt = kraw_table()
    def transform(terms, degree):
        return sum((mult*(kt[degree][w]+(kt[degree][256-w] if w != 128 else 0))*v
                    for prefix, mult in terms for n, v in x.items() if n.startswith(prefix+'_')
                    for w in [int(n.split('_')[1])]), F(0))
    moments = []
    for name, terms, dim in [('Qdual', [('q', 1)], 123), ('Pdual', [('q', 1), ('h', 255)], 131),
                              ('HQdual', [('r', 1), ('s', 255)], 93)]:
        for j in ((20, 22, 24, 26, 28) if name == 'HQdual' else (30, 32, 34, 36)):
            count = transform(terms, j)/(1 << dim)
            assert count >= 0
            moments.append(dict(code=name,weight=j,primal_count=base.encode(count),
                log2_count_diagnostic=-margin(count) if count else None,
                hypothetical_zero_excludes_primal=count != 0))
    # All half-spectrum counts of the existing feasible witness; not a real code.
    terms = {w: (x.get(f'q_{min(w,256-w)}', F(0))+31*x.get(f'h_{min(w,256-w)}', F(0)))*v for w, v in co.items()}
    full_primal_score = sum(terms.values(), F(0))
    files = [Path(__file__), Path(base.__file__), q1path, fullpath,
             *[folder/n for n in ('model.json','scales.json','h_38.sol','objective.json','audit.json')],
             *cap_paths, base.BCH/'code/verify_scaled_rational_solution.py', base.BCH/'code/audit_bch_closure_envelope.py']
    return dict(status='EXACT_RELAXATION_AUDIT_AND_CONDITIONAL_CAP_SENSITIVITY_NOT_NEW_SPECTRUM_FACTS',
        message_exponent=28,configuration='t64_s20',rows_replayed=checked,
        baseline_upper=base.encode(upper+higher),baseline_margin_bits=margin(upper+higher),
        rest_upper=base.encode(rest),higher_upper=base.encode(higher),rest_fraction=float(rest/upper),
        best_margin_if_all_three_low_shells_removed=margin(rest+higher),
        fixed_coefficient_reoptimization_gain_ceiling_bits=margin(full_primal_score)-margin(upper),
        low_shell_only_reoptimization_gain_ceiling_bits=margin(lower)-margin(upper),
        exact_feasible_score=base.encode(full_primal_score),
        primal_shells=[dict(weight=w,count=base.encode(primal_counts[w]),
            log2_count_diagnostic=-margin(primal_counts[w]) if primal_counts[w] else None,
            low_shell_score_fraction=float(primal_counts[w]*value[w]/lower),
            coefficient_ratio_to_old=float(value[w]/cost[w])) for w in cost],
        baseline_caps={str(w): int(v) for w, v in caps.items()},
        sensitivity=rows,moment_screen=moments,
        source_sha256={p.relative_to(base.ROOT).as_posix():base.sha(p) for p in files})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = build()
    if args.verify:
        assert result == base.read(args.output)
    else:
        base.write_new(args.output, result)
    for key in ('baseline_margin_bits','fixed_coefficient_reoptimization_gain_ceiling_bits',
                'best_margin_if_all_three_low_shells_removed','primal_shells','moment_screen'):
        print(key, result[key], flush=True)
    for row in result['sensitivity']:
        print(row['tightened_weights'],row['cap_reduction_bits'],round(row['conditional_margin_bits'],6),flush=True)
