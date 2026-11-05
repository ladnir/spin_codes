"""Bounded rank-one moment diagnostics for the eighteen cheapest pair survivors.

Selection preserves a local two-input minimum of six and adds nine XORs.
Neither criterion proves a whole-code distance bound. The complement screen
computes moment LOWER bounds, not upper bounds or a certificate.
"""
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import hd_lambda_multi_transvection as core
import hd_lambda_multi_transvection_pairs as pairs

DIRECTORY = Path(__file__).resolve().parent
PAIR_RECEIPT = DIRECTORY / 'hd_lambda_multi_transvection_pairs_complete.json'
PAIR_SHA256 = 'b6380e17da2effa4954dd8a9d15ccbd9b7a8dd9ea14834131860867211b4829c'


def run(max_seconds=60):
    if type(max_seconds) is not int or not 1 <= max_seconds <= 120:
        raise ValueError('Use a bounded budget from1 to120 seconds')
    pins = pairs.dependencies()
    raw = PAIR_RECEIPT.read_bytes()
    if sha256(raw).hexdigest() != PAIR_SHA256:
        raise ValueError('Observed complete pair receipt changed')
    saved = json.loads(raw)
    if (not saved['complete'] or saved['considered_pairs'] != 9730
            or saved['source_sha256'] != pins):
        raise ValueError('Complete frozen pair screen required')
    candidates = sorted(core.configuration(row['supports'])
        for row in saved['local_minimum6_survivors']
        if row['circuit']['extra_xors_over_Hlambda'] == 9)
    if len(candidates) != 18 or len(set(candidates)) != 18:
        raise ValueError('Expected all eighteen cheapest survivors')
    for candidate in candidates:
        if (pairs.exact_pair_minimum(candidate)['minimum_output_weight'] != 6
                or core.primitive_counts(candidate)['extra_xors_over_Hlambda'] != 9):
            raise ValueError('Candidate selection recheck failed')
    started, results = perf_counter(), []
    # Include both controls, but never transfer their proof receipts.
    plan = [(), ((0, 1, 2, 3),)] + candidates
    for candidate in plan:
        remaining = max_seconds - (perf_counter() - started)
        if remaining <= 0:
            break
        row = core.complement_screen(candidate, max_seconds=min(30, remaining))
        results.append(row)
        print(json.dumps(dict(completed=len(results), expected=len(plan),
            supports=candidate, complete=row['complete_all4368'],
            maximum_projective_class=row['maximum_projective_class'],
            lower_I4=str(row['maximum_zero_plus_rank1_moment_parts']['4']))), flush=True)
        if not row['complete_all4368']:
            break
    if pairs.dependencies() != pins or sha256(PAIR_RECEIPT.read_bytes()).hexdigest() != PAIR_SHA256:
        raise ValueError('Inputs changed during diagnostic')
    ranked = sorted((row for row in results[2:] if row['complete_all4368']),
        key=lambda row: Fraction(row['maximum_zero_plus_rank1_moment_parts']['4']))
    return dict(mode='cheapest-six-minimum-two-stage-rank1-diagnostic',
        expected_candidates=18, expected_controls=2, results=results,
        complete=len(results) == 20 and all(row['complete_all4368'] for row in results),
        ranked_completed_candidates=[dict(supports=row['supports'],
            maximum_projective_class=row['maximum_projective_class'],
            lower_I4=row['maximum_zero_plus_rank1_moment_parts']['4']) for row in ranked],
        pair_receipt_sha256=PAIR_SHA256, source_sha256=pins,
        controller_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        all_bounds_are_lower=True, whole_code_certificate=False,
        full_support_tables=0, regional_provider_runs=0, native_benchmarks=0,
        elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-seconds', type=int, default=60)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = run(args.max_seconds)
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items()
        if k not in ('results', 'ranked_completed_candidates', 'source_sha256')}))
