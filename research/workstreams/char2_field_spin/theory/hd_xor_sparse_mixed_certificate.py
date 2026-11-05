"""Fresh directed replay of the full r17 (one single-row,16 pair-row) cohort.

The input artifact supplies only literal Chernoff and potential witnesses.
All class counts, positive support-polynomial tails, ordered regional
operators, resource envelopes and normalization factors are recomputed.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb
from pathlib import Path
from time import perf_counter

from flint import arb, arb_poly, ctx

from hd_shifted_coverage import rational_ceiling
from hd_xor_geometry_strata import MultiplicityCapsArb
from hd_xor_joint_columns import categories
from hd_xor_pair_route import arb_pair_regions
from hd_xor_paired_joint_gate import ordered_tables
from hd_xor_sparse_mixed_refine import MU_BINS, geometry_count, metadata, source_counts
from hd_xor_multigroup import Q, LIMIT


def sparse_count_arb(caps, n1, n2, n=32, k=16, q=256):
    """Positive aggregate-support tail, with one fixed row set per group."""
    if any(type(value) is not int for value in (n1, n2, n, k, q)):
        raise ValueError('Integral code/cohort dimensions required')
    if min(n1, n2) < 0 or not 0 < k <= n or q < 2:
        raise ValueError('Invalid code/cohort dimensions')
    if set(caps) != set(categories(n1+n2)):
        raise ValueError('Exactly one cap per group support category required')
    active_rows = n1+2*n2
    coefficients = [arb(0) for _ in range(active_rows+1)]
    for one in range(n1+1):
        for single in range(n2+1):
            for double in range(n2-single+1):
                degree = one+single+2*double
                multiplicity = comb(n1, one)*comb(n2, single)*comb(n2-single, double)*2**single
                coefficients[degree] += multiplicity*caps[one+single, double]*(q-1)**degree
    expanded = arb_poly(coefficients)**n
    floor = (n-k+1)*active_rows
    return sum((expanded[degree] for degree in range(floor, n*active_rows+1)), arb(0))/(q-1)**((n-k)*active_rows)


def checked_witnesses(proposal):
    """Reject changed scope, incomplete partitions and duplicate classes."""
    for key, value in metadata().items():
        if type(proposal.get(key)) is not type(value) or proposal[key] != value:
            raise ValueError(f'Wrong sparse mixed proposal field: {key}')
    screen = proposal.get('screen', {})
    if screen.get('complete_cohort') is not True:
        raise ValueError('The proposed source partition must be complete')
    expected = {(source, bounds) for source in source_counts() for bounds in MU_BINS}
    checked = {}
    for row in screen.get('entries', []):
        source, bounds = tuple(row['source_matching_counts']), tuple(row['mu_range'])
        key = source, bounds
        if (any(type(value) is not int for value in source+bounds)
                or key not in expected or key in checked):
            raise ValueError('Invalid or duplicated source/multiplicity class')
        z, logp = F(row['mark']), F(row['log2_potential'])
        if not 0 < z < 1 or logp > 0:
            raise ValueError('Invalid literal witness')
        checked[key] = (z, logp)
    if set(checked) != expected:
        raise ValueError('The full81 by4 class partition is required')
    return checked


def run(proposal, precision, *, selected_sources=None, progress=None):
    if type(precision) is not int or precision < 128:
        raise ValueError('At least128 bits of directed precision required')
    witnesses = checked_witnesses(proposal)
    selected = set(source_counts()) if selected_sources is None else set(map(tuple, selected_sources))
    if not selected or not selected <= set(source_counts()):
        raise ValueError('Select a nonempty set of feasible source classes')
    started, entries = perf_counter(), []
    with ctx.workprec(precision):
        marks = {z for (source, _), (z, _) in witnesses.items() if source in selected}
        raw = {z: arb_pair_regions(z, 17) for z in marks}
        total = arb(0)
        for source in sorted(selected):
            evaluators = {}
            for bounds in MU_BINS:
                z, logp = witnesses[source, bounds]
                potential = arb(2)**(arb(logp.numerator)/logp.denominator)
                if (z, logp) not in evaluators:
                    evaluators[z, logp] = MultiplicityCapsArb(
                        ordered_tables(raw[z], *source, directed=True), potential, 17,
                        labelled_regions=True)
                caps = evaluators[z, logp](bounds[1])
                count = sparse_count_arb(caps, 1, 16)
                geometries = geometry_count(*bounds)
                value = (source_counts()[source]*17*8*geometries*count
                    *(arb(z.numerator)/z.denominator)**(-LIMIT)/(Q-1)/potential)
                if not value > 0 or not value.is_finite():
                    raise ArithmeticError('Invalid directed class contribution')
                row = dict(source_matching_counts=list(source), mu_range=list(bounds),
                    mark=str(z), log2_potential=str(logp), source_supports=source_counts()[source],
                    cohort_colorings=17, one_row_choices=8, two_row_geometries=geometries,
                    failure_log2_upper=str(rational_ceiling(value.log()/arb(2).log())))
                entries.append(row)
                total += value
                if progress:
                    progress(dict(completed=len(entries), total=4*len(selected), **row))
        ceiling = rational_ceiling(total.log()/arb(2).log())
    result = metadata()
    result.update(directed=True, entries=entries, precision_bits=precision,
        complete_cohort=selected == set(source_counts()), source_classes=len(selected),
        all_source_classes=len(source_counts()), fitted_new_witnesses=0,
        proposal_requires_fresh_directed_replay=False,
        failure_log2_upper=str(ceiling), elapsed_seconds=perf_counter()-started)
    return result


def verify_receipts(proposal, receipts):
    """Check both complete replays and independently sum rounded cells at512."""
    if len(receipts) != 2:
        raise ValueError('Two independently produced receipts required')
    witnesses = checked_witnesses(proposal)
    checked = []
    for precision, receipt in zip((256, 384), receipts):
        expected = metadata()
        expected.update(directed=True, precision_bits=precision, complete_cohort=True,
            source_classes=81, all_source_classes=81, fitted_new_witnesses=0,
            proposal_requires_fresh_directed_replay=False)
        for key, value in expected.items():
            if type(receipt.get(key)) is not type(value) or receipt[key] != value:
                raise ValueError(f'Wrong directed receipt field: {key}')
        rows = {}
        for row in receipt['entries']:
            source, bounds = tuple(row['source_matching_counts']), tuple(row['mu_range'])
            key = source, bounds
            if (any(type(value) is not int for value in source+bounds)
                    or key not in witnesses or key in rows):
                raise ValueError('Invalid or duplicate receipt class')
            if (F(row['mark']), F(row['log2_potential'])) != witnesses[key]:
                raise ValueError('Receipt witness differs from fixed proposal')
            counts = dict(source_supports=source_counts()[source], cohort_colorings=17,
                one_row_choices=8, two_row_geometries=geometry_count(*bounds))
            if any(type(row.get(name)) is not int or row[name] != value for name, value in counts.items()):
                raise ValueError('Incorrect once-only class multiplicity')
            if not isinstance(row.get('failure_log2_upper'), str):
                raise ValueError('Expected an exact rational logarithmic ceiling')
            rows[key] = F(row['failure_log2_upper'])
        if rows.keys() != witnesses.keys():
            raise ValueError('Incomplete directed class partition')
        checked.append(rows)
    if checked[0] != checked[1] or receipts[0]['failure_log2_upper'] != receipts[1]['failure_log2_upper']:
        raise ValueError('Independent precision replays disagree')
    with ctx.workprec(512):
        total = sum((arb(2)**(arb(value.numerator)/value.denominator)
                     for value in checked[0].values()), arb(0))
        ceiling = rational_ceiling(total.log()/arb(2).log())
    raw = F(receipts[0]['failure_log2_upper'])
    if ceiling < raw:
        raise ArithmeticError('Rounded-cell union lies below raw directed aggregate')
    return dict(construction=metadata()['construction'], cohort_counts=[1, 16, 0],
        active_groups=17, complete_cohort=True, verified_cells=len(checked[0]),
        precision_bits=[256, 384], independent_union_precision_bits=512,
        all_per_cell_ceilings_agree=True, raw_aggregate_log2_upper=str(raw),
        rounded_cell_union_log2_upper=str(ceiling), selected_cohort_reaches_40_bits=ceiling < -40,
        whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--witnesses', type=Path, required=True)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--source-pair', type=int, nargs=2, action='append')
    parser.add_argument('--verify', type=Path, nargs=2)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite directed receipts')
    data = args.witnesses.read_bytes()
    proposal = json.loads(data)
    if args.verify:
        contents = [path.read_bytes() for path in args.verify]
        result = verify_receipts(proposal, list(map(json.loads, contents)))
        result['input_sha256'] = [sha256(content).hexdigest() for content in contents]
    else:
        result = run(proposal, args.precision, selected_sources=args.source_pair,
            progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['witness_sha256'] = sha256(data).hexdigest()
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'entries'}), flush=True)
