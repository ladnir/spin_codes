"""Directed replay of a complete weight interval using fixed alias atlases.

The caller fixes exact z, alpha, coefficient markers, and one positive
state potential. Each directed atlas is built once at the requested
precision. Its selected residue entries are normalized separately for
every target histogram before taking entrywise minima across markers.

For each integer weight, run_sum_arb bounds all indexed scalar supports.
The band receipt sums these per-weight probability bounds. No weight is
interpolated, omitted, or inferred from its neighbors. The default resource
budget is33 weights. A caller can explicitly increase it up to the natural
8192-weight domain. A band is not a whole-code distance certificate.
"""
from fractions import Fraction as F
from pathlib import Path

from flint import arb, ctx

from hd_grouped_packet_screen import LIMIT
from hd_shifted_alias_certificate import rational
from hd_shifted_cover_backend import profile_count
from hd_shifted_coverage import RegionalAliases, rational_ceiling
from hd_shifted_good_pair_certificate import component_power
from hd_shifted_outer_density import (DEFAULT_OUTER_MEASURE, TOEPLITZ_OUTER_MEASURE,
    checked_outer_measure, log_density_arb, bound_description)
from hd_shifted_run_potential import run_sum_arb
from hd_shifted_swap_certificate import SWAP_CONSTRUCTION, swap_unmarked
from hd_shifted_swap_potential import full_counts, positive_vector
from hd_shifted_unoriented_certificate import unsigned_components


MAX_WEIGHTS = 33


def singleton_numerators(atlas, residues):
    """Cache singleton alias queries without constructing general range arrays.

    The stored binary64 upper endpoint and Arb scaling factor are exactly
    those used by RegionalAliases.numerator. Each distinct residue is scaled
    once per atlas; histogram normalization remains a separate later step.
    """
    upper, factors = atlas.upper, atlas.factors
    if upper.ndim != 3 or upper.shape[0] != 4 or len(factors) != 4:
        raise ValueError('A two-dimensional regional alias atlas is required')
    return {(singleton, full): tuple(arb(float(upper[index, singleton, full]))*factor
                                    for index, factor in enumerate(factors))
            for singleton, full in residues}


def checked_inputs(first, last, z, alpha, markers, grid, precision, bits, *, max_weights=MAX_WEIGHTS):
    if (not isinstance(first, int) or not isinstance(last, int)
            or type(max_weights) is not int or not 1 <= max_weights <= 8192
            or not 1 <= first <= last <= 8192 or last-first+1 > max_weights
            or isinstance(z, float) or isinstance(alpha, float)
            or not 0 < F(z) < 1 or not 0 < F(alpha) <= 1
            or not isinstance(grid, int) or grid < 64 or grid & (grid-1)
            or not isinstance(precision, int) or precision < 128
            or not isinstance(bits, int) or bits < 16):
        raise ValueError('Exact admissible parameters and an interval within the explicit weight budget are required')
    rows = []
    for row in markers:
        if len(row) != 2 or any(isinstance(value, float) for value in row):
            raise ValueError('Exact two-entry coefficient markers are required')
        row = tuple(map(F, row))
        if min(row) <= 0:
            raise ValueError('Coefficient markers must be positive')
        if row not in rows:
            rows.append(row)
    if not rows or len(rows) > 100:
        raise ValueError('A nonempty bank of at most100 distinct marker pairs is required')
    return F(z), F(alpha), tuple(rows)


def directed_band_regions(first, last, z, alpha, markers, *, grid=128,
                          precision=256, bits=60, max_weights=MAX_WEIGHTS, progress=None):
    """Return all normalized regional upper matrices, streaming each atlas once."""
    z, alpha, markers = checked_inputs(first, last, z, alpha, markers, grid, precision, bits,
                                       max_weights=max_weights)
    with ctx.workprec(precision):
        single, paired = swap_unmarked(z, alpha, precision, bits)
        histograms = [(weight, full, weight-2*full, 4096-weight+full)
                      for weight in range(first, last+1) for full in full_counts(4096, weight)]
        needed = {4096} | {value for _, full, singleton, empty in histograms
                           for value in (full, singleton, empty)}
        log_two = arb(2).log()
        factorial = {n: (arb(n)+1).lgamma()/log_two for n in needed}
        targets = [(weight, full, singleton,
                    factorial[4096]-factorial[empty]-factorial[singleton]-factorial[full],
                    (singleton % grid, full % grid))
                   for weight, full, singleton, empty in histograms]
        residues = set(target[-1] for target in targets)
        best = {weight: {} for weight in range(first, last+1)}
        for index, row in enumerate(markers):
            pieces = unsigned_components(single, paired, row)
            powered, scaling_vector, rho = component_power(*pieces, grid, dimensions=2)
            atlas = RegionalAliases.from_power(powered, scaling_vector, rho)
            numerators = singleton_numerators(atlas, residues)
            log_single, log_full = (rational(value).log()/log_two for value in row)
            for weight, full, singleton, log_mult, residue in targets:
                numerator = numerators[residue]
                scale = arb(2)**(-log_mult-singleton*log_single-full*log_full)
                values = tuple((entry*scale).upper() for entry in numerator)
                if len(values) != 4 or any(not value >= 0 or not value.is_finite() for value in values):
                    raise ArithmeticError('Invalid normalized regional upper matrix')
                previous = best[weight].get(full)
                best[weight][full] = (values if previous is None
                                      else tuple(min(a, b) for a, b in zip(previous, values)))
            if progress:
                progress(dict(stage='directed-band-atlas', precision=precision,
                              completed=index+1, total=len(markers),
                              weight_range=[first, last], normalized_histograms=len(targets)))
    return best


def band_receipt(first, last, z, alpha, markers, vector, *, grid=128,
                 precision=256, bits=60, max_weights=MAX_WEIGHTS,
                 outer_measure=DEFAULT_OUTER_MEASURE, progress=None):
    """Replay every weight and sum its positive upper probability bound."""
    z, alpha, markers = checked_inputs(first, last, z, alpha, markers, grid, precision, bits,
                                       max_weights=max_weights)
    outer_measure = checked_outer_measure(outer_measure)
    vector = positive_vector(vector)
    regions = directed_band_regions(first, last, z, alpha, markers, grid=grid,
                                     precision=precision, bits=bits, max_weights=max_weights, progress=progress)
    rows = []
    with ctx.workprec(precision):
        log_two = arb(2).log()
        log = lambda value: value.log()/log_two
        density = log_density_arb(outer_measure)
        constant = -rational(alpha)*(log(arb(255))+LIMIT*log(rational(z)))
        slope = rational(alpha)*density
        union = arb(0)
        for weight in range(first, last+1):
            # run_sum_arb includes v[0]/min(v), exactly once for this event.
            support_sum = run_sum_arb(4096, weight, regions[weight], vector)
            probability = (support_sum*arb(2)**(constant+slope*weight)).upper()
            if not probability > 0 or not probability.is_finite():
                raise ArithmeticError('Invalid directed per-weight probability bound')
            bound = log(probability)
            rows.append(dict(active_groups=weight, regional_histograms=len(regions[weight]),
                             log2_upper_rational=str(rational_ceiling(bound)),
                             log2_upper_enclosure=str(bound)))
            union += probability
            if progress and ((weight-first+1) % 8 == 0 or weight == last):
                progress(dict(stage='directed-band-weights', precision=precision,
                              completed=weight-first+1, total=last-first+1,
                              latest_log2_upper_rational=rows[-1]['log2_upper_rational']))
        union_log = log(union)
        ceiling = rational_ceiling(union_log)
    ranges = ((first, last), (0, 4096), (0, 4096))
    result = dict(ranges=ranges, weight_range=[first, last], results=rows,
                weights_replayed=len(rows), max_weights=max_weights,
                every_integer_weight_replayed=True, atlas_builds=len(markers),
                integer_profile_triples=profile_count(ranges),
                log2_upper_rational=str(ceiling), log2_upper_enclosure=str(union_log),
                directed=True, precision=precision, whole_code_certificate=False,
                covers_all_supports_at_each_weight=True, covers_all_even_odd_splits=True,
                construction=SWAP_CONSTRUCTION, pair_swaps=True,
                support_scope='scalar-code active input groups',
                source_dimension=1048576, output_length=2097152, outer_groups=8192,
                field_order=256, forbidden_weight=LIMIT,
                support_count_method='common-potential-run-convolution',
                coefficient_method='fixed-normalized-directed-alias-bank',
                witness=dict(z=str(z), alpha=str(alpha),
                             markers=[list(map(str, row)) for row in markers],
                             potential=list(map(str, vector)), grid=grid, precision=precision,
                             packet_mantissa_bits=bits,
                             regional_minimum='after-histogram-normalization'))
    if outer_measure != DEFAULT_OUTER_MEASURE:
        metadata = dict(outer_measure=outer_measure,
                        outer_measure_bound=bound_description(outer_measure))
        result.update(metadata)
        result['witness'].update(metadata)
    return result


def replay_saved_band(artifact, *, precision=256, max_weights=MAX_WEIGHTS,
                      outer_measure=None, progress=None):
    """Read only exact proposals from a FLOAT artifact, never its bound values."""
    witness = artifact.get('fixed_witness', {})
    if (artifact.get('construction') != SWAP_CONSTRUCTION
            or witness.get('construction') != SWAP_CONSTRUCTION
            or witness.get('pair_swaps') is not True
            or witness.get('orientation_aggregation') != 'independent-fair-swaps'):
        raise ValueError('An explicit independently swapped fixed-atlas proposal is required')
    source = witness.get('source_exact_witness')
    if not source:
        raise ValueError('Exact saved witness literals are required')
    saved_measure = checked_outer_measure(witness.get('outer_measure', DEFAULT_OUTER_MEASURE))
    measure = saved_measure if outer_measure is None else checked_outer_measure(outer_measure)
    first, last = artifact['requested_weight_range']
    markers = tuple(tuple(map(F, row.split())) for row in source['markers'])
    return band_receipt(first, last, F(source['z']), F(source['alpha']), markers,
                        (F(1), F(source['potential'])), grid=witness.get('grid', 128),
                        precision=precision, max_weights=max_weights, outer_measure=measure, progress=progress)


if __name__ == '__main__':
    import argparse, json
    from hashlib import sha256
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--proposal', type=Path, required=True)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--max-weights', type=int, default=MAX_WEIGHTS,
                        help='Explicit replay resource budget; default33, at most8192')
    parser.add_argument('--outer-measure', choices=(DEFAULT_OUTER_MEASURE, TOEPLITZ_OUTER_MEASURE),
                        help='Explicitly override the proposal outer bound; omitted uses its saved mode')
    parser.add_argument('--output', type=Path, help='Exclusive-create a new directed receipt artifact')
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error('Output exists; frozen artifacts are never overwritten')
    raw = args.proposal.read_bytes()
    result = replay_saved_band(json.loads(raw), precision=args.precision, max_weights=args.max_weights,
                               outer_measure=args.outer_measure,
                               progress=lambda row: print(json.dumps(row), flush=True))
    result['proposal_sha256'] = sha256(raw).hexdigest()
    if args.output:
        with args.output.open('x', encoding='utf-8') as target:
            json.dump(result, target, separators=(',', ':'), allow_nan=False)
            target.write('\n')
        print(json.dumps(dict(path=str(args.output.resolve()),
                              **{key: value for key, value in result.items()
                                 if key not in ('results', 'witness', 'log2_upper_enclosure')})), flush=True)
    else:
        print(json.dumps(result), flush=True)
