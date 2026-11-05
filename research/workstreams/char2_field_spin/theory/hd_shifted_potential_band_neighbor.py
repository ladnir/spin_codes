"""Screen neighboring weights with one fixed bank, optionally fitting it once centrally."""
from fractions import Fraction as F
from math import log2
from pathlib import Path

from hd_shifted_potential_band_screen import FixedAtlasScreen
from hd_shifted_outer_density import (DEFAULT_OUTER_MEASURE, TOEPLITZ_OUTER_MEASURE,
    checked_outer_measure)
from hd_shifted_swap_screen import CONSTRUCTION


def checked_window(first, last, max_weights=33):
    """Reject an oversized workload before fitting or constructing any atlas."""
    if type(max_weights) is not int or not 1 <= max_weights <= 8192:
        raise ValueError('max_weights must be an integer from 1 through 8192')
    if (type(first) is not int or type(last) is not int
            or not 1 <= first <= last <= 8192
            or last-first+1 > max_weights):
        raise ValueError('The nonzero weight interval must fit the explicit max_weights budget')
    return first, last, max_weights


def source_witness(artifact):
    """Validate the ensemble and return exact literals, never saved bounds."""
    witness = artifact.get('fixed_witness', {})
    if (artifact.get('construction') != CONSTRUCTION
            or witness.get('construction') != CONSTRUCTION
            or witness.get('pair_swaps') is not True
            or witness.get('orientation_aggregation') != 'independent-fair-swaps'):
        raise ValueError('An explicit swapped-route fixed atlas artifact is required')
    exact = witness.get('source_exact_witness')
    if not exact:
        raise ValueError('Exact saved witness literals are required')
    return exact, witness.get('grid', 128)


def artifact_outer_measure(artifact, override=None):
    saved = checked_outer_measure(artifact.get('fixed_witness', {}).get('outer_measure', DEFAULT_OUTER_MEASURE))
    return saved if override is None else checked_outer_measure(override)


def fixed_screen_from_artifact(artifact, *, outer_measure=None):
    """Use exact source literals; ignore the source band's FLOAT bound values."""
    exact, grid = source_witness(artifact)
    odds = tuple(tuple(log2(F(value)) for value in row.split()) for row in exact['markers'])
    screen = FixedAtlasScreen(F(exact['z']), F(exact['alpha']), odds,
                              log2(F(exact['potential'])), grid=grid,
                              outer_measure=artifact_outer_measure(artifact, outer_measure))
    screen.source_witness = dict(exact)
    return screen


def fit_one_central_bank(artifact, weight, *, outer_measure=None, marker_step=None):
    """Scale log(z) with density and fit one central bank; retain alpha and v."""
    exact, grid = source_witness(artifact)
    if not isinstance(weight, int) or not 2 <= weight <= 4096:
        raise ValueError('An explicit medium-weight central target is required')
    base_weight = exact['active_groups']
    maximum = weight//2
    if marker_step is not None and (type(marker_step) is not int or marker_step < 1):
        raise ValueError('marker_step must be a positive integer when supplied')
    base_centers = ({round(maximum*index/16) for index in range(17)} if marker_step is None
                    else set(range(0, maximum+1, marker_step)) | {maximum})
    centers = tuple(sorted(base_centers
                           | {maximum-step for step in (1, 2, 4, 8, 16) if maximum >= step}))
    if len(centers) > 100:
        raise ValueError('The explicit marker step generates more than 100 centers')
    screen = FixedAtlasScreen.fit_central_once(weight, float(F(exact['z']))**(weight/base_weight),
        F(exact['alpha']), centers, log2(F(exact['potential'])), grid=grid,
        outer_measure=artifact_outer_measure(artifact, outer_measure))
    # Fitting uses FLOAT proposals, but these two parameters were not retuned.
    # Preserve their source rationals rather than exactifying a round trip.
    screen.source_witness.update(alpha=exact['alpha'], potential=exact['potential'])
    screen.source_witness.update(density_scaled_from_weight=base_weight,
                                  central_z_rule='log(z_new)=log(z_source)*new_weight/source_weight',
                                  alpha_and_potential='unchanged-from-source')
    if marker_step is not None:
        screen.source_witness['central_marker_step'] = marker_step
    return screen


if __name__ == '__main__':
    import argparse, json
    from hashlib import sha256
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--first', type=int, required=True)
    parser.add_argument('--last', type=int, required=True)
    parser.add_argument('--threshold', type=float, default=-80.)
    parser.add_argument('--max-weights', type=int, default=33,
                        help='Explicit work budget (default 33; at most 8192); never inferred from an artifact')
    parser.add_argument('--central-weight', type=int,
                        help='Fit one density-scaled central bank; never fit individual neighboring weights')
    parser.add_argument('--marker-step', type=int,
                        help='Explicit full-pair spacing for the central bank; at most100 markers')
    parser.add_argument('--outer-measure', choices=(DEFAULT_OUTER_MEASURE, TOEPLITZ_OUTER_MEASURE),
                        help='Explicitly override the source outer bound; omitted preserves its saved mode')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output exists; frozen artifacts are never overwritten')
    if args.marker_step is not None and args.central_weight is None:
        parser.error('--marker-step requires --central-weight')
    try:
        checked_window(args.first, args.last, args.max_weights)
    except ValueError as error:
        parser.error(str(error))
    raw = args.source.read_bytes()
    artifact = json.loads(raw)
    screen = (fixed_screen_from_artifact(artifact, outer_measure=args.outer_measure)
              if args.central_weight is None else fit_one_central_bank(
                  artifact, args.central_weight, outer_measure=args.outer_measure,
                  marker_step=args.marker_step))
    result = screen.band(args.first, args.last, threshold=args.threshold,
                         max_weights=args.max_weights)
    result['source_artifact_sha256'] = sha256(raw).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, separators=(',', ':'), allow_nan=False)
        target.write('\n')
    print(json.dumps(dict(path=str(args.output.resolve()),
                          **{key: value for key, value in result.items()
                             if key not in ('results', 'fixed_witness')})), flush=True)
