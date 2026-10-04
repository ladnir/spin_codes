"""Outward dense coverage for the fixed K20 small-outer ensemble.

Fix q active groups and let T_j be the routing-conditioned two-step
fractional operator from replay_fractional.fractional_blocks. For a marker
probability 0<p<1, define M_p=sum_j C(32,j)p^j(1-p)^(32-j)T_j.
The regional polynomial has nonnegative matrix coefficients. Consequently,
its conditional q-slot operator R_q is entrywise at most M_p^256/c_q,
where c_q=C(8192,q)p^q(1-p)^(8192-q). Monotonicity of nonnegative matrix
products gives R_q^64 <= M_p^16384/c_q^64. The resulting union bound is

  C(8192,q) beta^(alpha*q) exp(alpha*tilt*209715)
      * e_zero M_p^16384 1 / c_q^64.

The group-subset union remains outside alpha. Thinning precedes the
fractional power, and the full matrix power retains state across regions.
At p=1 only q=8192 is covered. Every computed scalar is rounded outward.
The receipt covers its explicit interval only, never the whole code.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import json
from math import comb
from pathlib import Path
import sys
from time import monotonic

from flint import arb, arb_mat, ctx

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'k18_margin60'))
import replay_fractional as replay

SCHEMA = 'k20-small-outer-fractional-dense-complete-1'
K = 1 << 20
GROUPS = K // 128
REGIONS = 64
CUTOFF = (2 * K) // 10
BETA = F(1 << 256, 65535**8)


def marker_mixture(blocks, probability, *, matrix=arb_mat,
                   rational=replay.aq, rounding=replay.placement.rounded):
    """Mix complete fixed operators; zero and one are exact endpoints."""
    probability = F(probability)
    if not blocks or not 0 <= probability <= 1:
        raise ValueError('nonempty family and probability in [0,1] required')
    size, degree = blocks[0].nrows(), len(blocks) - 1
    if not size or any(block.nrows() != size or block.ncols() != size for block in blocks):
        raise ValueError('matching nonempty square matrices required')
    total = matrix(size, size)
    for j, block in enumerate(blocks):
        weight = F(comb(degree, j)) * probability**j * (1-probability)**(degree-j)
        if weight:
            total += block * rational(weight)
    return rounding(total)


def coefficient_log_bound(q, probability, log_moment, log_binomial, tilt, alpha,
                          *, groups=GROUPS, regions=REGIONS, cutoff=CUTOFF, beta=BETA):
    """Return an Arb enclosure, retaining the entire binomial logarithm."""
    probability, tilt, alpha, beta = map(F, (probability, tilt, alpha, beta))
    if (type(q) is not int or type(groups) is not int or type(regions) is not int
            or not 1 <= q <= groups or regions < 1 or cutoff < 0
            or tilt <= 0 or not 0 < alpha <= 1 or beta <= 0
            or not 0 < probability <= 1 or (probability == 1 and q != groups)):
        raise ValueError('valid geometry, occupancy and positive conditioning mass required')
    base = log_moment + replay.aq(alpha * tilt) * cutoff
    log_beta = replay.aq(alpha) * replay.aq(beta).log()
    if probability == 1:
        return base + groups * log_beta
    log_p = replay.aq(probability).log()
    log_complement = replay.aq(1-probability).log()
    return (base - regions*groups*log_complement + (1-regions)*log_binomial
            + q*(log_beta - regions*(log_p-log_complement)))


def _logs(first, last):
    """Exact consecutive binomials avoid repeated large integer products."""
    value = comb(GROUPS, first)
    result = {}
    for q in range(first, last+1):
        result[q] = arb(value).log()
        value = value * (GROUPS-q) // (q+1)
    return result


def _merge_pins(pins, incoming):
    for path, digest in incoming.items():
        if path in pins and pins[path] != digest:
            raise ArithmeticError(f'conflicting source pin: {path}')
        pins[path] = digest


def export_locals(output_prefix, tilts=('0.15', '0.25', '0.5', '1', '2.1972246')):
    """Export larger-tilt physical families from one fresh authenticated census."""
    import tune
    tilts, output_prefix = tuple(map(F, tilts)), Path(output_prefix)
    if not tilts or len(set(tilts)) != len(tilts) or min(tilts) <= 0:
        raise ValueError('distinct positive physical tilts required')
    paths = [output_prefix.with_name(output_prefix.name+'-t'+str(tilt).replace('/', '_')+'.json')
             for tilt in tilts]
    if any(path.exists() for path in paths):
        raise ValueError('fresh physical export paths required')
    start = monotonic()
    ctx.prec = 192
    print('Preparing one fresh certified s22 census for dense physical families', flush=True)
    data, record = tune.certified.prepare(22, birth_density='capped')
    tune.certified.authenticate(data, record)
    if record['map_sha256'] != replay.MAP:
        raise ArithmeticError('exact retained s22 map required')
    model = tune.ActivityModel(data, record, activity=F(1, 2), precision=192)
    pins = dict(model.sources)
    for path in (Path(__file__), Path(tune.__file__), Path(tune.dense.__file__)):
        _merge_pins(pins, {str(path.resolve()): replay.digest(path)})
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    print(f'Fresh census prepared in {monotonic()-start:.3f}s', flush=True)
    for path, tilt in zip(paths, tilts):
        model.physical(tilt)
        replay.authenticate(pins)
        tune.export_physical(path, model, tilt, pins)
        replay.authenticate(pins)
        print(f'Tilt {tilt}: elapsed {monotonic()-start:.3f}s; SHA256 {replay.digest(path)}', flush=True)
    return paths


def run(locals, output, *, minimum=513, maximum=GROUPS, alphas=('0.7', '0.8', '1'),
        markers=None, precision=256):
    """Check each integer in the interval with authenticated physical inputs."""
    locals, output = tuple(map(Path, locals)), Path(output)
    alphas = tuple(map(F, alphas))
    markers = tuple(F(j, 128) for j in range(1, 129)) if markers is None else tuple(map(F, markers))
    if (not locals or output.exists() or type(minimum) is not int or type(maximum) is not int
            or not 3 <= minimum <= maximum <= GROUPS or precision < 192
            or not alphas or any(not 0 < a <= 1 for a in alphas)
            or not markers or any(not 0 < p <= 1 for p in markers)
            or (minimum < GROUPS and not any(p < 1 for p in markers))):
        raise ValueError('fresh output, authenticated inputs, interval and valid witnesses required')
    start = monotonic()
    ctx.prec = precision
    pins = {}
    for path in (Path(__file__), Path(replay.__file__), Path(replay.placement.__file__),
                 HERE.parent/'k18_margin60'/'SLACK_REVIEW.md'):
        pins[str(path.resolve())] = replay.digest(path)
    logs = _logs(minimum, maximum)
    best = dict.fromkeys(logs)
    choices = dict.fromkeys(logs)
    trials = []
    map_record = None
    for local in locals:
        physical, tilt, source = replay.load_local(local)
        if tilt <= 0:
            raise ValueError('positive physical tilt required')
        if map_record is not None and source['map_record'] != map_record:
            raise ValueError('identical physical map records required')
        map_record = source['map_record']
        _merge_pins(pins, source['source_sha256'])
        _merge_pins(pins, {str(local.resolve()): replay.digest(local)})
        for alpha in alphas:
            trial_start = monotonic()
            blocks = replay.fractional_blocks(physical, alpha)
            for p in markers:
                if p == 1 and maximum != GROUPS:
                    continue
                mixed = marker_mixture(blocks, p)
                power = mixed ** (REGIONS * GROUPS // 32)
                moment = sum((power[0, j] for j in range(power.ncols())), arb(0))
                if not moment.is_finite() or not moment > 0:
                    raise ArithmeticError('finite positive all-terminal-state moment required')
                # Cache the affine dependence on q for this common moment.
                log_moment = moment.log()
                base = log_moment + replay.aq(alpha * tilt) * CUTOFF
                log_beta = replay.aq(alpha) * replay.aq(BETA).log()
                if p != 1:
                    log_complement = replay.aq(1-p).log()
                    base -= REGIONS*GROUPS*log_complement
                    slope = log_beta - REGIONS*(replay.aq(p).log()-log_complement)
                for q in ((GROUPS,) if p == 1 else logs):
                    value = (base + GROUPS*log_beta if p == 1 else
                             base + (1-REGIONS)*logs[q] + q*slope).upper()
                    if not value.is_finite():
                        raise ArithmeticError('finite logarithmic endpoint required')
                    if best[q] is None or value < best[q]:
                        best[q] = value
                        choices[q] = dict(local_receipt=str(local.resolve()), tilt=str(tilt),
                                          alpha=str(alpha), marker_probability=str(p))
            if any(value is None for value in best.values()):
                raise ArithmeticError('every requested integer needs a valid witness')
            replay.authenticate(pins)
            trial = dict(local_receipt=str(local.resolve()), tilt=str(tilt), alpha=str(alpha),
                         worst_margin_bits_display=float(-max(best.values())/arb(2).log()),
                         elapsed_seconds=monotonic()-trial_start)
            trials.append(trial)
            print(json.dumps(trial), flush=True)
    endpoints = {str(q): replay.endpoint(value.exp()) for q, value in best.items()}
    union = sum(map(replay.fraction, endpoints.values()), F(0))
    replay.authenticate(pins)
    if ctx.prec != precision:
        raise ArithmeticError('precision changed during replay')
    result = dict(schema=SCHEMA, whole_code_certificate=False, all_occupancies_covered=False,
        K=K, N=2*K, threshold=CUTOFF, map_record=map_record, state_bits=22,
        groups=GROUPS, regions=REGIONS, slots_per_region=GROUPS, group_dimension=128,
        group_output_bits=256, packet_bits=4, physical_t=64, physical_steps_per_macro=2,
        macros_per_region=GROUPS//32, zero_initial_state=True, continuous_state=True,
        final_flush=False, distance='1/10', beta=str(BETA), precision=precision,
        occupancy_covered=list(logs), occupancy_interval=[minimum, maximum],
        every_requested_occupancy_checked=True, occupancy_uppers=endpoints,
        occupancy_choices={str(q): value for q, value in choices.items()},
        exact_union=replay.encode_exact(union), margin_bits_display=replay.bits(union),
        local_census_recomputed_by_this_replay=False,
        local_receipts=[str(path.resolve()) for path in locals],
        alphas=list(map(str, alphas)), marker_probabilities=list(map(str, markers)),
        source_sha256=pins, trials=trials, elapsed_seconds=monotonic()-start,
        scope='Listed-occupancy setup-failure bound for independent ideal outer GL16, '
              'uniform routing, and inner GL22 updates; no seeded guarantee.')
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({key: result[key] for key in
          ('occupancy_interval', 'margin_bits_display', 'elapsed_seconds')}), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--local', type=Path, nargs='+')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--export-prefix', type=Path)
    parser.add_argument('--export-tilts', nargs='+', default=['.15', '.25', '.5', '1', '2.1972246'])
    parser.add_argument('--min-q', type=int, default=513)
    parser.add_argument('--max-q', type=int, default=GROUPS)
    parser.add_argument('--alphas', nargs='+', default=['.7', '.8', '1'])
    parser.add_argument('--markers', nargs='+')
    parser.add_argument('--precision', type=int, default=256)
    args = parser.parse_args()
    if args.export_prefix:
        if args.local or args.output:
            parser.error('physical export mode does not accept local or output')
        export_locals(args.export_prefix, args.export_tilts)
    else:
        if not args.local or not args.output:
            parser.error('dense replay mode requires local and output')
        run(args.local, args.output, minimum=args.min_q, maximum=args.max_q,
            alphas=args.alphas, markers=args.markers, precision=args.precision)
