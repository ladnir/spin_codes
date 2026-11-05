"""Isolated width-aware FLOAT alias bank, default full-RS[128,48].

The frozen width/state caps supply local matrices. Packet coefficient
shapes, regional pair powers, and Perron restoration depend explicitly on
width. The frozen FFT component transform, histogram normalization, and
all-profile common-potential engine are reused without module mutation.
This is a selected-weight screen, not an encoder or distance certificate.
"""

import ast
from fractions import Fraction as F
from hashlib import sha256
from math import isfinite, log2
from pathlib import Path

import numpy as np

import hd_field128_wide_rs_state_screen as local
import hd_field128_lift_alias_screen as alias
import hd_field128_lift_alias_cover as cover


old = local.old
DEFAULT_WIDTH, DEFAULT_STATE = 128, 48


class AliasCoefficients:
    """Width-dependent histograms; compose good-good matrices before alpha."""

    def __init__(self, z, alpha, *, q=old.FIELD_ORDER, width=DEFAULT_WIDTH,
                 state=DEFAULT_STATE, joint=True, max_seconds=60., _budget=None):
        local.construction(width, state)
        old.checked_q(q, construction=True)
        if not isfinite(z) or not isfinite(alpha) or not 0 < z < 1 or not 0 < alpha <= 1:
            raise ValueError('Expected finite valid mark and fractional power')
        budget = _budget or old.Budget(max_seconds)
        budget.check()
        self.z, self.alpha, self.q = float(z), float(alpha), q
        self.width, self.state, self.joint = width, state, joint
        pairs = width//2
        logs = local.local_logs(self.z, q, width, state, joint)
        self.single = np.full((4, pairs+1, pairs+1), old.NEG)
        for s, f, count in local.geometry(width, width=width):
            self.single[:, s, f] = log2(count)+self.alpha*logs[s+2*f].ravel()
        good = local.geometry(width, width=width, good=True)
        s = np.array([row[0] for row in good], dtype=int)
        f = np.array([row[1] for row in good], dtype=int)
        # Width128 good-type counts can exceed uint64; take logs before conversion.
        weights = np.array([log2(row[2]) for row in good])
        raw = logs[s+2*f]
        composed = np.logaddexp2(raw[:, None, :, 0, None]+raw[None, :, None, 0, :],
                                 raw[:, None, :, 1, None]+raw[None, :, None, 1, :])
        degrees = ((s[:, None]+s[None, :])*(width+1)+f[:, None]+f[None, :]).ravel()
        values = weights[:, None, None]+weights[None, :, None]+self.alpha*composed.reshape((width+1, width+1, 4))
        result = np.full((4, 3*(width+1)), old.NEG)
        for entry in range(4):
            budget.check()
            np.logaddexp2.at(result[entry], degrees, values[:, :, entry].ravel())
        self.composed = result.reshape((4, 3, width+1))
        self.packet = self
        budget.check()

    def pair_logs(self, odds, *, max_seconds=60., _budget=None):
        odds = alias.checked_markers((odds,))[0]
        budget = _budget or old.Budget(max_seconds)
        budget.check()
        pairs = self.width//2
        tilt = np.arange(pairs+1)[:, None]*odds[0]+np.arange(pairs+1)[None, :]*odds[1]
        marked = self.single+tilt[None, :, :]
        good = np.logaddexp2.reduce(marked[:, :2, :].reshape((4, -1)), axis=1).reshape((2, 2))
        dirty = np.logaddexp2.reduce(marked[:, 2:, :].reshape((4, -1)), axis=1).reshape((2, 2))
        tilt2 = np.arange(3)[:, None]*odds[0]+np.arange(self.width+1)[None, :]*odds[1]
        result = np.logaddexp2.reduce((self.composed+tilt2[None, :, :]).reshape((4, -1)), axis=1).reshape((2, 2))
        for left, right in ((good, dirty), (dirty, good), (dirty, dirty)):
            result = np.logaddexp2(result, np.asarray(old.multiply(left, right)))
        if np.any(np.isnan(result)) or np.any(np.isposinf(result)):
            raise ArithmeticError('Invalid positive width-aware paired logarithm')
        budget.check()
        return result


def alias_table(z, alpha, odds, *, q=old.FIELD_ORDER, width=DEFAULT_WIDTH, state=DEFAULT_STATE,
                joint=True, grid=128, cushion=1e-12, cancellation_floor=1e-10,
                max_seconds=60., _coefficients=None, _budget=None):
    """Use32 two-packet powers at width128,64 powers at width64."""
    local.construction(width, state)
    old.checked_q(q, construction=True)
    odds = alias.checked_markers((odds,))[0]
    if (type(grid) is not int or grid not in (64, 128, 256)
            or not isfinite(cushion) or not 0 < cushion <= 1e-8
            or not isfinite(cancellation_floor) or not 100*cushion <= cancellation_floor <= 1e-4):
        raise ValueError('Use grid64/128/256 and positive cushion/floor with floor>=100*cushion')
    budget = _budget or old.Budget(max_seconds)
    coefficients = _coefficients or AliasCoefficients(z, alpha, q=q, width=width, state=state,
                                                      joint=joint, _budget=budget)
    if ((coefficients.z, coefficients.alpha, coefficients.q, coefficients.width,
         coefficients.state, coefficients.joint) != (float(z), float(alpha), q, width, state, joint)):
        raise ValueError('Cached coefficient parameters differ from this witness')
    route = local.routing(width)
    steps = route['composed_pairs_per_region']
    squarings = steps.bit_length()-1
    if steps != 1 << squarings:
        raise ValueError('Regional packet-pair count must be a power of two')
    pair = coefficients.pair_logs(odds, _budget=budget)
    # This rho scales a TWO-packet matrix, so its restoration exponent is steps.
    log_rho, log_ratio = alias.perron_logs(pair)
    log_vector = np.array((0., log_ratio))
    similarity = (log_vector[None, :]-log_vector[:, None]).ravel()
    good, dirty = coefficients.single.copy(), coefficients.single.copy()
    good[:, 2:, :], dirty[:, :2, :] = old.NEG, old.NEG
    transforms, underflows = [], 0
    for polynomial, scale in ((good, log_rho/2), (dirty, log_rho/2),
                              (coefficients.composed, log_rho)):
        budget.check()
        frequency, count = alias._frequency(polynomial, odds, similarity-scale, grid)
        transforms.append(frequency)
        underflows += count
    good, dirty, composed = transforms
    values = (composed+alias._matrix_product(good, dirty)+alias._matrix_product(dirty, good)
              +alias._matrix_product(dirty, dirty))
    del transforms, good, dirty, composed
    row_sums = values[:, 0, 0].real.reshape((2, 2)).sum(axis=1)
    row_error = float(np.max(abs(row_sums-1)))
    if not np.all(np.isfinite(values)) or row_error > 1e-9:
        raise ArithmeticError('Normalized positive row sums failed their FLOAT check')
    for unused in range(squarings):
        budget.check()
        a, b, c, d = values
        trace, bc = a+d, b*c
        values = np.stack((a*a+bc, b*trace, c*trace, d*d+bc))
    aliases = np.fft.ifftn(values, axes=(1, 2))
    if not np.all(np.isfinite(aliases)):
        raise ArithmeticError('Nonfinite inverse FFT')
    imaginary, minimum = float(np.max(abs(aliases.imag))), float(np.min(aliases.real))
    if imaginary > 1e-10 or minimum < -1e-10:
        raise ArithmeticError('FLOAT alias cancellation diagnostic failed')
    reliable = aliases.real >= cancellation_floor
    logs = np.log2(np.maximum(aliases.real, 0)+cushion)
    logs += steps*log_rho-similarity[:, None, None]
    numerator = np.array(old.power(pair, steps)).ravel()
    budget.check()
    return alias.AliasTable(logs, reliable, numerator, dict(
        grid=grid, inner_width=width, state_dimension=state,
        packets_per_region=route['packets_per_region'], paired_region_power=steps,
        matrix_squarings=squarings, perron_pair_restoration_power=steps,
        perron_packet_restoration_power=route['packets_per_region'],
        log2_pair_scale=log_rho, log2_packet_scale=log_rho/2, log2_similarity_ratio=log_ratio,
        normalized_row_sum_error=row_error, maximum_imaginary_alias=imaginary,
        minimum_real_alias=minimum, cancellation_floor=cancellation_floor,
        diagnostic_cushion=cushion, reliable_entries=int(reliable.sum()),
        below_floor_entries=int(reliable.size-reliable.sum()), normalized_input_underflow_terms=underflows,
        underflow_interpretation='FLOAT terms below representation; not exact zeros',
        directed=False, whole_code_certificate=False))


class AliasBankScreen(alias.AliasBankScreen):
    """Frozen4096-pair histogram and common-potential engine, new atlas."""

    def __init__(self, z, alpha, marker_log_odds, log_ratio, *, width=DEFAULT_WIDTH,
                 state=DEFAULT_STATE, joint=True, q=old.FIELD_ORDER, distance=F(12, 25),
                 outer_measure='invertible-toeplitz', grid=128, floor_policy='positive-fallback',
                 max_seconds=60., max_table_bytes=384 << 20, _coefficients=None, _budget=None):
        self.construction = local.construction(width, state)
        self.width, self.state, self.joint = width, state, joint
        self.bits = old.checked_q(q, construction=True)
        self.q, self.distance = q, F(distance)
        if not 0 < self.distance < F(1, 2):
            raise ValueError('Expected a distance strictly between zero and one half')
        self.z, self.alpha, self.log_ratio = float(z), float(alpha), float(log_ratio)
        if not isfinite(self.log_ratio) or abs(self.log_ratio) > 4096:
            raise ValueError('Expected a finite common log-potential in[-4096,4096]')
        self.odds = alias.checked_markers(marker_log_odds)
        if floor_policy not in ('positive-fallback', 'legacy-cushion'):
            raise ValueError('Unknown cancellation-floor policy')
        if (type(grid) is not int or grid not in (64, 128, 256)
                or not 1 <= max_table_bytes <= 512 << 20
                or len(self.odds)*4*grid*grid*9 > max_table_bytes):
            raise ValueError('Alias bank exceeds its explicit table memory budget')
        self.grid, self.floor_policy, self.outer_measure = grid, floor_policy, outer_measure
        self.density = old.outer_log_density(q, outer_measure)
        self.forbidden = int(self.distance*old.OUTPUT_LENGTH)
        budget = _budget or old.Budget(max_seconds)
        coefficients = _coefficients or AliasCoefficients(z, alpha, q=q, width=width, state=state,
                                                          joint=joint, _budget=budget)
        self.constant = -self.alpha*(self.bits+old.log1mexp2(-self.bits)+self.forbidden*log2(self.z))
        self.tables = tuple(alias_table(self.z, self.alpha, odds, q=q, width=width, state=state,
                                       joint=joint, grid=grid, _coefficients=coefficients, _budget=budget)
                            for odds in self.odds)
        self.numerators = np.array([table.positive_numerator.reshape((2, 2)) for table in self.tables])
        self.regional_diagnostics = None
        budget.check()

    def metadata(self):
        pairs = self.width//2
        return dict(**local.metadata(self.width, self.state), joint_rank_return_cap=self.joint,
                    single_histogram_terms=(pairs+1)*(pairs+2)//2, good_geometry_types=self.width+1,
                    raw_good_pair_products=(self.width+1)**2, maximum_packet_occupancy=self.width,
                    maximum_paired_occupancy=2*self.width, regional_even_odd_pair_power=old.REGIONS//2,
                    directed=False, whole_code_certificate=False, arithmetic='FLOAT')

    def evaluate(self, weight, **kwargs):
        return dict(super().evaluate(weight, **kwargs), **self.metadata(),
                    scope='selected active-group weight only; no encoder, timing, or global-distance claim')

    def witness(self):
        return dict(super().witness(), **self.metadata())

    def box(self, ranges, **kwargs):
        return dict(super().box(ranges, **kwargs), **self.metadata())


def proposal_from_artifact(artifact, provenance, *, width=DEFAULT_WIDTH, state=DEFAULT_STATE,
                           distance=F(12, 25), r=None, old_source=False, z=None, alpha=None):
    name = local.construction(width, state)
    exact = cover.witness_from_artifact(artifact, old=old_source)
    mark = F(str(float(exact['z']))) if z is None else cover.base.exact(z)
    exponent = F(str(float(exact['alpha']))) if alpha is None else cover.base.exact(alpha)
    distance = F(distance)
    if not 0 < mark < 1 or not 0 < exponent <= 1 or not 0 < distance < F(1, 2):
        raise ValueError('Inadmissible new-construction proposal')
    calibrated = {key: artifact[key] for key in ('center', 'selected_weight') if key in artifact}
    if r is not None and (type(r) is not int or not 1 <= r <= old.GROUPS
                          or any(type(value) is not int or value != r for value in calibrated.values())):
        raise ValueError('Requested selected weight must match explicit source marker calibration')
    return dict(schema='field128-wide-rs-state-alias-witness-only-v1', construction=name,
                state_dimension=state, inner_width=width, field_order=old.FIELD_ORDER, field_bits=128,
                target_distance=str(distance), selected_weight=r, pair_swaps=True,
                full_field_random_coefficients=True, source=provenance,
                source_role='exact proposal literals only', old_numerical_bounds_reused=False,
                source_weight_calibration=calibrated, marker_calibrations=artifact.get('marker_calibrations', []),
                status='new-construction untested witness proposal', directed=False, whole_code_certificate=False,
                state_potential_ratio=str(exact['ratio']),
                witness=dict(z=str(mark), alpha=str(exponent), grid=exact['grid'],
                             markers=[list(map(str, pair)) for pair in exact['markers']]))


def log_rational(value):
    value = F(value)
    if value <= 0:
        raise ValueError('Expected a positive exact marker or potential')
    return log2(value.numerator)-log2(value.denominator)


def source_pins():
    directory = Path(__file__).resolve().parent
    pending, seen = [Path(__file__).name], {}
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        raw = (directory/name).read_bytes()
        seen[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw)):
            modules = ([item.name for item in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0]+'.py'
                if (directory/target).is_file() and target not in seen:
                    pending.append(target)
    return dict(sorted(seen.items()))


def main():
    import argparse, json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--old-source', action='store_true')
    parser.add_argument('--width', type=int, choices=local.WIDTHS, default=DEFAULT_WIDTH)
    parser.add_argument('--state', type=int, default=DEFAULT_STATE)
    parser.add_argument('--distance', default='12/25')
    parser.add_argument('--r', type=int, required=True)
    parser.add_argument('--z')
    parser.add_argument('--alpha')
    parser.add_argument('--grid', type=int)
    parser.add_argument('--log-ratio', type=float)
    parser.add_argument('--optimize-potential', action='store_true')
    parser.add_argument('--max-evaluations', type=int, default=40)
    parser.add_argument('--max-seconds', type=float, default=120.)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output must be a new file')
    budget = old.Budget(args.max_seconds)
    pins = source_pins()
    artifact, provenance = cover.read_json(args.source)
    proposal = proposal_from_artifact(artifact, provenance, width=args.width, state=args.state,
                                      distance=F(args.distance), r=args.r, old_source=args.old_source,
                                      z=args.z, alpha=args.alpha)
    w = proposal['witness']
    markers = tuple(tuple(log_rational(x) for x in pair) for pair in w['markers'])
    ratio = log_rational(proposal['state_potential_ratio']) if args.log_ratio is None else args.log_ratio
    bank = AliasBankScreen(float(F(w['z'])), float(F(w['alpha'])), markers, ratio,
                           width=args.width, state=args.state, distance=F(args.distance),
                           grid=w['grid'] if args.grid is None else args.grid, _budget=budget)
    result = bank.evaluate(args.r, method='runs', optimize_potential=args.optimize_potential,
                            max_evaluations=args.max_evaluations, _budget=budget)
    witness = bank.witness()
    witness.update(log2_potential_ratio=result['log2_potential_ratio'],
                   log2_potential_ratio_exact=str(F(result['log2_potential_ratio'])))
    if source_pins() != pins:
        raise RuntimeError('A pinned alias helper changed during execution')
    result.update(proposal=proposal, witness=witness, source_sha256=pins)
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ('proposal', 'witness', 'source_sha256')}))


if __name__ == '__main__':
    main()
