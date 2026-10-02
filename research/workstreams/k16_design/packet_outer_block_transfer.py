"""Authenticated tail transfer from the closed S22/K18 uniform-input bound.

See OUTER_BLOCK_TRANSFER.md for the zero-state coset argument. This is a
derived q>=3 bound, not a fresh replay or a whole-code certificate. It does
not transfer the old exact-shell q1/q2 components. No state census is run.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from math import comb
from pathlib import Path

from flint import arb, ctx

SCHEMA = 'rs-larger-outer-block-transfer-tail-1'
BASE_SHA256 = '20cd124704ff83d7ac0e098f5b80ce97e59aabf48f547bd118b1b3b957d5b60b'
MAP_SHA256 = 'e2467b9d5de8403ccd8263f0f55c73df128a8a2dd4c4cc405bc3a8e057fc3581'
L = 2048
BASE_K = 262144
BASE_BETA = Q(1 << 256, ((1 << 16)-1)**8)


def _digest(raw):
    return hashlib.sha256(raw).hexdigest()


def _canonical(record):
    return json.dumps(record, sort_keys=True, separators=(',', ':')).encode()


def _dyadic(pair):
    if (not isinstance(pair, (tuple, list)) or len(pair) != 2 or
            any(type(x) is not int for x in pair) or pair[0] <= 0 or abs(pair[1]) > 10_000_000):
        raise ValueError('positive finite integer dyadic endpoint required')
    return pair[0], pair[1]


def _fraction(pair):
    mantissa, exponent = _dyadic(pair)
    return Q(mantissa << max(exponent, 0), 1 << max(-exponent, 0))


def _sum_dyadics(pairs):
    values = [_dyadic(pair) for pair in pairs]
    exponent = min(e for _, e in values)
    return sum(m << (e-exponent) for m, e in values), exponent


def _up(value):
    if not value.is_finite() or not value > 0:
        raise ArithmeticError('finite positive enclosure required')
    return list(map(int, value.upper().man_exp()))


def _aq(value):
    value = Q(value)
    return arb(value.numerator) / value.denominator


def _margin(pair):
    m, e = _dyadic(pair)
    return str(-(arb(m).log()/arb(2).log()+e))


def beta_ratio(multiplier):
    """Exact beta_new / beta_old**m for the two specified larger outers."""
    if type(multiplier) is not int or multiplier not in (2, 4):
        raise ValueError('multiplier must be 2 or 4')
    return Q(65535, 65537)**(4*multiplier)


def authenticate(path):
    """Authenticate this specific audited source, its components and choices.

    The SHA anchor intentionally rejects other/reformatted receipts. A different
    source requires a separate audit, not a caller-controlled trust override.
    """
    path = Path(path).resolve()
    raw = path.read_bytes()
    if _digest(raw) != BASE_SHA256:
        raise ValueError('the audited S22/K18 whole-receipt hash is required')
    whole = json.loads(raw)
    expected = dict(schema='rs16-state-ladder-complete-replay-1', K=BASE_K,
        N=2*BASE_K, groups=L, regions=64, group_dimension=128,
        group_output_bits=256, distance='1/10', threshold=52428,
        physical_t=64, state_bits=22, macros_per_region=64,
        inner_distribution='uniform_gl', inner_group='GL(22,2)',
        physical_updates_independent=True, outer_symbol_group='GL(16,2)',
        zero_initial_state=True, final_flush=False, precision=256,
        state_continuity='retained_between_every_physical_step_and_region',
        terminal='sum of all mass-envelope coordinates; no flush',
        birth_density='capped', return_denominator=(1 << 22)-1,
        all_occupancies_covered=True, whole_code_certificate=True,
        fresh_computation=True, fresh_replay=True, target_met=True)
    if any(whole.get(k) != v for k, v in expected.items()):
        raise ValueError('source geometry, state, setup or completion metadata mismatch')
    if whole['map_record'].get('map_sha256') != MAP_SHA256:
        raise ValueError('source map identity mismatch')
    sources = whole.get('source_sha256', {})
    if len(sources) != 103:
        raise ValueError('complete audited source snapshot required')
    for filename, digest in sources.items():
        if _digest(Path(filename).read_bytes()) != digest:
            raise ValueError(f'proof source changed: {filename}')
    components = []
    for item in whole['components']:
        component = json.loads(Path(item['path']).read_text())
        if (item.get('hash_scope') != 'canonical record JSON' or
                _digest(_canonical(component)) != item['sha256'] or
                component['schema'] != item['schema'] or
                component['occupancy_covered'] != item['occupancy_covered']):
            raise ValueError('component identity or occupancy mismatch')
        if (component.get('map_record') != whole['map_record'] or
                component.get('source_sha256') != sources or
                any(component.get(k) != v for k, v in expected.items()
                    if k not in ('schema', 'all_occupancies_covered', 'whole_code_certificate',
                                 'fresh_replay', 'target_met', 'birth_density', 'return_denominator'))):
            raise ValueError('component scope mismatch')
        components.append(component)
    endpoints, choices = {}, {}
    for q in range(3, L+1):
        index = whole['selected_components'].get(str(q))
        if type(index) is not int or not 0 <= index < len(components):
            raise ValueError('selected component index required')
        record = components[index]
        if (record.get('schema') != 'rs16-state-ladder-tail-1' or
                record.get('method') not in ('exact', 'fugacity') or
                record.get('beta') != str(BASE_BETA) or
                record.get('every_requested_occupancy_checked') is not True or
                record.get('every_shell_checked') is not True or
                record.get('count_kind') != 'exact_expected_shells' or
                record.get('count_sha256') != whole['count_sha256'] or
                q not in record['occupancy_covered']):
            raise ValueError('selected tail must be an audited uniform-majorant bound')
        endpoint = record['occupancy_uppers'][str(q)]
        if endpoint != whole['occupancy_uppers'][str(q)]:
            raise ValueError('selected endpoint differs from whole receipt')
        _dyadic(endpoint)
        witness = record['occupancy_choices'][str(q)]
        tilt = Q(witness['tilt'])
        if not tilt > 0 or str(tilt) not in record['tilts']:
            raise ValueError('selected rational tilt missing from computed trials')
        if record['method'] == 'fugacity':
            marker = Q(witness['marker_probability'])
            if (not 0 < marker <= 1 or str(marker) not in record['marker_probabilities'] or
                    (marker == 1 and q != L)):
                raise ValueError('invalid conditioned-iid witness')
        endpoints[q] = tuple(endpoint)
        choices[q] = dict(witness, method=record['method'], source_component=index)
    return whole, endpoints, choices


def transfer_endpoint(endpoint, *, q, tilt, multiplier):
    """Outward bound exp(lambda*d) ratio**q U**m / choose(L,q)**(m-1).

    The caller must supply an authenticated uniform-majorant endpoint. This
    arithmetic primitive does not itself turn an arbitrary number into proof.
    """
    ratio = beta_ratio(multiplier)
    if type(q) is not int or not 3 <= q <= L or Q(tilt) <= 0:
        raise ValueError('tail occupancy 3..2048 and positive rational tilt required')
    mantissa, exponent = _dyadic(endpoint)
    old = arb(mantissa) * arb(2)**exponent
    cutoff_difference = (2*BASE_K*multiplier)//10-multiplier*((2*BASE_K)//10)
    value = ((_aq(tilt)*cutoff_difference).exp() * _aq(ratio)**q * old**multiplier /
             arb(comb(L, q))**(multiplier-1))
    return _up(value)


def run(path, *, multiplier, precision=256, output=None):
    """Return a derived tail receipt. Missing q1/q2 are never silently filled."""
    beta_ratio(multiplier)
    if type(precision) is not int or precision < 192:
        raise ValueError('at least 192-bit outward arithmetic required')
    if output is not None and Path(output).exists():
        raise ValueError('output must be fresh')
    whole, endpoints, choices = authenticate(path)
    previous_precision = ctx.prec
    try:
        ctx.prec = precision
        transferred = {q: transfer_endpoint(endpoint, q=q, tilt=choices[q]['tilt'],
                                            multiplier=multiplier)
                       for q, endpoint in endpoints.items()}
        total_m, total_e = _sum_dyadics(transferred.values())
        upper = _up(arb(total_m)*arb(2)**total_e)
        result = dict(schema=SCHEMA, K=BASE_K*multiplier, N=2*BASE_K*multiplier,
            groups=L, regions=64*multiplier, group_dimension=128*multiplier,
            group_output_bits=256*multiplier, physical_t=64, state_bits=22,
            map_sha256=MAP_SHA256, inner_distribution='uniform_gl',
            outer=dict(parallel_rows=4, base_field_bits=8,
                       rs_n=8*multiplier, rs_k=4*multiplier, symbol_map='GL(32,2)'),
            independent_outer_symbol_maps=True, independent_regional_shuffles=True,
            zero_initial_state=True, final_flush=False,
            state_continuity='retained_between_every_physical_step_and_region',
            distance='1/10', threshold=(2*BASE_K*multiplier)//10,
            multiplier=multiplier, beta_ratio=str(beta_ratio(multiplier)),
            source=dict(path=str(Path(path).resolve()), sha256=BASE_SHA256,
                        source_pin_count=len(whole['source_sha256'])),
            derivation='uniform-input coset dominance and conditional block induction',
            proof_note_sha256=_digest(Path(__file__).with_name('OUTER_BLOCK_TRANSFER.md').read_bytes()),
            helper_sha256=_digest(Path(__file__).read_bytes()), precision=precision,
            occupancy_covered=list(range(3, L+1)), occupancy_missing=[1, 2],
            occupancy_uppers={str(q): v for q, v in transferred.items()},
            source_choices={str(q): value for q, value in choices.items()},
            tail_upper=upper, tail_margin_bits=_margin(upper),
            whole_code_certificate=False, fresh_computation=False,
            scope='Derived tail bound for the specified ideal larger-outer ensemble. '
                  'Requires separate q1 and q2 bounds; no implementation or seeded guarantee.')
        # Source identities must remain stable during even this small replay.
        authenticate(path)
    finally:
        ctx.prec = previous_precision
    if output is not None:
        with Path(output).open('x') as handle:
            json.dump(result, handle, indent=2)
            handle.write('\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source')
    parser.add_argument('--multiplier', type=int, choices=(2, 4), required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    result = run(args.source, multiplier=args.multiplier, output=args.output)
    print(json.dumps({k: result[k] for k in ('K', 'tail_margin_bits', 'occupancy_missing')}))
