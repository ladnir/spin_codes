"""Exact bounded K20 transfer to eight GF16 RS[8,4] rows per group.

Eight independent nonzero GF(2^32) scalars randomize the aligned columns.
Four-bit packet routing and the authenticated s22 inner remain unchanged.
This uses saved endpoints, not a local census or encoder benchmark.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import importlib.util
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('k20_rs8_transfer_helpers', HERE / 'field128_transfer.py')
helpers = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = helpers
spec.loader.exec_module(helpers)
HELPER_HASH = '45d1fd4dbc4ecec20ea8ffac2fff26038fe044160582d6b0aa746e5ccc6b8153'
Q = 1 << 32
GROUPS = 8192


def transfer(reference, sparse, new_counts, old_counts, exact_tail):
    for counts in (new_counts, old_counts):
        if len(counts) != 65 or sum(counts) != (1 << 128) - 1:
            raise ArithmeticError('invalid expected group measure')
    if [v for v, a in enumerate(new_counts) if a] != list(range(5, 65)):
        raise ArithmeticError('RS[8,4] packet support must range from five to 64')
    missing = [v for v in range(65) if new_counts[v] and not old_counts[v]]
    if missing != list(range(5, 9)):
        raise ArithmeticError('unexpected missing reference supports')
    # Exclude the event that some local group has any word in these shells.
    # Products below use unconditional truncated expected measures, not a
    # new setup distribution conditioned on avoiding this event.
    rare = GROUPS * sum((new_counts[v] for v in missing), F())
    ratios = {v: new_counts[v] / old_counts[v] for v in range(9, 65)}
    witness = max(ratios, key=ratios.__getitem__)
    ratio = ratios[witness]
    if any(new_counts[v] > ratio * old_counts[v] for v in range(9, 65)):
        raise ArithmeticError('regular shell domination failed')
    probabilities = list(map(helpers.dyadic, sparse['support_probability_uppers']))
    if len(probabilities) != 65 or any(not 0 <= p <= 1 for p in probabilities):
        raise ArithmeticError('invalid sparse support endpoints')
    q1 = GROUPS * sum((a * p for a, p in zip(new_counts, probabilities)), F())
    q2 = ratio**2 * helpers.dyadic(reference['occupancy_uppers']['2'])
    # For h active MDS symbols, at most (Q-1)^(h-4) codewords have one
    # specified support. Randomization divides by (Q-1)^h. Thus every
    # atom is <= (Q-1)^-4, improving the old uniform-input envelope.
    beta_new = F(1 << 256, (Q - 1)**4)
    beta_old = F(1 << 256, 65535**8)
    if not beta_new < beta_old:
        raise ArithmeticError('new uniform envelope must improve the old one')
    values = dict(rare=rare, q1=q1, q2_regular=q2, tail=exact_tail)
    endpoints = {name: helpers.ceil_dyadic(value) for name, value in values.items()}
    upper = helpers.ceil_dyadic(sum(map(helpers.dyadic, endpoints.values()), F()))
    return dict(exact_expected_shells=list(map(str, new_counts)),
        excluded_packet_supports=missing, regular_shell_ratio=str(ratio),
        regular_ratio_witness=witness, beta_new=str(beta_new), beta_old=str(beta_old),
        component_uppers=endpoints,
        component_margin_bits_display={name: helpers.margin(value) for name, value in values.items()},
        union_upper=upper, margin_bits_display=helpers.margin(helpers.dyadic(upper)),
        target_40_met=helpers.dyadic(upper) < F(1, 1 << 40))


def run():
    if helpers.digest(helpers.__file__) != HELPER_HASH:
        raise ValueError('frozen arithmetic helper identity changed')
    if (helpers.digest(helpers.REFERENCE) != helpers.REFERENCE_HASH or
            helpers.digest(helpers.Q1) != helpers.Q1_HASH):
        raise ValueError('retained certificate identity changed')
    audit = helpers.load_module('k20_rs8_old_audit', HERE.parent / 'k20_small_outer_screen/verify_whole.py')
    authentication = audit.verify(helpers.REFERENCE)
    reference = json.loads(helpers.REFERENCE.read_text(encoding='utf-8'))
    sparse = json.loads(helpers.Q1.read_text(encoding='utf-8'))
    if not any(Path(c['path']).resolve() == helpers.Q1.resolve() for c in reference['components']):
        raise ValueError('sparse receipt is not authenticated by complete receipt')
    outer_path = HERE.parent / 'k16_design/rs_outer.py'
    inherited = {str(Path(p).resolve()): h for p, h in reference['source_sha256'].items()}
    if inherited.get(str(outer_path.resolve())) != helpers.digest(outer_path):
        raise ValueError('RS enumerator is not authenticated')
    outer = helpers.load_module('k20_rs8_outer', outer_path)
    new_counts = outer.expected_group_support_counts(n=8, k=4, packet_bits=4, packets_per_symbol=8)
    old_counts = outer.expected_group_support_counts(n=16, k=8, packet_bits=4, packets_per_symbol=4)
    tail = audit.sum_dyadics(pair for q, pair in reference['occupancy_uppers'].items() if int(q) >= 3)
    if tail != audit.exact_union(reference['exact_tail_union']):
        raise ArithmeticError('authenticated tail union changed')
    result = transfer(reference, sparse, new_counts, old_counts, tail)
    result.update(schema='k20-rs8-gf32-s22-transfer-1', K=1048576, N=2097152,
        distance='1/10', minimum_distance=209716, bad_cutoff=209715, groups=GROUPS,
        group_dimension=128, group_output_bits=256, packet_bits=4, regions=64,
        outer='eight parallel GF16 RS[8,4] words; one aligned 32-bit column per symbol position',
        randomizer='independent uniform nonzero GF(2^32) multiplier in each of eight columns per group',
        counting_alphabet='GF(2^32) MDS[8,4]',
        inner_and_route='unchanged exact certified22 maps and independent uniform four-bit-packet routing',
        map_sha256=reference['map_sha256'], continuous_state=True,
        zero_initial_state=True, final_flush=False,
        reference=dict(path=str(helpers.REFERENCE), sha256=helpers.REFERENCE_HASH),
        sparse_reference=dict(path=str(helpers.Q1), sha256=helpers.Q1_HASH),
        reference_audit=authentication,
        proof_method='global rare-outer event; exact q1 fold; unconditional restricted-measure q2 comparison; unchanged no-worse uniform-envelope tail',
        all_occupancies_covered=True, local_operators_recomputed=False,
        arithmetic='exact rationals; outward integer dyadic conversion; exact acceptance comparison',
        source_sha256={str(Path(path).resolve()): helpers.digest(path) for path in
                       (__file__, helpers.__file__, audit.__file__, outer_path)},
        scope='Proof-only ideal independent-setup transfer; no implementation or timing claim. Independent mathematical review required before promotion.')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('fresh output required')
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({name: result[name] for name in
                     ('margin_bits_display', 'component_margin_bits_display',
                      'regular_ratio_witness', 'target_40_met')}, indent=2))
