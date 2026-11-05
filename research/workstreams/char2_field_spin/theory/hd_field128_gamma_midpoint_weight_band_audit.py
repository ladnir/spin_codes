"""Independent count, partition, provenance, and saved-ceiling band audit.

This auditor does not reconstruct coefficients or regional matrices. It
accepts numerical cell caps conditionally on the reviewed, byte-pinned
producer, checks their stored enclosures, and recomputes their union.
"""

from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb
from pathlib import Path

from flint import arb, ctx

import hd_field128_gamma_midpoint_weight_bands as producer


HERE = Path(__file__).resolve().parent
base, tail = producer.base, producer.tail
PRODUCER_CLOSURE = '5564792fea3c6738aad902bb571e6028acc3cc48ffdfea3db65154a548096e33'
CONSTRUCTION = 'toeplitz-fullfield-full-rs64-30-gamma16-feedback-alternating-swapped-pairs-k20-rate-half-v1'
COUNT_METHOD = 'sum of raw midpoint rectangle caps over explicit original weights'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def source_pins():
    pins = producer.source_pins()
    require(sha256(canonical(pins)).hexdigest() == PRODUCER_CLOSURE,
            'The reviewed producer closure changed')
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def check_source_pins(pins):
    for name, digest in pins.items():
        require(sha256((HERE/name).read_bytes()).hexdigest() == digest,
                'Pinned audit source changed: '+name)


class RectangleCounts:
    """Independent five-prefix formula; one weight retained at a time.

    The square prefix adds twice each lower-triangle row minus its diagonal.
    A rectangular prefix extends the smaller square by complete rows. This
    does not call either frozen midpoint count or rectangle implementation.
    """

    def __init__(self, pairs, weight, budget):
        require(type(pairs) is int and pairs > 0 and type(weight) is int
                and 0 <= weight <= 2*pairs, 'Invalid count parameters')
        self.pairs, self.r = pairs, min(weight, 2*pairs-weight)
        self.shift = max(0, weight-pairs)
        low, high = self.r//2, (self.r+1)//2
        self.last, self.denominator = low, (pairs-low)*low*(pairs-high)*high
        v0, v1, a0, a1, square = ([0] for unused in range(5))
        left, right = comb(pairs-low, high), comb(pairs-high, low)
        for i in range(low+1):
            if not i % 64:
                budget.check()
            v0.append(v0[-1]+right)
            v1.append(v1[-1]+i*right)
            scale = (self.r-2*i+1)*left
            aa = scale*low*high*(pairs-self.r+i)
            bb = scale*(low*high-pairs*i)
            a0.append(a0[-1]+aa)
            a1.append(a1[-1]+bb)
            square.append(square[-1]+2*(aa*v0[-1]+bb*v1[-1])-right*(aa+i*bb))
            if i < low:
                numerator, denominator = (high-i)*(low-i), (pairs-self.r+i+1)*(i+1)
                left, rem0 = divmod(left*numerator, denominator)
                right, rem1 = divmod(right*numerator, denominator)
                require(not rem0 and not rem1, 'Nonintegral midpoint array recurrence')
        self.v0, self.v1, self.a0, self.a1, self.square = map(tuple, (v0, v1, a0, a1, square))

    def _prefix(self, x, y):
        x, y = min(x, self.last), min(y, self.last)
        if min(x, y) < 0:
            return 0
        if x < y:
            x, y = y, x
        return (self.square[y+1]+(self.a0[x+1]-self.a0[y+1])*self.v0[y+1]
                +(self.a1[x+1]-self.a1[y+1])*self.v1[y+1])

    def rectangle(self, al, ah, bl, bh):
        require(all(type(x) is int for x in (al, ah, bl, bh)), 'Noninteger count rectangle')
        al, ah = max(0, al-self.shift), min(self.last, ah-self.shift)
        bl, bh = max(0, bl-self.shift), min(self.last, bh-self.shift)
        if al > ah or bl > bh:
            return 0
        if self.r <= 1:
            return 1 if self.r == 0 else 2*self.pairs
        numerator = self._prefix(ah, bh)-self._prefix(al-1, bh)-self._prefix(ah, bl-1)+self._prefix(al-1, bl-1)
        value, remainder = divmod(self.pairs*numerator, self.denominator)
        require(not remainder and value >= 0, 'Nonintegral or negative midpoint rectangle')
        return value+(2 if self.r == self.pairs and al == bl == 0 else 0)


def expected_roots(pairs, weights):
    roots = []
    for eps in (0, 1):
        matching = [r for r in weights if r % 2 == eps]
        if matching:
            limit = max(min(r//2, pairs-r//2-eps) for r in matching)
            roots.append(dict(parity=eps, u_deficit_ranges=[[matching[0]//2, matching[-1]//2],
                                                           [0, limit], [0, limit]]))
    return roots


def _tree(record, roots):
    """Second, iterative partition check; no producer traversal is reused."""
    nodes = record['nodes']
    require(isinstance(nodes, list) and nodes, 'Missing full tree')
    identities = [row.get('id') for row in nodes]
    require(all(type(i) is int and i >= 0 for i in identities) and len(set(identities)) == len(nodes),
            'Invalid or duplicate tree identities')
    by_id = dict(zip(identities, nodes))
    pending = [(i, root, None) for i, root in reversed(list(enumerate(roots)))]
    seen, leaves = set(), []
    while pending:
        identity, expected, parent = pending.pop()
        require(identity not in seen and identity in by_id, 'Missing, shared, or cyclic child')
        seen.add(identity)
        row = by_id[identity]
        require(row['cell'] == expected and row['parent'] == parent, 'Tree geometry or parent mismatch')
        if row['status'] == 'accepted':
            leaves.append(identity)
        else:
            require(row['status'] == 'split', 'Incomplete tree cannot certify a band')
            axis, cut = row['split']['axis'], row['split']['integer_cut']
            require(type(axis) is int and axis in (0, 1, 2) and type(cut) is int, 'Invalid split')
            ranges = expected['u_deficit_ranges']
            require(ranges[axis][0] <= cut < ranges[axis][1], 'Nonpartitioning integer split')
            children = row['children']
            require(isinstance(children, list) and len(children) == 2
                    and all(type(i) is int for i in children), 'Invalid child list')
            for side in (1, 0):
                changed = [list(pair) for pair in ranges]
                changed[axis] = ([cut+1, ranges[axis][1]] if side else [ranges[axis][0], cut])
                pending.append((children[side], dict(parity=expected['parity'], u_deficit_ranges=changed), identity))
    require(seen == set(by_id), 'Unreachable tree node')
    require(record['accepted_leaves'] == leaves and record['pending_frontier'] == [], 'Leaf list mismatch')
    require(record['tree_nodes'] == len(nodes) == record['evaluated_nodes'], 'Incomplete node evaluation metadata')
    return by_id, leaves


def _enclosure(value, ceiling, label):
    enclosure = arb(value)
    require(enclosure.is_finite() and base.rational(F(ceiling)) >= enclosure.upper(), label+' understates enclosure')


def audit_numeric(record, pairs, budget):
    """Generic tiny-test core. Only audit_summary applies production guards."""
    first, last = record['requested_weight_range']
    weights = list(producer.checked_weights(pairs, first, last))
    context = record['context']
    digest = sha256(canonical(context)).hexdigest()
    require(record['context_sha256'] == digest, 'Context hash mismatch')
    roots = expected_roots(pairs, weights)
    require(record['requested_weights'] == context['requested_weights'] == weights, 'Omitted original integer weight')
    require(context['requested_weight_range'] == [first, last], 'Context range mismatch')
    require(record['root_cells'] == context['root_cells'] == roots, 'Incomplete parity roots')
    require(context['parity_weights'] == {str(eps): [r for r in weights if r % 2 == eps]
            for eps in (0, 1) if any(r % 2 == eps for r in weights)}, 'Parity weight list mismatch')
    histograms = sum(r//2-max(0, r-pairs)+1 for r in weights)
    require(context['regional_histograms_prepared'] == histograms, 'Incomplete original histogram preparation metadata')
    for key in ('complete_requested_band', 'covers_all_supports_in_band', 'covers_all_matching_profiles',
                'every_integer_weight_prepared', 'no_weight_interpolation'):
        require(record.get(key) is True, 'Missing completion claim: '+key)
    require(record['status'] == 'complete' and record['stop_reason'] is None, 'Partial producer result')
    by_id, leaves = _tree(record, roots)
    counts = dict.fromkeys(by_id, 0)
    # Stream one count table; never retain O(B*r) big integers or a profile-pair table.
    for r in weights:
        budget.check()
        table, u, eps = RectangleCounts(pairs, r, budget), r//2, r % 2
        limit = min(u, pairs-u-eps)
        for identity, row in by_id.items():
            (ul, uh), (dl, dh), (el, eh) = row['cell']['u_deficit_ranges']
            if row['cell']['parity'] == eps and ul <= u <= uh and dl <= limit and el <= limit:
                counts[identity] += table.rectangle(u-dh, u-dl, u-eh, u-el)
    failure_bits = record['cell_failure_bits']
    require(type(failure_bits) is int and 0 <= failure_bits <= 4096, 'Invalid cell budget')
    with ctx.workprec(384):
        for identity, row in by_id.items():
            budget.check()
            receipt = row['receipt']
            require(receipt['cell'] == row['cell'] and receipt['context_sha256'] == digest, 'Cell receipt linkage mismatch')
            require(receipt['complete_cell'] is True and receipt['status'] == 'complete', 'Missing completed node receipt')
            require(type(receipt['count_upper']) is int and receipt['count_upper'] == counts[identity], 'Midpoint cell count mismatch')
            require(receipt['count_method'] == COUNT_METHOD and receipt['failure_bits'] == failure_bits, 'Changed count method/budget')
            if receipt['zero'] is True:
                require(receipt['log2_upper_rational'] is None and receipt['log2_upper_enclosure'] is None,
                        'Zero cell has logarithmic bound')
                # A positive-count zero transfer is trusted only as a producer result;
                # no matrix data is embedded from which to rederive it here.
                passing = True
            else:
                require(receipt['zero'] is False and counts[identity] > 0, 'Invalid nonzero cell')
                _enclosure(receipt['log2_upper_enclosure'], receipt['log2_upper_rational'], 'Cell ceiling')
                passing = F(receipt['log2_upper_rational']) <= -failure_bits
            require(receipt['meets_cell_budget'] is passing, 'Cell pass flag mismatch')
            require(passing == (row['status'] == 'accepted'), 'Wrong split/accepted status')
        total = sum((arb(2)**base.rational(F(by_id[i]['receipt']['log2_upper_rational']))
                     for i in leaves if not by_id[i]['receipt']['zero']), arb(0))
        zero = total.is_zero()
        require(record['zero'] is zero, 'Band zero flag mismatch')
        value = None if zero else tail.log2arb(total)
        union = None if zero else str(base.rational_ceiling(value))
        if zero:
            require(record['log2_upper_rational'] is None and record['log2_upper_enclosure'] is None, 'Zero band has a log')
        else:
            require(F(record['log2_upper_rational']) >= F(union), 'Saved band ceiling understates fresh384 union')
            _enclosure(record['log2_upper_enclosure'], record['log2_upper_rational'], 'Band ceiling')
        target = record['requested_failure_bits']
        require(type(target) is int and 0 <= target <= 4096, 'Invalid band target')
        require(record['meets_requested_band_target'] is (zero or F(record['log2_upper_rational']) <= -target),
                'Band target flag mismatch')
    return dict(weights=weights, roots=roots, cells_checked=len(by_id), accepted_cells=len(leaves),
                normalized_histograms=histograms, zero=zero, union=union, union_enclosure=None if zero else str(value),
                positive_count_zero_cells=sum(row['receipt']['zero'] and counts[i] > 0 for i, row in by_id.items()))


EXPECTED = dict(schema='gamma-midpoint-weight-band-context-v1', construction=CONSTRUCTION,
    state_dimension=30, inner_width=64, inner_emission='single-RS[64,30]', state_nodes='128..157',
    field_bits=128, field_order=1 << 128, fixed_map_field_order=256, source_dimension=1 << 20,
    output_length=1 << 21, outer_groups=8192, pairs=4096, pair_slots_per_region=4096,
    regions=256, repetitions=256, paired_region_power=128, packets_per_region=128,
    composed_pairs_per_region=64, target_distance='2/5', forbidden_weight=838860,
    outer_measure='invertible-toeplitz', projective_divisor=str((1 << 128)-1),
    pair_swaps=True, initial_state='zero', final_state='unrestricted', state_persists_across_regions=True,
    all_pair_before_alpha=True, counts_outside_alpha=True, full_field_random_coefficients=True,
    gamma_extension_degree=16, gamma_base_field_order=256, CA_zero=True,
    coefficient_method='gamma-all-pair-before-alpha-directed-interval-radix2-DFT',
    no_coefficient_truncation=True, normalized_every_original_weight=True,
    no_weight_interpolation=True, count_only_complement_above_pairs=True,
    count_pointwise_minimum_used=False, whole_code_certificate=False, directed=False,
    complete_requested_band=False, no_profile_pair_table=True,
    scaling='upward 2^(bank.slope*r/256) times each normalized regional entry; bank.constant once per event',
    count_method='exact rectangle sums of raw midpoint caps, summed once over explicit weights',
    max_domain='exact per-cell U range; both deficit intervals marginally nonempty; componentwise maxima of normalized scaled matrices')


def _metadata(record):
    context = record['context']
    for key, value in EXPECTED.items():
        require(context.get(key) == value and type(context.get(key)) is type(value), 'Wrong construction/context field: '+key)
    for key, value in dict(schema='gamma-midpoint-weight-band-cover-v1', directed=True,
                          numerical_backend='actual-all-pair-bank', whole_code_certificate=False,
                          batch_source_verified_at_entry=True, batch_source_verified_at_exit=True).items():
        require(record.get(key) == value and type(record.get(key)) is type(value), 'Wrong producer result field: '+key)
    require(context['source_sha256'] == producer.source_pins(), 'Producer closure mismatch')
    require(context['bank_source_sha256'] == producer.provider.source_pins(), 'Bank closure mismatch')
    witness = context['witness']
    require(0 < F(witness['z']) < 1 and 0 < F(witness['alpha']) <= 1, 'Invalid exact z/alpha')
    for key, value in dict(state_dimension=30, inner_width=64, gamma_extension_degree=16,
                          gamma_base_field_order=256, all_pair_before_alpha=True, joint_rank_return_cap=True).items():
        require(witness.get(key) == value, 'Witness construction mismatch: '+key)
    require(type(witness['precision']) is int and 128 <= witness['precision'] <= 4096,
            'Invalid bank precision')
    require(witness['grid'] in (64, 128, 256) and type(witness['packet_mantissa_bits']) is int
            and 32 <= witness['packet_mantissa_bits'] <= 512, 'Invalid numerical witness parameters')
    markers = witness['markers']
    require(isinstance(markers, list) and 1 <= len(markers) <= 128
            and all(isinstance(pair, list) and len(pair) == 2 and all(F(x) > 0 for x in pair) for pair in markers),
            'Invalid positive marker bank')


class Inputs:
    def __init__(self):
        self.pins = {}

    def load(self, path, expected=None):
        path = Path(path).resolve()
        raw = path.read_bytes()
        digest = sha256(raw).hexdigest()
        require(expected is None or expected == digest, 'Input byte pin mismatch: '+str(path))
        previous = self.pins.get(str(path))
        pin = dict(path=str(path), sha256=digest)
        require(previous is None or previous == pin, 'Input changed during audit: '+str(path))
        self.pins[str(path)] = pin
        return json.loads(raw), pin

    def check(self):
        for path, pin in self.pins.items():
            require(sha256(Path(path).read_bytes()).hexdigest() == pin['sha256'], 'Pinned input changed: '+path)


def _provenance(value, inputs, sources):
    """Rehash explicit path/SHA links and source maps, without inheriting bounds."""
    if isinstance(value, dict):
        if isinstance(value.get('path'), str) and isinstance(value.get('sha256'), str):
            path = Path(value['path'])
            path = path if path.is_absolute() else HERE/path
            if str(path.resolve()) not in inputs.pins:
                linked, unused = inputs.load(path, value['sha256'])
                _provenance(linked, inputs, sources)
            else:
                require(inputs.pins[str(path.resolve())]['sha256'] == value['sha256'], 'Conflicting provenance byte pins')
        if isinstance(value.get('source_sha256'), dict):
            for name, digest in value['source_sha256'].items():
                require(name not in sources or sources[name] == digest, 'Conflicting helper source pin: '+name)
                require(sha256((HERE/name).read_bytes()).hexdigest() == digest, 'Provenance helper pin changed: '+name)
                sources[name] = digest
        for child in value.values():
            _provenance(child, inputs, sources)
    elif isinstance(value, list):
        for child in value:
            _provenance(child, inputs, sources)


def audit_summary(path, *, max_seconds=120.):
    """Audit a complete saved band plus its context and all node checkpoints."""
    budget, pins, inputs = base.Budget(max_seconds), source_pins(), Inputs()
    check_source_pins(pins)
    record, summary_pin = inputs.load(path)
    _metadata(record)
    context = record['context']
    sources = dict(pins)
    _provenance(context, inputs, sources)
    provenance = context.get('witness_source')
    require(isinstance(provenance, dict) and {'path', 'sha256'} <= provenance.keys(), 'Missing exact witness provenance')
    proposal_path = Path(provenance['path'])
    proposal, unused = inputs.load(proposal_path if proposal_path.is_absolute() else HERE/proposal_path, provenance['sha256'])
    for key in ('z', 'alpha', 'markers'):
        require(proposal['witness'][key] == context['witness'][key], 'Proposal/bank exact literal mismatch: '+key)
    directory = Path(summary_pin['path']).parent
    saved_context, context_pin = inputs.load(directory/'context.json')
    require(saved_context == dict(context=context, context_sha256=record['context_sha256'], directed=False,
                                  complete_requested_band=False, whole_code_certificate=False), 'Context artifact mismatch')
    checkpoint_count = record['checkpoints']
    require(type(checkpoint_count) is int and checkpoint_count == len(record['nodes'])+1, 'Missing node checkpoints')
    replay = {}
    for sequence in range(checkpoint_count):
        budget.check()
        checkpoint, unused = inputs.load(directory/f'node_{sequence:06d}.json')
        require(checkpoint['sequence'] == sequence and checkpoint['context_sha256'] == record['context_sha256'], 'Checkpoint linkage mismatch')
        for key in ('complete_requested_band', 'covers_all_supports_in_band', 'directed', 'batch_source_verified_at_exit', 'whole_code_certificate'):
            require(checkpoint[key] is False, 'Checkpoint overclaims completion')
        require(checkpoint['log2_upper_rational'] is None, 'Partial checkpoint has band bound')
        for row in checkpoint['changed_nodes']:
            replay[row['id']] = row
    require([replay[i] for i in sorted(replay)] == record['nodes'], 'Checkpoint replay differs from final tree')
    numeric = audit_numeric(record, 4096, budget)
    check_source_pins(sources)
    inputs.check()
    return dict(schema='gamma-midpoint-weight-band-independent-audit-v1', audit_result='PASS',
        input_summary=summary_pin, context_file=context_pin, context_sha256=record['context_sha256'],
        construction=CONSTRUCTION, state_dimension=30, field_order=1 << 128, target_distance='2/5',
        requested_weight_range=record['requested_weight_range'], requested_weights=numeric['weights'],
        complete_requested_band=True, covers_all_matching_profiles=True, integer_partition_checked=True,
        every_original_weight_preparation_metadata_checked=True, normalized_histograms=numeric['normalized_histograms'],
        counts_recomputed=True, count_algorithm='independent five-prefix square/rectangle formula with exact array recurrence; one original weight at a time',
        cells_checked=numeric['cells_checked'], accepted_cells=numeric['accepted_cells'], checkpoints_checked=checkpoint_count,
        cell_log_enclosures_checked=True, positive_count_zero_cells=numeric['positive_count_zero_cells'],
        directed_union_precision=384, saved_band_log2_upper_rational=record['log2_upper_rational'],
        reaggregated_cell_log2_upper_rational=numeric['union'], reaggregated_cell_log2_upper_enclosure=numeric['union_enclosure'],
        log2_upper_rational=record['log2_upper_rational'], zero=numeric['zero'],
        source_sha256=pins, checked_source_sha256=dict(sorted(sources.items())), checked_file_pins=list(inputs.pins.values()),
        source_pin_conflicts=False, source_and_input_bytes_checked_at_exit=True,
        coefficient_replay_in_this_audit=False, normalized_regional_replay=False, cell_transfer_recomputed=False,
        trusted_numerical_producer='reviewed and byte-pinned complete band producer; stored cell caps are not authenticated by this audit',
        no_weight_interpolation=True, whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('summary', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-seconds', type=float, default=120.)
    args = parser.parse_args()
    result = audit_summary(args.summary, max_seconds=args.max_seconds)
    print(json.dumps(producer.profiles.write_new(args.output, result)))
