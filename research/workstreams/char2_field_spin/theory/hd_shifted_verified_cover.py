"""Compose externally pinned, previously observed directed replay artifacts.

This is a trusted-replay-cache checker, not an arithmetic proof verifier.
The expected manifest SHA-256 is a trust anchor supplied outside the manifest.
It authorizes the manifest's byte pins and its record of observed fresh
256/384-bit runs. Hashes and cross-precision agreement do not prove arithmetic.

The checker independently checks exact witness identity, receipt scope,
disjoint geometry, exhaustive coverage, and a newly directed probability sum.
It never accepts a saved whole_code_certificate or a saved coverage summary.
An incomplete manifest always produces whole_code_certificate=False.

Manifest version1 has construction, trust_basis=TRUST_BASIS, and entries.
Each entry has id, kind, weight_range, observed_fresh_dual_replay=True,
and receipts={"256": artifact_ref, "384": artifact_ref}. An artifact_ref
has path and sha256, optionally format=json/jsonl-last/sparse-jsonl and
select (a list of dictionary keys or list indices). Paths are manifest-relative.
Potential/plan entries also have proposal=artifact_ref. Sparse batches
have z, alpha, maximum; sparse families use the frozen family_witness.
Potential entries explicitly name outer_measure; other paths retain the
inherited density. Optional implementation_pins pin reviewed source bytes.

Kinds: potential-band, sparse-family, sparse-batch, swap-plan,
unsigned-plan, terminal, checkpoint. Plan caches require accepted leaves and witnesses,
as emitted by the old cover ledger with --details or a programmatic export.
The assigned plan band can be narrower than its original root; whole leaf
mass is retained. Summary-only plan caches are deliberately insufficient.

A checkpoint entry additionally has checkpoint_kind=sparse/tail/dense.
Tail/dense entries have plan_bands (label,kind,weight_range,proposal) and
replay_source_pin for hd_shifted_cover_ledger.py. Their literal replay inputs
may be reconstructed from pinned proposals using that pinned implementation;
this is recorded, not misrepresented as receipt-stored witness metadata.

--fresh ignores cached bounds and invokes the existing replay_swap_cover
with threshold=-40, separately at both precisions. It is potentially costly.
"""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from flint import arb, ctx

from hd_shifted_alias_certificate import rational
from hd_shifted_cover_backend import CONSTRUCTION, profile_count, validate_box
from hd_shifted_cover_ledger import Box, Ledger
from hd_shifted_coverage import rational_ceiling
from hd_shifted_outer_density import DEFAULT_OUTER_MEASURE, checked_outer_measure, bound_description
from hd_shifted_swap_certificate import SWAP_CONSTRUCTION
from hd_shifted_swap_ledger import SwapLedger, UNSIGNED_METHODS
import hd_shifted_cover_ledger as replay_implementation


TRUST_BASIS = 'externally-pinned-observed-fresh-256-and-384-bit-replays'
# expected_plan_witness duplicates the conversion in this reviewed version.
# A future replay implementation needs an explicit compatibility review.
REPLAY_RECIPE_SOURCE_SHA256 = '17c86b79dc0559fc16fd035f485303c09332a547df798a4fa1e77c41dff902fd'
KINDS = frozenset(('potential-band', 'sparse-family', 'sparse-batch',
                   'swap-plan', 'unsigned-plan', 'terminal', 'checkpoint'))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def exact(value):
    require(type(value) in (str, int), 'An exact rational literal is required')
    return F(value)


def pinned_bytes(path, expected):
    require(isinstance(expected, str) and len(expected) == 64
            and all(c in '0123456789abcdef' for c in expected), 'Expected an external lowercase SHA-256 pin')
    raw = Path(path).read_bytes()
    require(sha256(raw).hexdigest() == expected, f'Byte hash mismatch: {path}')
    return raw


def strict_json(raw):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, f'Duplicate JSON key: {key}')
            result[key] = value
        return result
    def invalid(value):
        raise ValueError(f'Nonfinite JSON constant: {value}')
    return json.loads(raw, object_pairs_hook=pairs, parse_constant=invalid)


def read_artifact(base, reference):
    raw = pinned_bytes(base / reference['path'], reference['sha256'])
    form = reference.get('format', 'json')
    require(form in ('json', 'jsonl-last', 'sparse-jsonl'), 'Unknown artifact format')
    if form == 'json':
        result = strict_json(raw)
    else:
        rows = [strict_json(line) for line in raw.splitlines() if line.strip()]
        require(bool(rows), 'Empty replay log')
        result = rows[-1]
        if form == 'sparse-jsonl':
            require('weight_receipts' not in result, 'Ambiguous sparse log summary')
            result = dict(result, weight_receipts=rows[:-1])
    for key in reference.get('select', []):
        require(type(key) in (str, int), 'Invalid artifact selector')
        if isinstance(key, int):
            require(key >= 0, 'Negative artifact selector')
        result = result[key]
    require(isinstance(result, dict), 'Selected artifact must be an object')
    return result


def weight_range(value):
    require(isinstance(value, (tuple, list)) and len(value) == 2
            and all(type(x) is int for x in value)
            and 1 <= value[0] <= value[1] <= 8192, 'Invalid assigned weight range')
    return tuple(value)


def band_box(interval):
    return Box((weight_range(interval), (0, 4096), (0, 4096)))


def validate_manifest(manifest):
    require(type(manifest.get('schema_version')) is int and manifest['schema_version'] == 1 and manifest.get('construction') == SWAP_CONSTRUCTION
            and manifest.get('trust_basis') == TRUST_BASIS, 'Wrong manifest version, construction, or trust model')
    if 'target' in manifest:
        expected = dict(source_dimension=1048576, output_length=2097152, rate='1/2', field_order=256,
                        forbidden_weight=734003, minimum_margin_bits=40, group=4, inner_rank=4)
        require(all(manifest['target'].get(key) == value for key, value in expected.items()),
                'Manifest target differs from the established construction and40-bit test')
    entries = manifest['entries']
    require(isinstance(entries, list), 'Manifest entries must be a list')
    ids, assigned = set(), []
    for entry in entries:
        require(isinstance(entry.get('id'), str) and entry['id'] and entry['id'] not in ids, 'Unique entry IDs required')
        ids.add(entry['id'])
        require(entry.get('kind') in KINDS and entry.get('observed_fresh_dual_replay') is True,
                'Unsupported kind or missing observed-replay authorization')
        require(set(entry.get('receipts', {})) == {'256', '384'}, 'Both precision artifacts are required')
        box = band_box(entry['weight_range'])
        require(all(box.intersection(other) is None for other in assigned), 'Assigned weight bands overlap')
        assigned.append(box)
        if entry['kind'] == 'terminal':
            require(box == band_box((7649, 8192)), 'The established terminal range is fixed')
        if entry['kind'] == 'sparse-family':
            require(box.ranges[0][0] == 1 and box.ranges[0][1] <= 64, 'Sparse family is scoped to1..64')
        if entry['kind'] == 'sparse-batch':
            require(type(entry.get('maximum')) is int and box.ranges[0][1] <= entry['maximum'] <= 256
                    and 0 < exact(entry['z']) < 1 and 0 < exact(entry['alpha']) <= 1, 'Invalid sparse batch')
        if entry['kind'] == 'checkpoint':
            require(entry.get('checkpoint_kind') in ('sparse', 'tail', 'dense'), 'Unknown checkpoint kind')
        checked_outer_measure(entry.get('outer_measure', DEFAULT_OUTER_MEASURE))
        require(entry['kind'] == 'potential-band' or entry.get('outer_measure', DEFAULT_OUTER_MEASURE) == DEFAULT_OUTER_MEASURE,
                'Only the established potential path supports an explicit sharper density')
    return entries


def check_receipt(receipt, construction=SWAP_CONSTRUCTION):
    require(receipt.get('construction') == construction and receipt.get('field_order') == 256
            and receipt.get('forbidden_weight') == 734003 and receipt.get('directed') is True
            and receipt.get('covers_all_even_odd_splits') is True
            and receipt.get('whole_code_certificate') is False, 'Receipt construction, target, or scope mismatch')
    ranges = validate_box(receipt['ranges'])
    require(all(type(n) is int for row in receipt['ranges'] for n in row), 'Integer profile intervals required')
    if 'integer_profile_triples' in receipt:
        require(receipt['integer_profile_triples'] == profile_count(ranges), 'Incorrect declared profile count')
    exact(receipt['log2_upper_rational'])
    return ranges


def check_dimensions(receipt):
    require(receipt.get('pair_swaps') is True and receipt.get('source_dimension') == 1048576
            and receipt.get('output_length') == 2097152 and receipt.get('outer_groups') == 8192,
            'Wrong scalar-code dimensions or pair-swap scope')


def log_union(values, precision):
    with ctx.workprec(precision):
        total = sum((arb(2)**rational(value) for value in values), arb(0))
        require(total > 0 and total.is_finite(), 'Empty or nonfinite probability union')
        return str(rational_ceiling(total.log()/arb(2).log()))


def expected_potential_witness(proposal, entry, precision):
    fixed = proposal.get('fixed_witness', {})
    require(proposal.get('construction') == SWAP_CONSTRUCTION and fixed.get('construction') == SWAP_CONSTRUCTION
            and fixed.get('pair_swaps') is True and fixed.get('orientation_aggregation') == 'independent-fair-swaps',
            'Potential proposal is not the independent-fair-swap construction')
    require(weight_range(proposal['requested_weight_range']) == weight_range(entry['weight_range']), 'Proposal range differs from assignment')
    source = fixed['source_exact_witness']
    z, alpha, potential = (exact(source[key]) for key in ('z', 'alpha', 'potential'))
    require(0 < z < 1 and 0 < alpha <= 1 and potential > 0, 'Invalid exact potential witness')
    markers = []
    for row in source['markers']:
        row = tuple(map(exact, row.split()))
        require(len(row) == 2 and min(row) > 0, 'Invalid exact coefficient marker')
        if row not in markers:
            markers.append(row)
    require(1 <= len(markers) <= 100, 'Invalid marker bank size')
    mode = checked_outer_measure(entry['outer_measure'])
    require(checked_outer_measure(fixed.get('outer_measure', DEFAULT_OUTER_MEASURE)) == mode,
            'Saved outer measure differs from explicitly assigned measure')
    result = dict(z=str(z), alpha=str(alpha), markers=[list(map(str, row)) for row in markers],
                  potential=['1', str(potential)], grid=fixed.get('grid', 128), precision=precision,
                  packet_mantissa_bits=60, regional_minimum='after-histogram-normalization')
    grid = result['grid']
    require(type(grid) is int and grid >= 64 and not grid & (grid-1), 'Invalid replay grid')
    if mode != DEFAULT_OUTER_MEASURE:
        result.update(outer_measure=mode, outer_measure_bound=bound_description(mode))
    return result


def check_weight_cache(receipt, entry, proposal, precision):
    """Validate every requested weight; never interpolate or trust a union field."""
    box = band_box(entry['weight_range'])
    require(Box(check_receipt(receipt)) == box, 'Receipt range differs from assignment')
    first, last = box.ranges[0]
    potential = entry['kind'] == 'potential-band'
    if potential:
        check_dimensions(receipt)
        require(receipt.get('precision') == precision and receipt.get('every_integer_weight_replayed') is True
                and receipt.get('covers_all_supports_at_each_weight') is True
                and receipt.get('support_scope') == 'scalar-code active input groups'
                and receipt.get('support_count_method') == 'common-potential-run-convolution'
                and receipt.get('coefficient_method') == 'fixed-normalized-directed-alias-bank', 'Incomplete potential receipt scope')
        require(receipt.get('witness') == expected_potential_witness(proposal, entry, precision), 'Exact potential witness identity mismatch')
        if 'proposal_sha256' not in receipt:
            require(entry.get('legacy_missing_proposal_pin') is True,
                    'Legacy receipt requires explicit authorization for its absent proposal hash')
            # The independently pinned receipt still has to match every exact
            # reconstructed input above. This exception never permits a
            # present-but-wrong hash and does not authenticate arithmetic.
        else:
            require(receipt['proposal_sha256'] == entry['proposal']['sha256'], 'Receipt refers to a different proposal')
        require(receipt.get('outer_measure', DEFAULT_OUTER_MEASURE) == entry['outer_measure'], 'Receipt density mode mismatch')
        rows = receipt['results']
        require(receipt.get('weights_replayed') == last-first+1, 'Wrong replay count')
    else:
        from hd_shifted_swap_sparse import family_witness
        require(receipt.get('pair_swaps') is True
                and receipt.get('scope') == 'all indexed nonempty support sets in the stated weight range', 'Incomplete sparse family scope')
        rows = receipt['weight_receipts']
    require(len(rows) == last-first+1 and all(type(row.get('active_groups')) is int for row in rows)
            and [row.get('active_groups') for row in rows] == list(range(first, last+1)),
            'Missing, duplicate, or reordered integer weights')
    ceilings = []
    for r, row in zip(range(first, last+1), rows):
        if potential:
            require(row.get('regional_histograms') == r//2-max(0, r-4096)+1, 'Wrong regional histogram count')
        else:
            require(Box(check_receipt(row)) == band_box((r, r)) and row.get('covers_all_supports_at_weight') is True,
                    'Sparse row does not cover its whole weight')
            check_dimensions(row)
            require(row.get('coefficient_method') == 'carry-free-Arb-truncated-bivariate', 'Unknown sparse coefficient method')
            if entry['kind'] == 'sparse-family':
                z, alpha = family_witness(r)
                expected = dict(z=str(z), alpha=str(alpha), precision=precision, packet_mantissa_bits=60)
            else:
                maximum = entry['maximum']
                expected = dict(z=str(exact(entry['z'])), alpha=str(exact(entry['alpha'])), precision=precision,
                                packet_mantissa_bits=60, coefficient_maximum=maximum,
                                coefficient_marker=str(F(maximum, 4096)))
            require(row.get('witness') == expected, 'Exact sparse witness identity mismatch')
        ceilings.append(exact(row['log2_upper_rational']))
    if entry['kind'] == 'sparse-batch':
        check_dimensions(receipt)
        require(receipt.get('covers_all_supports_in_range') is True and receipt.get('witness') == expected,
                'Incomplete sparse batch identity')
    # Potential receipts sum unrounded probabilities, whereas sparse receipts
    # sum row ceilings. Recompute both from row ceilings; never reject a valid
    # potential receipt merely because that conservative sum is slightly larger.
    normalized = deepcopy(receipt)
    normalized['log2_upper_rational'] = log_union(ceilings, precision)
    normalized['cached_union_not_used'] = receipt['log2_upper_rational']
    return normalized, tuple(ceilings)


def expected_plan_witness(proposal, index, precision, swapped):
    w = proposal['witnesses'][index]
    if swapped:
        require(w.get('pair_swaps') is True and w.get('orientation_aggregation') == 'independent-fair-swaps',
                'Plan witness is not independently swapped')
    markers = [[str(F(float(2**value)).limit_denominator(10**9)) for value in row]
               for row in w['marker_log_odds']]
    grid = w.get('grid')
    method = 'positive' if grid is None else 'positive-alias'
    require(w.get('coefficient_method', method) == method, 'Inconsistent plan coefficient method')
    result = dict(index=index, z=str(F(w['mark'])), alpha=str(F(w['alpha'])), markers=markers,
                  precision=precision, packet_mantissa_bits=60, grid=grid,
                  uniform_aliases=bool(w.get('uniform_aliases', False)),
                  run_markers=None if w.get('run_markers') is None else list(map(lambda x: str(F(x)), w['run_markers'])),
                  coefficient_method=method)
    if swapped:
        result.update(pair_swaps=True, orientation_aggregation='independent-fair-swaps')
    return result


def check_plan_cache(cache, entry, proposal, precision):
    """Reconstruct source geometry from leaves, not from saved completeness flags."""
    swapped = entry['kind'] == 'swap-plan'
    construction = SWAP_CONSTRUCTION if swapped else CONSTRUCTION
    require((proposal.get('construction') if swapped else proposal.get('construction', construction)) == construction
            and proposal.get('profile_coordinates') == ['r', 'a', 'b']
            and proposal.get('all_even_odd_splits') is True, 'Wrong unsigned plan identity')
    require(cache.get('proposal_sha256') == entry['proposal']['sha256'], 'Plan cache proposal hash mismatch')
    proposed = {(Box(row['ranges']), row['witness_id']) for row in proposal['accepted_boxes']}
    require(len(proposed) == len(proposal['accepted_boxes']), 'Duplicate proposed leaf')
    require(cache.get('failed') == [] and cache.get('all_proposed_leaves_replayed') is True, 'Plan replay has failed or unobserved leaves')
    leaves = cache['accepted']
    actual = {(Box(row['ranges']), row['witness_index']) for row in leaves}
    require(actual == proposed and len(actual) == len(leaves), 'Cached leaves differ from replayed proposal')
    indices = sorted({index for _, index in proposed})
    witnesses = {row['index']: row for row in cache['witnesses']}
    require(len(witnesses) == len(cache['witnesses']) and set(witnesses) == set(indices), 'Incomplete plan witnesses')
    for index in indices:
        require(witnesses[index] == expected_plan_witness(proposal, index, precision, swapped), 'Exact plan witness identity mismatch')
    source = SwapLedger(band_box(proposal['requested_weight_range']))
    signatures = []
    for row in leaves:
        ranges = check_receipt(row, construction)
        require(row.get('support_count_method') in UNSIGNED_METHODS, 'Unestablished plan support-count method')
        if swapped:
            source.accept(row)
        else:
            source.accept_inherited(row)
        signatures.append((ranges, row['witness_index'], exact(row['log2_upper_rational'])))
    summary = source.band_summary(*entry['weight_range'], precision)
    require(summary['complete_local_cover'], 'Plan does not cover its assigned band')
    return source, tuple(sorted(signatures))


def leaf_signature(row):
    return (tuple(map(tuple, row['ranges'])), row.get('witness_index'),
            row.get('support_count_method'), exact(row['log2_upper_rational']))


def check_checkpoint(cache, entry, base, precision):
    """Inspect exported leaves; input metadata reconstruction does no arithmetic."""
    require(cache.get('construction') == SWAP_CONSTRUCTION and cache.get('precision') == precision
            and cache.get('kind') == entry['checkpoint_kind'] and cache.get('fresh_arithmetic') is True
            and cache.get('whole_code_certificate') is False, 'Checkpoint identity mismatch')
    root = band_box(entry['weight_range'])
    source = SwapLedger(root)
    signatures, reconstructed = [], []
    if entry['checkpoint_kind'] == 'sparse':
        from hd_shifted_swap_sparse import BATCH_SCHEDULE, BATCH_SCHEDULE_256
        schedule = {(a, b): (z, alpha, maximum) for maximum, rows in
                    ((128, BATCH_SCHEDULE), (256, BATCH_SCHEDULE_256)) for a, b, z, alpha in rows}
        for row in cache['accepted']:
            interval = tuple(row['ranges'][0])
            child = dict(weight_range=interval, kind='sparse-family')
            if interval != (1, 64):
                require(interval in schedule, 'Unknown frozen sparse checkpoint interval')
                z, alpha, maximum = schedule[interval]
                child.update(kind='sparse-batch', z=str(z), alpha=str(alpha), maximum=maximum)
            normalized, ceilings = check_weight_cache(row, child, None, precision)
            source.accept(normalized)
            signatures.append((interval, ceilings))
    else:
        pinned_bytes(base/entry['replay_source_pin']['path'], entry['replay_source_pin']['sha256'])
        require(Path(entry['replay_source_pin']['path']).name == 'hd_shifted_cover_ledger.py',
                'Pin the implementation that reconstructs replay inputs')
        actual_replay_hash = sha256(Path(replay_implementation.__file__).read_bytes()).hexdigest()
        require(entry['replay_source_pin']['sha256'] == actual_replay_hash == REPLAY_RECIPE_SOURCE_SHA256,
                'Replay source is not the reviewed implementation of the reconstructed recipe')
        expected_leaves, plans, bands = {}, {}, []
        for spec in entry['plan_bands']:
            require(spec['kind'] in ('swap-plan', 'unsigned-plan') and spec['label'] not in plans,
                    'Unique named established plan bands required')
            band = band_box(spec['weight_range'])
            require(root.intersection(band) == band and all(band.intersection(other) is None for other in bands),
                    'Invalid or overlapping checkpoint plan assignment')
            bands.append(band)
            proposal = read_artifact(base, spec['proposal'])
            swapped = spec['kind'] == 'swap-plan'
            construction = SWAP_CONSTRUCTION if swapped else CONSTRUCTION
            require((proposal.get('construction') if swapped else proposal.get('construction', construction)) == construction
                    and proposal.get('profile_coordinates') == ['r', 'a', 'b']
                    and proposal.get('all_even_odd_splits') is True, 'Checkpoint proposal construction mismatch')
            require(cache.get('proposal_sha256', {}).get(spec['label']) == spec['proposal']['sha256'],
                    'Checkpoint proposal pin mismatch')
            require(band_box(proposal['requested_weight_range']).intersection(band) == band,
                    'Assignment extends beyond original proposal root')
            plans[spec['label']] = spec
            for row in proposal['accepted_boxes']:
                piece = Box(row['ranges']).intersection(band)
                if piece is None or not profile_count(piece.ranges):
                    continue
                key = (piece.ranges, row['witness_id'])
                require(key not in expected_leaves, 'Duplicate checkpoint proposed leaf')
                expected_leaves[key] = spec['label']
            indices = sorted({row['witness_id'] for row in proposal['accepted_boxes']})
            witnesses = [expected_plan_witness(proposal, index, precision, swapped) for index in indices]
            saved_sources = [item for item in cache.get('sources', []) if item.get('label') == spec['label']]
            require(len(saved_sources) <= 1, 'Duplicate checkpoint source label')
            if saved_sources:
                saved = saved_sources[0].get('replay', {}).get('witnesses')
                if saved is not None:
                    require(sorted(saved, key=lambda x: x['index']) == witnesses, 'Saved checkpoint witnesses disagree with reconstruction')
            reconstructed.append(dict(label=spec['label'], witnesses=witnesses))
        observed = set()
        selected = cache.get('source_accepted', cache['accepted'])
        kept = []
        for original in selected:
            piece = Box(original['ranges']).intersection(root)
            if piece is None or not profile_count(piece.ranges):
                continue
            row = deepcopy(original)
            row['ranges'] = piece.ranges
            row['integer_profile_triples'] = profile_count(piece.ranges)
            check_receipt(row)
            if row.get('method') == 'inverse-image-terminal-range':
                from hd_shifted_cover_ledger import terminal_receipt
                fresh = terminal_receipt(precision)
                require(piece == band_box((7649, 8192)) and row.get('inherited_from') == CONSTRUCTION
                        and row.get('inheritance') == 'route-independent-inverse-image-bound', 'Invalid checkpoint terminal transfer')
                require(exact(row['log2_upper_rational']) == exact(fresh['log2_upper_rational']), 'Terminal arithmetic mismatch')
            else:
                key = (piece.ranges, row['witness_index'])
                require(key in expected_leaves and key not in observed, 'Unexpected or duplicate checkpoint leaf')
                label = expected_leaves[key]
                require(row.get('source_label') == label, 'Checkpoint leaf source identity mismatch')
                require(row.get('support_count_method') in UNSIGNED_METHODS, 'Unknown checkpoint support method')
                if plans[label]['kind'] == 'unsigned-plan':
                    require(row.get('inherited_from') == CONSTRUCTION
                            and row.get('inheritance') == 'independent-swap-unsigned-coefficient-domination',
                            'Unsigned checkpoint is missing its established transfer')
                else:
                    require('inherited_from' not in row, 'Swap plan must not be relabeled as inherited')
                observed.add(key)
            source.accept(row)
            kept.append(row)
        require(observed == set(expected_leaves), 'Checkpoint omitted proposed leaves')
        # A second import may change labels but must not alter its retained
        # event ranges, methods, or whole-leaf bounds.
        require(sorted(map(leaf_signature, kept)) == sorted(map(leaf_signature, cache['accepted'])),
                'Exported checkpoint leaves differ from source leaves')
        signatures = sorted(map(leaf_signature, kept))
    require(source.summary(precision)['complete_local_cover'], 'Checkpoint contains unresolved feasible profiles')
    return source, tuple(signatures), reconstructed


def _verify_entries(manifest, base, *, precision=384):
    """Internal entry point; callers must authenticate manifest bytes first."""
    entries = validate_manifest(manifest)
    require(type(precision) is int and precision >= 128, 'Invalid aggregation precision')
    for reference in manifest.get('implementation_pins', []):
        pinned_bytes(base / reference['path'], reference['sha256'])
    ledger, records = SwapLedger(), []
    for entry in entries:
        kind = entry['kind']
        proposal = read_artifact(base, entry['proposal']) if kind in ('potential-band', 'unsigned-plan', 'swap-plan') else None
        dual = []
        for bits in (256, 384):
            cached = read_artifact(base, entry['receipts'][str(bits)])
            if kind in ('potential-band', 'sparse-family', 'sparse-batch'):
                dual.append(check_weight_cache(cached, entry, proposal, bits))
            elif kind in ('unsigned-plan', 'swap-plan'):
                dual.append(check_plan_cache(cached, entry, proposal, bits))
            elif kind == 'checkpoint':
                dual.append(check_checkpoint(cached, entry, base, bits))
            else:
                require(Box(check_receipt(cached, CONSTRUCTION)) == band_box((7649, 8192))
                        and cached.get('method') == 'inverse-image-terminal-range'
                        and cached.get('precision') == bits, 'Invalid terminal receipt')
                # The terminal calculation is tiny; verify its arithmetic anew.
                from hd_shifted_cover_ledger import terminal_receipt
                require(cached == json.loads(json.dumps(terminal_receipt(bits))), 'Terminal receipt differs from fresh computation')
                dual.append((cached, exact(cached['log2_upper_rational'])))
        require(dual[0][1] == dual[1][1], '256/384-bit rational ceilings disagree')
        selected = dual[1][0]
        if kind in ('unsigned-plan', 'swap-plan', 'checkpoint'):
            imported = ledger.import_band(selected, *entry['weight_range'], label=entry['id'])
        elif kind == 'terminal':
            ledger.accept_inherited(selected)
            imported = ledger.band_summary(*entry['weight_range'], precision)
        else:
            ledger.accept(selected)
            imported = ledger.band_summary(*entry['weight_range'], precision)
        records.append(dict(id=entry['id'], kind=kind, weight_range=entry['weight_range'],
                            witness_metadata_reconstructed=kind == 'checkpoint' and entry['checkpoint_kind'] != 'sparse',
                            legacy_missing_proposal_pin_authorized=entry.get('legacy_missing_proposal_pin') is True,
                            dual_ceiling_agreement=True, **{k: imported[k] for k in
                            ('accepted_profile_triples', 'accepted_union_log2_upper', 'complete_local_cover')}))
    summary = ledger.summary(precision)
    # Weight-level gaps are recomputed independently of all saved scope flags.
    covered = {r for entry in entries for r in range(entry['weight_range'][0], entry['weight_range'][1]+1)}
    missing = sorted(set(range(1, 8193))-covered)
    require(summary['complete_local_cover'] == (not missing), 'Weight and profile coverage disagree')
    return dict(summary=summary, entries=records, missing_weights=missing,
                verification_mode='trusted-replay-cache-with-fresh-directed-aggregation',
                trust_basis=TRUST_BASIS, arithmetic_validity_from_hashes=False,
                cached_arithmetic_requires_trusted_observed_replays=True,
                directed_aggregation=True, aggregation_precision=precision,
                whole_code_certificate=bool(not missing and summary['whole_code_certificate']))


def load_manifest(path, expected_sha256):
    path = Path(path).resolve()
    manifest = strict_json(pinned_bytes(path, expected_sha256))
    validate_manifest(manifest)
    return path.parent, manifest


def verify_manifest(path, expected_sha256, *, precision=384):
    base, manifest = load_manifest(path, expected_sha256)
    return dict(_verify_entries(manifest, base, precision=precision), manifest_sha256=expected_sha256)


def fresh_manifest(path, expected_sha256, *, precision=256, progress=None):
    """Recompute from pinned proposals/schedules through the established API."""
    from hd_shifted_swap_ledger import replay_swap_cover
    base, manifest = load_manifest(path, expected_sha256)
    unsigned, swapped, potential, batches = [], [], [], []
    sparse_last = 0
    for entry in manifest['entries']:
        kind, interval = entry['kind'], entry['weight_range']
        if kind in ('potential-band', 'unsigned-plan', 'swap-plan'):
            proposal = read_artifact(base, entry['proposal'])
            if kind == 'potential-band':
                expected_potential_witness(proposal, entry, precision)
                potential.append((entry['id'], proposal))
            else:
                (swapped if kind == 'swap-plan' else unsigned).append((entry['id'], proposal, *interval))
        elif kind == 'sparse-family':
            sparse_last = interval[1]
        elif kind == 'sparse-batch':
            batches.append((*interval, exact(entry['z']), exact(entry['alpha']), entry['maximum']))
        elif kind == 'checkpoint':
            if entry['checkpoint_kind'] == 'sparse':
                from hd_shifted_swap_sparse import BATCH_SCHEDULE, BATCH_SCHEDULE_256
                require(tuple(interval) == (1, 256), 'Fresh sparse checkpoint is the frozen1..256 schedule')
                sparse_last = 64
                batches.extend((a, b, z, alpha, maximum) for maximum, rows in
                               ((128, BATCH_SCHEDULE), (256, BATCH_SCHEDULE_256)) for a, b, z, alpha in rows)
            else:
                for spec in entry['plan_bands']:
                    proposal = read_artifact(base, spec['proposal'])
                    (swapped if spec['kind'] == 'swap-plan' else unsigned).append(
                        (spec['label'], proposal, *spec['weight_range']))
    ledger, result = replay_swap_cover(unsigned, swapped, sparse_last=sparse_last, sparse_batches=batches,
        potential_bands=potential, potential_max_weights=max((b['requested_weight_range'][1]-b['requested_weight_range'][0]+1
                                                            for _, b in potential), default=33),
        precision=precision, threshold=-40, progress=progress)
    assigned = {r for entry in manifest['entries'] for r in range(entry['weight_range'][0], entry['weight_range'][1]+1)}
    missing = sorted(set(range(1, 8193))-assigned)
    if missing:
        # The replay API may freshly add a tail/endpoint absent from the
        # manifest; do not turn an incomplete manifest into a complete claim.
        result['summary']['whole_code_certificate'] = False
    # The existing API also supplies its fixed terminal and possibly endpoint.
    return dict(result, manifest_sha256=expected_sha256, verification_mode='fresh-directed-replay',
                cached_bounds_used=False, precision=precision,
                missing_manifest_weights=missing,
                accepted=ledger.accepted, unresolved=[box.ranges for box in ledger.unresolved])


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--expected-manifest-sha256', required=True,
                        help='Reviewed trust anchor supplied separately, never read from the manifest')
    parser.add_argument('--precision', type=int, default=384, help='Fresh aggregation precision')
    parser.add_argument('--fresh', action='store_true', help='Costly: independently recompute both256 and384 from proposals')
    parser.add_argument('--output', type=Path, help='Exclusive-create a new JSON verification receipt')
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error('Output exists; frozen verification receipts are never overwritten')
    if args.fresh:
        result = {str(bits): fresh_manifest(args.manifest, args.expected_manifest_sha256, precision=bits,
                    progress=lambda row: print(json.dumps(row), flush=True)) for bits in (256, 384)}
    else:
        result = verify_manifest(args.manifest, args.expected_manifest_sha256, precision=args.precision)
    if args.output:
        with args.output.open('x', encoding='utf-8') as target:
            json.dump(result, target, separators=(',', ':'), allow_nan=False)
            target.write('\n')
        print(json.dumps(dict(path=str(args.output.resolve()),
                              sha256=sha256(args.output.read_bytes()).hexdigest(),
                              verification_mode='fresh-directed-replay' if args.fresh else result['verification_mode'])),
              flush=True)
    else:
        print(json.dumps(result, allow_nan=False), flush=True)
