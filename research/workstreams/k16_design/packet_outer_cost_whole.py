"""Fresh complete proof replay for an explicit larger outer and t64 inner.

A proposal file contributes only positive rational tilts. Every selected
occupancy is recomputed with freshly prepared maps and outward Arb arithmetic.
The exact without-replacement prefix and conditioned-iid suffix are disjoint.
Their dyadic endpoints are added exactly before one final upward rounding.

Each group consists of parallel rate-one-half RS rows with common evaluation
points. An aligned symbol tuple is independently randomized by an invertible
binary map whose action on every fixed nonzero tuple is uniform nonzero.
Uniform GL matrices and multiplication by an independent uniform nonzero
field element both satisfy this condition. For every fixed message, symbol
images are independent across positions and groups. No joint independence
between the images of different messages is used by the first moment.

For r parallel [n,k] RS rows over GF(2^b), the aligned alphabet has Q=2^(br)
elements. A shortened code supported on h coordinates has a common information
set of size h-(n-k). Projection injects every fully nonzero aligned word into
nonzero tuples on that set, giving at most (Q-1)^(h-(n-k)) such words.
Uniform nonzero symbol images therefore give the same pointwise majorant
beta=Q^n/(Q-1)^(n-k) as an MDS code over GF(Q). The actual zero word is omitted;
the dominating fully uniform input deliberately includes artificial zero.

The 512-to-1024 numerical geometry uses the larger of two valid beta values:
four GF256 RS[32,16] rows with32-bit symbol maps, or sixteen GF16 RS[16,8]
rows with64-bit symbol maps. The latter beta is smaller. Thus one uniform
majorant replay bounds both explicitly declared ensembles. It does not
claim that their sampled code distributions are equal.

The t64 inner emits x+A*a and updates a to M*a+C*x. Its recorded fixed maps
come from the selected quadratic-prefix construction. Per-step M is independent
and sends each fixed nonzero state uniformly over the nonzero states. The
retained uniform-GL kernel requires only this action law for a fixed message.
State starts at zero, persists through every step and region, and is not
flushed. Results concern the ideal setup ensemble, not every seed.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
from time import monotonic

from flint import arb, ctx
import packet_q1 as q1
from packet_outer_cost_search import Model, prepare, geometry
import packet_outer_cost_dense as dense
from rs_uniform_envelope import UniformInputEnvelope


DEFAULT_DENSE_TILTS = tuple(map(Q, ('.03125', '.0625', '.09375', '.125', '.15625',
    '.1875', '.25', '.375', '.5', '.75', '1', '1.5', '2.1972246')))
DEFAULT_MARKERS = tuple(Q(j, 64) for j in range(1, 65))
REQUIRED_SOURCES = ('packet_outer_cost_whole.py', 'packet_outer_cost_dense.py',
    'packet_outer_cost_search.py', 'rs_uniform_envelope.py', 'packet_regional_power.py',
    'packet_uniform_tail.py', 'packet_rs_k20_dense.py', 'packet_rs_state_sparse.py')
PROOF_NOTE = Path(__file__).resolve().with_name('TRANSITIVE_SYMBOL_MAPS.md')


def outer_options(group_dimension):
    """Check base-field RS feasibility and domination for each declared outer."""
    if type(group_dimension) is not int or group_dimension not in (256, 512):
        raise ValueError('outer dimension must be256 or512')
    declarations = ((16, 8, 8, 4), (16, 8, 4, 8)) if group_dimension == 256 else (
        (32, 16, 8, 4), (16, 8, 4, 16))
    records = []
    envelopes = []
    for n, k, field_bits, rows in declarations:
        symbol_bits = field_bits*rows
        if n > 1 << field_bits or symbol_bits % 4:
            raise ArithmeticError('declared parallel base-field evaluation code is infeasible')
        env = UniformInputEnvelope(n, k, 4, symbol_bits//4)
        if env.message_bits != group_dimension:
            raise ArithmeticError('parallel outer has the wrong input dimension')
        envelopes.append(env)
        records.append(dict(n=n, k=k, base_field_bits=field_bits, parallel_rows=rows,
            symbol_bits=symbol_bits, evaluation_points='common n distinct base-field elements',
            symbol_maps='independent invertible maps with uniform nonzero fixed-input images',
            supported_symbol_maps=['uniform GL overF2',
                f'uniform nonzero GF(2^{symbol_bits}) multiplication',
                f'binary adjoints of uniform nonzero GF(2^{symbol_bits}) multiplication'],
            beta=str(env.beta)))
    reference = envelopes[0]
    if any(env.regions != reference.regions or env.beta > reference.beta for env in envelopes):
        raise ArithmeticError('all declared outers must be dominated at the same geometry')
    return reference, records


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':')).encode()


def _digest(raw):
    return hashlib.sha256(raw).hexdigest()


def _check_sources(sources):
    here = Path(__file__).resolve().parent
    required = {str((here/name).resolve()) for name in REQUIRED_SOURCES}
    if not required <= sources.keys():
        raise ValueError('whole replay is missing required mathematical source pins')
    for filename, digest in sources.items():
        if _digest(Path(filename).read_bytes()) != digest:
            raise RuntimeError(f'mathematical source changed: {filename}')


def exact_union(endpoints, precision=256):
    """Add positive dyadics exactly, then enclose only the final conversion."""
    endpoints = tuple(tuple(pair) for pair in endpoints)
    if (not endpoints or any(len(pair) != 2 or any(type(x) is not int for x in pair)
            or pair[0] <= 0 or abs(pair[1]) > 10_000_000 for pair in endpoints)):
        raise ValueError('nonempty positive finite dyadic endpoints required')
    if type(precision) is not int or precision < 256:
        raise ValueError('at least256-bit final union precision required')
    exponent = min(e for _, e in endpoints)
    mantissa = sum(m << (e-exponent) for m, e in endpoints)
    old = ctx.prec
    try:
        ctx.prec = precision
        upper = q1.kernel_t64.up(arb(mantissa)*arb(2)**exponent)
        return q1.endpoint(upper), str(-upper.log()/arb(2).log()), bool(upper < arb(2)**-40)
    finally:
        ctx.prec = old


def _proposal(path, *, K, bits, reference, prefix_end):
    raw = Path(path).read_bytes()
    data = json.loads(raw)
    if (data.get('proposal_only') is not True or data.get('whole_code_certificate') is not False
            or data.get('K') != K or data.get('state_bits') != bits
            or data.get('envelope') != reference.metadata()):
        raise ValueError('matching explicit outer, K, state size and proposal-only record required')
    best = data.get('best')
    if not isinstance(best, dict) or any(str(q) not in best for q in range(1, prefix_end+1)):
        raise ValueError('every exact-prefix occupancy needs a tilt proposal')
    groups = defaultdict(list)
    for q in range(1, prefix_end+1):
        tilt = Q(best[str(q)]['tilt'])
        if tilt <= 0:
            raise ValueError('positive rational tilt proposals required')
        groups[tilt].append(q)
    # No supplied margin, map profile, endpoint, or source hash is trusted.
    return dict(groups), dict(path=str(Path(path).resolve()), sha256=_digest(raw),
        role='positive rational tilt proposals only; all numeric endpoints recomputed')


def _write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2)
        handle.write('\n')


def replay(proposal, output, *, K=1 << 20, bits=16, group_dimension=512,
        prefix_end=221, precision=256, dense_tilts=DEFAULT_DENSE_TILTS,
        markers=DEFAULT_MARKERS, dense_proposal=None, progress=None):
    reference, options = outer_options(group_dimension)
    meta = geometry(K, reference)
    if type(prefix_end) is not int or not 1 <= prefix_end < meta['groups']:
        raise ValueError('proper nonempty exact prefix required')
    if type(precision) is not int or precision < 256:
        raise ValueError('whole replay requires at least256-bit interval precision')
    paths = [Path(output), Path(str(output)+'.exact.json'), Path(str(output)+'.dense.json')]
    if any(path.exists() for path in paths):
        raise ValueError('all three output files must be fresh')
    grouped, proposal_source = _proposal(proposal, K=K, bits=bits,
        reference=reference, prefix_end=prefix_end)
    dense_proposal_source = None
    if dense_proposal is not None:
        raw = Path(dense_proposal).read_bytes()
        proposal_record = json.loads(raw)
        if (proposal_record.get('K') != K or proposal_record.get('state_bits') != bits
                or proposal_record.get('beta') != str(reference.beta)
                or proposal_record.get('groups') != meta['groups']
                or proposal_record.get('regions') != meta['regions']
                or proposal_record.get('whole_code_certificate') is not False):
            raise ValueError('matching dense tilt/marker proposal required')
        dense_tilts = tuple(map(Q, proposal_record['tilts']))
        markers = tuple(map(Q, proposal_record['marker_probabilities']))
        dense_proposal_source = dict(path=str(Path(dense_proposal).resolve()),
            sha256=_digest(raw), role='rational tilts and marker probabilities only; endpoints discarded')
    started = monotonic()
    proof_note = dict(path=str(PROOF_NOTE), sha256=_digest(PROOF_NOTE.read_bytes()))
    model = Model(*prepare(bits), precision=precision)
    _check_sources(model.sources)
    values, choices, trials = {}, {}, []
    for tilt in sorted(grouped):
        result = model.replay(K=K, envelope=reference, occupancies=grouped[tilt], tilt=tilt)
        for q, entry in result['values'].items():
            values[q], choices[q] = entry['upper'], str(tilt)
        trials.append(dict(tilt=str(tilt), occupancies=grouped[tilt],
            canonical_component_sha256=_digest(_canonical(result))))
        if progress is not None:
            progress(dict(stage='exact', tilt=str(tilt), completed=len(values),
                prefix_end=prefix_end, elapsed=monotonic()-started))
    exact = dict(meta, schema='larger-outer-exact-prefix-component-1',
        state_bits=bits, physical_t=64, precision=precision,
        map_record=model.map_record, source_sha256=model.sources,
        proof_note_source=proof_note,
        occupancy_covered=list(range(1, prefix_end+1)), occupancy_uppers=values,
        occupancy_choices=choices, trials=trials, fresh_computation=True,
        whole_code_certificate=False, proposal_source=proposal_source)
    _write(paths[1], exact)

    def dense_progress(tilt, best, selected, elapsed):
        if progress is not None:
            progress(dict(stage='dense', tilt=str(tilt), elapsed=monotonic()-started))

    tail = dense.run(model.data, model.map_record, K=K, envelope=reference,
        occupancies=range(prefix_end+1, meta['groups']+1), tilts=dense_tilts,
        marker_probabilities=markers, precision=precision, progress=dense_progress)
    if tail['source_sha256'] != model.sources or tail['map_record'] != model.map_record:
        raise ValueError('exact and dense components differ in maps or mathematical sources')
    _write(paths[2], tail)
    values = dict(values)
    if set(values) & set(tail['occupancy_uppers']):
        raise ValueError('exact prefix and dense suffix overlap')
    values.update(tail['occupancy_uppers'])
    if set(values) != {str(q) for q in range(1, meta['groups']+1)}:
        raise ValueError('every nonzero occupancy must be covered exactly once')
    upper, margin, passes = exact_union(values.values(), precision)
    whole = dict(meta, schema='larger-outer-fresh-whole-1', outer_options=options,
        state_bits=bits, physical_t=64, precision=precision,
        inner_distribution='independent invertible maps with uniform nonzero fixed-input state images',
        routing='independent uniform within-group packet permutations and within-region group permutations',
        map_record=model.map_record, source_sha256=model.sources,
        proof_note_source=proof_note,
        occupancy_covered=[1, meta['groups']], all_occupancies_covered=True,
        occupancy_uppers=values, union_upper=upper, margin_bits=margin,
        target_margin_bits=40, target_minimum_distance=meta['cutoff']+1,
        target_met=passes, whole_code_certificate=passes, fresh_replay=True,
        components=[dict(path=str(path.resolve()), sha256=_digest(path.read_bytes())) for path in paths[1:]],
        proposal_source=proposal_source, dense_proposal_source=dense_proposal_source,
        elapsed_seconds=monotonic()-started,
        arithmetic='Exact sum of disjoint positive dyadic endpoints, followed by one upward Arb rounding.',
        scope='Complete first-moment bound for either explicitly listed ideal outer ensemble and '
              'the recorded inner maps. It is not a guarantee for every sampled setup or seed.')
    model.authenticate()
    _check_sources(model.sources)
    if _digest(PROOF_NOTE.read_bytes()) != proof_note['sha256']:
        raise RuntimeError('the transitive symbol-map proof note changed during replay')
    _write(paths[0], whole)
    return whole


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('proposal')
    parser.add_argument('--output', required=True)
    parser.add_argument('--K', type=int, default=1 << 20)
    parser.add_argument('--bits', type=int, default=16)
    parser.add_argument('--outer', type=int, choices=(256, 512), default=512)
    parser.add_argument('--prefix-end', type=int, default=221)
    parser.add_argument('--dense-proposal')
    args = parser.parse_args()
    result = replay(args.proposal, args.output, K=args.K, bits=args.bits,
        group_dimension=args.outer, prefix_end=args.prefix_end,
        dense_proposal=args.dense_proposal,
        progress=lambda value: print(json.dumps(value), flush=True))
    print(json.dumps({key: result[key] for key in ('K', 'state_bits', 'margin_bits', 'whole_code_certificate')}), flush=True)
