"""Fresh finite-length components for the unchanged RS16 packet ensemble.

For a positive message length K divisible by 4096, let L=K/128. Each of L
independent groups encodes 128 bits into 64 four-bit packets. The outer uses
four GF16 RS[16,8] rows and independent uniform GL16 maps at symbol positions.
Independent uniform group and regional permutations route packets into 64
regions, each containing L/32 ordered 128-bit proof macros.

The selected t64/s16 inner starts at zero and retains its state across every
physical step and region. Each physical step has a fresh independent uniform
GL16 update; the final state is discarded without a flush. This module varies
only K. Its setup distribution and selected expansion/feedback maps do not vary.

Q1 and Q2 reuse exact support-placement averages with expected RS shell counts.
For q>=3, the expected output measure of one active group is dominated by beta
times uniform 256-bit input, where beta=2^256/(2^16-1)^8. This comparison keeps
the active label even for its artificial zero output. Thus a region has
J~Bin(q,15/16) active packets in uniformly chosen distinct slots.

Let R_j(lambda) be the ordered regional average of local transition envelopes,
and R(q,lambda)=E[R_J(lambda)]. At positive tilt lambda the q contribution is
bounded by C(L,q)*beta^q*exp(lambda*floor(2*K/10))*e_zero*R(q,lambda)^64*1.
All terminal coordinates count mass. We minimize these bounds separately for
each requested q. A component covers only its explicit occupancy list; it is
not a whole-code certificate or a guarantee for a particular seeded setup.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from fractions import Fraction as Q
import hashlib
import json
from math import comb
from pathlib import Path
from time import monotonic

from flint import arb, ctx
import packet_q1 as q1
import packet_q2 as q2
import packet_rs_k20_sparse as algebra
import packet_uniform_tail as uniform

TAIL_SCHEMA = 'rs16-length-tail-1'
DEFAULT_TILTS = ('.00032', '.00128', '.00512', '.02048', '.08192')
exact_outer = algebra.exact_outer
count_hash = algebra.count_hash
positive_endpoint = algebra.positive_endpoint
endpoint_value = algebra.endpoint_value
checked_options = algebra.checked_options


def geometry(K):
    """Return the exact half-rate geometry; each region has complete macros."""
    if type(K) is not int or K <= 0 or K % 4096:
        raise ValueError('K must be a positive integer multiple of4096')
    return q1.Geometry(K // 128, 64, 128)


def authenticate_record(map_record):
    """Check selected-map identity without accepting saved spectra as premises.

    No census runs here. A numerical caller must freshly prepare the selected
    maps; a replay coordinator must compare its complete fresh map record.
    """
    path = q1.kernel_t64.SELECTED_MAP.resolve()
    raw = path.read_bytes()
    declared = json.loads(raw)
    rows = [int(value, 16) for value in declared['generator_rows_hex']]
    columns = [sum(((row >> p) & 1) << j for j, row in enumerate(rows)) for p in range(64)]
    source = dict(path=str(path), sha256=hashlib.sha256(raw).hexdigest())
    macro = dict(physical_steps=2, physical_step_bits=64, macro_step_bits=128,
        physical_windows=16, macro_windows=32, state_continuity='retained_between_halves',
        initial_state='zero', flush=False)
    identity = dict(bits=16, width=64, expansion_rows=list(map(hex, rows)),
                    feedback_columns=columns)
    digest = hashlib.sha256(json.dumps(identity, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    if (declared.get('step_bits') != 64 or declared.get('state_bits') != 16
            or len(rows) != 16 or declared.get('columns') != columns
            or not isinstance(map_record, dict)
            or map_record.get('schema') != 's16-selected-t64-fixed-maps-1'
            or map_record.get('t') != 64 or map_record.get('s') != 16
            or map_record.get('expansion_rows_hex') != list(map(hex, rows))
            or map_record.get('feedback_columns') != columns
            or map_record.get('source') != source or map_record.get('macro') != macro
            or map_record.get('distribution') != 'uniform_gl'
            or map_record.get('sampling') != 'independent uniform GL16 for every physical t64 step'
            or map_record.get('map_sha256') != digest
            or map_record.get('whole_code_certificate') is not False):
        raise ValueError('exact selected t64/s16 map record and independent uniform GL16 scope required')
    return digest


def source_snapshot():
    """Pin every loaded local proof module, this wrapper, and the fixed map."""
    sources = q1.source_snapshot()
    path = Path(__file__).resolve()
    sources[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return sources


def prepare(precision=192, data=None, map_record=None):
    """Prepare once, or authenticate shared data from the current fresh prepare.

    Shared data must come from kernel_t64.prepare(birth_density='capped').
    This interface does not load or accept saved transition matrices.
    """
    if type(precision) is not int or precision < 128:
        raise ValueError('integer precision >=128 required')
    if (data is None) != (map_record is None):
        raise ValueError('shared prepared data and map record must be supplied together')
    ctx.prec = precision
    if data is None:
        data, map_record = q1.kernel_t64.prepare(birth_density='capped')
    digest = authenticate_record(map_record)
    physical = q1.kernel_t64.authenticate(data)
    if (data.get('bits') != 16 or data.get('physical_step_bits') != 64
            or data.get('windows') != 32 or data.get('macro_step_bits') != 128
            or data.get('birth_density') != 'capped' or data.get('map_sha256') != digest
            or physical.get('bits') != 16 or physical.get('birth_density') != 'capped'
            or list(map(int, physical['columns'])) != map_record['feedback_columns']
            or [hex(int(physical['map_images'][1 << j])) for j in range(16)]
                != map_record['expansion_rows_hex']):
        raise ValueError('fresh capped-density data for the unchanged selected16 map required')
    return data, map_record, source_snapshot()


def premises():
    return dict(outer='four parallel GF16 RS[16,8] words',
        label_mixing='independent uniform GL16 per aligned four-nibble symbol',
        independent_setups_between_groups=True, exact_expected_shell_counts=True,
        independent_uniform_group_packet_shuffles=True,
        independent_uniform_regional_shuffles=True)


def _finish(record, sources, output):
    if sources != source_snapshot():
        raise RuntimeError('loaded proof source changed during length computation')
    if record.get('whole_code_certificate') is not False:
        raise ArithmeticError('a length component cannot claim a whole-code certificate')
    record['source_sha256'] = sources
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open('x', encoding='utf-8') as stream:
            stream.write(json.dumps(record, indent=2) + '\n')
    return record


def run_low(occupancy, *, K, tilts=DEFAULT_TILTS, precision=192, output=None,
            data=None, map_record=None):
    """Fresh q1 or q2 bound at the requested length, with all supports covered."""
    geom = geometry(K)
    if type(occupancy) is not int or occupancy not in (1, 2):
        raise ValueError('low occupancy must be one or two')
    tilts, output = checked_options(tilts, precision, output)
    _, counts = exact_outer()
    data, map_record, sources = prepare(precision, data, map_record)
    kwargs = dict(tilts=tilts, precision=precision, data=data,
                  map_record=map_record, metadata=premises(), output=None)
    if occupancy == 1:
        record = q1.evaluate_q1(counts, group_count=geom.group_count,
            regions=geom.regions, epochs_per_region=geom.macros_per_region,
            group_dimension=geom.group_dimension, count_kind='shells', **kwargs)
    else:
        record = q2.evaluate_q2(counts, geometry=geom, **kwargs)
    if (record.get('schema') != f'finite-packet-q{occupancy}-diagnostic-1'
            or record.get('geometry') != asdict(geom) or record.get('K') != K
            or record.get('N') != geom.N or record.get('threshold') != geom.N // 10
            or record.get('distance') != '1/10' or record.get('map_record') != map_record
            or record.get('zero_initial_state') is not True or record.get('final_flush') is not False
            or record.get('occupancy_covered') != [occupancy]
            or record.get('count_sha256') != count_hash(counts, cumulative=occupancy == 1)
            or record.get('source_sha256') != sources or ctx.prec != precision):
        raise ArithmeticError('underlying sparse evaluator returned inconsistent length metadata')
    endpoint_value(record[f'q{occupancy}_upper'])
    record.update(outer='rs16', length_component=f'q{occupancy}', fresh_computation=True,
                  macros_per_region=geom.macros_per_region)
    return _finish(record, sources, output)


def occupancy_upper(regional, *, K, occupancy, beta, tilt):
    """Evaluate the full regional expression with upward-rounded arithmetic."""
    geom = geometry(K)
    if type(occupancy) is not int or not 3 <= occupancy <= geom.group_count:
        raise ValueError('tail occupancy must lie in3..L')
    if Q(beta) <= 0 or Q(tilt) <= 0:
        raise ValueError('positive beta and tilt required')
    matrix = uniform.regional_uniform(regional, occupancy) ** geom.regions
    moment = sum((matrix[0, j] for j in range(matrix.ncols())), arb(0))
    upper = q1.kernel_t64.up(comb(geom.group_count, occupancy) *
        q1.kernel_t64.aq(Q(beta)) ** occupancy *
        (q1.kernel_t64.aq(Q(tilt)) * (geom.N // 10)).exp() * moment)
    positive_endpoint(upper)
    return upper


def run_tail(*, K, occupancies, tilts=DEFAULT_TILTS, precision=192, output=None,
             data=None, map_record=None):
    """Fresh bounds for an explicit, possibly nonconsecutive occupancy list."""
    geom = geometry(K)
    tilts, output = checked_options(tilts, precision, output)
    try:
        requested = tuple(occupancies)
    except TypeError as error:
        raise ValueError('an explicit nonempty occupancy sequence is required') from error
    if (not requested or any(type(q) is not int or not 3 <= q <= geom.group_count for q in requested)
            or len(set(requested)) != len(requested)):
        raise ValueError('distinct tail occupancies in3..L required')
    requested = sorted(requested)
    start = monotonic()
    beta, counts = exact_outer()
    data, map_record, sources = prepare(precision, data, map_record)
    best, choices = dict.fromkeys(requested), dict.fromkeys(requested)
    record = dict(schema=TAIL_SCHEMA, outer='rs16', length_component='tail', K=K, N=geom.N,
        geometry=asdict(geom), threshold=geom.N // 10, distance='1/10',
        occupancy_covered=requested, whole_code_certificate=False, fresh_computation=True,
        physical_t=64, state_bits=16, physical_steps=geom.N // 64,
        macro_t=128, macro_steps=geom.N // 128, physical_steps_per_macro=2,
        macros_per_region=geom.macros_per_region, zero_initial_state=True, final_flush=False,
        state_continuity='retained_between_every_physical_step_and_region',
        terminal='sum of all mass-envelope coordinates; no flush',
        regional_placement='exact without-replacement average of local upper envelopes',
        outer_comparison='beta times uniform full 256-bit group output; artificial zero retains active label',
        beta=str(beta), every_shell_checked=True, count_sha256=count_hash(counts),
        count_premises=premises(), precision=precision, tilts=list(tilts),
        source_sha256=sources, map_record=map_record, trials=[],
        scope='Expected-message first moment over ideal independent uniform GL16 outer '
              'symbol maps, group packet shuffles, regional shuffles, and physical GL16 '
              'updates. Only the explicit occupancy list is covered; no whole-code or seed certificate.')
    for tilt in tilts:
        trial_start = monotonic()
        # Activity1/2 selects a valid local density envelope. Packet activity
        # in the outer comparison remains15/16, with nonzero uniform labels.
        local = q1.kernel_t64.local_operators(data, Q(tilt), activity=Q(1, 2))
        regional = q1.placement(local, epochs=geom.macros_per_region,
            windows=geom.macro_windows, rounding=q1.rounded, maximum_groups=max(requested))
        for q in requested:
            upper = occupancy_upper(regional, K=K, occupancy=q, beta=beta, tilt=tilt)
            if best[q] is None or upper < best[q]:
                best[q], choices[q] = upper, tilt
        total = q1.kernel_t64.up(sum(best.values(), arb(0)))
        if ctx.prec != precision:
            raise ArithmeticError('precision changed during length computation')
        margin = str(-total.log() / arb(2).log())
        record['trials'].append(dict(tilt=tilt, union_upper=positive_endpoint(total),
            margin_bits=margin, elapsed_seconds=monotonic() - trial_start))
        record.update(union_upper=positive_endpoint(total), margin_bits=margin,
            occupancy_uppers={str(q): positive_endpoint(v) for q, v in best.items()},
            occupancy_choices={str(q): value for q, value in choices.items()},
            elapsed_seconds=monotonic() - start)
        if sources != source_snapshot():
            raise RuntimeError('loaded proof source changed during length-tail computation')
        print(f'RS16 K={K} q={requested} tilt={tilt} margin_bits={margin}', flush=True)
    return _finish(record, sources, output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('component', choices=('q1', 'q2', 'tail'))
    parser.add_argument('--K', type=int, required=True)
    parser.add_argument('--occupancies', type=int, nargs='+')
    parser.add_argument('--tilts', nargs='+', default=list(DEFAULT_TILTS))
    parser.add_argument('--precision', type=int, default=192)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    kwargs = dict(K=args.K, tilts=args.tilts, precision=args.precision, output=args.output)
    if args.component == 'tail':
        if args.occupancies is None:
            parser.error('tail requires an explicit --occupancies list')
        run_tail(occupancies=args.occupancies, **kwargs)
    else:
        if args.occupancies is not None:
            parser.error('q1/q2 reject --occupancies')
        run_low(int(args.component[1]), **kwargs)


if __name__ == '__main__':
    main()
