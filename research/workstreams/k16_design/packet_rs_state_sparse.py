"""Fresh q1/q2 bounds for RS16 groups and an explicit uniform-GL_s inner.

Each outer group encodes 128 bits as four GF16 RS[16,8] words. Independent
uniform GL(16,2) symbol maps give exact expected counts by four-bit support.
Independent uniform group shuffles make that support uniform conditional
on its size. Regional permutations place active packets in distinct slots.

The caller supplies freshly prepared physical maps A and C and their record.
For state a and a 64-bit input x, a physical step emits x+A*a, then updates
a to M*a+C*x. Each M is independent uniform GL(s,2), where s comes from the
prepared data. The state starts at zero, persists across all steps and
regions, and is discarded at the end. The record contains the derived maps;
its source entry identifies their base source, not necessarily a file that
directly declares the derived s-dimensional map.

Two physical steps form each 128-bit proof macro. Regional placement averages
the local upper-envelope operators over exact without-replacement positions.
Q1 averages over one uniform support. Q2 averages over two independent
supports and retains their shared regional routing. These ordered products
never reset the inner state. All terminal envelope coordinates count mass.

For each support, Chernoff bounds at distinct positive tilts may be minimized
before folding the exact expected outer counts. Q1 multiplies by L. Q2
multiplies by C(L,2), with distinct support pairs counted in both orders.
The cutoff is floor(N/10). Results cover only the selected occupancy, under
the ideal independent-randomness ensemble. They do not certify the whole
code, a particular seed, or every realized setup.
"""
from __future__ import annotations

from dataclasses import asdict
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
from time import monotonic

from flint import arb, ctx
import packet_q1 as q1
import packet_q2 as q2
from rs_outer import expected_group_support_counts


DEFAULT_GEOMETRY = q1.Geometry(8192, 64, 128)


def _base_sources(map_record):
    """Authenticate declared base files without attributing the derived map to them."""
    entries = [map_record.get('source')]
    if 'base_source' in map_record:
        entries.append(map_record['base_source'])
    if 'constructor_source' in map_record:
        entries.append(map_record['constructor_source'])
    result = {}
    for entry in entries:
        if (not isinstance(entry, dict) or not isinstance(entry.get('path'), str)
                or not isinstance(entry.get('sha256'), str)):
            raise ValueError('map record requires a base source with path and SHA256')
        path = Path(entry['path']).resolve()
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != entry['sha256']:
            raise ValueError(f'map base source changed: {path}')
        result[str(path)] = digest
    return result


def validated(data, map_record):
    """Check map/data identity and explicit uniform GL(s,2) scope; return s.

    This authenticates a fresh preparation; it does not regenerate its large
    exact spectra. The caller must prepare data from the recorded derived maps.
    """
    if not isinstance(data, dict) or not isinstance(map_record, dict):
        raise ValueError('prepared two-step data and complete derived-map record required')
    bits = data.get('bits')
    if type(bits) is not int or not 1 <= bits <= 24:
        raise ValueError('explicit supported positive state dimension required')
    if ('full_state_census' in map_record and map_record['full_state_census'] is not True
            or 'state_count' in map_record and map_record['state_count'] != 1 << bits):
        raise ValueError('a declared full-state census must be complete and dimension-matched')
    if (map_record.get('s') != bits or map_record.get('t') != 64
            or map_record.get('distribution') != 'uniform_gl'
            or map_record.get('map_sha256') != data.get('map_sha256')
            or data.get('physical_step_bits') != 64 or data.get('macro_windows') != 32
            or data.get('macro_step_bits') != 128):
        raise ValueError('record must match the actual t64/GL_s maps and macro geometry')
    physical = q1.kernel_t64.authenticate(data)
    if physical['bits'] != bits or physical['distribution'] != 'uniform_gl':
        raise ValueError('independent uniform GL_s physical updates required')
    rows = map_record.get('expansion_rows_hex')
    columns = map_record.get('feedback_columns')
    if (not isinstance(rows, list) or len(rows) != bits or
            not isinstance(columns, list) or len(columns) != 64):
        raise ValueError('complete derived expansion rows and feedback columns required')
    try:
        rows = [int(value, 16) for value in rows]
    except (TypeError, ValueError) as error:
        raise ValueError('derived expansion rows must be hexadecimal strings') from error
    actual_rows = [int(physical['map_images'][1 << j]) for j in range(bits)]
    if (rows != actual_rows or columns != list(map(int, physical['columns']))
            or any(not 0 <= row < 1 << 64 for row in rows)
            or any(type(column) is not int or not 0 <= column < 1 << bits for column in columns)):
        raise ValueError('derived-map record does not match prepared physical maps')
    _base_sources(map_record)
    return bits


def source_snapshot(map_record):
    """Pin loaded proof sources, this helper, and every declared map base file."""
    sources = q1.source_snapshot()
    path = Path(__file__).resolve()
    sources[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    sources.update(_base_sources(map_record))
    return sources


def shell_hash(counts):
    return hashlib.sha256(json.dumps(list(map(str, counts))).encode()).hexdigest()


def _options(occupancy, geometry, tilts, precision, output):
    if type(occupancy) is not int or occupancy not in (1, 2):
        raise ValueError('occupancy must be one or two')
    if (not isinstance(geometry, q1.Geometry) or geometry.regions != 64
            or geometry.group_dimension != 128 or geometry.group_count < occupancy):
        raise ValueError('RS16 128-to-256 group geometry required')
    if type(precision) is not int or precision < 128:
        raise ValueError('integer precision >=128 required')
    if isinstance(tilts, (str, bytes)):
        raise ValueError('a sequence of positive distinct rational tilts is required')
    try:
        values = tuple(Q(t) for t in tilts)
    except (TypeError, ValueError, ZeroDivisionError, OverflowError) as error:
        raise ValueError('finite rational tilts required') from error
    if not values or min(values) <= 0 or len(set(values)) != len(values):
        raise ValueError('positive distinct rational tilts required')
    output = Path(output) if output is not None else None
    if output is not None and output.exists():
        raise ValueError('output must be fresh')
    return tuple(map(str, values)), output


def _endpoint(value):
    if not value.is_finite() or not value > 0:
        raise ArithmeticError('positive finite outward endpoint required')
    return q1.endpoint(value)


def run_low(occupancy, *, data, map_record, geometry=DEFAULT_GEOMETRY,
            tilts, precision=192, output=None):
    """Freshly evaluate all q1 or q2 supports with the caller's actual inner maps."""
    tilts, output = _options(occupancy, geometry, tilts, precision, output)
    bits = validated(data, map_record)
    ctx.prec = precision
    start = monotonic()
    counts = expected_group_support_counts(n=16, k=8, packet_bits=4, packets_per_symbol=4)
    if len(counts) != 65 or counts[0] or sum(counts) != (1 << 128) - 1:
        raise ArithmeticError('complete exact RS16 nonzero-message shells required')
    sources = source_snapshot(map_record)
    regions, threshold = geometry.regions, geometry.N // 10
    best = ([arb(1) for _ in range(regions + 1)] if occupancy == 1 else
            [[arb(1) for _ in range(v + 1)] for v in range(regions + 1)])
    choices = ([None for _ in best] if occupancy == 1 else [[None for _ in row] for row in best])
    record = dict(schema=f'finite-packet-rs-state-q{occupancy}-1', outer='rs16',
        K=geometry.K, N=geometry.N, geometry=asdict(geometry), groups=geometry.group_count,
        regions=regions, group_dimension=128, group_output_bits=256,
        distance='1/10', threshold=threshold, physical_t=64, state_bits=bits,
        physical_steps=geometry.N // 64, macro_t=128, macro_steps=geometry.N // 128,
        macros_per_region=geometry.macros_per_region, physical_steps_per_macro=2,
        inner_distribution='uniform_gl', inner_group=f'GL({bits},2)',
        physical_updates_independent=True, outer_symbol_group='GL(16,2)',
        zero_initial_state=True, final_flush=False,
        state_continuity='retained_between_every_physical_step_and_region',
        terminal='sum of all mass-envelope coordinates; no flush',
        occupancy_covered=[occupancy], all_supports_covered=True,
        all_two_group_support_pairs_covered=occupancy == 2,
        whole_code_certificate=False, fresh_computation=True,
        count_kind='exact_expected_shells', count_sha256=shell_hash(counts),
        count_premises=dict(outer='four parallel GF16 RS[16,8] words',
            label_mixing='independent uniform GL16 per aligned four-nibble symbol',
            independent_setups_between_groups=True, exact_expected_shell_counts=True,
            independent_uniform_group_packet_shuffles=True,
            independent_uniform_regional_shuffles=True),
        precision=precision, tilts=list(tilts), map_record=map_record,
        map_source_role='declared base source; the map record specifies the actual derived maps',
        source_sha256=sources, local_selection_activity='1/2',
        local_selection_activity_role='valid birth-row proposal, not the packet input law',
        regional_placement='exact without-replacement average of local upper envelopes',
        scope=f'Expected-message first moment over independent ideal uniform outer GL16 '
              f'symbol maps, group packet shuffles, regional shuffles, and inner GL({bits},2) '
              f'updates. The actual t64/s{bits} maps are recorded, not inferred from a '
              f'base file. Only occupancy {occupancy} is covered; no whole-code or seed certificate.',
        trials=[])
    if occupancy == 2:
        record['support_pair_storage'] = 'row v contains u=0..v; off-diagonal count products doubled'
    for tilt in tilts:
        trial_start = monotonic()
        local = q1.kernel_t64.local_operators(data, Q(tilt), activity=Q(1, 2))
        regional = q1.placement(local, epochs=geometry.macros_per_region,
            windows=geometry.macro_windows, rounding=q1.rounded, maximum_groups=occupancy)
        factor = (q1.kernel_t64.aq(Q(tilt)) * threshold).exp()
        if occupancy == 1:
            moments = q1.support_moments(regional[0], regional[1], regions)
            for support, moment in enumerate(moments):
                candidate = q1.kernel_t64.up(factor * moment)
                if candidate < best[support]:
                    best[support], choices[support] = candidate, tilt
            upper = q1.kernel_t64.up(geometry.group_count * sum(
                (q1.kernel_t64.aq(c) * probability for c, probability in zip(counts, best)), arb(0)))
            record.update(support_choices=choices[:],
                support_probability_uppers=[_endpoint(v) for v in best])
        else:
            moments = q2.pair_support_moments(regional, regions)
            for v, row in enumerate(moments):
                for u, moment in enumerate(row):
                    candidate = q1.kernel_t64.up(factor * moment)
                    if candidate < best[v][u]:
                        best[v][u], choices[v][u] = candidate, tilt
            upper = q2.fold_shell_pairs(counts, best, geometry.group_count)
            record.update(support_pair_choices=[row[:] for row in choices],
                support_pair_probability_uppers=[[_endpoint(v) for v in row] for row in best])
        if ctx.prec != precision:
            raise ArithmeticError('precision changed during sparse-state computation')
        endpoint = _endpoint(upper)
        margin = str(-upper.log() / arb(2).log())
        record['trials'].append(dict(tilt=tilt, upper=endpoint, margin_bits=margin,
            elapsed_seconds=monotonic() - trial_start))
        record.update({f'q{occupancy}_upper': endpoint}, margin_bits=margin,
                      elapsed_seconds=monotonic() - start)
        if sources != source_snapshot(map_record):
            raise RuntimeError('loaded mathematical or map-base source changed during computation')
        print(f'RS16 t64/s{bits} q={occupancy} K={geometry.K} '
              f'tilt={tilt} margin_bits={margin}', flush=True)
    validated(data, map_record)
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open('x', encoding='utf-8') as stream:
            stream.write(json.dumps(record, indent=2) + '\n')
    return record
