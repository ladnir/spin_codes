"""Fresh finite-length occupancy-one diagnostic for the current packet code.

This is not a whole-code certificate. The ordered support-placement average
is exact conditional on support. Local state transitions and BCH/GL32 counts
are authenticated upper envelopes, not exact failure probabilities.
No retained numerical bound is reused and no retained proof source is edited.
"""
import argparse
from dataclasses import dataclass, asdict
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys
from time import monotonic

LOCALITY = Path(__file__).resolve().parents[1] / 'permutation_locality'
PACKED = LOCALITY / 'packed_mixing'
CLOSURE = PACKED / 's16_closure'
for directory in (LOCALITY, LOCALITY / 'gf16_packets', PACKED, CLOSURE):
    sys.path.insert(0, str(directory))

from flint import arb, ctx
import kernel_t64
from occupancy_model import placement
from occupancy_memory import rounded
from single_group import support_moments


@dataclass(frozen=True)
class Geometry:
    group_count: int
    regions: int
    group_dimension: int
    macro_windows: int = 32
    packet_bits: int = 4

    def __post_init__(self):
        if any(type(x) is not int or x < 1 for x in asdict(self).values()):
            raise ValueError('positive integer geometry required')
        if (self.macro_windows != 32 or self.packet_bits != 4 or
                self.group_count % self.macro_windows or
                self.group_dimension != 2 * self.regions):
            raise ValueError('half-rate groups of four-bit packets and complete128-bit macro steps required')
        if self.N != 2 * self.K:
            raise ValueError('half-rate geometry required')

    @property
    def K(self):
        return self.group_count * self.group_dimension

    @property
    def N(self):
        return self.group_count * self.regions * self.packet_bits

    @property
    def macros_per_region(self):
        return self.group_count // self.macro_windows


def source_snapshot():
    """Pin loaded local Python sources and the declared fixed map, not packages."""
    paths = {Path(__file__).resolve(), kernel_t64.SELECTED_MAP.resolve()}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename:
            path = Path(filename).resolve()
            if path.suffix == '.py' and (LOCALITY in path.parents or path.parent == Path(__file__).resolve().parent):
                paths.add(path)
    return {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(paths)}


def endpoint(value):
    return [int(v) for v in kernel_t64.up(value).man_exp()]


def fold_cdf(counts, weights):
    """Abel fold against a suffix majorant; CDF differences are not shell caps."""
    tails = list(weights)
    for i in range(len(tails) - 2, -1, -1):
        tails[i] = max(tails[i], tails[i + 1])
    return kernel_t64.up(sum((kernel_t64.aq(counts[i] - counts[i - 1]) * tails[i]
                             for i in range(1, len(counts))), arb(0)))


def evaluate_q1(counts, *, group_count, regions, epochs_per_region,
                group_dimension, count_kind, tilts, precision=192,
                data=None, map_record=None, metadata=None, output=None,
                canonical_caps=None, canonical_law=None):
    """Evaluate q=1 with explicit finite geometry and nonzero-message counts.

    count_kind is 'shells' for exact expected shell counts, or 'cdf' for
    cumulative upper caps. Both must exclude the zero message. The input law
    must supply independent uniform nonzero four-bit labels conditional on
    support, with an independent uniform column shuffle and regional shuffle.
    """
    if precision < 128 or any(type(x) is not int or x < 1 for x in
            (group_count, regions, epochs_per_region, group_dimension)):
        raise ValueError('positive finite geometry and precision >=128 required')
    geometry = Geometry(group_count, regions, group_dimension)
    if epochs_per_region != geometry.macros_per_region:
        raise ValueError('complete macro steps and half-rate group dimension required')
    if not tilts or any(Q(t) <= 0 for t in tilts) or len(set(map(Q, tilts))) != len(tilts):
        raise ValueError('distinct positive rational tilts required')
    if output is not None and output.exists():
        raise ValueError('output must be fresh')
    start = monotonic()
    counts = tuple(map(Q, counts))
    if len(counts) != regions + 1 or counts[0] != 0 or any(x < 0 for x in counts):
        raise ValueError('one nonzero-message support count per support0..regions required')
    shell_counts = counts if count_kind == 'shells' else None
    if count_kind == 'shells':
        total, cumulative = Q(0), []
        for value in counts:
            total += value
            cumulative.append(total)
        counts = tuple(cumulative)
    elif count_kind != 'cdf':
        raise ValueError('count_kind must explicitly be shells or cdf')
    if any(a > b for a, b in zip(counts, counts[1:])) or counts[-1] != (1 << group_dimension) - 1:
        raise ValueError('complete monotone count CDF with exact nonzero-message total required')
    k, n, groups, windows = geometry.K, geometry.N, group_count, geometry.macro_windows
    assert n == 2 * k and regions * epochs_per_region * 128 == n
    threshold = n // 10
    print(f'DIAGNOSTIC K={k} N={n} L={groups} cutoff={threshold} '
          f'regions={regions} macros_per_region={epochs_per_region}', flush=True)
    if data is None:
        data, map_record = kernel_t64.prepare(birth_density='capped')
    sources = source_snapshot()
    ctx.prec = precision
    best, choices = [arb(1)] * (regions + 1), [None] * (regions + 1)
    if (canonical_caps is None) != (canonical_law is None):
        raise ValueError('canonical caps and block law must be supplied together')
    if canonical_caps is not None:
        from sparse_hill_q1 import canonical_terms, support_laws
        laws = support_laws(canonical_law, len(canonical_caps) - 1)
        sources = source_snapshot()
    record = dict(schema='finite-packet-q1-diagnostic-1', K=k, N=n,
        geometry=asdict(geometry), group_output_bits=regions * 4,
        distance='1/10', threshold=threshold, groups=groups, regions=regions,
        physical_t=64, state_bits=16, physical_steps=n // 64,
        macro_t=128, macro_steps=n // 128, macros_per_region=epochs_per_region,
        physical_steps_per_macro=2, zero_initial_state=True, final_flush=False,
        occupancy_covered=[1], occupancies_not_covered=[2, groups],
        whole_code_certificate=False, precision=precision, tilts=tilts,
        count_kind=count_kind, group_dimension=group_dimension,
        count_sha256=hashlib.sha256(json.dumps(list(map(str, counts))).encode()).hexdigest(),
        count_premises=metadata, support_min=next(i for i, c in enumerate(counts) if c),
        map_record=map_record, source_sha256=sources, trials=[],
        scope='Explicit supplied outer expected support counts/caps; '
              'independent uniform shared-column/regional shuffles and '
              'uniform nonzero packet labels, selected t64/s16 '
              'maps and independent ideal uniform GL16 updates. '
              'Conditional placement is exact; count and transition envelopes '
              'are upper bounds. Not a seed-specific or whole-code certificate.')
    for tilt in tilts:
        ctx.prec = precision
        local = kernel_t64.local_operators(data, Q(tilt), activity=Q(1, 2))
        exact = placement(local, epochs=epochs_per_region, windows=windows,
                          rounding=rounded, maximum_groups=1)
        moments = support_moments(exact[0], exact[1], regions)
        factor = (kernel_t64.aq(Q(tilt)) * threshold).exp()
        for support, moment in enumerate(moments):
            candidate = kernel_t64.up(factor * moment)
            if candidate < best[support]:
                best[support], choices[support] = candidate, tilt
        cdf = kernel_t64.up(groups * fold_cdf(counts, best))
        upper, terms, largest = cdf, [], []
        shell_upper = None
        if shell_counts is not None:
            shell_upper = kernel_t64.up(groups * sum((kernel_t64.aq(c) * w
                for c, w in zip(shell_counts, best)), arb(0)))
            upper = min(upper, shell_upper)
        canonical = None
        if canonical_caps is not None:
            canonical, terms = canonical_terms(canonical_caps, laws, best, groups)
            upper = min(cdf, canonical)
            largest = sorted(range(1, len(terms) + 1), key=lambda h: float(terms[h - 1]), reverse=True)[:5]
        if ctx.prec != precision or not upper.is_finite() or not upper > 0:
            raise ArithmeticError('positive finite diagnostic endpoint required')
        margin = str(-upper.log() / arb(2).log())
        record['trials'].append(dict(tilt=tilt, margin_bits=margin,
            cdf_upper=endpoint(cdf), canonical_upper=endpoint(canonical) if canonical is not None else None,
            exact_shell_upper=endpoint(shell_upper) if shell_upper is not None else None,
            q1_upper=endpoint(upper), dominant_canonical_H=largest))
        record.update(q1_upper=endpoint(upper), margin_bits=margin,
            q1_below_2_minus_40=bool(upper < arb(2) ** -40),
            support_choices=choices[:],
            support_probability_uppers=[endpoint(v) for v in best],
            canonical_H_terms=[endpoint(v) for v in terms],
            elapsed_seconds=monotonic() - start)
        if sources != source_snapshot():
            raise RuntimeError('loaded local mathematical source changed during diagnostic')
        if output is not None:
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(json.dumps(record, indent=2) + '\n')
        print(f'tilt={tilt} q1_margin_bits={margin} dominant_H={largest}', flush=True)
    print('Only q=1 is covered. q=2..L and the whole-code union remain unchecked.', flush=True)
    return record


def run(exponent, tilts, precision, output=None, outer='bch'):
    if exponent not in (16, 20):
        raise ValueError('explicit K16/K20 diagnostic required')
    if outer in ('random-systematic', 'rs'):
        if outer == 'random-systematic':
            import random_systematic
            counts = random_systematic.expected_group_support_counts(128, 4)
            regions, dimension = 256, 512
            premises = dict(outer='four independent random systematic binary[256,128] codes',
                setup='independent uniform binary128x128 parity matrix per constituent',
                label_mixing='independent uniform GL4 per column', exact_expected_shell_counts=True)
        else:
            import rs_outer
            counts = rs_outer.expected_group_support_counts()
            regions, dimension = 64, 128
            premises = dict(outer='four parallel GF256 RS[8,4] words',
                label_mixing='independent uniform GL32 per aligned four-byte symbol',
                exact_expected_shell_counts=True)
        groups = (1 << exponent) // dimension
        return evaluate_q1(counts, group_count=groups, regions=regions,
            epochs_per_region=groups // 32, group_dimension=dimension, count_kind='shells',
            tilts=tilts, precision=precision, metadata=premises, output=output)
    if outer != 'bch':
        raise ValueError('unknown outer diagnostic')
    from outer_hill_intersection import authenticated_bch_cdf
    from local_models import full_block
    print('Authenticating local BCH/GL32 count premises', flush=True)
    counts, premises = authenticated_bch_cdf()
    groups = (1 << exponent) // 512
    return evaluate_q1(counts, group_count=groups, regions=256,
        epochs_per_region=groups // 32, group_dimension=512, count_kind='cdf',
        tilts=tilts, precision=precision, metadata=premises, output=output,
        canonical_caps=premises['canonical_cdf'], canonical_law=full_block(8))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exponent', type=int, choices=(16, 20), default=16)
    parser.add_argument('--outer', choices=('bch', 'random-systematic', 'rs'), default='bch')
    parser.add_argument('--tilts', nargs='+', default=['.00256', '.00512', '.01024', '.0256', '.1024', '.256'])
    parser.add_argument('--precision', type=int, default=192)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    run(args.exponent, args.tilts, args.precision, args.output, args.outer)


if __name__ == '__main__':
    main()
