"""Fresh conditioned-iid components for explicitly chosen larger MDS outers.

This uses the retained physical kernel and marker-conditioning argument at
the caller's actual t64 state size. It never calls an entry point that
silently fixes s=16 or the 128-to-256 outer. The explicit envelope determines
the number of groups, routing regions, beta, and output length.

Within each region, mark its L slots independently with probability p.
Conditioning on exactly q markers has probability C(L,q)p^q(1-p)^(L-q).
Assign independent uniform four-bit values to the marked slots, including
zero. Conditional on this event, this is the uniform-input majorant for q
active outer groups. The independent conditioning events across R regions
cost the inverse R-th power of that probability. State continues between
regions. All local row choices are complete valid positive envelopes.

Only listed occupancy contributions are covered. A whole certificate must
independently validate the outer construction, cover every other occupancy,
and add the outward endpoints with the same map and setup distributions.
"""
from fractions import Fraction as Q
from math import comb
from time import monotonic

from flint import arb, ctx
import packet_q1 as q1
import packet_rs_state_sparse as sparse
import packet_rs_k20_dense as dense
from packet_outer_cost_search import geometry, _occupancies


def run(data, map_record, *, K, envelope, occupancies, tilts,
        marker_probabilities, precision=192, progress=None):
    meta = geometry(K, envelope)
    qs = _occupancies(occupancies, meta['groups'])
    tilts, markers = tuple(map(Q, tilts)), tuple(map(Q, marker_probabilities))
    if (not tilts or len(set(tilts)) != len(tilts) or min(tilts) <= 0 or
            not markers or len(set(markers)) != len(markers) or min(markers) <= 0 or
            max(markers) > 1 or qs[0] < meta['groups'] and not any(p < 1 for p in markers)):
        raise ValueError('distinct positive tilts and feasible marker probabilities required')
    if type(precision) is not int or precision < 192:
        raise ValueError('at least 192-bit interval precision required')
    bits = sparse.validated(data, map_record)
    sources = sparse.source_snapshot(map_record)
    finite = q1.Geometry(meta['groups'], meta['regions'], meta['outer_dimension'])
    best, choices = dict.fromkeys(qs), dict.fromkeys(qs)
    old = ctx.prec
    started = monotonic()
    try:
        ctx.prec = precision
        log_binomial = {q: arb(comb(meta['groups'], q)).log() for q in qs}
        for tilt in tilts:
            candidates = dense.physical_candidates(data, tilt)
            for p in markers:
                matrix, rows = dense.select_physical(candidates, Q(15, 16)*p)
                power = matrix**(meta['N']//64)
                moment = sum((power[0, j] for j in range(power.ncols())), arb(0))
                if not moment.is_finite() or not moment > 0:
                    raise ArithmeticError('finite positive iid moment required')
                logged = moment.log()
                for q in qs:
                    if p == 1 and q != meta['groups']:
                        continue
                    bound = dense.coefficient_log_bound(finite_geometry=finite,
                        q=q, tilt=tilt, beta=envelope.beta, marker_probability=p,
                        log_moment=logged, log_binomial=log_binomial[q])
                    upper = q1.kernel_t64.up(bound.exp())
                    if not upper.is_finite() or not upper > 0:
                        raise ArithmeticError('finite positive outward occupancy endpoint required')
                    if best[q] is None or upper < best[q]:
                        best[q] = upper
                        choices[q] = dict(tilt=str(tilt), marker_probability=str(p),
                            physical_birth_rows=list(rows))
            sparse.validated(data, map_record)
            if sources != sparse.source_snapshot(map_record):
                raise RuntimeError('loaded mathematical sources changed during dense component')
            if progress is not None:
                progress(tilt, best, choices, monotonic()-started)
        if any(value is None for value in best.values()):
            raise ArithmeticError('some listed occupancy received no valid bound')
        return dict(meta, schema='larger-outer-conditioned-iid-component-1',
            physical_t=64, state_bits=bits, precision=precision, fresh_computation=True,
            map_record=map_record, source_sha256=sources, occupancy_covered=list(qs),
            occupancy_uppers={str(q): q1.endpoint(value) for q, value in best.items()},
            occupancy_choices={str(q): value for q, value in choices.items()},
            occupancy_margin_bits={str(q): str(-v.log()/arb(2).log()) for q, v in best.items()},
            tilts=list(map(str, tilts)), marker_probabilities=list(map(str, markers)),
            whole_code_certificate=False, elapsed_seconds=monotonic()-started,
            scope='Only listed occupancies of the explicit ideal ensemble are bounded. '
                  'Full coverage and an exact union sum are required for a whole certificate.')
    finally:
        ctx.prec = old
