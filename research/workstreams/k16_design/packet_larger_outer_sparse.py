"""Fresh uniform-majorant q1/q2 bounds for the larger GF256 RS outers.

The inner is the retained t64/s22 construction. Each group consists of four
parallel GF256 RS[n,n/2] words with independent uniform GL32 symbol maps.
An independent packet shuffle in each group and independent regional
permutations follow. Setup matrices are sampled independently for each
physical inner step. State starts at zero, persists, and is not flushed.

The n=16 and n=32 outers have dimensions 256 and 512 respectively. With
2048 groups they give K=2^19 and K=2^20. We dominate active-group output
measures by beta times fully uniform bits, including artificial zero words.
This comparison works at q1/q2 too, although it can be looser than exact
support counts. No old sparse numerical bound is reused.
"""
from fractions import Fraction as Q
from math import comb
from flint import arb, ctx
import packet_inner_quadratic_extension as maps
import packet_q1 as q1
from packet_regional_power import placement_power
from packet_uniform_tail import regional_uniform
from rs_uniform_envelope import UniformInputEnvelope


def outer(multiplier):
    if type(multiplier) is not int or multiplier not in (2, 4):
        raise ValueError('outer multiplier must be two or four')
    return UniformInputEnvelope(8*multiplier, 4*multiplier, 4, 8)


def occupancy_bound(regional, *, q, groups, regions, beta, tilt, cutoff):
    """Positive Arb composition; all state coordinates retain terminal mass."""
    if type(q) is not int or not 1 <= q <= groups:
        raise ValueError('positive feasible occupancy required')
    aq, up = q1.kernel_t64.aq, q1.kernel_t64.up
    matrix = regional_uniform(regional, q)**regions
    moment = sum((matrix[0, j] for j in range(matrix.ncols())), arb(0))
    result = up(comb(groups, q)*aq(beta)**q*(aq(Q(tilt))*cutoff).exp()*moment)
    if not result.is_finite() or not result > 0:
        raise ArithmeticError('positive finite outward endpoint required')
    return result


def run(data, map_record, *, multipliers=(2, 4), tilts, precision=256, progress=None):
    """Return a separate partial receipt for each explicitly selected outer."""
    multipliers, tilts = tuple(multipliers), tuple(map(Q, tilts))
    if not multipliers or len(set(multipliers)) != len(multipliers):
        raise ValueError('nonempty distinct outer multipliers required')
    envelopes = {m: outer(m) for m in multipliers}
    if not tilts or len(set(tilts)) != len(tilts) or min(tilts) <= 0:
        raise ValueError('distinct positive rational tilts required')
    if type(precision) is not int or precision < 192:
        raise ValueError('at least192 bits of interval precision required')
    maps.authenticate(data, map_record)
    if data['bits'] != 22 or data['birth_density'] != 'capped':
        raise ValueError('retained t64/s22 capped map preparation required')
    sources = q1.source_snapshot()
    best = {m: {1: None, 2: None} for m in multipliers}
    choices = {m: {} for m in multipliers}
    previous = ctx.prec
    try:
        ctx.prec = precision
        for tilt in tilts:
            local = q1.kernel_t64.local_operators(data, tilt, activity=Q(1, 2))
            regional = placement_power(local, epochs=64, windows=32, maximum_groups=2)
            for m, envelope in envelopes.items():
                N = 2048*envelope.output_bits
                for q in (1, 2):
                    value = occupancy_bound(regional, q=q, groups=2048,
                        regions=envelope.regions, beta=envelope.beta, tilt=tilt, cutoff=N//10)
                    if best[m][q] is None or value < best[m][q]:
                        best[m][q], choices[m][q] = value, str(tilt)
            if progress is not None:
                progress(str(tilt), {m: {q: str(-v.log()/arb(2).log())
                    for q, v in values.items()} for m, values in best.items()})
            maps.authenticate(data, map_record)
            if sources != q1.source_snapshot():
                raise RuntimeError('proof sources changed during sparse comparison')
        receipts = {}
        for m, envelope in envelopes.items():
            total = q1.kernel_t64.up(sum(best[m].values(), arb(0)))
            receipts[m] = dict(schema='larger-rs-outer-uniform-sparse-1',
                multiplier=m, K=2048*envelope.message_bits, N=2048*envelope.output_bits,
                groups=2048, regions=envelope.regions, outer_n=envelope.n, outer_k=envelope.k,
                parallel_rows=4, base_field_size=256, symbol_mixer='independent uniform GL32',
                outer_dimension=envelope.message_bits, outer_length=envelope.output_bits,
                beta=str(envelope.beta), distance='1/10', cutoff=2048*envelope.output_bits//10,
                occupancy_covered=[1, 2], whole_code_certificate=False,
                fresh_computation=True, precision=precision, physical_t=64, state_bits=22,
                zero_initial_state=True, final_flush=False, continuous_state=True,
                map_record=map_record, source_sha256=sources, tilts=list(map(str, tilts)),
                occupancy_choices={str(q): t for q, t in choices[m].items()},
                occupancy_uppers={str(q): q1.endpoint(v) for q, v in best[m].items()},
                occupancy_margin_bits={str(q): str(-v.log()/arb(2).log()) for q,v in best[m].items()},
                union_upper=q1.endpoint(total), margin_bits=str(-total.log()/arb(2).log()),
                scope='Only q1/q2 first-moment contributions for the stated ideal ensemble. '
                      'Not a whole-code certificate, seed guarantee, or performance measurement.')
        return receipts
    finally:
        ctx.prec = previous
