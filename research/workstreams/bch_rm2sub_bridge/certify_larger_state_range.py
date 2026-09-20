"""Outward larger-state occupancy ranges, retaining exact all-one bands."""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from flint import arb, ctx
import bridge as base
import larger_state_maps as maps
import activation_bridge as q1
import polynomial_regions as poly
import scaled_adaptive as scaled
import positive_line_hull as hull
from general_batch_certificate import density_cost
from certify_full_exactcaps_range import checked_caps
from audit_bch_q1_full_arb import rational


def costs_for(bands, ps, caps):
    assert len(bands) == len(ps)
    result = []
    for band, p in zip(bands, ps):
        if p == 1:
            assert tuple(band) == (256,) and caps[256] == 1
            result.append(F(1))
        else:
            assert 0 < p < 1
            result.append(density_cost(band, p, caps))
    return result


def run(screen_name, anchor, lower, upper, tag, verify=False):
    assert Path(screen_name).name == screen_name and tag.isidentifier()
    assert 1 <= lower <= upper <= 8192
    screen_path = base.HERE/'generated'/screen_name
    output = base.HERE/'generated'/f'larger_range_{tag}_outward.json'
    old = base.read(output) if verify else None
    if not verify:
        assert not output.exists()
    screen = base.read(screen_path)
    caps = checked_caps(screen)
    bands = screen['bands']
    assert sorted(w for band in bands for w in band) == list(base.WEIGHTS)
    name = screen['configuration']
    t, s, spectrum, kernel = maps.load(name)
    witnesses = [row for row in screen['rows'] if row['occupation'] == anchor]
    assert len(witnesses) == 1
    witness = witnesses[0]
    if old:
        assert old['occupancy_range'] == [lower, upper] and old['configuration'] == name
        for source, digest in old['local_sha256'].items():
            assert base.sha(base.HERE/source) == digest
    ctx.prec = 512 if verify else 256
    tilt = witness['witness_tenth']
    lam = (arb(tilt)/10).exp()
    correction = (base.CUTOFF*lam).exp()
    region = poly.regions(t, s, spectrum, kernel, (-lam).exp(), upper)
    ps = [base.decode(p) for p in witness['p']]
    costs = costs_for(bands, ps, caps)
    roots = [((arb(c.numerator)/c.denominator).log()/256).exp().upper() for c in costs]
    probabilities = [arb(p.numerator)/p.denominator for p in ps]
    left = np.nextafter(np.array([float((r*(1-p)).upper()) for r,p in zip(roots,probabilities)]), np.inf)
    right = np.nextafter(np.array([float((r*p).upper()) for r,p in zip(roots,probabilities)]), np.inf)
    keep = hull.indices(left, right)
    mantissas, exponents = scaled.initial(region)
    rows = []
    total = F(0)
    for q, matrix, exponent in scaled.matrices(mantissas, exponents, left[keep], right[keep]):
        if q < lower:
            continue
        value = tuple(arb(float(v)) for v in matrix.flat)
        for _ in range(8):
            value = q1.positive_mul(value, value)
        bound = math.comb(8192,q)*len(bands)**q*rational((sum(value[:3],arb(0))*correction*arb(2)**(256*exponent)).upper())
        assert bound > 0
        # Keep at most 80 certified bits per occupancy. The target needs
        # only 40; this avoids storing enormous dyadic denominators.
        numerator, denominator = bound.numerator, bound.denominator
        power = numerator.bit_length()-denominator.bit_length()
        candidate = F(1 << power) if power >= 0 else F(1, 1 << (-power))
        if candidate < bound:
            power += 1
        power = max(-80, power)
        compact = F(1 << power) if power >= 0 else F(1, 1 << (-power))
        assert compact >= bound
        total += compact
        if old:
            assert bound <= base.decode(old['rows'][q-lower]['upper'])
        rows.append(dict(occupation=q, upper=base.encode(compact),
                         margin_bits_diagnostic=math.log2(denominator)-math.log2(numerator)))
    if old:
        assert len(rows) == len(old['rows'])
        # Compare the saved compact bounds, rather than requiring identical
        # 256/512-bit rounding decisions at a dyadic boundary.
        assert sum((base.decode(row['upper']) for row in old['rows']), F(0)) == base.decode(old['range_upper'])
        base.write_new(output.with_name(f'larger_range_{tag}_replay.json'), dict(
            status='LARGER_STATE_512_BIT_RANGE_REPLAY_PASSED', rows_checked=len(rows),
            producer_sha256=base.sha(output), verifier_sha256=base.sha(Path(__file__))))
        print('512-bit range replay passed', lower, upper, flush=True)
        return
    dependencies = [Path(__file__), Path(poly.__file__), Path(scaled.__file__), Path(hull.__file__), Path(q1.__file__),
        base.HERE/'general_occupancy.py', base.HERE/'tightened_occupancy.py', base.HERE/'general_batch_certificate.py',
        base.HERE/'certify_full_exactcaps_range.py', screen_path]
    base.write_new(output, dict(status='LARGER_STATE_OUTWARD_RANGE', configuration=name,
        occupancy_range=[lower,upper], anchor=anchor, rows=rows, range_upper=base.encode(total),
        range_below_2_to_minus_60=total<F(1,1<<60), all_occupations_certified=False,
        parameters=dict(message_bits=1<<20, output_bits=1<<21, step_bits=t, state_bits=s, distance_cutoff=base.CUTOFF),
        local_sha256={**screen['local_sha256'], **{str(p.relative_to(base.HERE)):base.sha(p) for p in dependencies}}))
    print('Range', lower, upper, 'margin', math.log2(total.denominator)-math.log2(total.numerator),
          'passing rows', sum(base.decode(r['upper']) < F(1,1<<60) for r in rows), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--screen', required=True)
    parser.add_argument('--anchor', type=int, required=True)
    parser.add_argument('--lower', type=int, required=True)
    parser.add_argument('--upper', type=int, required=True)
    parser.add_argument('--tag', required=True)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    run(args.screen, args.anchor, args.lower, args.upper, args.tag, args.verify)
