"""Larger-state Q1 producer plus independent 512-bit coefficient replay."""
import argparse
import math
from pathlib import Path
from fractions import Fraction as F
from flint import arb, ctx
import bridge as base
import larger_state_maps as maps
import activation_bridge as q1
import verify_bridge as replay
from audit_bch_q1_full_arb import rational


def run(name, verify=False):
    screen_path = base.HERE/'generated'/f'larger_{name}_q1_screen.json'
    output = base.HERE/'generated'/f'larger_{name}_q1_outward.json'
    old = base.read(output) if verify else None
    if not verify:
        assert not output.exists()
    t, s, spectrum, _ = maps.load(name)
    screen = base.read(screen_path)
    assert screen['configuration'] == name
    for source, digest in screen['local_sha256'].items():
        assert base.sha(base.HERE/source) == digest
    if old:
        for source, digest in old['local_sha256'].items():
            assert base.sha(base.HERE/source) == digest
        for source, digest in old['outer_sha256'].items():
            assert base.sha(base.BCH/source) == digest
        assert old['configuration'] == name
    groups = {}
    for w in base.WEIGHTS:
        tilt = screen['coefficient_rows'][str(w)]['witness_tenth']
        assert isinstance(tilt, int) and -120 <= tilt <= 0
        groups.setdefault(tilt, []).append(w)
    ctx.prec = 512 if verify else 256
    coefficients = {}
    for tilt, weights in sorted(groups.items()):
        if verify:
            values = replay.coefficients(t, s, spectrum, tilt)
        else:
            lam = (arb(tilt)/10).exp()
            values = q1.positive_coefficients(*q1.positive_regions(t, s, spectrum, (-lam).exp(), arb, base.ROWS//t), arb, 256)
            correction = base.ROWS*(base.CUTOFF*lam).exp()
            values = [v*correction for v in values]
        for w in weights:
            coefficients[w] = min(F(base.ROWS), rational(values[w].upper()))
            assert coefficients[w] > 0
            if old:
                assert coefficients[w] <= base.decode(old['coefficient_upper'][str(w)])
        print(name, 'replayed' if verify else 'certified', 'tilt', tilt, flush=True)
    upper, factor, rest = base.bch_bound(coefficients)
    assert upper < F(1, 1 << 40)
    if old:
        assert upper <= base.decode(old['Q1_upper'])
        stored = {int(w):base.decode(v) for w,v in old['coefficient_upper'].items()}
        assert base.bch_bound(stored) == tuple(base.decode(old[k]) for k in ('Q1_upper','dual_domination_factor','remaining_shell_upper'))
        base.write_new(output.with_name(f'larger_{name}_q1_replay.json'), dict(
            status='INDEPENDENT_512_BIT_Q1_REPLAY_PASSED', configuration=name, coefficients_checked=len(coefficients),
            producer_sha256=base.sha(output), verifier_sha256=base.sha(Path(__file__)),
            independent_coefficient_source_sha256=base.sha(Path(replay.__file__))))
        print(name, 'all 92 independent coefficients passed', flush=True)
        return
    margin = math.log2(upper.denominator)-math.log2(upper.numerator)
    dependencies = [Path(__file__), Path(q1.__file__), Path(replay.__file__), screen_path]
    base.write_new(output, dict(status='LARGER_STATE_OUTWARD_Q1_BOUND', configuration=name,
        parameters=dict(message_bits=1<<20, output_bits=1<<21, step_bits=t, state_bits=s, distance_cutoff=base.CUTOFF),
        coefficient_upper={str(w):base.encode(v) for w,v in coefficients.items()},
        Q1_upper=base.encode(upper), dual_domination_factor=base.encode(factor), remaining_shell_upper=base.encode(rest),
        margin_bits_diagnostic=margin, all_occupations_certified=False,
        local_sha256={**screen['local_sha256'], **{str(p.relative_to(base.HERE)):base.sha(p) for p in dependencies}},
        outer_sha256={f'generated/shift_rank_oa29_joint/{f}.json':base.sha(base.BCH/f'generated/shift_rank_oa29_joint/{f}.json') for f in ('audit','objective')}))
    print(name, 'outward Q1 margin', margin, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--configuration', choices=maps.NAMES, required=True)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    run(args.configuration, args.verify)
