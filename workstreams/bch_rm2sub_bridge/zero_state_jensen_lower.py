"""Sharp actual first-moment lower via parity-conditioned single marginals.

All occupied rows have weight80. After row permutation a row is a uniform
weight80 slice. Conditioning one column leaves uniform weight79/80 slices
on255 columns, which gives an explicit conditional parity mixing bound.
Jensen is applied only after conditioning on the actual all-even good event.
"""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path
from flint import arb, arb_poly, ctx
import bridge as base
import tightened_occupancy as tight
from certify_tail_parity_mixing import kraw_values


def run(verify=False):
    output = base.HERE/'generated/zero_state_jensen_lower.json'
    old = base.read(output) if verify else None
    if old:
        for name, digest in old['local_sha256'].items():
            assert base.sha(base.HERE/name) == digest
    else:
        assert not output.exists()
    ctx.prec = 512 if verify else 256
    n, q, lo, hi = 256, 2620, 64, 1536
    base.load_map(tight.NAME)
    family_path = base.HERE/'generated/weight80_family.json'
    family = base.read(family_path)
    for name, digest in family['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    gaps = [max(F(abs(v), math.comb(255, j))
                for j, v in enumerate(kraw_values(255, w)) if 0 < j < 255)
            for w in [79, 80]]
    assert gaps == [F(97, 255), F(19, 51)]
    epsilon = ((1 << 255)-2)*max(gaps)**q
    assert epsilon < F(1, 1 << 3300)
    inflation = 1+arb(epsilon.numerator)/epsilon.denominator
    spectrum = tight.kernel_spectrum()
    kernel = arb_poly([spectrum.get(j, 0) for j in range(129)])**64
    mass = (arb(11)/16)**q
    outside = arb(0)
    cost = arb(0)
    for j in range(q+1):
        if j % 2 == 0:
            if lo <= j <= hi:
                beta = kernel[j]/math.comb(8192, j)
                assert beta > 0 and beta <= 1
                cost -= mass*beta.log()
            else:
                outside += mass
        if j < q:
            mass *= arb(q-j)*5/((j+1)*11)
    good_relative = (1-2*n*inflation*outside).lower()
    assert good_relative > arb(999)/1000
    # P(all-even good) >=2^-255*good_relative. The numerator of the
    # conditional expected negative log-kernel product is <=
    # 256*2^-254*(1+epsilon)*cost.
    expected_cost = (2*n*inflation*cost/good_relative).upper()
    probability = (arb(2)**-255*good_relative*(-expected_cost).exp()).lower()
    from audit_bch_q1_full_arb import rational
    lower = rational(probability)
    mean = math.comb(8192, q)*family['shell_size_lower']**q*lower
    previous = base.decode(family['first_moment_lower'])
    assert mean > previous
    if old:
        assert lower >= base.decode(old['zero_state_probability_lower'])
        assert mean >= base.decode(old['first_moment_lower'])
        print('512-bit Jensen first-moment replay passed', flush=True)
        return
    result = dict(status='ACTUAL_PARITY_CONDITIONED_JENSEN_FIRST_MOMENT_LOWER',
                  configuration=tight.NAME, occupation=q, row_weight=80,
                  eligible_even_region_weights=[lo, hi],
                  conditional_parity_error_upper=base.encode(epsilon),
                  conditional_slice_character_gaps=[base.encode(g) for g in gaps],
                  good_parity_relative_probability_lower=base.encode(rational(good_relative)),
                  conditional_expected_negative_log_kernel_upper=base.encode(rational(expected_cost)),
                  zero_state_probability_lower=base.encode(lower), first_moment_lower=base.encode(mean),
                  setup_failure_probability_lower_bound=False,
                  local_sha256={str(p.relative_to(base.HERE)):base.sha(p)
                                for p in [Path(__file__), Path(tight.__file__), family_path,
                                          base.HERE/'certify_tail_parity_mixing.py',
                                          base.HERE/'inputs/t128_s15_b_kernel_spectrum.json']})
    base.write_new(output, result)
    log2 = lambda value: math.log2(value.numerator)-math.log2(value.denominator)
    print('Actual Jensen probability lower log2', log2(lower), 'mean lower log2', log2(mean),
          'gain over previous weight80 mean', log2(mean/previous), 'bits', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    run(parser.parse_args().verify)
