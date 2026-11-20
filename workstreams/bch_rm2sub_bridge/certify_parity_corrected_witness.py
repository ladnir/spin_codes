"""Actual parity-corrected message family; the encoder is unchanged.

Three reserved BCH rows solve column parity. The rank-failure bound is
exact, not inferred from the sample right inverse. All other rows retain
the original independent setup permutations.
"""
import argparse
import math
import random
from fractions import Fraction as F
from pathlib import Path
from flint import arb, ctx
import bridge as base
import christoffel_caps as christoffel
import parity_corrector as correction


def run(verify=False):
    path = base.HERE/'generated/parity_corrected_witness.json'
    old = base.read(path) if verify else None
    if old:
        for name, digest in old['local_sha256'].items():
            assert base.sha(base.HERE/name) == digest
        for name, digest in old['outer_sha256'].items():
            assert base.sha(base.BCH/name) == digest
    else:
        assert not path.exists()
    christoffel.build(verify=True)
    # C has d>=38, so C^perp has OA strength37. C^perp is even, contains1,
    # and has d>=30 from the retained rank certificate checked above.
    caps = {w:min(math.comb(256, w), math.floor(christoffel.shell_bound(256, 128, 18, w)))
            for w in range(30, 129, 2)}
    failure = sum((F(cap**3, math.comb(256, w)**2)*(1 if w == 128 else 2)
                   for w, cap in caps.items()), F(0))/2
    assert failure < F(1, 1 << 100)
    basis = correction.fixed_c_basis()
    rng = random.Random(83742109)
    permutations = [rng.sample(range(256), 256) for _ in range(3)]
    rank, section = correction.right_inverse(basis, permutations)
    assert rank == 255 and len(section) == 255
    ctx.prec = 512 if verify else 256
    q, low, high = 2610, 740, 897
    mass, tail = (arb(11)/16)**q, arb(0)
    for j in range(q+1):
        if not low <= j <= high:
            tail += mass
        if j < q:
            mass *= arb(q-j)*5/((j+1)*11)
    retained = (1-256*tail).lower()
    normalized_mean = ((1-arb(failure.numerator)/failure.denominator)*retained).lower()
    assert normalized_mean > arb(3)/4
    assert 80*q+3*256 == 209568 < 209716
    from audit_bch_q1_full_arb import rational
    if old:
        assert failure == base.decode(old['corrector_rank_failure_upper'])
        assert rational(normalized_mean) >= base.decode(old['weighted_mean_over_core_family_size_lower'])
        assert section == old['sample_right_inverse_coefficients']
        print('Exact rank and512-bit parity-corrected witness replay passed', flush=True)
        return
    local = [Path(__file__), Path(correction.__file__), Path(christoffel.__file__),
             base.HERE/'generated/christoffel_oa29_caps.json', base.HERE/'generated/weight80_family.json']
    outer = [base.BCH/'code/bch_quotient.py', base.BCH/'code/affine_wambach.py',
             base.BCH/'generated/bch256_shift_rank_q30_refined.json']
    base.write_new(path, dict(status='ACTUAL_PARITY_CORRECTED_INVERSE_WEIGHT_WITNESS', configuration='t128_s15',
        core_occupied_rows=q, core_row_weight=80, available_core_positions=8189,
        reserved_correction_positions=[8189, 8190, 8191], total_occupation_range=[2610, 2613],
        output_weight_upper=209568, core_region_window=[low, high], actual_even_region_window=[740, 900],
        dual_shell_caps=[dict(weight=w, cap=cap) for w, cap in caps.items()],
        corrector_rank_failure_upper=base.encode(failure), core_window_probability_lower=base.encode(rational(retained)),
        weighted_mean_over_core_family_size_lower=base.encode(rational(normalized_mean)),
        core_family_size_formula='binom(8189,2610)*|T80|^2610',
        sample_seed=83742109, sample_corrector_rank=rank, sample_right_inverse_coefficients=section,
        sample_used_as_probability_evidence=False, same_encoder_and_setup_distribution=True,
        full_second_moment_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in local},
        outer_sha256={str(p.relative_to(base.BCH)):base.sha(p) for p in outer}))
    print('Corrector rank failure <2^-100; normalized weighted mean >=', float(normalized_mean),
          'output <=209568; actual sample right inverse checked on255 parity targets', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    run(parser.parse_args().verify)
