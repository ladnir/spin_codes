"""A lower bound on the current comparison, NOT on actual setup failure.

One artificial-law composition and one multinomial count vector force
every globally covering log-convex envelope comparison above2^23000.
This remains true with a different envelope for each total overlap.
"""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path
from flint import arb, arb_poly, ctx
import bridge as base
import tightened_occupancy as tight


def run(verify=False):
    output = base.HERE/'generated/convex_route_obstruction.json'
    old = base.read(output) if verify else None
    if old:
        for name, digest in old['local_sha256'].items():
            assert base.sha(base.HERE/name) == digest
    else:
        assert not output.exists()
    law_path = base.HERE/'generated/johnson_convex_law.json'
    law = base.read(law_path)
    for name, digest in law['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    probabilities = {row['overlap']:base.decode(row['probability']) for row in law['rows']}
    composition = [(80, 1592), (61, 1000), (60, 18)]
    assert sum(count for _, count in composition) == 2610
    total = sum(a*count for a, count in composition)
    assert total == 189440 == 256*740
    assert all(probabilities[a] > 0 for a, _ in composition)
    ctx.prec = 512 if verify else 256
    support = -arb(math.comb(8189, 2610)).log()
    row_mass = arb(math.factorial(2610)).log()
    for a, count in composition:
        p = probabilities[a]
        row_mass += count*(arb(p.numerator)/p.denominator).log()-arb(math.factorial(count)).log()
    multinomial_atom = arb(math.factorial(total)).log()-256*arb(math.factorial(740)).log()-total*arb(256).log()
    spectrum = tight.kernel_spectrum()
    kernel = arb_poly([spectrum.get(j, 0) for j in range(129)])**64
    beta = kernel[740]/math.comb(8192, 740)
    assert beta > 0
    bound = ((support+row_mass+multinomial_atom-256*beta.log())/arb(2).log()).lower()
    from audit_bch_q1_full_arb import rational
    exact = rational(bound)
    assert exact > 23000
    if old:
        assert exact >= base.decode(old['comparison_lower_log2'])
        print('512-bit comparison-obstruction replay passed', flush=True)
        return
    rounded = F(math.floor(exact*(1 << 24)), 1 << 24)
    base.write_new(output, dict(status='RIGOROUS_OBSTRUCTION_TO_CURRENT_CONVEX_COMPARISON_ONLY',
        comparison_lower_log2=base.encode(rounded), lower_exceeds2pow23000=True,
        common_support_rows=2610, overlap_composition=[dict(overlap=a, rows=count) for a, count in composition],
        total_overlap=total, multinomial_region_counts=[740]*256,
        forced_kernel_ratio='R(740,740,740)=1/beta740',
        applies_to_per_total_envelopes=True, actual_failure_probability_lower_bound=False,
        actual_target_refuted=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in
                      [Path(__file__), law_path, Path(tight.__file__), base.HERE/'inputs/t128_s15_b_kernel_spectrum.json']}))
    print('Current convex-comparison lower log2 >', float(rounded),
          '; this does NOT lower-bound actual setup failure', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    run(parser.parse_args().verify)
