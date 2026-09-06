"""Greatest convex-order law under all retained row-distance atom caps.

Fill every atom outside[5,35] to its cap; use partial masses at5 and35.
A chord of any convex function through these two points proves maximality.
"""
import argparse
from fractions import Fraction as F
from pathlib import Path
import bridge as base
import verify_weight80_dependence as replay


def run(verify=False):
    output = base.HERE/'generated/overlap_extremal_law.json'
    old = base.read(output) if verify else None
    if old:
        for name, digest in old['local_sha256'].items():
            assert base.sha(base.HERE/name) == digest
    else:
        assert not output.exists()
    replay.run()
    source = base.HERE/'generated/weight80_fourth.json'
    saved = base.read(source)
    caps = {80-row['difference_weight']//2:base.decode(row['probability_upper'])
            for row in saved['rows']}
    low, high = 5, 35
    law = {a:(cap if a < low or a > high else F(0)) for a, cap in caps.items()}
    remaining_mass = 1-sum(law.values())
    remaining_mean = 25-sum(a*p for a, p in law.items())
    law[high] = (remaining_mean-low*remaining_mass)/(high-low)
    law[low] = remaining_mass-law[high]
    assert all(0 <= p <= caps[a] for a, p in law.items())
    assert sum(law.values()) == 1 and sum(a*p for a, p in law.items()) == 25
    variance = sum((a-25)**2*p for a, p in law.items())
    assert variance == base.decode(saved['centered_inner_product_second_moment_upper'])
    # Stop-loss functions span all discrete convex functions modulo affine
    # functions. Check the chord dual and saturation for every hinge.
    stop_losses = []
    for t in range(81):
        g = lambda a:max(a-t, 0)
        slope = F(g(high)-g(low), high-low)
        intercept = g(low)-slope*low
        upper = intercept+25*slope
        for a, cap in caps.items():
            residual = F(g(a))-intercept-slope*a
            if a < low or a > high:
                assert residual >= 0 and law[a] == cap
            elif low < a < high:
                assert residual <= 0 and law[a] == 0
            else:
                assert residual == 0
            upper += cap*max(residual, F(0))
        primal = sum(p*g(a) for a, p in law.items())
        assert upper == primal
        stop_losses.append(base.encode(primal))
    result = dict(status='EXACT_GREATEST_CONVEX_LAW_UNDER_ALL_ROW_ATOM_CAPS',
        actual_overlap_mean=25, boundary_atoms=[low, high], variance=base.encode(variance),
        rows=[dict(overlap=a, probability=base.encode(p), atom_cap=base.encode(caps[a]))
              for a, p in sorted(law.items())], stop_loss_expectations=stop_losses,
        convex_order_only=True, full_second_moment_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p)
                      for p in [Path(__file__), Path(replay.__file__), source]})
    if old:
        assert result == old
        print('Exact extremal row-overlap law replay passed', flush=True)
    else:
        base.write_new(output, result)
        print('Greatest convex law: boundary masses', float(law[5]), float(law[35]),
              'variance', float(variance), 'all81 stop-loss duals tight', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    run(parser.parse_args().verify)
