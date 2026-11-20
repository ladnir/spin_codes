"""Read-only exact replay: no optimizer and no floating-point arithmetic."""
import math
from fractions import Fraction as F
import bridge as base
import christoffel_caps as christoffel
import screen_exponential_modes as cap_source
from certify_tail_parity_mixing import kraw_values


def run():
    caps, _ = cap_source.latest_caps()
    saved = base.read(base.HERE/'generated/weight80_fourth.json')
    for name, digest in saved['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    family = base.read(base.HERE/'generated/weight80_family.json')
    a = family['shell_size_lower']
    assert saved['occupation'] == 2620 and saved['row_weight'] == 80
    assert [row['difference_weight'] for row in saved['rows']] == [0]+list(range(38, 161, 2))
    alpha, slope = base.decode(saved['dual_intercept']), base.decode(saved['dual_slope'])
    second = alpha+110*slope
    for row in saved['rows']:
        c = row['difference_weight']
        upper = base.decode(row['probability_upper'])
        if c == 0:
            assert upper == F(1, a)
        else:
            i, j = c//2, 80-c//2
            kernel = sum((F(christoffel.kraw(c, r, i)**2 *
                            christoffel.kraw(256-c, s, j)**2,
                            math.comb(c, r)*math.comb(256-c, s))
                          for r in range(15) for s in range(15-r)), F(0))
            exact_cap = min(caps[80], math.comb(c, i)*math.comb(256-c, j),
                            math.floor(F(1 << 128)/kernel))
            assert row['intersection_upper'] >= exact_cap
            assert upper == min(F(1), F(caps[c]*row['intersection_upper'], a*a))
        hinge = base.decode(row['dual_hinge'])
        assert hinge == max(F(0), (55-c//2)**2-alpha-slope*c)
        second += upper*hinge
    assert second == base.decode(saved['centered_inner_product_second_moment_upper'])
    p = F(5, 16)
    v = p*(1-p)
    pair_mean = F(80*79, 256*255)
    elementary = pair_mean*(p-pair_mean)
    tau = min(elementary, (second-F(256*256, 255)*v*v)/(256*255))
    assert elementary == base.decode(saved['elementary_variance_upper'])
    assert tau == base.decode(saved['off_diagonal_covariance_variance_upper'])
    assert tau/(2620*v*v) == base.decode(saved['normalized_shared_edge_mixed_fourth_upper'])
    assert tau/(2620*v*v) < F(1, 36000)
    characters = kraw_values(256, 80)
    single = max(F(abs(characters[j]), math.comb(256, j)) for j in range(1, 256))
    assert single == F(3, 8)
    far_tables = [kraw_values(256, c) for c in range(82, 161, 2)]
    far = max(F(abs(values[j]), math.comb(256, j))
              for values in far_tables for j in range(1, 256))
    assert far == F(23, 64)
    close = base.decode(family['close_pair_probability_upper'])
    paired = close+(1-close)*far
    assert paired >= single*single
    nontrivial = (1 << 256)-2
    epsilon = nontrivial*single**2620+F(nontrivial**2, 4)*paired**2620
    assert epsilon < F(1, 1 << 3300)
    box = base.read(base.HERE/'generated/pair_type_central_box_outward.json')
    for name, digest in box['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    assert len(box['rows']) == 125 and not box['all_below_eleven_tenths']
    assert all(base.decode(row['ratio_upper']) < F(1101, 1000) for row in box['rows'])
    print('Exact rational replay passed: all 63 difference caps, fourth dual, '
          'normalized bound <1/36000, weight80 parity relative error <2^-3300, '
          'and 125-type <1.101 ledger', flush=True)


if __name__ == '__main__':
    run()
