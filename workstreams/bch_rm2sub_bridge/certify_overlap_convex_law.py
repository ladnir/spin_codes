"""Exact convex-order majorant for one shared T80 row's overlap.

For independent uniform U,V in T80, a=|supp(U) intersect supp(V)| has
mean25. Replace a<=39 by endpoints0,39, then increase each a>39 atom
to its certified cap while preserving the mean. This increases every
convex expectation. The result is NOT a stochastic-order upper bound.
"""
import argparse
from fractions import Fraction as F
from pathlib import Path
import bridge as base
import screen_exponential_modes as caps_source


def build(verify=False):
    output = base.HERE/'generated/overlap_convex_law.json'
    old = base.read(output) if verify else None
    if old:
        for name, digest in old['local_sha256'].items():
            assert base.sha(base.HERE/name) == digest
    else:
        assert not output.exists()
    source = base.HERE/'generated/weight80_family.json'
    family = base.read(source)
    for name, digest in family['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    caps, sources = caps_source.latest_caps()
    size = family['shell_size_lower']
    law = {80:F(1, size)}
    for row in family['rows']:
        difference = row['difference_weight']
        a = 80-difference//2
        assert 40 <= a <= 61
        law[a] = F(caps[difference]*row['intersection_upper'], size*size)
    tail_mass = sum(law.values())
    assert tail_mass == base.decode(family['close_pair_probability_upper'])
    tail_mean = sum(a*p for a, p in law.items())
    law[39] = (25-tail_mean)/39
    law[0] = 1-tail_mass-law[39]
    assert all(p >= 0 for p in law.values())
    assert sum(law.values()) == 1 and sum(a*p for a, p in law.items()) == 25
    result = dict(status='EXACT_CONVEX_ORDER_MAJORANT_FOR_T80_ROW_OVERLAP',
        actual_overlap_mean=25, endpoint=39, actual_overlap_support=[0, 61, 80],
        support_description='All integers0..61, plus80; positive differences have weight>=38.',
        rows=[dict(overlap=a, probability=base.encode(p)) for a, p in sorted(law.items())],
        convex_order_only=True, tail_mass_upper=base.encode(tail_mass),
        variance=base.encode(sum((a-25)**2*p for a, p in law.items())),
        full_second_moment_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p)
                      for p in [Path(__file__), source, Path(caps_source.__file__)]+sources})
    if old:
        assert result == old
        print('Exact row-overlap convex-law replay passed', flush=True)
    else:
        base.write_new(output, result)
        print('Row-overlap convex majorant: mean25, variance', float(base.decode(result['variance'])),
              'tail mass', float(tail_mass), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    build(parser.parse_args().verify)
