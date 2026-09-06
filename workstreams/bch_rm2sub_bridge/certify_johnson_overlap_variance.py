"""Exact rational Hahn-positivity dual for actual T80 overlap variance."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import bridge as base
from johnson_overlap import hahn


def run(verify=False):
    output = base.HERE/'generated/johnson_overlap_variance_attempt01.json'
    old = base.read(output) if verify else None
    if old:
        for name, digest in old['local_sha256'].items():
            assert base.sha(base.HERE/name) == digest
    else:
        assert not output.exists()
    source = base.HERE/'generated/johnson_overlap_screen.json'
    screen = base.read(source)
    for name, digest in screen['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    caps_source = base.HERE/'generated/weight80_fourth.json'
    saved = base.read(caps_source)
    for name, digest in saved['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    caps = {80-row['difference_weight']//2:base.decode(row['probability_upper']) for row in saved['rows']}
    proposal = screen['rows'][-1]
    assert proposal['maximum_degree'] == 80 and proposal['success']
    intercept, slope = map(F.from_float, proposal['equality_dual'])
    multipliers = [max(F(0), F.from_float(v)) for v in proposal['positivity_multipliers']]
    norms = list(map(base.decode, screen['normalizers']))
    assert len(multipliers) == len(norms) == 79 and min(norms) > 0
    hinges = []
    upper = intercept+25*slope
    for a in sorted(caps):
        affine = intercept+slope*a-sum((mult*hahn(256, 80, j, 80-a)/norm
                                         for j, mult, norm in zip(range(2, 81), multipliers, norms)), F(0))
        hinge = max(F(0), F((a-25)**2)-affine)
        assert affine+hinge >= (a-25)**2
        upper += caps[a]*hinge
        hinges.append(dict(overlap=a, hinge=base.encode(hinge)))
    print('Exact dual diagnostic', float(upper), flush=True)
    assert 0 < upper
    result = dict(status='EXACT_HAHN_POSITIVITY_T80_OVERLAP_VARIANCE_UPPER',
        variance_upper=base.encode(upper), exact_mean=25, upper_below48=upper<48,
        dual_intercept=base.encode(intercept), dual_slope=base.encode(slope),
        positivity_multipliers=[base.encode(v) for v in multipliers],
        normalizers=screen['normalizers'], hinges=hinges,
        premise='Normalized Hahn pair averages are nonnegative on any constant-weight set; mean25 uses T80 translation symmetry.',
        mathematical_source='https://arxiv.org/html/2405.07666v2#S4.SS2',
        full_second_moment_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in
                      [Path(__file__), base.HERE/'johnson_overlap.py', source, caps_source]})
    if old:
        assert result == old
        print('Exact Hahn-positivity variance replay passed', flush=True)
    else:
        base.write_new(output, result)
        print('Exact actual row-overlap variance upper', float(upper), 'desired <48:', upper<48, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    run(parser.parse_args().verify)
