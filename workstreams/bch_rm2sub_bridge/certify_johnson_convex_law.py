"""A certified convex upper law from Hahn-constrained stop-loss bounds."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import bridge as base
import johnson_moment_lp as moments


def lower_convex_hull(values):
    vertices = []
    for t, value in enumerate(values):
        while len(vertices) >= 2:
            a, b = vertices[-2:]
            if (values[b]-values[a])/(b-a) < (value-values[b])/(t-b):
                break
            vertices.pop()
        vertices.append(t)
    hull = [None]*len(values)
    for a, b in zip(vertices, vertices[1:]):
        for t in range(a, b+1):
            hull[t] = (values[a]*(b-t)+values[b]*(t-a))/(b-a)
    return vertices, hull


def run(verify=False):
    output = base.HERE/'generated/johnson_convex_law.json'
    old = base.read(output) if verify else None
    if old:
        for name, digest in old['local_sha256'].items():
            assert base.sha(base.HERE/name) == digest
    else:
        assert not output.exists()
    problem = moments.MomentProblem(24)
    ext_path = base.HERE/'generated/overlap_extremal_law.json'
    ext = base.read(ext_path)
    for name, digest in ext['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    ext_calls = list(map(base.decode, ext['stop_loss_expectations']))
    uppers, proofs = [F(25)], []
    for t in range(1, 80):
        values = [F(max(a-t, 0)) for a in problem.atoms]
        if old:
            proof = old['stop_loss_proofs'][t-1]
            candidates = [problem.replay(values, row) for row in proof['candidates']]
        else:
            candidates_rows = [problem.solve(values, degree) for degree in [10, 16, 20, 24]]
            proof = dict(threshold=t, candidates=candidates_rows)
            candidates = [base.decode(row['upper']) for row in candidates_rows]
        assert proof['threshold'] == t
        upper = min([ext_calls[t]]+candidates)
        assert upper >= max(F(0), F(25-t))
        uppers.append(upper)
        proofs.append(proof)
        if t % 10 == 0:
            print('Johnson stop-loss', t, float(upper), 'old', float(ext_calls[t]), flush=True)
    uppers.append(F(0))
    vertices, calls = lower_convex_hull(uppers)
    assert calls[0] == 25 and calls[80] == 0
    assert all(c <= u for c, u in zip(calls, uppers))
    differences = [calls[t+1]-calls[t] for t in range(80)]
    assert all(-1 <= d <= 0 for d in differences)
    assert all(a <= b for a, b in zip(differences, differences[1:]))
    law = [1+differences[0]]+[differences[a]-differences[a-1] for a in range(1, 80)]+[-differences[-1]]
    assert min(law) >= 0 and sum(law) == 1 and sum(a*p for a, p in enumerate(law)) == 25
    assert all(sum(max(a-t, 0)*p for a, p in enumerate(law)) == calls[t] for t in range(81))
    variance = sum((a-25)**2*p for a, p in enumerate(law))
    result = dict(status='EXACT_HAHN_STOP_LOSS_CONVEX_UPPER_LAW',
        rows=[dict(overlap=a, probability=base.encode(p)) for a, p in enumerate(law) if p],
        actual_mean=25, variance=base.encode(variance), hull_vertices=vertices,
        stop_loss_upper=[base.encode(v) for v in uppers], stop_loss_expectations=[base.encode(v) for v in calls],
        stop_loss_proofs=proofs, convex_order_only=True, full_second_moment_certified=False,
        mathematical_source='https://arxiv.org/html/2405.07666v2#S4.SS2',
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in
                      [Path(__file__), Path(moments.__file__), base.HERE/'johnson_overlap.py',
                       base.HERE/'generated/weight80_fourth.json', ext_path]})
    if old:
        assert result == old
        print('Exact Hahn convex-law replay passed', flush=True)
    else:
        base.write_new(output, result)
        print('Certified Hahn convex law variance', float(variance), 'hull vertices', vertices, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    run(parser.parse_args().verify)
