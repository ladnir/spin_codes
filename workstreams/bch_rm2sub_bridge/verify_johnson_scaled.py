"""Independent exact replay of the cap-scaled Hahn variance duals."""
from fractions import Fraction as F
import bridge as base
from johnson_moment_lp import MomentProblem


def run():
    saved = base.read(base.HERE/'generated/johnson_scaled_variance.json')
    for name, digest in saved['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    problem = MomentProblem(80)
    assert problem.atoms == saved['atoms']
    assert problem.norms == list(map(base.decode, saved['normalizers']))
    values = [F((a-25)**2) for a in problem.atoms]
    bounds = []
    for row in saved['rows']:
        assert row['optimizer_success']
        intercept, slope = map(base.decode, row['equality_dual'])
        upper = problem.exact_bound(values, row['maximum_degree'], intercept, slope,
                                    list(map(base.decode, row['positivity_multipliers'])))
        assert upper == base.decode(row['variance_upper'])
        bounds.append(upper)
    assert min(bounds) < 44
    print('All7 exact scaled-Hahn duals replayed; best variance upper', float(min(bounds)), flush=True)


if __name__ == '__main__':
    run()
