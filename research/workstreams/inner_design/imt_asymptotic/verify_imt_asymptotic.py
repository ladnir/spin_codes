"""Assemble IMT numerical evidence; analytic proof review remains separate.

Reconstructs the sparse polynomials and small outward inequalities. For the
658 dense boxes, authenticates the retained 512-bit replay and checks its
geometry/enclosure; use certify_imt_dense.py for a fresh dense replay.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

from flint import ctx
import screen
import certify_imt_dense as dense
import certify_imt_sparse as sparse_cert
import certify_imt_one_two as one_two

sys.path.insert(0, str(screen.ROOT/'artifact'))
import reproduce


def authenticate(record):
    for name, digest in record['source_sha256'].items():
        path = reproduce.checked_path(name)
        if screen.sha(path) != digest:
            raise ValueError(f'Changed certificate dependency: {name}')


def check_fixed(saved, segments):
    num = dense.number
    sigma, tau, v = (F(saved[k]) for k in ('sigma', 'tau', 'v'))
    assert (sigma, tau, v) == (F(127,250), F(133,125), F(3,1600))
    r, p0 = F(1,524287), F(262144,524287)
    assert saved['impulse_upper'] == [['0','1'],[str(r),'1']]
    assert F(saved['p0']) == p0
    assert F(saved['delta']) == F(11,100)
    assert F(saved['block_constant']) == F(39,4)
    d = (1-sigma)/sigma
    critical = F(4,3)/tau-1/d
    assert 0 < critical < 1
    spacing = max(dense.upper(-num(tau*y)+num(F(4,3))*(1+num(d*y)).log())
                  for y in (F(0), critical, F(1)))
    expected = {(s.index, x) for s in segments for x in (s.lower, s.upper)}
    assert len(saved['checks']) == len(expected) == 62
    assert {(row['segment'],F(row['x'])) for row in saved['checks']} == expected
    values = []
    for row in saved['checks']:
        segment = segments[row['segment']]
        x, u = F(row['x']), F(row['u'])
        assert u > 0
        norm = max(1+u*sigma*v, sigma+u*(r/v+sigma))
        entropy = -num(x)*num(x).log()-num(1-x)*num(1-x).log()
        val = (num(segment.slope*x+segment.intercept)-entropy-num(x)*num(u).log()
               +num(norm).log()+num(spacing)+num(F(11,100)*tau/p0)
               +num(2).log()/num(F(39,4)))
        upper = dense.upper(val)
        assert upper <= F(row['exponent_upper']) < 0
        values.append(upper)
    return max(values)


def main():
    if not __debug__:
        raise RuntimeError('Assertions are certificate checks; do not use -O.')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    reproduce.check_asymptotic()
    names = ['DENSE_OUTWARD_D1099_v3.json','DENSE_REPLAY_D1099_v3.json',
             'SPARSE_EXACT.json','SPARSE_EXACT_REPLAY.json',
             'FIXED_LIMIT.json','ONE_TWO_OUTWARD.json']
    records = [json.loads((screen.HERE/name).read_text()) for name in names]
    for record in records:
        authenticate(record)
    original, replay, sparse, sparse_replay, fixed, small = records
    model = screen.Model()
    inner = model.engine.identity()['inner']
    assert all(record['inner'] == inner for record in records[:4])
    assert (inner['t'],inner['s'],inner['transvection_rounds'],inner['feedback_name']) == (128,19,1,'weight5_seed0')
    assert all(0 < c < 1<<19 for c in inner['expansion_columns']+inner['feedback_columns'])
    majorant = screen.FROZEN/'golay_ba3_concave_majorant.json'
    segments = screen.load_segments(majorant)
    assert len(segments) == 31
    assert segments[0].lower == F(13,125) and segments[-1].upper == F(112,125)
    assert all(a.upper == b.lower for a,b in zip(segments,segments[1:]))
    for record in (original,replay):
        assert record['delta'] == '1099/10000'
        assert record['alpha_range'] == ['1/10000','1']
        assert record['row_density_range'] == ['13/125','112/125']
        assert record['boxes'] == len(record['leaves']) == 658
        dense.check_geometry(record['leaves'],segments)
    assert original['precision_bits'] == 256 and replay['precision_bits'] == 512
    for a,b in zip(original['leaves'],replay['leaves']):
        assert all(a[k] == b[k] for k in ('segment','alpha','row_density','rational_witness'))
        assert F(b['exponent_upper']) <= F(a['exponent_upper']) < 0
    dense_max = max(F(row['exponent_upper']) for row in original['leaves'])
    print('Dense replay authenticated; all 658 boxes cover the required domain.',flush=True)

    assert sparse['alpha_interval'] == ['0','1/10000']
    assert sparse['contraction'] == '1-96*alpha'
    assert sparse['beta'] == '(4/5)*alpha' and sparse['z'] == '1-(8/5)*alpha'
    assert sparse['strict_for_positive_alpha'] is True
    vector, _, residuals, dominance = sparse_cert.build(model)
    assert [[str(c) for c in p.coeffs()] for p in vector] == sparse['vector_coefficients']
    assert len(residuals) == 7 and len(dominance) == sparse['dominance_comparisons']
    checks = []
    for i, poly in enumerate(residuals):
        assert poly[0] == 0
        checks.append(dict(row=i,degree=poly.degree(),
                           coefficients_sha256=hashlib.sha256(str(poly).encode()).hexdigest(),
                           **sparse_cert.sign_certificate(poly,strict=True)))
    assert checks == sparse['row_checks'] == sparse_replay['row_checks']
    assert vector[0].degree() == 0 and vector[0][0] == 1
    minimum = min(min(F(str(p[0])), F(str(p[0]+p[1]))) for p in vector)
    assert all(p.degree() <= 1 for p in vector) and minimum >= F(1,2048)
    print('Sparse polynomial certificate reconstructed exactly; uniform prefactor <=2048.',flush=True)

    ctx.prec = 512
    num = dense.number
    likelihood = []
    for segment in segments:
        for x in (segment.lower,segment.upper):
            entropy = -num(x)*num(x).log()-num(1-x)*num(1-x).log()
            likelihood.append(dense.upper(num(segment.slope*x+segment.intercept)-entropy+num(2).log()/2))
    assert max(likelihood) < F(1281,100000)
    coefficient = (num(2).log()/2+num(F(5,8)).log()
                   +num(F(3,5)+F(11,100)*F(8,5))/num(1-F(8,5)*F(1,10000))
                   -num(F(96,128))+num(F(1281,100000))+num(2).log()/num(F(39,4)))
    coefficient_upper = dense.upper(coefficient)
    assert coefficient_upper <= F(sparse['conditional_structured_exponent_coefficient_upper'])
    conditioning_upper = dense.upper(num(512).log()/4096)
    assert conditioning_upper < F(2,1000)
    # For large b, epsilon_b+1/b+ln(2048)/(Qb) is uniformly below .004.
    uniform_upper = coefficient_upper+conditioning_upper+F(4,1000)
    assert uniform_upper < -F(6,1000)
    fixed_max = check_fixed(fixed,segments)
    assert small['delta'] == '11/100' and small['block_constant'] == '39/4'
    assert small['outer_likelihood_excess'] == '1281/100000'
    small_replay = one_two.evaluate(512)
    assert [r['q'] for r in small['checks']] == [1,2]
    for a,b in zip(small['checks'],small_replay):
        assert F(b['exponent_per_outer_bit_upper']) <= F(a['exponent_per_outer_bit_upper']) < 0
    print('Outer likelihood, fixed-Q inequalities, and sparse cutoff pass at 512 bits.',flush=True)

    paths = [screen.HERE/name for name in names]
    paths += [Path(__file__),screen.HERE/'ASYMPTOTIC_PROOF.md',Path(reproduce.__file__),
              *sorted(screen.HERE.glob('test_*.py'))]
    paths += [screen.ROOT/name for name,_ in reproduce.asymptotic_entries()]
    result = dict(status='IMT_D1099_COMPLETE_PROOF_DRAFT_NUMERICAL_ASSEMBLY_PASSED',
                  delta='1099/10000',block_constant='39/4',inner=inner,
                  dense_boxes=658,dense_maximum_upper=str(dense_max),
                  sparse_exact_rows=7,sparse_prefactor=2048,
                  sparse_coefficient_upper=str(coefficient_upper),
                  sparse_cutoff=4096,uniform_large_Q_exponent_upper=str(uniform_upper),
                  fixed_Q_ge_3_maximum_upper=str(fixed_max),
                  one_two_512=small_replay,outer_likelihood_excess_upper=str(max(likelihood)),
                  outer_manifest_entries=31,fresh_dense_replay=False,
                  analytic_proof_machine_checked=False,independent_analytic_review_pending=True,
                  paper_or_default_changed=False,
                  source_sha256={p.relative_to(screen.ROOT).as_posix():screen.sha(p) for p in paths})
    with args.output.open('x',encoding='utf-8') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print('PASS: numerical assembly. Analytic proof review is separate.',flush=True)


if __name__ == '__main__':
    main()
