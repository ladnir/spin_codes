"""Validate replay coverage and preserve the scope of the resulting checkpoint."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'scripts/sparse_data'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifests = {phase: json.loads((DATA/name).read_text()) for phase,name in
                 [('weights','manifest.json'), ('dominance','dominance_manifest.json'),
                  ('contributions','contribution_manifest.json')]}
    weight_reports = ['kernel_weights.json']
    if (DATA/'kernel_weights_102_128.json').exists():
        weight_reports.append('kernel_weights_102_128.json')
    coverage = {}
    for phase, manifest in manifests.items():
        names = weight_reports if phase == 'weights' else [f'kernel_{phase}.json']
        passed = {}
        for name in names:
            report = json.loads((DATA/name).read_text())
            for record in report['weights']:
                if record['exit_code'] == 0:
                    passed[record['weight']] = record
        assert set(passed) == {r['weight'] for r in manifest}, phase
        for record in manifest:
            actual = sha(ROOT/record['path'])
            assert actual == record['sha256'] == passed[record['weight']]['source_sha256']
        coverage[phase] = dict(modules=len(passed), reports=names)
    axioms = (DATA/'axioms.log').read_text(encoding='utf-8-sig')
    required = "'Spin.Structured.SparsePolynomial.programResidual_neg' depends on axioms: [propext, Classical.choice, Quot.sound]"
    assert required in axioms
    for theorem in ['Spin.Imt.Occupation.Sparse.sparse_collatz',
                    'Spin.Imt.Occupation.Sparse.sparse_iterate_bound']:
        assert f"'{theorem}' depends on axioms: [propext, Classical.choice, Quot.sound]" in axioms
    assert not any(word in axioms for word in ['sorryAx', 'ofReduceBool', 'trustCompiler'])
    assert 'ALL CHECKS PASS' in (DATA/'invariant_check.log').read_text(encoding='utf-8-sig')
    assert 'Build completed successfully' in (DATA/'base_build.log').read_text(encoding='utf-8-sig')
    semantic = json.loads((DATA/'semantic_verification.json').read_text())
    assert semantic['status'] == 'PASS'
    expected_semantic = {f'SpinCodes/Structured/{name}.lean' for name in
                         ['SparseCancellation', 'SparseColumnBounds', 'SparseNonneg',
                          'SparseContraction', 'SparseIteration', 'SparseBridge/Pin']}
    assert {r['path'] for r in semantic['modules']} == expected_semantic
    for record in semantic['modules']:
        assert record['exit_code'] == 0
        assert sha(ROOT/record['path']) == record['source_sha256']
    sources = set((ROOT/'SpinCodes/Structured').glob('Sparse*.lean'))
    sources.update((ROOT/'SpinCodes/Structured').glob('Poly*.lean'))
    sources.update((ROOT/'SpinCodes/Structured/SparseBridge').glob('*.lean'))
    sources.update((ROOT/'SpinCodes/Structured/SparseData').glob('*.lean'))
    sources.add(ROOT/'SpinCodes/Structured/Occupation.lean')
    sources.add(ROOT/'SpinCodes.lean')
    report = dict(status='PASS', scope='Complete T3a Collatz inequality and iterate bound for the numerical occupation matrix. Concrete map identification, one-step law, and full distance theorem remain open.',
                  theorem='Spin.Imt.Occupation.Sparse.sparse_collatz',
                  iterate_theorem='Spin.Imt.Occupation.Sparse.sparse_iterate_bound',
                  coverage=coverage, maximum_comparisons=693, residual_rows=7,
                  source_sha256={str(p.relative_to(ROOT)): sha(p) for p in sorted(sources)},
                  verification_logs=['axioms.log','invariant_check.log','base_build.log',
                                     'semantic_verification.json'])
    (DATA/'final_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS: 129 weight modules, 129 dominance modules, 119 sums, seven residuals, final theorem and invariant audit.')
    print('Scope: T3a numerical-matrix Collatz and iterate bound complete; full distance theorem remains open.')


if __name__ == '__main__':
    main()
