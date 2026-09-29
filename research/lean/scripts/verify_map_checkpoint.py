"""Audit the proved structural checkpoint without claiming full spectra."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'scripts/map_data'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    semantic = json.loads((DATA/'structural_verification.json').read_text())
    assert semantic['status'] == 'PASS'
    expected = {f'SpinCodes/Structured/{n}.lean' for n in
        ['PackedMapDefs','ConcreteMapData','PackedMap','ConcreteMaps','ConcreteMapPin',
         'MapSpectrumDefs','MapSpectrum','MapSpectrumSumDefs','MapSpectrumSum']}
    assert {r['path'] for r in semantic['modules']} == expected
    for record in semantic['modules']:
        assert record['exit_code'] == 0
        assert sha(ROOT/record['path']) == record['source_sha256']
    for name in ['kernel_basis.json','kernel_sample.json']:
        report = json.loads((DATA/name).read_text())
        assert report['status'] == 'PASS'
        for record in report['blocks']:
            assert record['exit_code'] == 0
            assert sha(ROOT/record['path']) == record['source_sha256']
        for dependency,digest in report['dependency_olean_sha256'].items():
            assert sha(ROOT/f'.lake/build/lib/lean/SpinCodes/Structured/{dependency}.olean') == digest
    axioms = (DATA/'axioms.log').read_text(encoding='utf-8-sig')
    for theorem in ['ConcreteMaps.Aset_injective','ConcreteMaps.Cset_surjective',
                    'ConcreteMaps.CtransposeSet_injective','ConcreteMaps.C_inputBasis_weight',
                    'ConcreteMaps.C_inputBasis_distinct','MapSpectrum.block_split_weights',
                    'MapSpectrum.histogram_blocks']:
        assert f"'Spin.Structured.{theorem}' depends on axioms: [propext, Classical.choice, Quot.sound]" in axioms
    assert not any(word in axioms for word in ['sorryAx','ofReduceBool','trustCompiler'])
    assert 'ALL CHECKS PASS' in (DATA/'invariant_check.log').read_text(encoding='utf-8-sig')
    paper = ROOT.parent/'paper/structured_imt_appendix.tex'
    result = dict(status='PASS',scope='Concrete map algebra, rank, and column properties. Full spectra, transfer-law instantiation, and distance theorem remain open.',
        paper_sha256=sha(paper), modules=semantic['modules'],
        finite_checks=13, spectrum_samples=['kernel_basis.json','kernel_sample.json'],
        verification_logs=['axioms.log','invariant_check.log'])
    (DATA/'checkpoint.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: concrete maps, rank and column properties, counting lemmas, standard axioms, and project invariants.')
    print('Full spectra and the concrete distance theorem remain open.')


if __name__ == '__main__':
    main()
