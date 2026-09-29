"""Emit the total-spectrum checks and semantic assembly after complete replay."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'scripts/map_data'


def main():
    manifest = json.loads((DATA/'manifest.json').read_text())
    report = json.loads((DATA/'kernel_blocks.json').read_text())
    assert report['status'] == 'PASS', 'The full block replay has not passed.'
    passed = {r['block']: r for r in report['blocks'] if r['exit_code'] == 0}
    assert set(passed) == set(range(512))
    assert len(manifest) == 512 and {r['block'] for r in manifest} == set(range(512))
    for dependency, digest in report['dependency_olean_sha256'].items():
        artifact = ROOT/f'.lake/build/lib/lean/SpinCodes/Structured/{dependency}.olean'
        assert hashlib.sha256(artifact.read_bytes()).hexdigest() == digest
    for record in manifest:
        assert record['size'] == 1024 and record['start'] == record['block']*1024
        assert hashlib.sha256((ROOT/record['path']).read_bytes()).hexdigest() == record['sha256'] == passed[record['block']]['source_sha256']
        artifact = (ROOT/'.lake/build/lib/lean'/record['path']).with_suffix('.olean')
        assert hashlib.sha256(artifact.read_bytes()).hexdigest() == passed[record['block']]['output_sha256']
    totals = json.loads((DATA/'spectra.json').read_text())
    pure = ['import SpinCodes.Structured.MapSpectrumSumDefs']
    pure += [f'import SpinCodes.Structured.MapSpectrumData.Block{i}' for i in range(512)]
    pure += ['namespace Spin.Structured.MapSpectrum',
             'set_option maxRecDepth 1000000', 'set_option maxHeartbeats 0']
    for name in ['a','ct']:
        pure += [f'def {name}Blocks : List (List Nat) := [' + ', '.join(f'{name}{i}' for i in range(512)) + ']',
                 f'def {name}Totals : List Nat := [' + ', '.join(map(str,totals[name])) + ']',
                 f'theorem {name}_lengths_checked : {name}Blocks.all (fun h => h.length == 129) = true := by decide +kernel',
                 f'theorem {name}_totals_checked : sumHistograms {name}Blocks = {name}Totals := by decide +kernel']
    pure += ['end Spin.Structured.MapSpectrum', '']
    (ROOT/'SpinCodes/Structured/MapSpectrumTotals.lean').write_text('\n'.join(pure),encoding='utf-8')
    proof = ['import SpinCodes.Structured.MapSpectrumTotals',
             'import SpinCodes.Structured.MapSpectrumSum',
             'import SpinCodes.Structured.ConcreteWeightSums',
             'namespace Spin.Structured.MapSpectrum',
             'set_option maxRecDepth 100000', 'set_option maxHeartbeats 0']
    for name,rows in [('a','aRows'),('ct','cTransposeRows')]:
        for i in range(512):
            proof += [f'theorem {name}{i}_checked : block ConcreteMaps.{rows} {i*1024} 1024 = {name}{i} :=',
                f'  (block_split_weights ConcreteMaps.{rows} 10 {i} (by decide) {name}Low {name}{i}High {name}{i}Weights',
                f'    {name}Low_checked {name}{i}_high_checked {name}{i}_weights_checked).trans {name}{i}_hist_checked']
        term = '(rfl : ([] : List (List Nat)) = [])'
        for i in reversed(range(512)):
            term = f'(congrArg₂ List.cons {name}{i}_checked {term})'
        proof += [f'theorem {name}_blocks_eq : (List.range 512).map (fun b => block ConcreteMaps.{rows} (b * 1024) 1024) = {name}Blocks := by',
                  f'  exact {term}',
                  f'theorem histogram_{name}_getD {{i : ℕ}} (hi : i < 129) :',
                  f'    (histogram ConcreteMaps.{rows} (List.range (2 ^ 19))).getD i 0 = {name}Totals.getD i 0 := by',
                  f'  have hb := histogram_blocks ConcreteMaps.{rows} 1024 512 hi',
                  f'  have he := congrArg (fun hs : List (List ℕ) => (hs.map fun h => h.getD i 0).sum) {name}_blocks_eq',
                  '  simp only [List.map_map, Function.comp_def] at he',
                  f'  have ht := getD_sumHistograms {name}Blocks',
                  f'    (by simpa only [List.all_eq_true, beq_iff_eq] using {name}_lengths_checked) hi',
                  f'  rw [{name}_totals_checked] at ht',
                  '  exact hb.trans (he.trans ht.symm)']
    proof += ['end Spin.Structured.MapSpectrum',
              'namespace Spin.Structured.ConcreteMaps',
              'open MapSpectrum PackedMap']
    for name,mapname,rows in [('a','A','aRows'),('ct','Ctranspose','cTransposeRows')]:
        setname = 'Aset' if mapname == 'A' else 'CtransposeSet'
        proof += [f'theorem {mapname}_weight_card {{i : ℕ}} (hi : i < 129) :',
                  f'    (Finset.univ.filter fun q : Fin (2 ^ 19) => weight 128 ({mapname} q) = i).card = {name}Totals.getD i 0 := by',
                  f'  have h := histogram_range_card {rows} (2 ^ 19) hi',
                  f'  rw [histogram_{name}_getD hi] at h',
                  f'  exact h.symm',
                  f'theorem {setname}_weightCounts {{i : ℕ}} (hi : i < 129) :',
                  f'    weightCounts {setname} i = {name}Totals.getD i 0 :=',
                  f'  ({setname}_weightCounts_packed i).trans ({mapname}_weight_card hi)']
    proof += ['end Spin.Structured.ConcreteMaps', '']
    (ROOT/'SpinCodes/Structured/MapSpectrumBridge.lean').write_text('\n'.join(proof),encoding='utf-8')
    print('Emitted total-spectrum checks and the exact cardinality bridge. Compile both files next.')


if __name__ == '__main__':
    main()
