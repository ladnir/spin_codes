"""Assemble all checked Fourier/cap rows; keep the concrete spectrum premise explicit."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'scripts/map_data'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    report = json.loads((DATA/'kernel_fiber_data.json').read_text())
    assert report['status'] == 'PASS', 'The full Fourier-row replay has not passed.'
    manifest = json.loads((DATA/'fiber_manifest.json').read_text())
    passed = {r['weight']:r for r in report['weights'] if r['exit_code'] == 0}
    assert len(manifest) == 129 and set(passed) == {r['weight'] for r in manifest} == set(range(129))
    for name,digest in report['dependency_olean_sha256'].items():
        assert sha(ROOT/f'.lake/build/lib/lean/SpinCodes/Structured/{name}.olean') == digest
    for record in manifest:
        result = passed[record['weight']]
        assert sha(ROOT/record['path']) == record['sha256'] == result['source_sha256']
        assert sha((ROOT/'.lake/build/lib/lean'/record['path']).with_suffix('.olean')) == result['output_sha256']
    candidates = json.loads((DATA/'fiber_candidates.json').read_text())['records']
    lines = ['import SpinCodes.Structured.FiberFrozen']
    lines += [f'import SpinCodes.Structured.FiberNumericsData.Weight{j}' for j in range(129)]
    lines += ['noncomputable section', 'namespace Spin.Structured.FiberNumerics',
              'set_option maxRecDepth 100000', 'set_option maxHeartbeats 0',
              'private theorem cases129 (P : Fin 129 → Prop)']
    # Keep elimination abstract so elaboration cannot unfold concrete finite sets.
    lines += [f'    (h{j} : P ⟨{j}, by decide⟩)' for j in range(129)]
    lines += ['    (j : Fin 129) : P j := by', '  fin_cases j']
    lines += [f'  · exact h{j}' for j in range(129)]
    lines += ['private theorem complement_of_count (j k cap : Nat) {q : Finset (Fin 19)}',
              '    (hq : q ≠ ∅) (hk : (ConcreteMaps.kernelLayer j).card = k)',
              '    (hc : polyChoose 128 j - k ≤ cap) :',
              '    (ConcreteMaps.syndromeFiber j q).card ≤ cap := by',
              '  apply ConcreteMaps.fiber_cap_complement j cap hq',
              '  rw [hk, ← polyChoose_eq]',
              '  exact hc']
    lines += [
              'def pairTotals : List Nat := ['+', '.join(str(r['pairs']) for r in candidates)+']',
              'variable (hs : ∀ w : Fin 129, weightCounts ConcreteMaps.CtransposeSet w = Data.spectrum.getD w 0)',
              'include hs',
              'theorem kernel_counts (j : Fin 129) :',
              '    (ConcreteMaps.kernelLayer j).card = Data.kernels.getD j 0 := by',
              '  apply cases129 (fun j => (ConcreteMaps.kernelLayer j).card = Data.kernels.getD j 0)',
              '    <;> clear j']
    for j in range(129):
        lines += [f'  · exact kernel_count_of_certificate Data.spectrum (j := {j}) (k := Data.kernels.getD {j} 0) hs Data.values{j}_checked Data.signed{j}_checked']
    lines += ['theorem pair_counts (j : Fin 129) :',
              '    (ConcreteMaps.equalSyndromePairs j).card = pairTotals.getD j 0 := by',
              '  apply cases129 (fun j => (ConcreteMaps.equalSyndromePairs j).card = pairTotals.getD j 0)',
              '    <;> clear j']
    for j in range(129):
        lines += [f'  · exact pair_count_of_certificate Data.spectrum (j := {j}) (V := pairTotals.getD {j} 0) hs Data.values{j}_checked Data.square{j}_checked']
    for record in manifest:
        j,method = record['weight'],record['method']
        lines += [f'private theorem fiber_cap_{j} {{q : Finset (Fin 19)}} (_hq : q ≠ ∅) :',
                  f'    (ConcreteMaps.syndromeFiber {j} q).card ≤ Data.caps.getD {j} 0 := by']
        if method == 'fourier':
            lines += [f'  apply fiber_cap_of_fourier_certificate Data.spectrum (j := {j}) (cap := Data.caps.getD {j} 0) hs Data.values{j}_checked',
                      f'  rw [Data.absolute{j}_checked]',f'  exact Data.cap{j}_checked']
        elif method == 'complement':
            lines += [f'  exact complement_of_count {j} _ _ _hq',
                      f'    (kernel_counts hs (⟨{j}, by decide⟩ : Fin 129)) Data.cap{j}_checked']
        elif method in ['packing','packing_compl']:
            lines += [f'  apply ConcreteMaps.fiber_cap_{method} {j} (Data.caps.getD {j} 0) q',
                      f'  simpa only [polyChoose_eq] using Data.cap{j}_checked']
        else:
            raise AssertionError(method)
    lines += ['theorem fiber_caps (j : Fin 129) {q : Finset (Fin 19)} (hq : q ≠ ∅) :',
              '    (ConcreteMaps.syndromeFiber j q).card ≤ Data.caps.getD j 0 := by',
              '  apply cases129 (fun j => (ConcreteMaps.syndromeFiber j q).card ≤ Data.caps.getD j 0)',
              '    <;> clear j']
    lines += [f'  · exact fiber_cap_{j} hs hq' for j in range(129)]
    lines += ['theorem kernel_counts_frozen (j : Fin 129) :',
              '    (ConcreteMaps.kernelLayer j).card = (SparsePolynomial.weightN j).kernel :=',
              '  (kernel_counts hs j).trans (kernel_frozen_eq j).symm',
              'theorem fiber_caps_frozen (j : Fin 129) {q : Finset (Fin 19)} (hq : q ≠ ∅) :',
              '    (ConcreteMaps.syndromeFiber j q).card ≤ (SparsePolynomial.weightN j).cap := by',
              '  rw [cap_frozen_eq]', '  exact fiber_caps hs j hq',
              'end Spin.Structured.FiberNumerics','']
    (ROOT/'SpinCodes/Structured/FiberNumericsAll.lean').write_text('\n'.join(lines),encoding='utf-8')
    print('Emitted FiberNumericsAll.lean. It still requires Lean compilation and the concrete spectrum proof.')


if __name__ == '__main__':
    main()
