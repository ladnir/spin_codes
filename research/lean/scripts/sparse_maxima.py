"""Bind the numerical maximum certificates to real-valued transfer moments."""
import argparse
import json
from pathlib import Path
import re
from sparse_bridge_data import lean_pattern

ROOT = Path(__file__).resolve().parents[1]


def moment_proof(j):
    path = ROOT/f'SpinCodes/Structured/SparseBridge/Dominance{j}.lean'
    text = path.read_text(encoding='utf-8')
    selected = int(re.search(r'selected_checked : .*momentMax = (\d+)', text)[1])
    w = [48,56,64,72,80][selected]
    ns = f'Dominance{j}'
    lines = [f'theorem maximum_{j} (x : ℝ) (h0 : 0 ≤ x) (h1 : x ≤ 1) :',
             f'    Spin.Imt.Occupation.Sparse.maximumMoment {j} (1 - x / 6250) ≤',
             f'      (arbitrary {j} Data.weight{j}).eval x := by',
             f'  have hs : (arbitrary {j} Data.weight{j}).eval x = {ns}.m{selected}.eval x := by',
             f'    change (moment {w} {j}).eval x = _',
             f'    exact RatPoly.checkEq_sound _ _ {ns}.m{selected}_checked x',
             f'  apply maximumMoment_le_polynomial {j} x h1 (arbitrary {j} Data.weight{j})']
    for i in range(5):
        lines += [f'  · exact (RatPoly.checkEq_sound _ _ {ns}.m{i}_checked x).le.trans',
                  f'      ((RatPoly.checkLE_sound _ _ {ns}.m{i}_le_checked h0 h1).trans_eq hs.symm)']
    return lines


def low_proof(j, model):
    text = (ROOT/f'SpinCodes/Structured/SparseBridge/Dominance{j}.lean').read_text(encoding='utf-8')
    selected = int(re.search(r'selected_low_checked : .*lowMax = (\d+)', text)[1])
    ps = model['low'][str(j)]['patterns']
    literal = '[' + ', '.join(lean_pattern(p) for p in ps) + ']'
    ns = f'Dominance{j}'
    lines = [f'theorem lowMaximum_{j} (x : ℝ) (h0 : 0 ≤ x) (h1 : x ≤ 1) :',
             f'    lowMaximum {j} Data.weight{j} (1 - x / 6250) ≤',
             f'      (selectedLowPolynomial {j} Data.weight{j}).eval x := by',
             f'  have he : selectedLowPolynomial {j} Data.weight{j} = pattern {j} {lean_pattern(ps[selected])} := rfl',
             f'  have hs : (selectedLowPolynomial {j} Data.weight{j}).eval x = {ns}.low{selected}.eval x := by',
             f'    rw [he]', f'    exact RatPoly.checkEq_sound _ _ {ns}.low{selected}_checked x',
             f'  have hn : 0 ≤ (selectedLowPolynomial {j} Data.weight{j}).eval x := by',
             f'    rw [he, eval_pattern]',
             f'    exact Spin.Imt.Occupation.Sparse.pattern_nonneg _ _ (by linarith)',
             f'  change ({literal}.map fun p => Spin.Imt.Occupation.Sparse.pattern {j} p (1 - x / 6250)).foldr max 0 ≤ _',
             f'  apply maximumPattern_le_polynomial {j} _ x _ hn',
             '  intro p hp',
             '  simp only [List.mem_cons, List.not_mem_nil, or_false] at hp',
             '  rcases hp with ' + ' | '.join('rfl' for _ in ps)]
    for i in range(len(ps)):
        lines += [f'  · exact (RatPoly.checkEq_sound _ _ {ns}.low{i}_checked x).le.trans',
                  f'      ((RatPoly.checkLE_sound _ _ {ns}.low{i}_le_checked h0 h1).trans_eq hs.symm)']
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--weight', type=int)
    args = parser.parse_args()
    weights = list(range(129)) if args.weight is None else [args.weight]
    model = json.loads((ROOT/'scripts/sparse_data/model.json').read_text())
    lines = ['import SpinCodes.Structured.SparseSelection']
    lines += [f'import SpinCodes.Structured.SparseBridge.Dominance{j}' for j in weights]
    lines += ['noncomputable section', 'namespace Spin.Structured.SparsePolynomial',
              'set_option maxRecDepth 100000', 'set_option maxHeartbeats 0']
    for j in weights:
        lines += moment_proof(j)
        if j in [1,2]:
            lines += low_proof(j, model)
    lines.append('end Spin.Structured.SparsePolynomial')
    file = 'SparseMaxima.lean' if args.weight is None else f'SparseBridge/Maxima{args.weight}.lean'
    (ROOT/'SpinCodes/Structured'/file).write_text('\n'.join(lines)+'\n',encoding='utf-8')
    if args.weight is None:
        data = ['import SpinCodes.Structured.SparseModelData',
                'namespace Spin.Structured.SparsePolynomial.Data', 'set_option maxRecDepth 100000',
                'def weights : List WeightData := [' + ', '.join(f'weight{j}' for j in weights) + ']',
                'def weight (j : Fin 129) : WeightData := weights.getD j weight0',
                'theorem weights_length : weights.length = 129 := by decide',
                'end Spin.Structured.SparsePolynomial.Data']
        (ROOT/'SpinCodes/Structured/SparseWeights.lean').write_text('\n'.join(data)+'\n',encoding='utf-8')
        all_proofs = ['import SpinCodes.Structured.SparseMaxima',
                      'import SpinCodes.Structured.SparseWeights',
                      'import SpinCodes.Structured.SparseDataChecks',
                      'noncomputable section', 'namespace Spin.Structured.SparsePolynomial',
                      'set_option maxRecDepth 100000', 'set_option maxHeartbeats 0',
                      'theorem maximum_all (x : ℝ) (h0 : 0 ≤ x) (h1 : x ≤ 1) (j : Fin 129) :',
                      '    Spin.Imt.Occupation.Sparse.maximumMoment j (1 - x / 6250) ≤',
                      '      (arbitrary j (Data.weight j)).eval x := by', '  fin_cases j']
        all_proofs += [f'  · exact maximum_{j} x h0 h1' for j in weights]
        all_proofs += ['theorem lowMaximum_all (x : ℝ) (h0 : 0 ≤ x) (h1 : x ≤ 1) (j : Fin 129) :',
                       '    lowMaximum j (Data.weight j) (1 - x / 6250) ≤',
                       '      (selectedLowPolynomial j (Data.weight j)).eval x := by', '  fin_cases j']
        for j in weights:
            if j in [1,2]:
                all_proofs += [f'  · exact lowMaximum_{j} x h0 h1']
            else:
                all_proofs += ['  · change (0 : ℝ) ≤ zero.eval x', '    norm_num [zero, RatPoly.eval_constant]']
        all_proofs += ['theorem data_all_valid (j : Fin 129) : checkData j (Data.weight j) = true := by',
                       '  fin_cases j']
        all_proofs += [f'  · exact Data.weight{j}_valid' for j in weights]
        all_proofs += ['end Spin.Structured.SparsePolynomial',
                       'namespace Spin.Imt.Occupation.Sparse',
                       'def numericalMatrix (β z : ℝ) : Transfer 5 :=',
                       '  matrix 128 count 524287 β (fun j => row j (Spin.Structured.SparsePolynomial.Data.weight j) z)',
                       'end Spin.Imt.Occupation.Sparse']
        (ROOT/'SpinCodes/Structured/SparseMaximumAll.lean').write_text('\n'.join(all_proofs)+'\n',encoding='utf-8')


if __name__ == '__main__':
    main()

