"""Emit the semantic composition of the checked product sums."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def emit_program():
    records = json.loads((ROOT/'scripts/sparse_data/contribution_manifest.json').read_text())
    chunks = [r['input_weights'] for r in records if r['row'] == 0]
    literal = '[' + ', '.join('[' + ', '.join(map(str, ns)) + ']' for ns in chunks) + ']'
    definitions = ['import SpinCodes.Structured.SparseWeights',
                   'namespace Spin.Structured.SparsePolynomial',
                   'set_option maxRecDepth 100000',
                   'def weightN (j : Nat) : WeightData := Data.weights.getD j Data.weight0',
                   f'def sourceChunks : List (List Nat) := {literal}',
                   'theorem sourceChunks_cover : sourceChunks.flatten = List.range 129 := by decide',
                   'end Spin.Structured.SparsePolynomial']
    (ROOT/'SpinCodes/Structured/SparseProgramDefs.lean').write_text('\n'.join(definitions)+'\n',encoding='utf-8')
    lines = ['import SpinCodes.Structured.SparseProgramDefs',
             'import SpinCodes.Structured.SparseResidualsSound',
             'noncomputable section', 'namespace Spin.Structured.SparsePolynomial',
             'set_option maxRecDepth 100000', 'set_option maxHeartbeats 0',
             'open Residuals',
             'def programResidual (i : ℕ) (x : ℝ) : ℝ :=',
             '  ((List.range 129).map (fun j => (contribution j (weightN j) i).eval x)).sum -',
             '    (1 - 96 * (x / 10000)) * (witness i).eval x',
             'theorem sum_chunks (chunks : List (List ℕ)) (f : ℕ → ℝ) :',
             '    (chunks.map (fun ns => (ns.map f).sum)).sum = (chunks.flatten.map f).sum := by',
             '  induction chunks with', '  | nil => simp',
             '  | cons ns chunks ih => simp [ih, List.sum_append]']
    for row in range(7):
        lines += [f'theorem programResidual_eq_{row} (x : ℝ) : programResidual {row} x = sourceResidual{row} x := by',
                  f'  have h := sum_chunks sourceChunks (fun j => (contribution j (weightN j) {row}).eval x)',
                  '  rw [sourceChunks_cover] at h',
                  '  rw [programResidual, ← h]', '  rfl']
    lines += ['theorem programResidual_neg (i : Fin 7) {x : ℝ} (h0 : 0 < x) (h1 : x ≤ 1) :',
              '    programResidual i x < 0 := by', '  fin_cases i']
    for row in range(7):
        lines += [f'  · rw [programResidual_eq_{row}]', f'    exact sourceResidual{row}_neg h0 h1']
    lines += ['end Spin.Structured.SparsePolynomial']
    (ROOT/'SpinCodes/Structured/SparseProgram.lean').write_text('\n'.join(lines)+'\n',encoding='utf-8')


def main():
    records = json.loads((ROOT/'scripts/sparse_data/contribution_manifest.json').read_text())
    lines = ['import SpinCodes.Structured.PolyPacked']
    lines += [f'import SpinCodes.Structured.SparseBridge.Row{r["row"]}Part{r["part"]}' for r in records]
    lines += ['noncomputable section', 'namespace Spin.Structured.SparsePolynomial',
              'set_option maxRecDepth 100000', 'set_option maxHeartbeats 0']
    for r in records:
        row, part, weights = r['row'], r['part'], r['input_weights']
        ns = f'Row{row}Part{part}'
        lines += [f'namespace {ns}', 'def sourceValue (x : ℝ) : ℝ := [' +
                  ', '.join(f'(contribution {j} Data.weight{j} {row}).eval x' for j in weights) + '].sum',
                  'theorem expected_eval (x : ℝ) : expected.eval x = sourceValue x := by']
        for j in weights:
            lines += [f'  have h{j} : Weight{j}.prob.eval x * Weight{j}.a{row}.eval x =',
                      f'      (contribution {j} Data.weight{j} {row}).eval x := by',
                      '    rw [contribution, RatPoly.eval_mul]',
                      '    exact congrArg₂ (· * ·)',
                      f'      (RatPoly.checkEq_sound _ _ Weight{j}.prob_checked x).symm',
                      f'      (RatPoly.checkEq_sound _ _ Weight{j}.a{row}_checked x).symm']
        lines += ['  rw [checkSumProducts_sound _ _ _ _ checked]',
                  '  simp only [terms, packedValue, sourceValue, List.sum_cons, List.sum_nil]']
        chain = 'rfl'
        for j in reversed(weights):
            chain = f'(congrArg₂ (· + ·) h{j} {chain})'
        lines += [f'  exact {chain}', f'end {ns}']
    lines += ['end Spin.Structured.SparsePolynomial']
    (ROOT/'SpinCodes/Structured/SparseContributionSound.lean').write_text('\n'.join(lines)+'\n',encoding='utf-8')

    model = json.loads((ROOT/'scripts/sparse_data/model.json').read_text())
    lines = ['import SpinCodes.Structured.SparseResiduals', 'import SpinCodes.Structured.SparseContributionSound',
             'import SpinCodes.Structured.PolyPacked',
             'noncomputable section', 'namespace Spin.Structured.SparsePolynomial.Residuals']
    for row in range(7):
        den = model['normalizations'][row]['denominator']
        lines += [f'theorem residual{row}_neg {{x : ℝ}} (h0 : 0 < x) (h1 : x ≤ 1) :',
                  f'    packedValue terms{row} x < 0 := by',
                  f'  rw [← checkSumProducts_sound _ _ _ _ row{row}_checked x]',
                  f'  change listEval (0 :: SparseData.row{row}) x / (({den} : ℕ) : ℝ) < 0',
                  '  rw [listEval, Int.cast_zero, zero_add]',
                  '  apply div_neg_of_neg_of_pos',
                  f'  · exact residual_neg_of_tailBound _ SparseData.row{row}_cert h0 h1',
                  '  · positivity']
        parts = [r['part'] for r in records if r['row'] == row]
        lines += [f'def sourceResidual{row} (x : ℝ) : ℝ :=',
                  '  [' + ', '.join(f'Row{row}Part{p}.sourceValue x' for p in parts) + '].sum -',
                  f'    (1 - 96 * (x / 10000)) * (witness {row}).eval x',
                  f'theorem sourceResidual{row}_neg {{x : ℝ}} (h0 : 0 < x) (h1 : x ≤ 1) :',
                  f'    sourceResidual{row} x < 0 := by',
                  '  have hc : negativeContraction.eval x = -(1 - 96 * (x / 10000)) := by',
                  '    norm_num [negativeContraction, RatPoly.eval, listEval]', '    ring',
                  f'  have hw : vector{row}.eval x = (witness {row}).eval x :=',
                  f'    (RatPoly.checkEq_sound _ _ vector{row}_checked x).symm',
                  f'  have he : packedValue terms{row} x = sourceResidual{row} x := by',
                  f'    simp only [sourceResidual{row}, terms{row}, packedValue, RatPoly.eval_constant,',
                  '      Int.cast_one, Nat.cast_one, div_one, mul_one, List.sum_cons, List.sum_nil,',
                  '      ' + ', '.join(f'Row{row}Part{p}.expected_eval' for p in parts) + ', hc, hw]',
                  '    ring',
                  f'  exact he ▸ residual{row}_neg h0 h1']
    lines += ['end Spin.Structured.SparsePolynomial.Residuals']
    (ROOT/'SpinCodes/Structured/SparseResidualsSound.lean').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    emit_program()


if __name__ == '__main__':
    main()

