"""Emit the final concrete map/count connection after both numerical replays pass."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'scripts/map_data'


def main():
    for report in ['kernel_blocks.json','kernel_fiber_data.json']:
        assert json.loads((DATA/report).read_text())['status'] == 'PASS', report
    for module in ['MapSpectrumBridge','FiberNumericsAll']:
        assert (ROOT/f'.lake/build/lib/lean/SpinCodes/Structured/{module}.olean').exists(), module
    lines = ['import SpinCodes.Structured.MapSpectrumBridge',
             'import SpinCodes.Structured.FiberNumericsAll',
             'import SpinCodes.Structured.ConcreteShells',
             '', '/-! Identification of the frozen numerical data with the actual table maps. -/',
             'noncomputable section', 'namespace Spin.Structured.ConcreteMaps',
             '', 'theorem transpose_spectrum (w : Fin 129) :',
             '    weightCounts CtransposeSet w = FiberNumerics.Data.spectrum.getD w 0 := by',
             '  rw [CtransposeSet_weightCounts w.isLt]', '  rfl',
             '', 'theorem actual_shell_counts (i : Fin 5) :',
             '    weightCounts Aset (shellWeight i) = shellCount i := by',
             '  fin_cases i']
    for w,count in [(48,5166),(56,110288),(64,293455),(72,110128),(80,5250)]:
        lines += [f'  · change weightCounts Aset {w} = {count}',
                  f'    exact Aset_weightCounts (i := {w}) (by decide)']
    lines += ['', 'def actualShellSystem : Spin.Imt.ShellSystem 19 5 := shellSystem actual_shell_counts',
              '', 'theorem kernel_card_frozen (j : Fin 129) :',
              '    (kernelLayer j).card = (SparsePolynomial.weightN j).kernel :=',
              '  FiberNumerics.kernel_counts_frozen transpose_spectrum j',
              '', 'theorem syndromeFiber_cap_frozen (j : Fin 129) {q : Finset (Fin 19)} (hq : q ≠ ∅) :',
              '    (syndromeFiber j q).card ≤ (SparsePolynomial.weightN j).cap :=',
              '  FiberNumerics.fiber_caps_frozen transpose_spectrum j hq',
              '', 'theorem pair_card_frozen (j : Fin 129) :',
              '    (equalSyndromePairs j).card = FiberNumerics.pairTotals.getD j 0 :=',
              '  FiberNumerics.pair_counts transpose_spectrum j',
              '', 'theorem actualShellSystem_card (i : Fin 5) :',
              '    ((actualShellSystem.shell i).card : ℝ) = Spin.Imt.Occupation.Sparse.count i :=',
              '  shellSystem_card actual_shell_counts i',
              '', 'end Spin.Structured.ConcreteMaps', '']
    (ROOT/'SpinCodes/Structured/ConcreteCounts.lean').write_text('\n'.join(lines),encoding='utf-8')
    print('Emitted ConcreteCounts.lean. Compile and audit it before claiming the concrete count connection.')


if __name__ == '__main__':
    main()
