import SpinCodes.UniformMixer.Mean

/-!
# The paper's 22-state expansion map

`rowMask` is the 22-row, 64-column binary matrix printed in the paper and
in `spin/src/certified22/productionMap.h`. Bit `j` of row `i` is the
coefficient of state coordinate `i` in output coordinate `j`.
-/

noncomputable section

namespace Spin.UniformMixer

def rowMask : Fin 22 → ℕ := ![
  0xffffffffffffffff, 0xaaaaaaaaaaaaaaaa, 0xcccccccccccccccc,
  0xf0f0f0f0f0f0f0f0, 0xff00ff00ff00ff00, 0xffff0000ffff0000,
  0xffffffff00000000, 0x3693af0a9cc6faa0, 0xd882d88272287228,
  0x666669960ff00000, 0xb874edded1e28448, 0xa5ffcc96cc965a00,
  0x8e2b4de88ed4b2e8, 0x2de14b78ddee4488, 0xcca5cc5a69009600,
  0x9556599a300cfcc0, 0x8888888888888888, 0xa0a0a0a0a0a0a0a0,
  0xaa00aa00aa00aa00, 0xaaaa0000aaaa0000, 0xaaaaaaaa00000000,
  0xcccc0000cccc0000]

def coordCoeff (j : Fin 64) (i : Fin 22) : ZMod 2 :=
  if (rowMask i).testBit j.val then 1 else 0

def coordMap (j : Fin 64) : State22 →ₗ[ZMod 2] ZMod 2 where
  toFun q := ∑ i : Fin 22, coordCoeff j i * q i
  map_add' q r := by
    simp [Pi.add_apply, mul_add, Finset.sum_add_distrib]
  map_smul' a q := by
    simp only [Pi.smul_apply, smul_eq_mul, Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro i _
    simp only [RingHom.id_apply]
    ring

theorem every_output_covered :
    ∀ j : Fin 64, ∃ i : Fin 22, (rowMask i).testBit j.val = true := by
  decide

def basisVector (i : Fin 22) : State22 :=
  fun k => if k = i then 1 else 0

theorem coordMap_basis (j : Fin 64) (i : Fin 22) :
    coordMap j (basisVector i) = coordCoeff j i := by
  simp [coordMap, basisVector, Finset.sum_ite_eq']

theorem coordMap_nonzero (j : Fin 64) :
    ∃ u : State22, coordMap j u = 1 := by
  obtain ⟨i, hi⟩ := every_output_covered j
  refine ⟨basisVector i, ?_⟩
  rw [coordMap_basis]
  simp [coordCoeff, hi]

theorem production_liveWeightAverage :
    liveWeightAverage coordMap = 64 * ((2 ^ 21 : ℝ) / (2 ^ 22 - 1)) :=
  state22_liveWeightAverage_formula coordMap coordMap_nonzero

end Spin.UniformMixer
