import SpinCodes.UniformMixer.ScalableDenseLeaf

/-! The hardest box reported by the external Arb search, replayed in Lean. -/

noncomputable section
set_option maxRecDepth 8192
set_option maxHeartbeats 1000000

namespace Spin.UniformMixer

open Spin.Numeric Spin.Structured.DenseOccupationFixed

private def hardP : QInput := ⟨49999308654867053, 50000000000000000⟩
private def hardY : QInput := ⟨6366695006060831, 12500000000000000⟩
private def hardZ : QInput :=
  ⟨123595507696772821266419713391, 1000000000000000000000000000000⟩
private def hardRadius : QInput :=
  ⟨93998315052091, 1000000000000000000000000000000⟩
private def hardV0 : QInput := ⟨1, 1⟩
private def hardV1 : QInput := ⟨2496174743187, 976562500000⟩
private def hardV2 : QInput := ⟨1013920728007387, 1000000000000000⟩

private def hardVector : Fin 3 → QInput := ![hardV0, hardV1, hardV2]
private def hardVectorFix : Fin 3 → Fix := fun i => (hardVector i).interval
private def hardVectorReal : Fin 3 → ℝ := fun i => (hardVector i).real

private theorem hard_transfer_check :
    transferCheck (Fix.mul hardP.interval hardY.interval)
      hardZ.interval hardRadius.interval hardVectorFix = true := by
  decide +kernel

#print axioms hard_transfer_check

private theorem hard_transfer :
    ∀ i, (paperSpectrumTransfer (hardP.real * hardY.real) hardZ.real).mulVec
      hardVectorReal i ≤ hardRadius.real * hardVectorReal i := by
  apply transferCheck_sound (v := hardVectorFix) (w := hardVectorReal)
    (Fix.mul_mem (hardP.interval_mem (by decide))
      (hardY.interval_mem (by decide)))
    (hardZ.interval_mem (by decide))
    (hardRadius.interval_mem (by decide))
  · norm_num [hardP, hardY, QInput.real]
  · norm_num [hardP, hardY, QInput.real]
  · norm_num [hardZ, QInput.real]
  · intro i
    simpa [hardVectorFix, hardVectorReal] using
      (hardVector i).interval_mem (by fin_cases i <;> decide)
  · exact hard_transfer_check

private def hardSlope : QInput := ⟨-1, 100⟩
private def hardIntercept : QInput :=
  ⟨1495524372286176733882761540124617467907,
    4253529586511730793292182592897102643200⟩
private def hardA0 : QInput := ⟨2097142001, 2097152000⟩
private def hardA1 : QInput := ⟨20971430009, 20971520000⟩
private def hardX0 : QInput :=
  ⟨4283433002144158131081961361172224591779859482217843,
    8507059173023461586584365185794205286400000000000000⟩
private def hardX1 : QInput :=
  ⟨17133734864128980751326149544381562483231436604095831,
    34028236692093846346337460743176821145600000000000000⟩

private theorem hard_vertex00 :
    scalableVertexCheck hardSlope hardIntercept hardP hardY hardRadius hardZ
      hardA0 hardX0 50 (-100000000000000000000000) = true := by
  decide +kernel
private theorem hard_vertex01 :
    scalableVertexCheck hardSlope hardIntercept hardP hardY hardRadius hardZ
      hardA0 hardX1 50 (-100000000000000000000000) = true := by
  decide +kernel
private theorem hard_vertex10 :
    scalableVertexCheck hardSlope hardIntercept hardP hardY hardRadius hardZ
      hardA1 hardX0 50 (-100000000000000000000000) = true := by
  decide +kernel
private theorem hard_vertex11 :
    scalableVertexCheck hardSlope hardIntercept hardP hardY hardRadius hardZ
      hardA1 hardX1 50 (-100000000000000000000000) = true := by
  decide +kernel

private theorem hard_corner (a x : QInput)
    (ha : 0 < a.den) (hx : 0 < x.den)
    (h : scalableVertexCheck hardSlope hardIntercept hardP hardY
      hardRadius hardZ a x 50 (-100000000000000000000000) = true) :
    scalableBoxExponent hardSlope.real hardIntercept.real hardP.real
      hardY.real hardRadius.real hardZ.real a.real x.real ≤
      -(1 / 10000000 : ℝ) := by
  have hh := scalableVertexCheck_sound hardSlope hardIntercept
    hardP hardY hardRadius hardZ a x
    (by decide) (by decide) ha hx
    (by decide) (by decide) (by decide) (by decide) h
  norm_num [Spin.Numeric.scale] at hh
  exact hh

private theorem hard_box {α x : ℝ}
    (hα : α ∈ Set.Icc hardA0.real hardA1.real)
    (hx : x ∈ Set.Icc hardX0.real hardX1.real) :
    scalableBoxExponent hardSlope.real hardIntercept.real hardP.real
      hardY.real hardRadius.real hardZ.real α x ≤
      -(1 / 10000000 : ℝ) := by
  have h := scalableBoxExponent_le_vertices
    hardSlope.real hardIntercept.real hardP.real hardY.real
    hardRadius.real hardZ.real
    (by norm_num [hardP, QInput.real])
    (by norm_num [hardP, QInput.real])
    (by norm_num [hardY, QInput.real])
    (by norm_num [hardY, QInput.real])
    (by norm_num [hardA0, QInput.real])
    (by norm_num [hardA1, QInput.real])
    (by norm_num [hardX0, QInput.real])
    (by norm_num [hardX1, QInput.real])
    hα hx
  apply h.trans
  apply max_le (max_le ?_ ?_) (max_le ?_ ?_)
  · exact hard_corner hardA0 hardX0 (by decide) (by decide) hard_vertex00
  · exact hard_corner hardA0 hardX1 (by decide) (by decide) hard_vertex01
  · exact hard_corner hardA1 hardX0 (by decide) (by decide) hard_vertex10
  · exact hard_corner hardA1 hardX1 (by decide) (by decide) hard_vertex11

private theorem hard_pointRate {α x : ℝ}
    (hα : α ∈ Set.Icc hardA0.real hardA1.real)
    (hx : x ∈ Set.Icc hardX0.real hardX1.real) :
    ScalableDensePointRate α x hardSlope.real hardIntercept.real
      (1 / 10000000) 1 := by
  have hα0 : 0 < α := lt_of_lt_of_le
    (by norm_num [hardA0, QInput.real]) hα.1
  have hα1 : α ≤ 1 := le_trans hα.2
    (by norm_num [hardA1, QInput.real])
  have hx0 : 0 < x := lt_of_lt_of_le
    (by norm_num [hardX0, QInput.real]) hx.1
  have hx1 : x ≤ 1 := le_trans hx.2
    (by norm_num [hardX1, QInput.real])
  have hones : ∀ i, (1 : ℝ) ≤ 1 * hardVectorReal i := by
    intro i
    fin_cases i <;> norm_num [hardVectorReal, hardVector,
      hardV0, hardV1, hardV2, QInput.real]
  have h := scalable_dense_witness_pointRate
    hα0 hα1 hx0 hx1
    (by norm_num [hardP, QInput.real])
    (by norm_num [hardP, QInput.real])
    (by norm_num [hardY, QInput.real])
    (by norm_num [hardY, QInput.real])
    (by norm_num [hardRadius, QInput.real])
    (by norm_num [hardZ, QInput.real])
    (by norm_num [hardZ, QInput.real])
    hardVectorReal (by norm_num : (0 : ℝ) ≤ 1)
    hones hard_transfer (hard_box hα hx)
  simpa [hardVectorReal, hardVector, hardV0, QInput.real] using h

private theorem hard_support :
    ((-1 / 100 : ℚ),
      (1495524372286176733882761540124617467907 /
        4253529586511730793292182592897102643200 : ℚ)) ∈
      Spin.Majorant.refined.supports := by
  exact Spin.Majorant.member20

private theorem hard_pointRate_via_leaf {α x : ℝ}
    (hα : α ∈ Set.Icc hardA0.real hardA1.real)
    (hx : x ∈ Set.Icc hardX0.real hardX1.real) :
    ScalableDensePointRate α x hardSlope.real hardIntercept.real
      (1 / 10000000) 1 := by
  have h := scalable_dense_checked_leaf hardSlope hardIntercept
    hardP hardY hardRadius hardZ hardA0 hardA1 hardX0 hardX1
    hardVector 1
    (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide)
    (by norm_num [hardP, QInput.real])
    (by norm_num [hardP, QInput.real])
    (by norm_num [hardY, QInput.real])
    (by norm_num [hardY, QInput.real])
    (by norm_num [hardRadius, QInput.real])
    (by norm_num [hardZ, QInput.real])
    (by norm_num [hardZ, QInput.real])
    (by norm_num [hardA0, QInput.real])
    (by norm_num [hardA1, QInput.real])
    (by norm_num [hardX0, QInput.real])
    (by norm_num [hardX1, QInput.real])
    (by intro i; fin_cases i <;> decide)
    (by norm_num : (0 : ℝ) ≤ 1)
    (by intro i; fin_cases i <;> norm_num
      [hardVector, hardV0, hardV1, hardV2, QInput.real])
    (by exact hard_transfer_check)
    hard_vertex00 hard_vertex01 hard_vertex10 hard_vertex11 hα hx
  simpa [hardVector, hardV0, QInput.real] using h

/-- The tightest external box reaches the selected-outer rate interface. -/
theorem scalable_dense_hard_pilot_coverage {α x : ℝ}
    (hα : α ∈ Set.Icc hardA0.real hardA1.real)
    (hx : x ∈ Set.Icc hardX0.real hardX1.real) :
    ∃ p : ℚ × ℚ, p ∈ Spin.Majorant.refined.supports ∧
      ScalableDensePointRate α x p.1 p.2 (1 / 10000000) 1 := by
  refine ⟨((-1 / 100 : ℚ),
      (1495524372286176733882761540124617467907 /
        4253529586511730793292182592897102643200 : ℚ)),
    hard_support, ?_⟩
  simpa [hardSlope, hardIntercept, QInput.real] using
    (hard_pointRate_via_leaf hα hx)

#print axioms hard_pointRate
#print axioms hard_pointRate_via_leaf
#print axioms scalable_dense_hard_pilot_coverage

#print axioms hard_transfer
#print axioms hard_vertex00
#print axioms hard_vertex01
#print axioms hard_vertex10
#print axioms hard_vertex11

end Spin.UniformMixer
