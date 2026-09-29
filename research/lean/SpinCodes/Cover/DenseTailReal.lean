/-
`eq:structured-ba-dense-tail`, stated over the paper's own intervals and
constant.

    g(α) + π(α,β) + π(β,x) < -7.68·10⁻⁸

at every feasible point with `α ∈ [1/125, 1]`, `β ∈ [0, 1]`,
`x ∈ [1/500, 13/125]`.

Hand-written, so that regenerating the parts cannot clobber it.  The proof is
arithmetic only: every endpoint is exact at `scale = 10^30`, since
`10^30/125 = 8·10^27`, `10^30/500 = 2·10^27` and `13·10^30/125 = 104·10^27`.
No second enclosure is involved.
-/
import SpinCodes.Cover.DenseTail

namespace Spin.Cover

open Spin.Numeric Spin.Numeric.Fix

private lemma scaleR_eq : ((scale : Int) : ℝ) = 10 ^ 30 := by
  unfold scale; norm_num

/-- **The dense-tail gap, over the paper's intervals.** -/
theorem denseTail_real {a b w : ℝ}
    (ha : a ∈ Set.Icc ((1 : ℝ) / 125) 1) (hb : b ∈ Set.Icc (0 : ℝ) 1)
    (hw : w ∈ Set.Icc ((1 : ℝ) / 500) (13 / 125))
    (hf : a / 2 < b ∧ b < 1 - a / 2 ∧ b / 2 < w ∧ w < 1 - b / 2
      ∧ 0 ≤ a ∧ a < 1) :
    gBA a + piBA a b + piBA b w < -(768 / 10 ^ 10) := by
  have hA : Mem rootA a := by
    refine ⟨?_, ?_⟩
    · have he : ((rootA.lo : Int) : ℝ) / ((scale : Int) : ℝ) = 1 / 125 := by
        rw [scaleR_eq, rootA]; norm_num
      rw [he]; exact ha.1
    · have he : ((rootA.hi : Int) : ℝ) / ((scale : Int) : ℝ) = 1 := by
        rw [scaleR_eq, rootA]; norm_num
      rw [he]; exact ha.2
  have hB : Mem rootB b := by
    refine ⟨?_, ?_⟩
    · have he : ((rootB.lo : Int) : ℝ) / ((scale : Int) : ℝ) = 0 := by
        rw [scaleR_eq, rootB]; norm_num
      rw [he]; exact hb.1
    · have he : ((rootB.hi : Int) : ℝ) / ((scale : Int) : ℝ) = 1 := by
        rw [scaleR_eq, rootB]; norm_num
      rw [he]; exact hb.2
  have hW : Mem rootW w := by
    refine ⟨?_, ?_⟩
    · have he : ((rootW.lo : Int) : ℝ) / ((scale : Int) : ℝ) = 1 / 500 := by
        rw [scaleR_eq, rootW]; norm_num
      rw [he]; exact hw.1
    · have he : ((rootW.hi : Int) : ℝ) / ((scale : Int) : ℝ) = 13 / 125 := by
        rw [scaleR_eq, rootW]; norm_num
      rw [he]; exact hw.2
  have h := denseTail a b w hA hB hW hf
  have he : ((thr : Int) : ℝ) / ((scale : Int) : ℝ) = -(768 / 10 ^ 10) := by
    rw [scaleR_eq, thr]; norm_num
  rwa [he] at h

end Spin.Cover
