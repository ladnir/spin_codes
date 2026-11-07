import SpinCodes.UniformMixer.QuadraticRankBlocks

/-!
# The finite rank histogram counts actual quadratic forms

The range checker enumerates exactly the fifteen-bit coefficient vectors.
This theorem is independent of the values of the histogram entries.
-/

namespace Spin.UniformMixer

theorem length_filter_range_eq_sum (p : Nat → Bool) (n : Nat) :
    ((List.range n).filter p).length =
      ∑ j ∈ Finset.range n, if p j then 1 else 0 := by
  induction n with
  | zero => simp
  | succ n ih =>
      simp only [List.range_succ, List.filter_append, List.length_append,
        List.filter_singleton, Finset.sum_range_succ]
      split_ifs <;> simp_all

theorem packedRadicalCount_eq_fastRange (radicalSize : Nat) :
    packedRadicalCount radicalSize =
      fastRadicalRangeCount 0 (2 ^ 15) radicalSize := by
  classical
  unfold packedRadicalCount
  rw [Finset.card_filter]
  simp_rw [packedRadicalSize_eq_fast]
  rw [Fin.sum_univ_eq_sum_range
    (fun j => if fastRadicalSize j = radicalSize then (1 : Nat) else 0)]
  have hlen := length_filter_range_eq_sum
    (fun j => decide (fastRadicalSize j = radicalSize)) (2 ^ 15)
  simp only [decide_eq_true_eq] at hlen
  rw [← hlen]
  rw [fastRadicalRangeCount_eq_block]
  simp [fastRadicalBlockCount]

end Spin.UniformMixer
