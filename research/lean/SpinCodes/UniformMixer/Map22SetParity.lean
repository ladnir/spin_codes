import Mathlib.Data.Finset.Interval
import Mathlib.Data.ZMod.Basic
import Mathlib.Tactic.NormNum

/-!
# Parity of a Boolean-lattice interval

The submasks between two six-bit masks form an interval in the lattice
of subsets. Its size is a power of two, so its parity is one exactly
when the endpoints coincide. This lemma is independent of the SPIN map.
-/

namespace Spin.UniformMixer

theorem finset_interval_parity {α : Type*} [DecidableEq α]
    (k j : Finset α) :
    (∑ _m ∈ Finset.Icc k j, (1 : ZMod 2)) =
      if k = j then 1 else 0 := by
  by_cases hsub : k ⊆ j
  · by_cases heq : k = j
    · subst j
      simp
    · have hlt : k.card < j.card :=
        Finset.card_lt_card (Finset.ssubset_iff_subset_ne.mpr ⟨hsub, heq⟩)
      have hpos : 0 < j.card - k.card := Nat.sub_pos_of_lt hlt
      cases hdiff : j.card - k.card with
      | zero => omega
      | succ d =>
          rw [Finset.sum_const, nsmul_eq_mul, mul_one,
            Finset.card_Icc_finset hsub, hdiff]
          have htwo : (2 : ZMod 2) = 0 := by decide
          simp [heq, pow_succ, htwo]
  · have hempty : Finset.Icc k j = ∅ := by
      ext m
      constructor
      · intro hm
        exact (hsub ((Finset.mem_Icc.mp hm).1.trans
          (Finset.mem_Icc.mp hm).2)).elim
      · intro hm
        have hf : False := by simp at hm
        exact hf.elim
    have hne : k ≠ j := by
      intro heq
      exact hsub heq.le
    simp [hempty, hne]

end Spin.UniformMixer
