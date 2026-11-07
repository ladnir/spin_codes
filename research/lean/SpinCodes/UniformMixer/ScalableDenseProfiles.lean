import SpinCodes.UniformMixer.ScalableDenseRateDefs
import SpinCodes.UniformMixer.ScalableFixedLargeProfileCount
import SpinCodes.Structured.DenseOccupationFamilyRateProfiles

/-! Apply a dense point-rate certificate to every selected outer profile. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteRoute Spin.Structured.ConcreteNativeFamily
  Spin.Structured.DenseOccupationFixed Finset

theorem ScalableDenseRates.profile_bound {η C : ℝ}
    (h : ScalableDenseRates η C) (hC : 0 ≤ C)
    (m : ℕ) (seed : NativeSeed m) (hg : nativeGood m seed)
    (S : Finset (Fin (Lsched m))) (hS : S.Nonempty)
    (hα : (S.card : ℝ) / Lsched m ∈ Set.Icc (1 / 10000) 1)
    (w : WeightProfile S (nativeBlocks m * 24))
    (d : ℕ)
    (hd : (d : ℝ) ≤ (11 / 100) *
      ((Lsched m : ℝ) * bsched m)) :
    (∑ x ∈ profileMessages seed S (profileWeights S w),
      (law (fun i => (rowSupports seed x i).map
        (finCongr (native_width m)).toEmbedding)).expect
        (scalableRegionFailure m d)) ≤
      outerCost (bsched m) ^ S.card *
        (((bsched m : ℝ) + 1) ^ S.card *
          ((Lsched m : ℝ) + 1) ^ bsched m) * C *
          Real.exp (-η * ((Lsched m : ℝ) * bsched m)) := by
  classical
  by_cases hempty : (profileMessages seed S (profileWeights S w)).Nonempty
  · obtain ⟨message, hm⟩ := hempty
    have hact : message ∈ activeMessages S := by
      have hf : message ∈ (activeMessages S).filter
          (fun x => messageProfile seed S x = w) := by
        rw [profile_fiber_eq]
        exact hm
      exact (Finset.mem_filter.mp hf).1
    have hk : 0 < nativeBlocks m := by
      have hb := bsched_pos m
      rw [← native_width m] at hb
      omega
    have hs := selected_active_rows hk seed hg S hact
    obtain ⟨x, hx, hmean, hden, hactive, hw0, hw1⟩ :=
      selected_message_parameters (Lsched_pos m) hk seed hg S hS hact
    obtain ⟨p, hp, hrate⟩ := h ((S.card : ℝ) / Lsched m) x hα hx
    let rows := rowSupports seed message
    let rowsB : Fin (Lsched m) → Finset (Fin (bsched m)) :=
      fun i => (rows i).map (finCongr (native_width m)).toEmbedding
    have hw0B : 0 < totalWeight rowsB := by
      simpa only [rowsB, rows, totalWeight, Finset.card_map,
        native_width] using hw0
    have hw1B : totalWeight rowsB < Lsched m * bsched m := by
      simpa only [rowsB, rows, totalWeight, Finset.card_map,
        native_width] using hw1
    have hdenB : density rowsB = ((S.card : ℝ) / Lsched m) * x := by
      simpa only [rowsB, rows, density, totalWeight, Finset.card_map,
        native_width] using hden
    have hweights : ∀ i ∈ S,
        (rowSupports seed message i).card = profileWeights S w i := by
      intro i hi
      have hh := (Finset.mem_filter.mp hm).2 i
      simp only [hi, ite_true] at hh
      simpa only [rowSupports, support_card] using hh.2
    rw [← profileMessages_congr seed S hweights,
      scalable_profile_failure_eq m seed S _ hs.1 d]
    have hpcount := selected_profile_support hk seed
      (weightWindow (nativeBlocks m * 24)) hg S
      (fun i => (rowSupports seed message i).card)
      hs.2.1 hs.2.2 hmean hp
    have hL : (Lsched m : ℝ) ≠ 0 := by
      exact_mod_cast (Lsched_pos m).ne'
    have he : ((nativeBlocks m * 24 : ℕ) : ℝ) * (S.card : ℝ) *
        ((p.1 : ℝ) * x + p.2) =
        ((Lsched m : ℝ) * bsched m) *
          ((S.card : ℝ) / Lsched m) *
          ((p.1 : ℝ) * x + p.2) := by
      rw [native_width]
      field_simp
    rw [he] at hpcount
    have hh := hrate.scaled m d rowsB hw0B hw1B hdenB hd
      (pow_pos (outerCost_pos (bsched_pos m)) _)
      (Finset.prod_nonneg (fun _ _ => Nat.cast_nonneg _))
      (by simpa only [native_width] using hpcount)
    have hactiveB : activeRows rowsB = S.card := by
      have heq : activeRows rowsB = activeRows rows := by
        simp [activeRows, rowsB, Finset.map_eq_empty]
      exact heq.trans (by simpa only [rows] using hactive)
    rw [scalable_route_probability_eq]
    simpa only [rowsB, rows, hactiveB, native_width] using hh
  · rw [Finset.not_nonempty_iff_eq_empty.mp hempty, Finset.sum_empty]
    have hroute : 0 ≤ ((bsched m : ℝ) + 1) ^ S.card *
        ((Lsched m : ℝ) + 1) ^ bsched m := by positivity
    have hout : 0 ≤ outerCost (bsched m) ^ S.card :=
      (pow_pos (outerCost_pos (bsched_pos m)) _).le
    positivity

#print axioms ScalableDenseRates.profile_bound

end Spin.UniformMixer
