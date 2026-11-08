import SpinCodes.Turbo.RS8ShuffledPlacement
import SpinCodes.Turbo.RS8SpectrumCount

/-!
# A complete-atom envelope for the revised RS8 outer

The finite tail transfer compares expected message counts at each complete
256-bit group word with a uniform-value reference. Before using that
comparison, we bound the number of RS8 messages having any one specified
column support. This uses the already checked shortening equivalence,
without importing numerical receipts or assuming an RS16 distribution.
-/

noncomputable section

namespace Spin.Turbo

open Finset

/-- A nonzero message with a fixed support lies in the shortened code
whose columns outside that support are zero. -/
theorem rs8_exact_support_le_shortened
    {F : Type*} [Field F] [Fintype F] [DecidableEq F]
    (points : Fin 8 → F) (s : Finset (Fin 8)) :
    (Finset.univ.filter (fun rows : RS8Group F =>
      rs8Support points rows = s)).card ≤
      Fintype.card {rows : RS8Group F //
        rs8ZeroOnSet points sᶜ rows} := by
  classical
  have hsub :
      (Finset.univ.filter (fun rows : RS8Group F =>
        rs8Support points rows = s)) ⊆
      (Finset.univ.filter (fun rows : RS8Group F =>
        rs8ZeroOnSet points sᶜ rows)) := by
    intro rows hr
    have hs := (Finset.mem_filter.mp hr).2
    apply Finset.mem_filter.mpr
    refine ⟨Finset.mem_univ _, ?_⟩
    exact (rs8SupportedIn_iff_zero_compl points s rows).mp
      (by intro j hj; rw [hs] at hj; exact hj)
  exact (Finset.card_le_card hsub).trans (by
    rw [Fintype.card_subtype])

/-- If `h≥4` columns may be nonzero, shortening leaves only `h-4`
free 32-bit columns. The exact-support fiber is no larger. -/
theorem rs8_exact_support_count_le
    (points : Fin 8 → Field16) (hpoints : Function.Injective points)
    (s : Finset (Fin 8)) (hs : 4 ≤ s.card) :
    (Finset.univ.filter (fun rows : RS8Group Field16 =>
      rs8Support points rows = s)).card ≤
        (2 ^ 32) ^ (s.card - 4) := by
  classical
  have hz : sᶜ.card ≤ 4 := by
    have h := Finset.card_le_card (Finset.subset_univ s)
    simp only [Finset.card_univ, Fintype.card_fin] at h
    rw [Finset.card_compl]
    simp only [Fintype.card_fin]
    omega
  calc
    _ ≤ Fintype.card {rows : RS8Group Field16 //
        rs8ZeroOnSet points sᶜ rows} :=
      rs8_exact_support_le_shortened points s
    _ = (Fintype.card Field16) ^ (8 * (4 - sᶜ.card)) :=
      rs8ZeroOnSet_card points hpoints sᶜ hz
    _ = (2 ^ 32) ^ (s.card - 4) := by
      have hcard : sᶜ.card = 8 - s.card := by
        rw [Finset.card_compl]
        simp
      rw [hcard, field16_card]
      have he : 8 * (4 - (8 - s.card)) = 8 * (s.card - 4) := by
        have hle := Finset.card_le_card (Finset.subset_univ s)
        simp at hle
        omega
      rw [he]
      norm_num
      conv_lhs => rw [pow_mul]
      norm_num

def packetColumnWordSupport22 (y : Fin 8 → PacketColumn) :
    Finset (Fin 8) :=
  Finset.univ.filter (fun j => y j ≠ 0)

/-- The complete word atom has one factor `1/(2^32-1)` per active
column, and vanishes unless the message and target column supports match. -/
theorem sampledRS8PacketColumnsLaw_apply_support
    (points : Fin 8 → Field16) (rows : RS8Group Field16)
    (y : Fin 8 → PacketColumn) :
    (sampledRS8PacketColumnsLaw points rows).p y =
      if rs8Support points rows = packetColumnWordSupport22 y then
        (1 / (2 ^ 32 - 1 : ℝ)) ^ (rs8Support points rows).card
      else 0 := by
  classical
  rw [sampledRS8PacketColumnsLaw_apply]
  by_cases hs : rs8Support points rows = packetColumnWordSupport22 y
  · rw [if_pos hs]
    have hterm (j : Fin 8) :
        (if j ∈ rs8Support points rows then
          (if y j = 0 then (0 : ℝ) else 1 / (2 ^ 32 - 1 : ℝ))
        else (if y j = 0 then 1 else 0)) =
          if j ∈ rs8Support points rows then
            1 / (2 ^ 32 - 1 : ℝ) else 1 := by
      have hj : (j ∈ rs8Support points rows) ↔ y j ≠ 0 := by
        rw [hs]
        simp [packetColumnWordSupport22]
      by_cases hsrc : j ∈ rs8Support points rows
      · simp [hsrc, hj.mp hsrc]
      · have hy : y j = 0 := by
          by_contra hn
          exact hsrc (hj.mpr hn)
        simp [hsrc, hy]
    simp_rw [hterm]
    rw [Finset.prod_ite_mem_eq, Finset.prod_const]
  · rw [if_neg hs]
    have hdiff : ∃ j : Fin 8,
        ¬ (j ∈ rs8Support points rows ↔
            j ∈ packetColumnWordSupport22 y) := by
      by_contra hn
      apply hs
      ext j
      exact not_exists.mp hn j |> not_not.mp
    obtain ⟨j, hj⟩ := hdiff
    apply Finset.prod_eq_zero (Finset.mem_univ j)
    by_cases hsrc : j ∈ rs8Support points rows
    · have ht : j ∉ packetColumnWordSupport22 y := by
        intro ht
        exact hj ⟨fun _ => ht, fun _ => hsrc⟩
      have hy : y j = 0 := by
        by_contra hn
        exact ht (Finset.mem_filter.mpr ⟨Finset.mem_univ _, hn⟩)
      simp [hsrc, hy]
    · have ht : j ∈ packetColumnWordSupport22 y := by
        by_contra hnot
        exact hj ⟨fun h => (hsrc h).elim, fun h => (hnot h).elim⟩
      have hy : y j ≠ 0 := (Finset.mem_filter.mp ht).2
      simp [hsrc, hy]

/-- The expected complete-word count is a support-fiber count divided
by the number of possible nonzero randomized column values. This is a
sum over messages under one sampled setup, not an independent-message
probability law. -/
theorem rs8_all_messages_atom_sum
    (points : Fin 8 → Field16) (y : Fin 8 → PacketColumn) :
    (∑ rows : RS8Group Field16,
      (sampledRS8PacketColumnsLaw points rows).p y) =
      ((Finset.univ.filter (fun rows : RS8Group Field16 =>
        rs8Support points rows = packetColumnWordSupport22 y)).card : ℝ) /
        (2 ^ 32 - 1 : ℝ) ^ (packetColumnWordSupport22 y).card := by
  classical
  simp_rw [sampledRS8PacketColumnsLaw_apply_support]
  rw [← Finset.sum_filter]
  have hterm (rows : RS8Group Field16)
      (hr : rows ∈ Finset.univ.filter (fun rows : RS8Group Field16 =>
        rs8Support points rows = packetColumnWordSupport22 y)) :
      (1 / (2 ^ 32 - 1 : ℝ)) ^ (rs8Support points rows).card =
        1 / (2 ^ 32 - 1 : ℝ) ^ (packetColumnWordSupport22 y).card := by
    rw [(Finset.mem_filter.mp hr).2]
    rw [one_div_pow]
  rw [Finset.sum_congr rfl hterm, Finset.sum_const, nsmul_eq_mul]
  ring

/-- Expected count over the active group messages at one complete
256-bit atom, before the independent row shuffle. -/
def rs8NonzeroAtomCount (points : Fin 8 → Field16)
    (y : Fin 8 → PacketColumn) : ℝ :=
  ∑ rows ∈ Spin.nonzeroMsgs (RS8Group Field16),
    (sampledRS8PacketColumnsLaw points rows).p y

theorem rs8NonzeroAtomCount_eq_support_count
    (points : Fin 8 → Field16) (y : Fin 8 → PacketColumn)
    (hy : 1 ≤ (packetColumnWordSupport22 y).card) :
    rs8NonzeroAtomCount points y =
      ((Finset.univ.filter (fun rows : RS8Group Field16 =>
        rs8Support points rows = packetColumnWordSupport22 y)).card : ℝ) /
        (2 ^ 32 - 1 : ℝ) ^ (packetColumnWordSupport22 y).card := by
  classical
  have hz : rs8Support points (0 : RS8Group Field16) = ∅ := by
    have hcol (j : Fin 8) :
        rs8Column points (0 : RS8Group Field16) j = 0 := by
      funext r
      change (0 : Polynomial Field16).eval (points j) = 0
      simp
    ext j
    simp [rs8Support, hcol]
  have hzero :
      (sampledRS8PacketColumnsLaw points (0 : RS8Group Field16)).p y = 0 := by
    rw [sampledRS8PacketColumnsLaw_apply_support, hz]
    have hne : packetColumnWordSupport22 y ≠ ∅ := by
      intro he
      rw [he] at hy
      simp at hy
    rw [if_neg (Ne.symm hne)]
  calc
    rs8NonzeroAtomCount points y =
        ∑ rows : RS8Group Field16,
          (sampledRS8PacketColumnsLaw points rows).p y := by
      unfold rs8NonzeroAtomCount Spin.nonzeroMsgs
      rw [Finset.sum_filter]
      apply Finset.sum_congr rfl
      intro rows _
      by_cases hr : rows = 0
      · subst rows
        simp [hzero]
      · simp [hr]
    _ = _ := rs8_all_messages_atom_sum points y

/-- The actual expected atom is bounded using the checked shortened-code
count. The numerator is the number of free 32-bit columns. -/
theorem rs8NonzeroAtomCount_le_shortened
    (points : Fin 8 → Field16) (hpoints : Function.Injective points)
    (y : Fin 8 → PacketColumn)
    (hy : 5 ≤ (packetColumnWordSupport22 y).card) :
    rs8NonzeroAtomCount points y ≤
      ((2 ^ 32 : ℕ) ^ ((packetColumnWordSupport22 y).card - 4) : ℝ) /
        (2 ^ 32 - 1 : ℝ) ^ (packetColumnWordSupport22 y).card := by
  rw [rs8NonzeroAtomCount_eq_support_count points y (by omega)]
  apply div_le_div_of_nonneg_right
    (by exact_mod_cast (rs8_exact_support_count_le points hpoints
      (packetColumnWordSupport22 y) (by omega)))
  positivity

theorem rs8NonzeroAtomCount_zero_below_five
    (points : Fin 8 → Field16) (hpoints : Function.Injective points)
    (y : Fin 8 → PacketColumn)
    (hy : (packetColumnWordSupport22 y).card < 5) :
    rs8NonzeroAtomCount points y = 0 := by
  classical
  unfold rs8NonzeroAtomCount
  apply Finset.sum_eq_zero
  intro rows hr
  rw [sampledRS8PacketColumnsLaw_apply_support]
  rw [if_neg]
  intro he
  have hrow : rows ≠ 0 := (Finset.mem_filter.mp hr).2
  have hfloor := rs8Support_ge_five points hpoints rows hrow
  have hc := congrArg Finset.card he
  omega

/-- The simpler shortening bound still improves the retained RS16
uniform-input factor for every possible nonzero support size. -/
theorem rs8_shortened_atom_le_old_factor
    (h : ℕ) (h5 : 5 ≤ h) (h8 : h ≤ 8) :
    ((2 ^ 32 : ℕ) ^ (h - 4) : ℝ) /
        (2 ^ 32 - 1 : ℝ) ^ h ≤
      1 / (65535 : ℝ) ^ 8 := by
  interval_cases h <;> norm_num

/-- Pointwise expected-message domination by the retained uniform
reference atom bound. This is before the row permutation and regional
route, so those independent pushforwards can use the same factor. -/
theorem rs8NonzeroAtomCount_le_old_uniform_atom
    (points : Fin 8 → Field16) (hpoints : Function.Injective points)
    (y : Fin 8 → PacketColumn) :
    rs8NonzeroAtomCount points y ≤ 1 / (65535 : ℝ) ^ 8 := by
  let h := (packetColumnWordSupport22 y).card
  have h8 : h ≤ 8 := by
    have hc := Finset.card_le_card (Finset.filter_subset
      (fun j : Fin 8 => y j ≠ 0) Finset.univ)
    simpa [h, packetColumnWordSupport22] using hc
  by_cases h5 : 5 ≤ h
  · exact (rs8NonzeroAtomCount_le_shortened points hpoints y h5).trans
      (rs8_shortened_atom_le_old_factor h h5 h8)
  · have hz : rs8NonzeroAtomCount points y = 0 :=
      rs8NonzeroAtomCount_zero_below_five points hpoints y (by omega)
    rw [hz]
    positivity

/-- A fixed row permutation acts bijectively on complete packet words. -/
def packetWordShuffleEquiv (σ : Equiv.Perm (Fin 64)) :
    PacketWord64 ≃ PacketWord64 where
  toFun := shufflePacketWord σ
  invFun := shufflePacketWord σ.symm
  left_inv := by intro x; funext i; simp [shufflePacketWord]
  right_inv := by intro x; funext i; simp [shufflePacketWord]

/-- Reordering the two independent samples exposes a uniform mixture of
bijective pushforwards of the original complete-word law. -/
theorem shuffledPacketLaw_eq_perm_mixture (P : FinPMF PacketWord64) :
    shuffledPacketLaw P =
      (FinPMF.uniform (Equiv.Perm (Fin 64))).bind
        (fun σ => P.map (packetWordShuffleEquiv σ)) := by
  apply FinPMF.eq_of_expect_eq
  intro f
  rw [shuffledPacketLaw, FinPMF.expect_map, FinPMF.expect_prod,
    FinPMF.expect_bind]
  simp only [FinPMF.expect_map]
  unfold FinPMF.expect
  change (∑ x : PacketWord64, P.p x *
      (∑ σ : Equiv.Perm (Fin 64),
        (FinPMF.uniform (Equiv.Perm (Fin 64))).p σ *
          f (shufflePacketWord σ x))) =
    ∑ σ : Equiv.Perm (Fin 64),
      (FinPMF.uniform (Equiv.Perm (Fin 64))).p σ *
        (∑ x : PacketWord64, P.p x * f (shufflePacketWord σ x))
  calc
    _ = ∑ x : PacketWord64, ∑ σ : Equiv.Perm (Fin 64),
          P.p x * ((FinPMF.uniform (Equiv.Perm (Fin 64))).p σ *
            f (shufflePacketWord σ x)) := by simp only [Finset.mul_sum]
    _ = ∑ σ : Equiv.Perm (Fin 64), ∑ x : PacketWord64,
          P.p x * ((FinPMF.uniform (Equiv.Perm (Fin 64))).p σ *
            f (shufflePacketWord σ x)) := Finset.sum_comm
    _ = _ := by
      apply Finset.sum_congr rfl
      intro σ _
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro x _
      ring

/-- The shuffled atom is the average of source atoms at the inverse
images of the sampled permutations. -/
theorem shuffledPacketLaw_apply_inverse (P : FinPMF PacketWord64)
    (y : PacketWord64) :
    (shuffledPacketLaw P).p y =
      (FinPMF.uniform (Equiv.Perm (Fin 64))).expect
        (fun σ => P.p (shufflePacketWord σ.symm y)) := by
  rw [shuffledPacketLaw_eq_perm_mixture, FinPMF.bind_p]
  simp_rw [FinPMF.map_equiv_apply]
  rfl

/-- The pre-shuffle packet-word count is the column-word count through
the checked flattening equivalence. -/
theorem rs8_raw_word_count_eq_atom_count
    (points : Fin 8 → Field16) (y : PacketWord64) :
    (∑ rows ∈ Spin.nonzeroMsgs (RS8Group Field16),
      (sampledRS8RowLaw points rows).p y) =
        rs8NonzeroAtomCount points (packetColumnsWordEquiv.symm y) := by
  unfold rs8NonzeroAtomCount
  apply Finset.sum_congr rfl
  intro rows _
  rw [sampledRS8RowLaw, FinPMF.map_equiv_apply]

/-- The pointwise expected-count bound survives the actual sampled row
permutation, while retaining all packet values. -/
theorem rs8ExpectedShuffledWordCount_le_old_uniform_atom
    (points : Fin 8 → Field16) (hpoints : Function.Injective points)
    (y : PacketWord64) :
    rs8ExpectedShuffledWordCount points y ≤
      1 / (65535 : ℝ) ^ 8 := by
  unfold rs8ExpectedShuffledWordCount
  simp_rw [shuffledPacketLaw_apply_inverse]
  rw [← FinPMF.expect_sum]
  simp_rw [rs8_raw_word_count_eq_atom_count]
  let U := FinPMF.uniform (Equiv.Perm (Fin 64))
  calc
    U.expect (fun σ => rs8NonzeroAtomCount points
        (packetColumnsWordEquiv.symm (shufflePacketWord σ.symm y))) ≤
        U.expect (fun _ => 1 / (65535 : ℝ) ^ 8) :=
      FinPMF.expect_mono U (fun σ =>
        rs8NonzeroAtomCount_le_old_uniform_atom points hpoints
          (packetColumnsWordEquiv.symm (shufflePacketWord σ.symm y)))
    _ = _ := FinPMF.expect_const U _

#print axioms rs8_exact_support_le_shortened
#print axioms rs8_exact_support_count_le
#print axioms sampledRS8PacketColumnsLaw_apply_support
#print axioms rs8_all_messages_atom_sum
#print axioms rs8NonzeroAtomCount_eq_support_count
#print axioms rs8NonzeroAtomCount_le_shortened
#print axioms rs8NonzeroAtomCount_le_old_uniform_atom
#print axioms rs8ExpectedShuffledWordCount_le_old_uniform_atom

end Spin.Turbo
