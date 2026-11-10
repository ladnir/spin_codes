import SpinCodes.Turbo.Occupation22
import SpinCodes.Turbo.PacketProductMoment

/-!
# Actual rare outer event for the Turbo sparse layers

Audit: `ConcreteFixedRoutedTwo.actual_two_fair_probability` retains two
marks under one regional permutation. `ScalableFixedSelected` conditions
its outer law; here we retain unconditional truncated measures, so no
setup independence is asserted after excluding the rare event.
-/
set_option maxHeartbeats 200000
noncomputable section
namespace Spin.Turbo
open Finset

/-- A sampled constituent contains a nonzero message of weight 5..8. -/
def turboSparseRareGroup22 (points : Fin 8 → Field16)
    (coeff : Fin 8 → {c : Field32 // c ≠ 0}) : Prop :=
  ∃ msg ∈ Spin.nonzeroMsgs (RS8Group Field16),
    ∃ v ∈ Finset.Icc 5 8, rs8SampledPacketWeight points msg coeff = v

noncomputable instance (points : Fin 8 → Field16) :
    DecidablePred (turboSparseRareGroup22 points) :=
  fun _ => Classical.propDecidable _

/-- Charge this event once for the whole sampled outer encoder. -/
def turboSparseRare22 (R : ℕ) (points : Fin 8 → Field16)
    (seed : TurboOuterSeed22 R) : Prop :=
  ∃ i : Fin (R * 16), turboSparseRareGroup22 points (seed i)

noncomputable instance (R : ℕ) (points : Fin 8 → Field16) :
    DecidablePred (turboSparseRare22 R points) :=
  fun _ => Classical.propDecidable _

private theorem finite_double_union_bound
    {Ω α β : Type*} [Fintype Ω] (P : FinPMF Ω)
    (s : Finset α) (t : Finset β) (E : α → β → Ω → Prop)
    [∀ a b, DecidablePred (E a b)]
    [∀ a, DecidablePred (fun ω => ∃ b ∈ t, E a b ω)]
    [DecidablePred (fun ω => ∃ a ∈ s, ∃ b ∈ t, E a b ω)] :
    P.prob (fun ω => ∃ a ∈ s, ∃ b ∈ t, E a b ω) ≤
      ∑ b ∈ t, ∑ a ∈ s, P.prob (E a b) := by
  calc
    _ ≤ ∑ a ∈ s, P.prob (fun ω => ∃ b ∈ t, E a b ω) :=
      P.prob_biUnion_le s (fun a ω => ∃ b ∈ t, E a b ω)
    _ ≤ ∑ a ∈ s, ∑ b ∈ t, P.prob (E a b) := by
      exact Finset.sum_le_sum fun a _ => P.prob_biUnion_le t (E a)
    _ = _ := Finset.sum_comm

theorem turboSparseRareGroup22_prob_le (points : Fin 8 → Field16) :
    (piPMF (fun _ : Fin 8 =>
      FinPMF.uniform {c : Field32 // c ≠ 0})).prob
        (turboSparseRareGroup22 points) ≤
      ∑ v ∈ Finset.Icc 5 8, rs8ExpectedPacketCount points v := by
  letI : ∀ msg : RS8Group Field16, DecidablePred (fun coeff =>
      ∃ v ∈ Finset.Icc 5 8, rs8SampledPacketWeight points msg coeff = v) :=
    fun _ _ => Classical.propDecidable _
  letI : DecidablePred (fun coeff =>
      ∃ msg ∈ Spin.nonzeroMsgs (RS8Group Field16),
        ∃ v ∈ Finset.Icc 5 8, rs8SampledPacketWeight points msg coeff = v) :=
    fun _ => Classical.propDecidable _
  let P := piPMF (fun _ : Fin 8 => FinPMF.uniform {c : Field32 // c ≠ 0})
  have he : P.prob (turboSparseRareGroup22 points) =
      P.prob (fun coeff => ∃ msg ∈ Spin.nonzeroMsgs (RS8Group Field16),
        ∃ v ∈ Finset.Icc 5 8, rs8SampledPacketWeight points msg coeff = v) :=
    P.prob_congr (s := turboSparseRareGroup22 points)
      (t := fun coeff => ∃ msg ∈ Spin.nonzeroMsgs (RS8Group Field16),
        ∃ v ∈ Finset.Icc 5 8, rs8SampledPacketWeight points msg coeff = v)
      (fun _ => Iff.rfl)
  have hb := finite_double_union_bound P
    (Spin.nonzeroMsgs (RS8Group Field16)) (Finset.Icc 5 8)
    (fun msg v coeff => rs8SampledPacketWeight points msg coeff = v)
  apply le_trans (he.le.trans hb)
  apply le_of_eq
  apply Finset.sum_congr rfl
  intro v _
  unfold rs8ExpectedPacketCount
  apply Finset.sum_congr rfl
  intro msg _
  exact P.prob_congr
    (s := fun coeff => rs8SampledPacketWeight points msg coeff = v)
    (t := fun coeff => rs8SampledPacketWeight points msg coeff = v)
    (fun _ => Iff.rfl)
/-- A coordinate of the actual independent outer law has its declared law. -/
theorem turboOuterLaw22_prob_coordinate (R : ℕ) (i : Fin (R * 16))
    (E : (Fin 8 → {c : Field32 // c ≠ 0}) → Prop) [DecidablePred E] :
    (turboOuterLaw22 R).prob (fun seed => E (seed i)) =
      (piPMF (fun _ : Fin 8 =>
        FinPMF.uniform {c : Field32 // c ≠ 0})).prob E := by
  let P := piPMF (fun _ : Fin 8 =>
    FinPMF.uniform {c : Field32 // c ≠ 0})
  let f := fun coeff => if E coeff then (1 : ℝ) else 0
  rw [FinPMF.prob_eq_expect_indicator, FinPMF.prob_eq_expect_indicator]
  change (piPMF (fun _ : Fin (R * 16) => P)).expect
    (fun seed => f (seed i)) = P.expect f
  have h := piPMF_expect_prod (fun _ : Fin (R * 16) => P)
    (fun j coeff => if j = i then f coeff else 1)
  have he (j : Fin (R * 16)) :
      P.expect (fun coeff => if j = i then f coeff else 1) =
        if j = i then P.expect f else 1 := by
    by_cases hj : j = i <;> simp [hj, FinPMF.expect_const]
  simpa only [Finset.prod_ite_eq', Finset.mem_univ, ite_true, he] using h

theorem turboSparseRare22_prob_le (R : ℕ) (points : Fin 8 → Field16) :
    (turboOuterLaw22 R).prob (turboSparseRare22 R points) ≤
      (R * 16 : ℕ) * ∑ v ∈ Finset.Icc 5 8, rs8ExpectedPacketCount points v := by
  calc
    _ = (turboOuterLaw22 R).prob (fun seed =>
        ∃ i ∈ (Finset.univ : Finset (Fin (R * 16))),
          turboSparseRareGroup22 points (seed i)) :=
      (turboOuterLaw22 R).prob_congr (s := turboSparseRare22 R points)
        (t := fun seed => ∃ i ∈ (Finset.univ : Finset (Fin (R * 16))),
          turboSparseRareGroup22 points (seed i))
        (fun _ => by simp only [turboSparseRare22, Finset.mem_univ, true_and])
    _ ≤ ∑ i : Fin (R * 16), (turboOuterLaw22 R).prob
        (fun seed => turboSparseRareGroup22 points (seed i)) :=
      (turboOuterLaw22 R).prob_biUnion_le Finset.univ
        (fun i seed => turboSparseRareGroup22 points (seed i))
    _ = ∑ _i : Fin (R * 16), (piPMF (fun _ : Fin 8 =>
        FinPMF.uniform {c : Field32 // c ≠ 0})).prob
          (turboSparseRareGroup22 points) := by
      apply Finset.sum_congr rfl
      intro i _
      exact turboOuterLaw22_prob_coordinate R i _
    _ ≤ ∑ _i : Fin (R * 16),
        ∑ v ∈ Finset.Icc 5 8, rs8ExpectedPacketCount points v := by
      exact Finset.sum_le_sum fun _ _ => turboSparseRareGroup22_prob_le points
    _ = _ := by simp
/-- All active messages in all actual groups have at least nine packets
on the complement, simultaneously over the shared outer setup. -/
theorem turboSparseRare22_compl_weight_ge_nine
    (R : ℕ) (points : Fin 8 → Field16) (hp : Function.Injective points)
    (seed : TurboOuterSeed22 R) (hg : ¬ turboSparseRare22 R points seed)
    (i : Fin (R * 16)) (msg : RS8Group Field16) (hm : msg ≠ 0) :
    9 ≤ rs8SampledPacketWeight points msg (seed i) := by
  have hr := rs8_sampled_packet_weight_range points hp msg hm (seed i)
  by_contra h
  apply hg
  refine ⟨i, msg, ?_, rs8SampledPacketWeight points msg (seed i), ?_, rfl⟩
  · simp [Spin.nonzeroMsgs, hm]
  · simp only [Finset.mem_Icc]
    omega

#print axioms turboOuterLaw22_prob_coordinate
#print axioms turboSparseRare22_prob_le
#print axioms turboSparseRareGroup22_prob_le
#print axioms turboSparseRare22_compl_weight_ge_nine
end Spin.Turbo
















