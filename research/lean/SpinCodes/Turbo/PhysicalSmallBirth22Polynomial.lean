import SpinCodes.Turbo.PhysicalSmallBirth22Census
import SpinCodes.Turbo.PhysicalCappedBirth22

noncomputable section
namespace Spin.Turbo
open Finset Spin.UniformMixer

/- The historical `stepKernel_empty` and `emptyStep_zero` keep an entering
zero state exactly zero. The exact-j=0 layer first has to be joined to that
empty input; a coarse nonzero birth budget is unnecessary here. -/
theorem physicalSmallBirth22_zero_word (x : PacketStepLayer22 0) : x.val = 0 := by
  have hc : (packetStepSupport22 x.val).card = 0 := x.property
  have hs : packetStepSupport22 x.val = ∅ := card_eq_zero.mp hc
  funext i
  by_contra hi
  have hm : i ∈ packetStepSupport22 x.val := mem_filter.mpr ⟨mem_univ _, hi⟩
  rw [hs] at hm
  simpa using hm

theorem physicalSmallBirth22_zero_expect (f : PacketStep16 → ℝ) :
    (packetStepLaw22 0).expect f = f 0 := by
  rw [packetStepLaw22, FinPMF.expect_map]
  have hf : (fun x : PacketStepLayer22 0 => f x.val) = fun _ => f 0 := by
    funext x
    rw [physicalSmallBirth22_zero_word x]
  rw [hf, FinPMF.expect_const]

theorem physicalSmallBirth22_zero_occupancy (z : ℝ) (q r : State22) :
    turboPhysicalKernel22 0 z q r = emptyStep z (coordinateWeight coordMap) q r := by
  unfold turboPhysicalKernel22
  rw [physicalSmallBirth22_zero_expect]
  change stepKernel z 0 q r = _
  rw [stepKernel_empty]

theorem physicalSmallBirthKernel22_zero (z : ℝ) (r : State22) :
    turboPhysicalKernel22 0 z 0 r = if r = 0 then 1 else 0 := by
  rw [physicalSmallBirth22_zero_occupancy]
  exact emptyStep_zero z (coordinateWeight coordMap) (coordinateWeight_zero coordMap) r

theorem physicalSmallBirth22_zero_cone (z : ℝ) :
    physicalCone22 physicalInitialEnvelope22 (fun r => turboPhysicalKernel22 0 z 0 r) := by
  have hr : (fun r => turboPhysicalKernel22 0 z 0 r) = physicalInitialMass22 := by
    funext r
    exact physicalSmallBirthKernel22_zero z r
  rw [hr]
  exact physicalInitialCone22
theorem physicalSmallBirth22_zero_exact_cone (z : ℝ) :
    physicalCone22 physicalInitialEnvelope22 (turboExactBirthMass22 0 z) := by
  simpa only [turboPhysicalKernel22_zero_entry] using physicalSmallBirth22_zero_cone z
private theorem uniformSmallBirthPolynomial22 {Ω : Type*} [Fintype Ω] [Nonempty Ω]
    (c : Ω → Fin 11) (w : Ω → Nat) (hw : ∀ x, w x < 9) (hist : List Nat)
    (hcount : ∀ b < 99,
      (univ.filter (fun x => 9*(c x).val+w x = b)).card = hist.getD b 0)
    (z : ℝ) (i : Fin 11) :
    (FinPMF.uniform Ω).expect (fun x => if c x = i then z ^ w x else 0) =
      (1 / (Fintype.card Ω : ℝ)) *
        ∑ t : Fin 9, (hist.getD (9*i.val+t.val) 0 : ℝ) * z ^ t.val := by
  letI : DecidableEq Ω := Classical.decEq Ω
  let X : Ω → Fin 11 × Fin 9 := fun x => (c x, ⟨w x, hw x⟩)
  let g : Fin 11 × Fin 9 → ℝ := fun a => if a.1 = i then z ^ a.2.val else 0
  have he : (FinPMF.uniform Ω).expect (fun x => if c x = i then z ^ w x else 0) =
      (FinPMF.uniform Ω).expect (fun x => g (X x)) := rfl
  rw [he, FinPMF.expect_comp _ X univ (by intro x; exact mem_univ _) g]
  simp only [Fintype.sum_prod_type]
  have hr (a : Fin 11) :
      (∑ t : Fin 9, (FinPMF.uniform Ω).prob (fun x => X x = (a,t)) * g (a,t)) =
        if a = i then ∑ t : Fin 9,
          (FinPMF.uniform Ω).prob (fun x => X x = (i,t)) * z ^ t.val else 0 := by
    by_cases h : a = i
    · subst a; simp [g]
    · simp [g, h]
  simp_rw [hr]
  simp only [sum_ite_eq', mem_univ, ite_true]
  rw [mul_sum]
  apply sum_congr rfl
  intro t _
  have hb : 9*i.val+t.val < 99 := by have hi := i.isLt; have ht := t.isLt; omega
  have hfilter :
      univ.filter (fun x => X x = (i,t)) =
        univ.filter (fun x => 9*(c x).val+w x = 9*i.val+t.val) := by
    ext x
    simp only [mem_filter, mem_univ, true_and]
    constructor
    · intro h
      have h1 : c x = i := congrArg Prod.fst h
      have h2 : w x = t.val := congrArg (fun a : Fin 11 × Fin 9 => a.2.val) h
      rw [h1, h2]
    · intro h
      have hwx := hw x
      have ht := t.isLt
      have hc : (c x).val = i.val := by omega
      have hw' : w x = t.val := by omega
      apply Prod.ext
      · exact Fin.ext hc
      · exact Fin.ext hw'
  rw [FinPMF.uniform_prob, hfilter, hcount _ hb]
  ring

set_option maxRecDepth 10000 in
private theorem smallBirthLabel22_weight_lt : ∀ v : Fin 15,
    smallBirthPacketWeight22 (smallBirthLabel22 v) < 5 := by decide +kernel

private theorem smallBirthOneLayer22_weight_lt (x : PacketStepLayer22 1) :
    stepWeight (packetStepInput22 x.val) 0 < 9 := by
  obtain ⟨a, rfl⟩ := smallBirthOneEquiv22.surjective x
  change stepWeight (packetStepInput22 (smallBirthOneWord22 a)) 0 < 9
  rw [smallBirthOneWord22_weight]
  exact (smallBirthLabel22_weight_lt a.2).trans (by decide)

private theorem smallBirthTwoLayer22_weight_lt (x : PacketStepLayer22 2) :
    stepWeight (packetStepInput22 x.val) 0 < 9 := by
  obtain ⟨a, rfl⟩ := smallBirthTwoEquiv22.surjective x
  change stepWeight (packetStepInput22 (smallBirthTwoWord22 a)) 0 < 9
  rw [smallBirthTwoWord22_weight]
  have h1 := smallBirthLabel22_weight_lt a.2.1
  have h2 := smallBirthLabel22_weight_lt a.2.2
  omega

/-- Exact-j physical birth class mass, before selected-slot thinning. -/
def physicalSmallBirthClassMass22 (j : Fin 17) (z : ℝ) (i : Fin 11) : ℝ :=
  (packetStepLaw22 j).expect (fun x =>
    if physicalBirthIndex22 (stepFeedback (packetStepInput22 x)) = i then
      z ^ stepWeight (packetStepInput22 x) 0 else 0)

theorem physicalSmallBirthClassMass22_nonneg (j : Fin 17) (z : ℝ)
    (hz : 0 ≤ z) (i : Fin 11) : 0 ≤ physicalSmallBirthClassMass22 j z i := by
  unfold physicalSmallBirthClassMass22
  apply FinPMF.expect_nonneg
  intro x
  split_ifs
  · exact pow_nonneg hz _
  · exact le_refl 0
theorem physicalSmallBirthClassMass22_empty (z : ℝ) (i : Fin 11) :
    physicalSmallBirthClassMass22 0 z i = physicalInitialEnvelope22 i := by
  unfold physicalSmallBirthClassMass22
  rw [physicalSmallBirth22_zero_expect]
  have hi : packetStepInput22 (0 : PacketStep16) = 0 := rfl
  rw [hi, stepFeedback_empty, stepWeight_empty, coordinateWeight_zero]
  have hz0 : physicalBirthIndex22 (0 : State22) = 0 := by simp [physicalBirthIndex22]
  rw [hz0]
  by_cases hi0 : i = 0
  · subst i
    simp [physicalInitialEnvelope22]
  · have hval : i.val ≠ 0 := fun h => hi0 (Fin.ext h)
    simp only [physicalInitialEnvelope22, Ne.symm hi0, hval, ite_false]
theorem physicalSmallBirthClassMass22_one (z : ℝ) (i : Fin 11) :
    physicalSmallBirthClassMass22 1 z i =
      (1 / 240 : ℝ) * ∑ t : Fin 9,
        (smallBirthOneHistogram22.getD (9*i.val+t.val) 0 : ℝ) * z ^ t.val := by
  rw [physicalSmallBirthClassMass22, packetStepLaw22, FinPMF.expect_map]
  have h := uniformSmallBirthPolynomial22
    (fun x : PacketStepLayer22 1 => physicalBirthIndex22 (stepFeedback (packetStepInput22 x.val)))
    (fun x => stepWeight (packetStepInput22 x.val) 0) smallBirthOneLayer22_weight_lt
    smallBirthOneHistogram22 smallBirthOneLayer22_count z i
  have hc : (Fintype.card (PacketStepLayer22 1) : ℝ) = 240 := by
    rw [packetStepLayer_card22]
    norm_num [Nat.choose]
  simpa only [hc] using h

theorem physicalSmallBirthClassMass22_two (z : ℝ) (i : Fin 11) :
    physicalSmallBirthClassMass22 2 z i =
      (1 / 27000 : ℝ) * ∑ t : Fin 9,
        (smallBirthTwoHistogram22.getD (9*i.val+t.val) 0 : ℝ) * z ^ t.val := by
  rw [physicalSmallBirthClassMass22, packetStepLaw22, FinPMF.expect_map]
  have h := uniformSmallBirthPolynomial22
    (fun x : PacketStepLayer22 2 => physicalBirthIndex22 (stepFeedback (packetStepInput22 x.val)))
    (fun x => stepWeight (packetStepInput22 x.val) 0) smallBirthTwoLayer22_weight_lt
    smallBirthTwoHistogram22 smallBirthTwoLayer22_count z i
  have hc : (Fintype.card (PacketStepLayer22 2) : ℝ) = 27000 := by
    rw [packetStepLayer_card22]
    norm_num [Nat.choose]
  simpa only [hc] using h

/-- This is the actual full-state zero row, grouped only at its terminal state. -/
theorem physicalSmallBirthClassMass22_eq_kernel (j : Fin 17) (z : ℝ) (i : Fin 11) :
    physicalSmallBirthClassMass22 j z i =
      ∑ r : State22, if physicalBirthIndex22 r = i then turboPhysicalKernel22 j z 0 r else 0 := by
  unfold physicalSmallBirthClassMass22 turboPhysicalKernel22
  have hi (r : State22) :
      (if physicalBirthIndex22 r = i then
        (packetStepLaw22 j).expect (fun x => stepKernel z (packetStepInput22 x) 0 r) else 0) =
      (packetStepLaw22 j).expect (fun x =>
        if physicalBirthIndex22 r = i then stepKernel z (packetStepInput22 x) 0 r else 0) := by
    by_cases h : physicalBirthIndex22 r = i
    · simp only [h, ite_true]
    · simp only [h, ite_false, FinPMF.expect_const]
  simp_rw [hi]
  rw [← FinPMF.expect_sum]
  apply congrArg
  funext x
  classical
  by_cases h : physicalBirthIndex22 (stepFeedback (packetStepInput22 x)) = i
  · rw [if_pos h]
    symm
    rw [sum_eq_single (stepFeedback (packetStepInput22 x))]
    · simp [stepKernel, update_zero, h]
    · intro r _ hr
      simp [stepKernel, update_zero, hr, Ne.symm hr]
    · intro hnot
      exact (hnot (mem_univ _)).elim
  · rw [if_neg h]
    symm
    apply sum_eq_zero
    intro r _
    by_cases hr : r = stepFeedback (packetStepInput22 x)
    · subst r
      simp [h]
    · simp [stepKernel, update_zero, hr, Ne.symm hr]

theorem physicalSmallBirthClassMass22_eq_exact (j : Fin 17) (z : ℝ) (i : Fin 11) :
    physicalSmallBirthClassMass22 j z i =
      ∑ r : State22, if physicalBirthIndex22 r = i then turboExactBirthMass22 j z r else 0 := by
  rw [physicalSmallBirthClassMass22_eq_kernel]
  simp_rw [turboPhysicalKernel22_zero_entry]
/- The historical selected-law `turboBirthParts22_piece` partitions one
measure by physicalBirthIndex22. The same partition here uses the actual
exact-j birth law; all numerical class-budget bounds remain explicit. -/
theorem physicalSmallBirth22_cone_of_class_bounds (j : Fin 17) (z : ℝ)
    (hz : 0 ≤ z) (v : Fin 11 → ℝ)
    (hv : ∀ i, physicalSmallBirthClassMass22 j z i ≤ v i) :
    physicalCone22 v (turboExactBirthMass22 j z) := by
  let parts : Fin 11 → State22 → ℝ := fun i r =>
    if physicalBirthIndex22 r = i then turboExactBirthMass22 j z r else 0
  refine ⟨parts, ?_, ?_⟩
  · intro r
    simp only [parts, Finset.sum_ite_eq, Finset.mem_univ, ite_true]
  · intro i
    refine ⟨?_, ?_, ?_, ?_⟩
    · intro r
      dsimp only [parts]
      split_ifs
      · exact turboExactBirthMass22_nonneg j z hz r
      · exact le_refl 0
    · intro r hr
      have hi : physicalBirthIndex22 r = i := by
        by_contra h
        exact hr (by simp only [parts, h, ite_false])
      rw [← hi]
      exact physicalBirthIndex22_allows r
    · change (∑ r : State22, if physicalBirthIndex22 r = i then turboExactBirthMass22 j z r else 0) ≤ v i
      rw [← physicalSmallBirthClassMass22_eq_exact]
      exact hv i
    · intro hi r
      have hi2 : i = 2 := Fin.ext hi
      have hindex : physicalBirthIndex22 r ≠ i := by
        rw [hi2]
        exact physicalBirthIndex22_ne_two r
      have hp : parts i r = 0 := by simp only [parts, hindex, ite_false]
      rw [hp]
      exact div_nonneg ((physicalSmallBirthClassMass22_nonneg j z hz i).trans (hv i)) (by norm_num)

def physicalSmallBirthOneBudget22 (z : ℝ) (i : Fin 11) : ℝ :=
  (1 / 240 : ℝ) * ∑ t : Fin 9,
    (smallBirthOneHistogram22.getD (9*i.val+t.val) 0 : ℝ) * z ^ t.val

def physicalSmallBirthTwoBudget22 (z : ℝ) (i : Fin 11) : ℝ :=
  (1 / 27000 : ℝ) * ∑ t : Fin 9,
    (smallBirthTwoHistogram22.getD (9*i.val+t.val) 0 : ℝ) * z ^ t.val

theorem physicalSmallBirth22_one_cone (z : ℝ) (hz : 0 ≤ z) :
    physicalCone22 (physicalSmallBirthOneBudget22 z) (turboExactBirthMass22 1 z) := by
  apply physicalSmallBirth22_cone_of_class_bounds 1 z hz
  intro i
  exact (physicalSmallBirthClassMass22_one z i).le

theorem physicalSmallBirth22_two_cone (z : ℝ) (hz : 0 ≤ z) :
    physicalCone22 (physicalSmallBirthTwoBudget22 z) (turboExactBirthMass22 2 z) := by
  apply physicalSmallBirth22_cone_of_class_bounds 2 z hz
  intro i
  exact (physicalSmallBirthClassMass22_two z i).le

theorem physicalSmallBirth22_one_cone_of_le (z : ℝ) (hz : 0 ≤ z)
    (v : Fin 11 → ℝ) (hv : ∀ i, physicalSmallBirthOneBudget22 z i ≤ v i) :
    physicalCone22 v (turboExactBirthMass22 1 z) := by
  apply physicalSmallBirth22_cone_of_class_bounds 1 z hz
  intro i
  exact (physicalSmallBirthClassMass22_one z i).le.trans (hv i)

theorem physicalSmallBirth22_two_cone_of_le (z : ℝ) (hz : 0 ≤ z)
    (v : Fin 11 → ℝ) (hv : ∀ i, physicalSmallBirthTwoBudget22 z i ≤ v i) :
    physicalCone22 v (turboExactBirthMass22 2 z) := by
  apply physicalSmallBirth22_cone_of_class_bounds 2 z hz
  intro i
  exact (physicalSmallBirthClassMass22_two z i).le.trans (hv i)
theorem physicalSmallBirthClassMass22_zero (j : Fin 17) (z : ℝ) :
    physicalSmallBirthClassMass22 j z 0 = turboPhysicalKernel22 j z 0 0 := by
  letI : DecidableEq State22 := fun a b => Classical.propDecidable (a = b)
  rw [physicalSmallBirthClassMass22_eq_kernel]
  have hsum : (∑ r : State22, if physicalBirthIndex22 r = 0 then turboPhysicalKernel22 j z 0 r else 0) =
      (if physicalBirthIndex22 0 = 0 then turboPhysicalKernel22 j z 0 0 else 0) := by
    apply sum_eq_single 0
    · intro r _ hr
      by_cases hi : physicalBirthIndex22 r = 0
      · have ha := physicalBirthIndex22_allows r
        rw [hi] at ha
        have hz : r = 0 := by simpa only [physicalClassAllows22, Fin.val_zero, ite_true] using ha
        exact (hr hz).elim
      · simp only [hi, ite_false]
    · intro h
      exact (h (mem_univ _)).elim
  simpa [physicalBirthIndex22] using hsum
theorem physicalSmallBirthZero22_one (z : ℝ) :
    turboPhysicalKernel22 1 z 0 0 = (1 / 240 : ℝ) *
      ∑ t : Fin 9, (smallBirthOneHistogram22.getD t.val 0 : ℝ) * z ^ t.val := by
  calc
    _ = physicalSmallBirthClassMass22 1 z 0 := (physicalSmallBirthClassMass22_zero 1 z).symm
    _ = _ := by simpa only [Fin.val_zero, Nat.mul_zero, Nat.zero_add] using
      (physicalSmallBirthClassMass22_one z 0)

theorem physicalSmallBirthZero22_two (z : ℝ) :
    turboPhysicalKernel22 2 z 0 0 = (1 / 27000 : ℝ) *
      ∑ t : Fin 9, (smallBirthTwoHistogram22.getD t.val 0 : ℝ) * z ^ t.val := by
  calc
    _ = physicalSmallBirthClassMass22 2 z 0 := (physicalSmallBirthClassMass22_zero 2 z).symm
    _ = _ := by simpa only [Fin.val_zero, Nat.mul_zero, Nat.zero_add] using
      (physicalSmallBirthClassMass22_two z 0)
#print axioms physicalSmallBirthClassMass22_one
#print axioms physicalSmallBirthClassMass22_two
#print axioms physicalSmallBirthZero22_two
#print axioms physicalSmallBirth22_zero_cone
#print axioms physicalSmallBirth22_two_cone_of_le
end Spin.Turbo















