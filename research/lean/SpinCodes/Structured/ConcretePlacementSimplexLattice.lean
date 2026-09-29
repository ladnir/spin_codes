import SpinCodes.Structured.ConcretePlacementSimplexDomain

/-! Exact bijection between sorted block subsets and ordered-site lattice points. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing ConcreteEncoder Filter MeasureTheory Bornology
open scoped Topology Pointwise

def normalizedSites {R a : Nat} (B : BlockSubset R a) (i : Fin a) : ℝ :=
  (orderedBlocks B i).val / (R : ℝ)

theorem integerGrid_iff {a R : Nat} (hR : 0 < R) (x : Fin a → ℝ) :
    x ∈ integerGrid a R ↔ ∀ i, ∃ z : ℤ, (z : ℝ) = (R : ℝ) * x i := by
  letI : NeZero R := ⟨Nat.ne_of_gt hR⟩
  unfold integerGrid
  rw [← Submodule.coe_pointwise_smul]
  convert! (BoxIntegral.unitPartition.mem_smul_span_iff (n := R) (v := x)) using 1 <;> simp

theorem normalizedSites_mem {R a : Nat} (hR : 0 < R) (B : BlockSubset R a) :
    normalizedSites B ∈ orderedSiteDomain a ∩ integerGrid a R := by
  have hRp : (0 : ℝ) < R := by exact_mod_cast hR
  constructor
  · constructor
    · intro i
      exact ⟨div_nonneg (Nat.cast_nonneg _) hRp.le,
        (div_lt_one hRp).mpr (by exact_mod_cast (orderedBlocks B i).isLt)⟩
    · intro i j hij
      exact div_lt_div_of_pos_right (by exact_mod_cast orderedBlocks_strictMono B hij) hRp
  · rw [integerGrid_iff hR]
    intro i
    refine ⟨(orderedBlocks B i).val, ?_⟩
    dsimp [normalizedSites]
    push_cast
    field_simp

theorem normalizedSites_injective {R a : Nat} (hR : 0 < R) :
    Function.Injective (normalizedSites (R := R) (a := a)) := by
  intro B C he
  have hRp : (R : ℝ) ≠ 0 := by exact_mod_cast Nat.ne_of_gt hR
  have hf : (orderedBlocks B : Fin a → Fin R) = orderedBlocks C := by
    funext i
    have hi := congrFun he i
    dsimp [normalizedSites] at hi
    have hv : ((orderedBlocks B i).val : ℝ) = (orderedBlocks C i).val := by
      exact (div_left_inj' hRp).mp hi
    exact Fin.ext (Nat.cast_injective hv)
  apply Subtype.ext
  rw [← orderedBlocks_image B, ← orderedBlocks_image C, hf]

theorem normalizedSites_surjective {R a : Nat} (hR : 0 < R)
    (x : Fin a → ℝ) (hx : x ∈ orderedSiteDomain a ∩ integerGrid a R) :
    ∃ B : BlockSubset R a, normalizedSites B = x := by
  have hRp : (0 : ℝ) < R := by exact_mod_cast hR
  obtain ⟨hd, hg⟩ := hx
  choose z hz using (integerGrid_iff hR x).mp hg
  have hz0 (i : Fin a) : 0 ≤ z i := by
    have h : (0 : ℝ) ≤ (z i : ℝ) := by rw [hz]; exact mul_nonneg hRp.le (hd.1 i).1
    exact_mod_cast h
  have hn (i : Fin a) : ((z i).toNat : ℝ) = (R : ℝ) * x i := by
    have h := Int.toNat_of_nonneg (hz0 i)
    have hc : ((z i).toNat : ℝ) = (z i : ℝ) := by exact_mod_cast h
    exact hc.trans (hz i)
  let f : Fin a → Fin R := fun i => ⟨(z i).toNat, by
    have h : ((z i).toNat : ℝ) < R := by rw [hn]; nlinarith [(hd.1 i).2]
    exact_mod_cast h⟩
  have hmono : StrictMono f := by
    intro i j hij
    change (z i).toNat < (z j).toNat
    have h : ((z i).toNat : ℝ) < (z j).toNat := by rw [hn, hn]; exact mul_lt_mul_of_pos_left (hd.2 hij) hRp
    exact_mod_cast h
  let B : BlockSubset R a := ⟨univ.image f, by rw [card_image_of_injective _ hmono.injective]; simp⟩
  have hf : f = orderedBlocks B :=
    Finset.orderEmbOfFin_unique B.property (fun i => mem_image.mpr ⟨i, mem_univ _, rfl⟩) hmono
  refine ⟨B, ?_⟩
  funext i
  unfold normalizedSites
  rw [← hf]
  change ((z i).toNat : ℝ) / R = x i
  rw [hn]
  field_simp

def orderedGridEquiv {R a : Nat} (hR : 0 < R) :
    BlockSubset R a ≃ ↑(orderedSiteDomain a ∩ integerGrid a R) :=
  Equiv.ofBijective (fun B => ⟨normalizedSites B, normalizedSites_mem hR B⟩) ⟨
    fun _ _ h => normalizedSites_injective hR (congrArg Subtype.val h), by
      intro x
      obtain ⟨B, hB⟩ := normalizedSites_surjective hR x.val x.property
      exact ⟨B, Subtype.ext hB⟩⟩

theorem orderedSite_sum_eq_lattice {R a : Nat} (hR : 0 < R) (F : (Fin a → ℝ) → ℝ) :
    (∑ B : BlockSubset R a, F (normalizedSites B)) =
      ∑' x : ↑(orderedSiteDomain a ∩ integerGrid a R), F x := by
  have h := (orderedGridEquiv (a := a) hR).tsum_eq
    (fun x : ↑(orderedSiteDomain a ∩ integerGrid a R) => F x.val)
  rw [tsum_fintype] at h
  exact h

theorem orderedSite_sum_tendsto {a : Nat} (F : (Fin a → ℝ) → ℝ) (hF : Continuous F) :
    Tendsto (fun R : Nat => (∑ B : BlockSubset R a, F (normalizedSites B)) / (R : ℝ)^a)
      atTop (𝓝 (∫ x in orderedSiteDomain a, F x)) := by
  apply (orderedSite_riemann_tendsto F hF).congr'
  filter_upwards [eventually_gt_atTop 0] with R hR
  rw [orderedSite_sum_eq_lattice hR]

end Spin.Structured.Placement




