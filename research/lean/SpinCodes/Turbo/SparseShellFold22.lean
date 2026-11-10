import SpinCodes.Turbo.SparseSupportCoefficient22

/-! Group the exact complete support choices by their packet cardinalities.
Only finite sum rearrangement is used; all Cv values, including 5..8, remain
present until a separately justified truncation is applied.

Proof engineering: every rearrangement is proved once over an abstract index
set and then instantiated. No broad simplification runs over the concrete
`Finset (Fin 64)` sums or over matrix applications of `List.prod`. -/
noncomputable section
namespace Spin.Turbo
open Finset SparseSupportCoefficient22

private theorem finite_sum_weighted_layers
    {X A M : Type*} [Fintype X] [DecidableEq A]
    [AddCommMonoid M] [Module ℝ M]
    (c : X → A) (t : Finset A) (hc : ∀ x, c x ∈ t) (w : A → ℝ) (F : X → M) :
    (∑ x, w (c x) • F x) = ∑ a ∈ t, w a • ∑ x, if c x = a then F x else 0 := by
  simp only [Finset.smul_sum, smul_ite, smul_zero]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro x _
  simp only [Finset.sum_ite_eq, hc x, ite_true]

private theorem support_card_mem (T : Finset (Fin 64)) : T.card ∈ Finset.range 65 := by
  have ht := Finset.card_le_card (Finset.subset_univ T)
  simp only [Finset.card_univ, Fintype.card_fin] at ht
  simp only [Finset.mem_range]
  omega

/-- Two-coordinate layering over an abstract `range n`. Instantiating the
product set at the concrete `range 65` makes the unifier evaluate
`List.range 65` and exceed the recursion limit, so `n` stays abstract here. -/
private theorem two_weighted_layers {X M : Type*} [Fintype X]
    [AddCommMonoid M] [Module ℝ M]
    (c₁ c₂ : X → ℕ) (n : ℕ) (h₁ : ∀ x, c₁ x ∈ Finset.range n)
    (h₂ : ∀ x, c₂ x ∈ Finset.range n) (w : ℕ → ℝ) (F : X → M) :
    (∑ x, (w (c₁ x) * w (c₂ x)) • F x) =
      ∑ u ∈ Finset.range n, ∑ v ∈ Finset.range n, (w u * w v) •
        ∑ x, if c₁ x = u ∧ c₂ x = v then F x else 0 := by
  have hc : ∀ x, (c₁ x, c₂ x) ∈ Finset.range n ×ˢ Finset.range n :=
    fun x => Finset.mk_mem_product (h₁ x) (h₂ x)
  refine (finite_sum_weighted_layers (fun x => (c₁ x, c₂ x))
    (Finset.range n ×ˢ Finset.range n) hc (fun uv : ℕ × ℕ => w uv.1 * w uv.2) F).trans ?_
  rw [Finset.sum_product]
  refine Finset.sum_congr rfl fun u _ => Finset.sum_congr rfl fun v _ => ?_
  exact congrArg (fun m : M => (w u * w v) • m)
    (Finset.sum_congr rfl fun x _ => if_congr Prod.mk_inj rfl rfl)

/-! Generic terminal-mass identities over abstract index sets. -/

private theorem mass_sum {X I : Type*} [Fintype I]
    (s : Finset X) (M : X → Matrix I I ℝ) (i : I) :
    (∑ k : I, (Finset.sum s M) i k) = ∑ x ∈ s, ∑ k : I, M x i k := by
  simp only [Matrix.sum_apply]
  exact Finset.sum_comm

private theorem mass_smul {I : Type*} [Fintype I] (a : ℝ) (M : Matrix I I ℝ) (i : I) :
    (∑ k : I, (a • M) i k) = a * ∑ k : I, M i k := by
  simp only [Matrix.smul_apply, smul_eq_mul, Finset.mul_sum]

private theorem mass_weighted_eq {X Y I : Type*} [Fintype I]
    (s : Finset X) (t : Finset Y) (a : X → ℝ) (b : Y → ℝ)
    (M : X → Matrix I I ℝ) (N : Y → Matrix I I ℝ) (i : I)
    (h : (∑ x ∈ s, a x • M x) = ∑ y ∈ t, b y • N y) :
    (∑ x ∈ s, a x * ∑ k : I, M x i k) = ∑ y ∈ t, b y * ∑ k : I, N y i k := by
  calc (∑ x ∈ s, a x * ∑ k : I, M x i k)
      = ∑ k : I, (∑ x ∈ s, a x • M x) i k := by
        rw [mass_sum]
        exact Finset.sum_congr rfl fun x _ => (mass_smul (a x) (M x) i).symm
    _ = ∑ k : I, (∑ y ∈ t, b y • N y) i k := by rw [h]
    _ = _ := by
        rw [mass_sum]
        exact Finset.sum_congr rfl fun y _ => mass_smul (b y) (N y) i

private theorem mass_weighted_eq₂ {X Y Z I : Type*} [Fintype I]
    (s : Finset X) (t : Finset Y) (t' : Finset Z) (a : X → ℝ) (b : Y → Z → ℝ)
    (M : X → Matrix I I ℝ) (N : Y → Z → Matrix I I ℝ) (i : I)
    (h : (∑ x ∈ s, a x • M x) = ∑ u ∈ t, ∑ v ∈ t', b u v • N u v) :
    (∑ x ∈ s, a x * ∑ k : I, M x i k) =
      ∑ u ∈ t, ∑ v ∈ t', b u v * ∑ k : I, N u v i k := by
  calc (∑ x ∈ s, a x * ∑ k : I, M x i k)
      = ∑ k : I, (∑ x ∈ s, a x • M x) i k := by
        rw [mass_sum]
        exact Finset.sum_congr rfl fun x _ => (mass_smul (a x) (M x) i).symm
    _ = ∑ k : I, (∑ u ∈ t, ∑ v ∈ t', b u v • N u v) i k := by rw [h]
    _ = _ := by
        rw [mass_sum]
        refine Finset.sum_congr rfl fun u _ => ?_
        rw [mass_sum]
        exact Finset.sum_congr rfl fun v _ => mass_smul (b u v) (N u v) i

/-! Matrix-valued support folds. -/

theorem sparseOneWeightedSupportFold22 {I : Type*} [Fintype I] [DecidableEq I]
    (w : ℕ → ℝ) (W : Bool → Matrix I I ℝ) :
    (∑ T : Finset (Fin 64), w T.card •
      (List.ofFn (fun region : Fin 64 => W (decide (region ∈ T)))).prod) =
      ∑ v ∈ Finset.range 65, w v •
        orderedMatrixCoefficient sparseBoolDegree22 W 64 v := by
  have h := finite_sum_weighted_layers Finset.card (Finset.range 65)
    support_card_mem w (fun T : Finset (Fin 64) =>
      (List.ofFn (fun region : Fin 64 => W (decide (region ∈ T)))).prod)
  refine h.trans (Finset.sum_congr rfl fun v _ => ?_)
  exact congrArg (fun M : Matrix I I ℝ => w v • M) (oneSupportCoefficient_eq_ordered W v)

theorem sparseTwoWeightedSupportFold22 {I : Type*} [Fintype I] [DecidableEq I]
    (w : ℕ → ℝ) (W : Bool × Bool → Matrix I I ℝ) :
    (∑ AB : Finset (Fin 64) × Finset (Fin 64), (w AB.1.card * w AB.2.card) •
      (List.ofFn (fun region : Fin 64 => W
        (decide (region ∈ AB.1), decide (region ∈ AB.2)))).prod) =
      ∑ u ∈ Finset.range 65, ∑ v ∈ Finset.range 65,
        (w u * w v) • sparseTwoShellCoefficient22 W u v := by
  refine (two_weighted_layers (fun AB : Finset (Fin 64) × Finset (Fin 64) => AB.1.card)
    (fun AB : Finset (Fin 64) × Finset (Fin 64) => AB.2.card) 65
    (fun AB => support_card_mem AB.1) (fun AB => support_card_mem AB.2) w
    (fun AB : Finset (Fin 64) × Finset (Fin 64) => (List.ofFn (fun region : Fin 64 => W
      (decide (region ∈ AB.1), decide (region ∈ AB.2)))).prod)).trans
    (Finset.sum_congr rfl fun u _ => Finset.sum_congr rfl fun v _ => ?_)
  exact congrArg (fun M : Matrix I I ℝ => (w u * w v) • M) (twoSupportCoefficient_eq_shell W u v)

/-- The kernel family may depend on the shell before the support sum. -/
theorem sparseOneWeightedSupportFold22_family {I : Type*} [Fintype I] [DecidableEq I]
    (w : ℕ → ℝ) (W : ℕ → Bool → Matrix I I ℝ) :
    (∑ T : Finset (Fin 64), w T.card •
      (List.ofFn (fun region : Fin 64 => W T.card (decide (region ∈ T)))).prod) =
      ∑ v ∈ Finset.range 65, w v •
        orderedMatrixCoefficient sparseBoolDegree22 (W v) 64 v := by
  have h := finite_sum_weighted_layers Finset.card (Finset.range 65)
    support_card_mem w (fun T : Finset (Fin 64) =>
      (List.ofFn (fun region : Fin 64 => W T.card (decide (region ∈ T)))).prod)
  refine h.trans (Finset.sum_congr rfl fun v _ => ?_)
  refine congrArg (fun M : Matrix I I ℝ => w v • M) ?_
  rw [← oneSupportCoefficient_eq_ordered]
  unfold oneSupportCoefficient
  refine Finset.sum_congr rfl fun T _ => ?_
  by_cases hT : T.card = v
  · subst hT
    rfl
  · rw [if_neg hT, if_neg hT]

/-! Terminal masses. -/

theorem sparseOneWeightedSupportFold22_family_mass
    {I : Type*} [Fintype I] [DecidableEq I]
    (initial : I) (w : ℕ → ℝ) (W : ℕ → Bool → Matrix I I ℝ) :
    (∑ T : Finset (Fin 64), w T.card * ∑ k : I,
      (List.ofFn (fun region : Fin 64 => W T.card (decide (region ∈ T)))).prod initial k) =
      ∑ v ∈ Finset.range 65, w v * ∑ k : I,
        orderedMatrixCoefficient sparseBoolDegree22 (W v) 64 v initial k :=
  mass_weighted_eq Finset.univ (Finset.range 65) (fun T : Finset (Fin 64) => w T.card) w
    (fun T : Finset (Fin 64) =>
      (List.ofFn (fun region : Fin 64 => W T.card (decide (region ∈ T)))).prod)
    (fun v => orderedMatrixCoefficient sparseBoolDegree22 (W v) 64 v) initial
    (sparseOneWeightedSupportFold22_family w W)

/-- Terminal mass of the same full coefficient identity. -/
theorem sparseOneWeightedSupportFold22_mass {I : Type*} [Fintype I] [DecidableEq I]
    (initial : I) (w : ℕ → ℝ) (W : Bool → Matrix I I ℝ) :
    (∑ T : Finset (Fin 64), w T.card * ∑ k : I,
      (List.ofFn (fun region : Fin 64 => W (decide (region ∈ T)))).prod initial k) =
      ∑ v ∈ Finset.range 65, w v * ∑ k : I,
        orderedMatrixCoefficient sparseBoolDegree22 W 64 v initial k :=
  mass_weighted_eq Finset.univ (Finset.range 65) (fun T : Finset (Fin 64) => w T.card) w
    (fun T : Finset (Fin 64) =>
      (List.ofFn (fun region : Fin 64 => W (decide (region ∈ T)))).prod)
    (fun v => orderedMatrixCoefficient sparseBoolDegree22 W 64 v) initial
    (sparseOneWeightedSupportFold22 w W)

theorem sparseTwoWeightedSupportFold22_mass {I : Type*} [Fintype I] [DecidableEq I]
    (initial : I) (w : ℕ → ℝ) (W : Bool × Bool → Matrix I I ℝ) :
    (∑ AB : Finset (Fin 64) × Finset (Fin 64), (w AB.1.card * w AB.2.card) *
      ∑ k : I, (List.ofFn (fun region : Fin 64 => W
        (decide (region ∈ AB.1), decide (region ∈ AB.2)))).prod initial k) =
      ∑ u ∈ Finset.range 65, ∑ v ∈ Finset.range 65,
        (w u * w v) * ∑ k : I, sparseTwoShellCoefficient22 W u v initial k :=
  mass_weighted_eq₂ Finset.univ (Finset.range 65) (Finset.range 65)
    (fun AB : Finset (Fin 64) × Finset (Fin 64) => w AB.1.card * w AB.2.card)
    (fun u v => w u * w v)
    (fun AB : Finset (Fin 64) × Finset (Fin 64) => (List.ofFn (fun region : Fin 64 => W
      (decide (region ∈ AB.1), decide (region ∈ AB.2)))).prod)
    (fun u v => sparseTwoShellCoefficient22 W u v) initial
    (sparseTwoWeightedSupportFold22 w W)

/-! Scalar support count and per-shell choice. -/

private theorem sparseOneSupportLayer22_count (v : ℕ) :
    (∑ T : Finset (Fin 64), if T.card = v then (1 : ℝ) else 0) =
      (Nat.choose 64 v : ℝ) := by
  have hs : (Finset.univ : Finset (Finset (Fin 64))).filter (fun T => T.card = v) =
      (Finset.univ : Finset (Fin 64)).powersetCard v := by
    ext T
    simp [Finset.mem_powersetCard]
  rw [← Finset.sum_filter, hs]
  simp

theorem sparseOneSupportScalarFold22 (w : ℕ → ℝ) :
    (∑ T : Finset (Fin 64), w T.card) =
      ∑ v ∈ Finset.range 65, w v * (Nat.choose 64 v : ℝ) := by
  have h := finite_sum_weighted_layers Finset.card (Finset.range 65)
    support_card_mem w (fun _ : Finset (Fin 64) => (1 : ℝ))
  calc (∑ T : Finset (Fin 64), w T.card)
      = ∑ T : Finset (Fin 64), w T.card • (1 : ℝ) :=
        Finset.sum_congr rfl fun T _ => (smul_eq_mul (w T.card) 1 ▸ mul_one (w T.card)).symm
    _ = _ := h.trans (Finset.sum_congr rfl fun v _ => by
        rw [smul_eq_mul]
        exact congrArg (fun r : ℝ => w v * r) (sparseOneSupportLayer22_count v))

private def chernoffWeight22 (w den : ℕ → ℝ) (useChernoff : ℕ → Bool) (v : ℕ) : ℝ :=
  if useChernoff v then w v / den v else 0

private def trivialWeight22 (w : ℕ → ℝ) (useChernoff : ℕ → Bool) (v : ℕ) : ℝ :=
  if useChernoff v then 0 else w v

/-- Choose the trivial probability bound or the Chernoff bound separately
for each shell. The trivial branch counts all its support choices. -/
theorem sparseOneWeightedSupportFold22_choice_mass
    {I : Type*} [Fintype I] [DecidableEq I]
    (initial : I) (w den : ℕ → ℝ) (useChernoff : ℕ → Bool)
    (W : ℕ → Bool → Matrix I I ℝ) :
    (∑ T : Finset (Fin 64), w T.card *
      (if useChernoff T.card then
        (∑ k : I, (List.ofFn (fun region : Fin 64 =>
          W T.card (decide (region ∈ T)))).prod initial k) / den T.card else 1)) =
      ∑ v ∈ Finset.range 65, if useChernoff v then
        (w v / den v) * ∑ k : I,
          orderedMatrixCoefficient sparseBoolDegree22 (W v) 64 v initial k
        else w v * (Nat.choose 64 v : ℝ) := by
  have hp (T : Finset (Fin 64)) : w T.card *
      (if useChernoff T.card then
        (∑ k : I, (List.ofFn (fun region : Fin 64 =>
          W T.card (decide (region ∈ T)))).prod initial k) / den T.card else 1) =
      chernoffWeight22 w den useChernoff T.card *
        (∑ k : I, (List.ofFn (fun region : Fin 64 =>
          W T.card (decide (region ∈ T)))).prod initial k) +
        trivialWeight22 w useChernoff T.card := by
    unfold chernoffWeight22 trivialWeight22
    split_ifs <;> ring
  calc _ = ∑ T : Finset (Fin 64), (chernoffWeight22 w den useChernoff T.card *
          (∑ k : I, (List.ofFn (fun region : Fin 64 =>
            W T.card (decide (region ∈ T)))).prod initial k) +
          trivialWeight22 w useChernoff T.card) :=
        Finset.sum_congr rfl fun T _ => hp T
    _ = (∑ T : Finset (Fin 64), chernoffWeight22 w den useChernoff T.card *
          (∑ k : I, (List.ofFn (fun region : Fin 64 =>
            W T.card (decide (region ∈ T)))).prod initial k)) +
          ∑ T : Finset (Fin 64), trivialWeight22 w useChernoff T.card :=
        Finset.sum_add_distrib
    _ = (∑ v ∈ Finset.range 65, chernoffWeight22 w den useChernoff v * ∑ k : I,
          orderedMatrixCoefficient sparseBoolDegree22 (W v) 64 v initial k) +
          ∑ v ∈ Finset.range 65, trivialWeight22 w useChernoff v * (Nat.choose 64 v : ℝ) := by
        rw [sparseOneWeightedSupportFold22_family_mass initial
            (chernoffWeight22 w den useChernoff) W,
          sparseOneSupportScalarFold22 (trivialWeight22 w useChernoff)]
    _ = _ := by
        rw [← Finset.sum_add_distrib]
        refine Finset.sum_congr rfl fun v _ => ?_
        unfold chernoffWeight22 trivialWeight22
        split_ifs <;> ring

#print axioms sparseOneWeightedSupportFold22
#print axioms sparseTwoWeightedSupportFold22
#print axioms sparseOneWeightedSupportFold22_family
#print axioms sparseOneWeightedSupportFold22_family_mass
#print axioms sparseOneWeightedSupportFold22_mass
#print axioms sparseTwoWeightedSupportFold22_mass
#print axioms sparseOneSupportScalarFold22
#print axioms sparseOneWeightedSupportFold22_choice_mass
end Spin.Turbo
