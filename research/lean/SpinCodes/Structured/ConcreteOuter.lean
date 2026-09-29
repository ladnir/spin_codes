import SpinCodes.Structured.ConcreteOuterGolay
import SpinCodes.Structured.ConcreteRouteTranspose
import SpinCodes.FiniteLaw

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Finset

abbrev LocalMessage (k : ℕ) := Fin k → (Fin 12 → Bool)
abbrev Seed (k : ℕ) := Equiv.Perm (Fin (k * 24)) × Equiv.Perm (Fin (k * 24))

/-- One pair of independent permutations is sampled for the shared constituent. -/
def seedLaw (k : ℕ) : FinPMF (Seed k) :=
  (FinPMF.uniform (Equiv.Perm (Fin (k * 24)))).prod
    (FinPMF.uniform (Equiv.Perm (Fin (k * 24))))

/-- The actual Golay direct sum followed by two permute--accumulate stages. -/
def encode {k : ℕ} (seed : Seed k) (x : LocalMessage k) : Fin (k * 24) → Bool :=
  accF ((accF (outerWord Golay.enc x ∘ seed.1)) ∘ seed.2)

/-- Bit support used by the concrete route. -/
def support {n : ℕ} (x : Fin n → Bool) : Finset (Fin n) := univ.filter (fun j => x j = true)

@[simp] theorem support_card {n : ℕ} (x : Fin n → Bool) : (support x).card = wtF x :=
  (wtF_eq_card x).symm

@[simp] theorem support_eq_empty {n : ℕ} (x : Fin n → Bool) : support x = ∅ ↔ x = 0 := by
  rw [← Finset.card_eq_zero, support_card, wtF_eq_zero]

theorem support_injective {n : ℕ} : Function.Injective (@support n) := by
  intro x y h
  funext j
  have hh := Finset.ext_iff.mp h j
  simp only [support, mem_filter, mem_univ, true_and] at hh
  cases hx : x j <;> cases hy : y j <;> simp_all

theorem golay_zero : Golay.enc (0 : Fin 12 → Bool) = 0 :=
  (wtF_eq_zero _).mp Golay.W_zero

theorem outerWord_injective {k : ℕ} : Function.Injective
    (outerWord Golay.enc : LocalMessage k → (Fin (k * 24) → Bool)) := by
  intro x y h
  funext i
  apply golay_injective
  funext j
  have hh := congrFun h (finProdFinEquiv (i,j))
  simpa only [outerWord, concatBlocks, Equiv.symm_apply_apply] using hh

theorem permute_injective {n : ℕ} (σ : Equiv.Perm (Fin n)) :
    Function.Injective (fun x : Fin n → Bool => x ∘ σ) := by
  intro x y h
  funext j
  simpa only [Function.comp_apply, Equiv.apply_symm_apply] using congrFun h (σ.symm j)

theorem encode_injective {k : ℕ} (seed : Seed k) : Function.Injective (encode seed) :=
  accF_injective.comp ((permute_injective seed.2).comp
    (accF_injective.comp ((permute_injective seed.1).comp outerWord_injective)))

@[simp] theorem encode_zero {k : ℕ} (seed : Seed k) : encode seed 0 = 0 := by
  have h : outerWord Golay.enc (0 : LocalMessage k) = 0 := by
    funext j
    simp only [outerWord, concatBlocks, Pi.zero_apply, golay_zero]
  rw [encode, h]
  change accF (accF (0 : Fin (k * 24) → Bool) ∘ seed.2) = 0
  rw [accF_zero]
  exact accF_zero _

@[simp] theorem encode_eq_zero {k : ℕ} (seed : Seed k) (x : LocalMessage k) :
    encode seed x = 0 ↔ x = 0 :=
  ⟨fun h => encode_injective seed (h.trans (encode_zero seed).symm), fun h => h ▸ encode_zero seed⟩

/-- Reuse the same sampled constituent at every outer position. -/
def rowSupports {L k : ℕ} (seed : Seed k) (messages : Fin L → LocalMessage k) :
    Fin L → Finset (Fin (k * 24)) := fun i => support (encode seed (messages i))

@[simp] theorem rowSupports_empty_iff {L k : ℕ} (seed : Seed k)
    (messages : Fin L → LocalMessage k) (i : Fin L) :
    rowSupports seed messages i = ∅ ↔ messages i = 0 := by
  simp only [rowSupports, support_eq_empty, encode_eq_zero]

/-- Concrete occupation counts nonzero message rows and equals routed occupation. -/
def occupation {L k : ℕ} (messages : Fin L → LocalMessage k) : ℕ :=
  (univ.filter (fun i => messages i ≠ 0)).card

theorem activeRows_eq {L k : ℕ} (seed : Seed k) (messages : Fin L → LocalMessage k) :
    ConcreteRoute.activeRows (rowSupports seed messages) = occupation messages := by
  simp only [ConcreteRoute.activeRows, occupation, ne_eq, rowSupports_empty_iff]

theorem occupation_pos {L k : ℕ} (messages : Fin L → LocalMessage k) (hne : messages ≠ 0) :
    0 < occupation messages := by
  apply Finset.card_pos.mpr
  obtain ⟨i, hi⟩ : ∃ i, messages i ≠ 0 := by
    by_contra h
    simp only [not_exists, not_not] at h
    exact hne (funext h)
  exact ⟨i, by simp [hi]⟩

theorem occupation_le {L k : ℕ} (messages : Fin L → LocalMessage k) :
    occupation messages ≤ L := by
  simpa only [occupation, card_univ, Fintype.card_fin] using
    Finset.card_le_card (filter_subset (fun i => messages i ≠ 0) univ)

theorem totalWeight_eq {L k : ℕ} (seed : Seed k) (messages : Fin L → LocalMessage k) :
    ConcreteRoute.totalWeight (rowSupports seed messages) = ∑ i, wtF (encode seed (messages i)) := by
  simp only [ConcreteRoute.totalWeight, rowSupports, support_card]

theorem rowSupports_injective {L k : ℕ} (seed : Seed k) :
    Function.Injective (rowSupports (L := L) seed) := by
  intro x y h
  funext i
  exact encode_injective seed (support_injective (congrFun h i))

/-- The number of nonzero local messages with the indicated realized output weight. -/
def spectrum {k : ℕ} (seed : Seed k) (w : ℕ) : ℕ :=
  (univ.filter (fun x : LocalMessage k => x ≠ 0 ∧ wtF (encode seed x) = w)).card

theorem spectrum_eq_sum {k : ℕ} (seed : Seed k) (w : ℕ) :
    (spectrum seed w : ℝ) =
      ∑ x ∈ univ.filter (fun x : LocalMessage k => x ≠ 0),
        if wtF (encode seed x) = w then (1 : ℝ) else 0 := by
  simp only [spectrum, Finset.card_filter, Nat.cast_sum, Nat.cast_ite, Nat.cast_one,
    Nat.cast_zero, Finset.sum_filter]
  apply Finset.sum_congr rfl
  intro x _
  split_ifs <;> simp_all

/-- The established Golay-BAA spectrum is the expectation of this concrete encoder. -/
theorem expected_spectrum (k w : ℕ) :
    (seedLaw k).expect (fun seed => (spectrum seed w : ℝ)) =
      ∑ a ∈ Icc 1 (k * 24),
        ((((1 : Polynomial ℕ) + Polynomial.C 759 * Polynomial.X ^ 8 +
          Polynomial.C 2576 * Polynomial.X ^ 12 + Polynomial.C 759 * Polynomial.X ^ 16 +
          Polynomial.X ^ 24) ^ k).coeff a : ℝ) *
          ∑ c ∈ range (k * 24 + 1), Pt (k * 24) a c * Pt (k * 24) c w := by
  simp_rw [spectrum_eq_sum]
  rw [FinPMF.expect_sum]
  simp_rw [← FinPMF.prob_eq_expect_indicator]
  simpa only [seedLaw, encode, wtF_accF] using Golay.ba_expected_spectrum k w

end Spin.Structured.ConcreteOuter

