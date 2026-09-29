/-
Binary Krawtchouk polynomials.

`app:imt-finite-transfers` gets its integer fiber bounds from

    k_j = 2^{-19} (n_j + ∑_w b_w K_j(w)),
    V_j = 2^{-19} (n_j² + ∑_w b_w K_j(w)²),

which are Fourier inversion and Parseval over `𝔽₂^19`.  Both rest on a single
combinatorial fact: the character sum over a weight layer,

    ∑_{|x| = j} (-1)^{|x ∩ v|},

depends on `v` only through `|v|`, and equals `K_j(|v|)`.  Mathlib has neither
Krawtchouk polynomials nor MacWilliams, so this is built from scratch.

Vectors are represented as subsets, so weight is `Finset.card` and the inner
product `⟨v,x⟩` is `|v ∩ x|` mod 2.
-/
import SpinCodes.Structured.Domination

namespace Spin

open Finset

/-- The binary Krawtchouk polynomial `K_j(w)` for length `n`. -/
def krawtchouk (n j w : ℕ) : ℤ :=
  ∑ h ∈ range (j + 1), (-1 : ℤ) ^ h * (w.choose h : ℤ) * ((n - w).choose (j - h) : ℤ)

/-- Splitting the weight-`j` layer by how much it meets `v`: the part inside
`v` and the part outside are chosen independently. -/
lemma card_meet_fiber {n : ℕ} (v : Finset (Fin n)) {j h : ℕ} (hh : h ≤ j) :
    ((((univ : Finset (Fin n)).powersetCard j)).filter
        (fun x => (x ∩ v).card = h)).card
      = v.card.choose h * (n - v.card).choose (j - h) := by
  classical
  have hcompl : (vᶜ : Finset (Fin n)).card = n - v.card := by
    rw [Finset.card_compl, Fintype.card_fin]
  have hbij : ((((univ : Finset (Fin n)).powersetCard j)).filter
        (fun x => (x ∩ v).card = h)).card
      = ((v.powersetCard h) ×ˢ (vᶜ.powersetCard (j - h))).card := by
    refine Finset.card_nbij' (fun x => (x ∩ v, x \ v)) (fun p => p.1 ∪ p.2)
      ?_ ?_ ?_ ?_
    · intro x hx
      simp only [Finset.mem_coe, Finset.mem_filter, Finset.mem_powersetCard] at hx
      obtain ⟨⟨-, hcard⟩, hmeet⟩ := hx
      have hsplit : (x ∩ v).card + (x \ v).card = x.card :=
        Finset.card_inter_add_card_sdiff x v
      simp only [Finset.mem_coe]
      refine Finset.mem_product.mpr ⟨?_, ?_⟩
      · exact Finset.mem_powersetCard.mpr ⟨Finset.inter_subset_right, hmeet⟩
      · refine Finset.mem_powersetCard.mpr ⟨?_, ?_⟩
        · intro a ha
          rw [Finset.mem_sdiff] at ha
          simpa using ha.2
        · show (x \ v).card = j - h
          omega
    · intro p hp
      simp only [Finset.mem_coe, Finset.mem_product, Finset.mem_powersetCard] at hp
      obtain ⟨⟨hp1, hc1⟩, ⟨hp2, hc2⟩⟩ := hp
      have hdisj : Disjoint p.1 p.2 := by
        refine Finset.disjoint_left.mpr fun a ha1 ha2 => ?_
        have := hp2 ha2
        rw [Finset.mem_compl] at this
        exact this (hp1 ha1)
      have hunion : (p.1 ∪ p.2).card = j := by
        rw [Finset.card_union_of_disjoint hdisj, hc1, hc2]; omega
      have hmeet : ((p.1 ∪ p.2) ∩ v) = p.1 := by
        rw [Finset.union_inter_distrib_right,
          Finset.inter_eq_left.mpr hp1]
        have : p.2 ∩ v = ∅ := by
          refine Finset.eq_empty_of_forall_notMem fun a ha => ?_
          rw [Finset.mem_inter] at ha
          have := hp2 ha.1
          rw [Finset.mem_compl] at this
          exact this ha.2
        rw [this, Finset.union_empty]
      simp only [Finset.mem_coe, Finset.mem_filter, Finset.mem_powersetCard]
      exact ⟨⟨Finset.subset_univ _, hunion⟩, by rw [hmeet, hc1]⟩
    · intro x _
      show (x ∩ v) ∪ (x \ v) = x
      rw [Finset.union_comm]
      exact Finset.sdiff_union_inter x v
    · intro p hp
      simp only [Finset.mem_coe, Finset.mem_product, Finset.mem_powersetCard] at hp
      obtain ⟨⟨hp1, -⟩, ⟨hp2, -⟩⟩ := hp
      have hp2v : ∀ a ∈ p.2, a ∉ v := by
        intro a ha
        have := hp2 ha
        rwa [Finset.mem_compl] at this
      have h1 : (p.1 ∪ p.2) ∩ v = p.1 := by
        ext a
        simp only [Finset.mem_inter, Finset.mem_union]
        exact ⟨fun ⟨hu, hv⟩ => hu.resolve_right (fun h2 => hp2v a h2 hv),
               fun ha => ⟨Or.inl ha, hp1 ha⟩⟩
      have h2 : (p.1 ∪ p.2) \ v = p.2 := by
        ext a
        simp only [Finset.mem_sdiff, Finset.mem_union]
        exact ⟨fun ⟨hu, hv⟩ => hu.resolve_left (fun h1' => hv (hp1 h1')),
               fun ha => ⟨Or.inr ha, hp2v a ha⟩⟩
      show ((p.1 ∪ p.2) ∩ v, (p.1 ∪ p.2) \ v) = p
      rw [h1, h2]
  rw [hbij, Finset.card_product, Finset.card_powersetCard, Finset.card_powersetCard,
    hcompl]

/-- **The character sum over a weight layer is a Krawtchouk value.**  In
particular it depends on `v` only through `|v|` — which is what makes the
MacWilliams-style identities well posed. -/
theorem krawtchouk_char_sum {n : ℕ} (v : Finset (Fin n)) (j : ℕ) :
    ∑ x ∈ (univ : Finset (Fin n)).powersetCard j, (-1 : ℤ) ^ (x ∩ v).card
      = krawtchouk n j v.card := by
  classical
  have hmaps : ∀ x ∈ (univ : Finset (Fin n)).powersetCard j,
      (x ∩ v).card ∈ range (j + 1) := by
    intro x hx
    rw [Finset.mem_powersetCard] at hx
    rw [Finset.mem_range]
    have := Finset.card_le_card (Finset.inter_subset_left (s₁ := x) (s₂ := v))
    omega
  rw [← Finset.sum_fiberwise_of_maps_to hmaps]
  refine Finset.sum_congr rfl fun h hh => ?_
  rw [Finset.mem_range] at hh
  have hconst : ∀ x ∈ ((univ : Finset (Fin n)).powersetCard j).filter
      (fun x => (x ∩ v).card = h), (-1 : ℤ) ^ (x ∩ v).card = (-1 : ℤ) ^ h := by
    intro x hx
    rw [(Finset.mem_filter.mp hx).2]
  rw [Finset.sum_congr rfl hconst, Finset.sum_const, nsmul_eq_mul,
    card_meet_fiber v (by omega : h ≤ j)]
  push_cast
  ring

/-- `K_j(0) = C(n, j)`: the empty character is trivial. -/
lemma krawtchouk_zero {n j : ℕ} (hj : j ≤ n) : krawtchouk n j 0 = (n.choose j : ℤ) := by
  unfold krawtchouk
  rw [Finset.sum_eq_single 0]
  · simp
  · intro h _ hne
    simp [Nat.choose_eq_zero_of_lt (Nat.pos_of_ne_zero hne)]
  · intro hcontra
    simp at hcontra

end Spin
