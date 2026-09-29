import SpinCodes.Structured.ConcreteFourierProduct

noncomputable section
namespace Spin.Structured.ConcreteFourier
open Finset

theorem overlap_monomial_mono {a b : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) (hab : a ≤ b)
    {w h k : ℕ} (hh : h ≤ k) (hk : k ≤ w) : a^(w-h)*b^h ≤ a^(w-k)*b^k := by
  have he : w-h = (w-k)+(k-h) := by omega
  rw [he, pow_add]
  calc
    _ ≤ a^(w-k)*b^(k-h)*b^h := by
      apply mul_le_mul_of_nonneg_right _ (pow_nonneg hb h)
      exact mul_le_mul_of_nonneg_left (pow_le_pow_left₀ ha hab _) (pow_nonneg ha _)
    _ = a^(w-k)*b^k := by rw [mul_assoc, ← pow_add, Nat.sub_add_cancel hh]

theorem overlap_monomial_le_endpoints {a b : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b)
    {w lo hi h : ℕ} (hlo : lo ≤ h) (hhi : h ≤ hi) (hiw : hi ≤ w) :
    a^(w-h)*b^h ≤ max (a^(w-lo)*b^lo) (a^(w-hi)*b^hi) := by
  rcases le_total a b with hab | hba
  · exact (overlap_monomial_mono ha hb hab hhi hiw).trans (le_max_right _ _)
  · have hh : w-h ≤ w-lo := by omega
    have hw : w-lo ≤ w := Nat.sub_le _ _
    have hc := overlap_monomial_mono hb ha hba hh hw
    have hh2 : w-(w-h) = h := by omega
    have hl2 : w-(w-lo) = lo := by omega
    rw [hh2, hl2] at hc
    have hc' : a^(w-h)*b^h ≤ a^(w-lo)*b^lo := by simpa only [mul_comm] using hc
    exact hc'.trans (le_max_left _ _)

theorem overlap_product {n : ℕ} (Y v : Finset (Fin n)) (a b : ℝ) :
    (∏ i ∈ v, if i ∈ Y then b else a) =
      a^(v.card-(v ∩ Y).card)*b^(v ∩ Y).card := by
  rw [prod_ite]
  have hy : v.filter (fun i => i ∈ Y) = v ∩ Y := by ext; simp
  have hn : v.filter (fun i => i ∉ Y) = v \ Y := by ext; simp
  rw [hy, hn]
  simp only [prod_const, card_sdiff]
  rw [inter_comm Y v]
  ring

def overlapCap (n d w : ℕ) (a b : ℝ) : ℝ :=
  max (a^(w-(d+w-n))*b^(d+w-n)) (a^(w-min d w)*b^(min d w))

/-- The Fourier overlap maximum is attained at an endpoint, including zero bases. -/
theorem overlap_product_le {n : ℕ} (Y v : Finset (Fin n)) {a b : ℝ}
    (ha : 0 ≤ a) (hb : 0 ≤ b) :
    (∏ i ∈ v, if i ∈ Y then b else a) ≤ overlapCap n Y.card v.card a b := by
  rw [overlap_product]
  apply overlap_monomial_le_endpoints ha hb
  · have h := Finset.card_union_add_card_inter Y v
    have hc : (Y ∪ v).card ≤ n := by simpa using Finset.card_le_univ (Y ∪ v)
    rw [inter_comm Y v] at h
    omega
  · exact le_min (card_le_card inter_subset_right) (card_le_card inter_subset_left)
  · exact min_le_right _ _

end Spin.Structured.ConcreteFourier
