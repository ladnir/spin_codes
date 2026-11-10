import SpinCodes.UniformMixer.Map22PacketProfileProduct

/-! Packet-profile shapes determined by the low-variable truth table. -/
noncomputable section
namespace Spin.UniformMixer
open Finset

def profileOddShape (w : ℕ) : Fin 5 → ℕ := ![0, 16-w, 0, w, 0]
def profileEvenShape (size w : ℕ) : Fin 5 → ℕ := ![size-w, 0, 16-size, 0, w]

private theorem bit_ne_zero_iff_one (b : ZMod 2) : b ≠ 0 ↔ b = 1 := by
  have h : ∀ b : ZMod 2, b ≠ 0 ↔ b = 1 := by decide
  exact h b

theorem profile_zero_count (Z : Finset (Fin 16)) (b : Fin 16 → ZMod 2) :
    (Z.filter (fun v => b v = 0)).card =
      Z.card - (Z.filter (fun v => b v = 1)).card := by
  have h := Finset.card_filter_add_card_filter_not (s := Z) (fun v => b v = 0)
  simp only [bit_ne_zero_iff_one] at h
  omega

theorem profileOfNibbles_odd_shape
    (f : Fin 16 × Fin 4 → ZMod 2) (b : Fin 16 → ZMod 2)
    (h : ∀ v, profileNibbleWordWeight f v = if b v = 0 then 1 else 3) :
    profileOfNibbles f = profileOddShape
      ((Finset.univ.filter (fun v => b v = 1)).card) := by
  have hz := profile_zero_count Finset.univ b
  simp only [Finset.card_univ, Fintype.card_fin] at hz
  funext r
  fin_cases r <;>
    simp [profileOfNibbles, h, profileOddShape, ite_eq_iff,
      bit_ne_zero_iff_one, hz]
  simpa only [← bit_ne_zero_iff_one, not_not] using hz

theorem profileOfNibbles_even_shape
    (f : Fin 16 × Fin 4 → ZMod 2) (b : Fin 16 → ZMod 2)
    (Z : Finset (Fin 16))
    (h : ∀ v, profileNibbleWordWeight f v =
      if v ∈ Z then (if b v = 0 then 0 else 4) else 2) :
    profileOfNibbles f = profileEvenShape Z.card
      ((Z.filter (fun v => b v = 1)).card) := by
  have hz := profile_zero_count Z b
  have hc : (Finset.univ.filter (fun v : Fin 16 => v ∉ Z)).card = 16-Z.card := by
    have he : Finset.univ.filter (fun v : Fin 16 => v ∉ Z) = Finset.univ \ Z := by
      ext v; simp
    rw [he, Finset.card_sdiff]
    simp
  have hfilter (P : Fin 16 → Prop) [DecidablePred P] :
      Finset.univ.filter (fun v => v ∈ Z ∧ P v) = Z.filter P := by
    ext v; simp
  funext r
  fin_cases r <;>
    simp [profileOfNibbles, h, profileEvenShape, ite_eq_iff,
      bit_ne_zero_iff_one, hfilter, hz, hc]
  simpa only [← bit_ne_zero_iff_one, not_not] using hz

/-- Product decomposition gives this profile without examining any state. -/
theorem profileSplit6_odd_shape (b c : ProfileAffine4) (q : ProfileQuad4) :
    profileOfNibbles (profileSplit6Eval (1,b,c,q)) =
      profileOddShape (profileQuad4Weight (q+profileAffineProduct b c)) := by
  have he := profileOfNibbles_odd_shape _ _ (profileOddNibbleWeight b c q)
  have hw : (univ.filter (fun v => profileQuad4Eval
      (q+profileAffineProduct b c) v = 1)).card =
      profileQuad4Weight (q+profileAffineProduct b c) := by
    simp only [profileQuad4Weight, Finset.card_filter]
  exact he.trans (congrArg profileOddShape hw)

#print axioms profileOfNibbles_odd_shape
#print axioms profileOfNibbles_even_shape
#print axioms profileSplit6_odd_shape
end Spin.UniformMixer
