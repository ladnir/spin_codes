import SpinCodes.Structured.PolyIdentity
import SpinCodes.Structured.PolyPackedDefs

namespace Spin.Structured

theorem coefficientMass_add (p q : List ℤ) :
    coefficientMass (polyAdd p q) ≤ coefficientMass p + coefficientMass q := by
  induction p generalizing q with
  | nil => simp [polyAdd, coefficientMass]
  | cons a p ih =>
    cases q with
    | nil => simp [polyAdd, coefficientMass]
    | cons b q =>
      have ht := ih q
      have ha := Int.natAbs_add_le a b
      simp only [polyAdd, coefficientMass]
      omega

theorem coefficientMass_scale (a : ℤ) (p : List ℤ) :
    coefficientMass (polyScale a p) = a.natAbs * coefficientMass p := by
  induction p with
  | nil => simp [polyScale, coefficientMass]
  | cons b p ih => simp [polyScale, coefficientMass, ih, Int.natAbs_mul]; ring

theorem coefficientMass_mul (p q : List ℤ) :
    coefficientMass (polyMul p q) ≤ coefficientMass p * coefficientMass q := by
  induction p with
  | nil => simp [polyMul, coefficientMass]
  | cons a p ih =>
    have hh := coefficientMass_add (polyScale a q) (0 :: polyMul p q)
    simp only [coefficientMass_scale, coefficientMass, Int.natAbs_zero, zero_add] at hh
    simp only [polyMul, coefficientMass]
    nlinarith

theorem encodePoly_add (p q : List ℤ) (b : ℤ) :
    encodePoly (polyAdd p q) b = encodePoly p b + encodePoly q b := by
  induction p generalizing q with
  | nil => simp [polyAdd, encodePoly]
  | cons a p ih =>
    cases q with
    | nil => simp [polyAdd, encodePoly]
    | cons c q => simp [polyAdd, encodePoly, ih]; ring

theorem encodePoly_scale (a : ℤ) (p : List ℤ) (b : ℤ) :
    encodePoly (polyScale a p) b = a * encodePoly p b := by
  induction p with
  | nil => simp [polyScale, encodePoly]
  | cons c p ih => simp [polyScale, encodePoly, ih]; ring

theorem encodePoly_mul (p q : List ℤ) (b : ℤ) :
    encodePoly (polyMul p q) b = encodePoly p b * encodePoly q b := by
  induction p with
  | nil => simp [polyMul, encodePoly]
  | cons a p ih =>
    simp [polyMul, encodePoly_add, encodePoly_scale, encodePoly, ih]
    ring

/-- A zero integer encoding with a sufficiently large radix is the zero
polynomial. The bound excludes every possible carry between coefficients. -/
theorem listEval_zero_of_encode (p : List ℤ) {B : ℕ}
    (hB : coefficientMass p < B) (he : encodePoly p B = 0) (x : ℝ) :
    listEval p x = 0 := by
  induction p with
  | nil => rfl
  | cons a p ih =>
    have haB : a.natAbs < B := by simp only [coefficientMass] at hB; omega
    have hd : (B : ℤ) ∣ a := by
      refine ⟨-encodePoly p B, ?_⟩
      simp only [encodePoly] at he
      linarith
    have ha : a = 0 := Int.eq_zero_of_dvd_of_natAbs_lt_natAbs hd (by simpa using haB)
    have hn : (B : ℤ) ≠ 0 := by exact_mod_cast (by omega : B ≠ 0)
    have ht : encodePoly p B = 0 := by
      simp only [encodePoly, ha, zero_add] at he
      exact (mul_eq_zero.mp he).resolve_left hn
    have hp : coefficientMass p < B := by simp only [coefficientMass] at hB; omega
    simp [listEval, ha, ih hp ht]

theorem packedCoefficients_mass (ts : List PackedTerm) :
    coefficientMass (packedCoefficients ts) ≤ packedMass ts := by
  induction ts with
  | nil => simp [packedCoefficients, coefficientMass, packedMass]
  | cons t ts ih =>
    have ha := coefficientMass_add t.coefficients (packedCoefficients ts)
    have hm := coefficientMass_mul t.left.num t.right.num
    have hscale := Nat.mul_le_mul_left t.factor hm
    simp only [PackedTerm.coefficients, coefficientMass_scale, Int.natAbs_natCast] at ha
    simp only [packedCoefficients, packedMass, PackedTerm.massBound, PackedTerm.coefficients]
    nlinarith

theorem encode_packedCoefficients (ts : List PackedTerm) (b : ℤ) :
    encodePoly (packedCoefficients ts) b = packedEncoded ts b := by
  induction ts with
  | nil => rfl
  | cons t ts ih =>
    simp only [packedCoefficients, packedEncoded, encodePoly_add, PackedTerm.coefficients,
      encodePoly_scale, encodePoly_mul, PackedTerm.encoded, ih]
    ring

noncomputable def packedValue : List PackedTerm → ℝ → ℝ
  | [], _ => 0
  | t :: ts, x => t.left.eval x * t.right.eval x + packedValue ts x

theorem packedCoefficients_eval (D : ℕ) (ts : List PackedTerm)
    (hv : ∀ t ∈ ts, t.valid D = true) (x : ℝ) :
    listEval (packedCoefficients ts) x = (D : ℝ) * packedValue ts x := by
  induction ts with
  | nil => simp [packedCoefficients, packedValue, listEval]
  | cons t ts ih =>
    have ht := hv t (List.mem_cons_self ..)
    simp only [PackedTerm.valid, Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq] at ht
    obtain ⟨⟨hl, hr⟩, hfactor⟩ := ht
    have hl' : (t.left.den : ℝ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hl)
    have hr' : (t.right.den : ℝ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hr)
    have hf : (t.left.den : ℝ) * t.right.den * t.factor = D := by exact_mod_cast hfactor
    have htval : listEval t.coefficients x = (D : ℝ) * (t.left.eval x * t.right.eval x) := by
      simp only [PackedTerm.coefficients, listEval_polyScale, listEval_polyMul,
        Int.cast_natCast, RatPoly.eval]
      rw [← hf]
      field_simp
    rw [packedCoefficients, listEval_polyAdd, htval,
      ih (fun t ht => hv t (List.mem_cons_of_mem _ ht))]
    simp only [packedValue]
    ring

theorem checkSumProducts_sound (D : ℕ) (ts : List PackedTerm) (r : RatPoly) (rf : ℕ)
    (h : checkSumProducts D ts r rf = true) (x : ℝ) : r.eval x = packedValue ts x := by
  simp only [checkSumProducts, Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq,
    List.all_eq_true] at h
  obtain ⟨⟨⟨⟨hD, hr⟩, hf⟩, hv⟩, he⟩ := h
  let c := polyAdd (packedCoefficients ts) (polyScale (-(rf : ℤ)) r.num)
  have hm : coefficientMass c < packedBase ts r rf := by
    have ha := coefficientMass_add (packedCoefficients ts) (polyScale (-(rf : ℤ)) r.num)
    have hb := packedCoefficients_mass ts
    simp only [coefficientMass_scale, Int.natAbs_neg, Int.natAbs_natCast] at ha
    dsimp [c, packedBase]
    omega
  have hec : encodePoly c (packedBase ts r rf) = 0 := by
    simp only [c, encodePoly_add, encode_packedCoefficients, encodePoly_scale]
    rw [he]
    ring
  have hz := listEval_zero_of_encode c hm hec x
  simp only [c, listEval_polyAdd, packedCoefficients_eval D ts hv,
    listEval_polyScale, Int.cast_neg, Int.cast_natCast] at hz
  have hf' : (r.den : ℝ) * rf = D := by exact_mod_cast hf
  have hr' : (r.den : ℝ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hr)
  have hD' : (0 : ℝ) < D := by exact_mod_cast hD
  have hh : (rf : ℝ) * listEval r.num x = (D : ℝ) * r.eval x := by
    rw [← hf', RatPoly.eval]
    field_simp
  rw [neg_mul, hh] at hz
  nlinarith

end Spin.Structured

