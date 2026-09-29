/-
Ingredients of `lem:structured-route-domination`.

The lemma bounds the structured route's law by an i.i.d. Bernoulli law at a
multiplicative cost `(b+1)^Q (L+1)^b`, and its proof rests on three facts:

1. conditioning `b` Bernoulli variables on their sum hitting a *mode* costs at
   most `b+1`;
2. the method-of-types lower bound
   `P[Bin(L,q) = k] ≥ exp(-L·D(k/L‖q)) / (L+1)`;
3. the Chernoff upper bound for a Poisson-binomial sum,
   `P[S = k] ≤ exp(-L·D(k/L‖q))`.

The first two both reduce to one observation — a mass function on `n` points
puts at least `1/n` on its largest atom — which is what this file starts from.
-/
import SpinCodes.Prob

set_option linter.unusedSectionVars false

namespace Spin

open Finset

variable {Ω : Type*} [Fintype Ω]

/-! ## The mode of a finite mass function

Both the conditioning cost and the method-of-types bound are this pigeonhole:
the atoms sum to one, so the largest is at least the average. -/

/-- **A mode carries at least `1/n` of the mass.** -/
theorem le_mode (P : FinPMF Ω) {m : Ω} (hm : ∀ ω, P.p ω ≤ P.p m) :
    1 / (Fintype.card Ω : ℝ) ≤ P.p m := by
  have hcard : (0 : ℝ) < (Fintype.card Ω : ℝ) := by
    have : 0 < Fintype.card Ω := Fintype.card_pos_iff.mpr ⟨m⟩
    exact_mod_cast this
  have hsum : (1 : ℝ) = ∑ ω, P.p ω := P.total.symm
  have hle : ∑ ω : Ω, P.p ω ≤ ∑ _ω : Ω, P.p m := Finset.sum_le_sum fun ω _ => hm ω
  rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul] at hle
  rw [div_le_iff₀ hcard]
  linarith [hsum, hle]

/-- The same, as the bound on the conditioning cost the lemma actually uses:
dividing by the mode's mass multiplies by at most `n`. -/
theorem inv_mode_le (P : FinPMF Ω) {m : Ω} (hm : ∀ ω, P.p ω ≤ P.p m)
    (hpos : 0 < P.p m) : 1 / P.p m ≤ (Fintype.card Ω : ℝ) := by
  have h := le_mode P hm
  have hcard : (0 : ℝ) < (Fintype.card Ω : ℝ) := by
    have : 0 < Fintype.card Ω := Fintype.card_pos_iff.mpr ⟨m⟩
    exact_mod_cast this
  rw [div_le_iff₀ hpos]
  rw [div_le_iff₀ hcard] at h
  linarith


/-! ## The binomial mass function

Carried as a `FinPMF` on `Fin (L+1)` so that `le_mode` applies to it directly.
Totality is the binomial theorem. -/

/-- `Bin(L, u)` as a mass function on `Fin (L+1)`. -/
noncomputable def binPMF (L : ℕ) (u : ℝ) (hu0 : 0 ≤ u) (hu1 : u ≤ 1) :
    FinPMF (Fin (L + 1)) where
  p := fun j => (L.choose (j : ℕ) : ℝ) * u ^ (j : ℕ) * (1 - u) ^ (L - (j : ℕ))
  nonneg := fun j => by
    have h1 : (0 : ℝ) ≤ 1 - u := by linarith
    positivity
  total := by
    rw [Fin.sum_univ_eq_sum_range
      (fun j => (L.choose j : ℝ) * u ^ j * (1 - u) ^ (L - j)) (L + 1)]
    have h := add_pow u (1 - u) L
    have hone : u + (1 - u) = 1 := by ring
    rw [hone, one_pow] at h
    refine Eq.trans ?_ h.symm
    exact Finset.sum_congr rfl fun j _ => by ring

@[simp] lemma binPMF_apply (L : ℕ) (u : ℝ) (hu0 : 0 ≤ u) (hu1 : u ≤ 1)
    (j : Fin (L + 1)) :
    (binPMF L u hu0 hu1).p j
      = (L.choose (j : ℕ) : ℝ) * u ^ (j : ℕ) * (1 - u) ^ (L - (j : ℕ)) := rfl

/-- **The method-of-types lower bound**, in the form the route lemma uses:
at `u = k/L` the binomial mass at `k` is at least `1/(L+1)`.

Dividing through by `u^k (1-u)^(L-k)` turns this into
`C(L,k) ≥ exp(L·H(k/L)) / (L+1)`, which is the usual statement. -/
theorem types_lower_bound {L : ℕ} {u : ℝ} (hu0 : 0 ≤ u) (hu1 : u ≤ 1)
    {k : Fin (L + 1)} (hmode : ∀ j : Fin (L + 1),
      (binPMF L u hu0 hu1).p j ≤ (binPMF L u hu0 hu1).p k) :
    1 / ((L : ℝ) + 1)
      ≤ (L.choose (k : ℕ) : ℝ) * u ^ (k : ℕ) * (1 - u) ^ (L - (k : ℕ)) := by
  have h := le_mode (binPMF L u hu0 hu1) hmode
  rwa [Fintype.card_fin, Nat.cast_add, Nat.cast_one] at h


/-! ## The mode of `Bin(L, k/L)`

The empirical frequency is its own most likely outcome.  Clearing denominators
by `L^L` leaves a statement purely about naturals,

    C(L,j) · k^j · (L-k)^(L-j)  ≤  C(L,k) · k^k · (L-k)^(L-k)

which needs neither the ratio test nor Stirling — only two factorial-versus-
power bounds, and a symmetry that handles `j ≥ k` by exchanging the roles of
`k` and `L-k`. -/

/-- `A! ≤ j! · A^(A-j)`: the factors `j+1, …, A` are each at most `A`. -/
theorem factorial_le_factorial_mul_pow {j A : ℕ} (h : j ≤ A) :
    Nat.factorial A ≤ Nat.factorial j * A ^ (A - j) := by
  obtain ⟨d, rfl⟩ := Nat.exists_eq_add_of_le h
  simp only [Nat.add_sub_cancel_left]
  induction d with
  | zero => simp
  | succ d ih =>
      have hstep : Nat.factorial (j + (d + 1)) = (j + d + 1) * Nat.factorial (j + d) := by
        rw [show j + (d + 1) = (j + d) + 1 by ring, Nat.factorial_succ]
      have hmono : (j + d) ^ d ≤ (j + (d + 1)) ^ d := Nat.pow_le_pow_left (by omega) d
      calc Nat.factorial (j + (d + 1)) = (j + d + 1) * Nat.factorial (j + d) := hstep
        _ ≤ (j + d + 1) * (Nat.factorial j * (j + d) ^ d) := Nat.mul_le_mul_left _ (ih (by omega))
        _ ≤ (j + d + 1) * (Nat.factorial j * (j + (d + 1)) ^ d) :=
            Nat.mul_le_mul_left _ (Nat.mul_le_mul_left _ hmono)
        _ = Nat.factorial j * (j + (d + 1)) ^ (d + 1) := by
            rw [show j + (d + 1) = j + d + 1 by ring]; ring

/-- The half of the mode argument that does the work, stated symmetrically in
the two "sides" `A` and `B` so that one lemma covers `j ≤ k` and `j ≥ k`. -/
theorem mode_half {A B j : ℕ} (h : j ≤ A) :
    Nat.factorial A * Nat.factorial B * (A ^ j * B ^ (A + B - j))
      ≤ Nat.factorial j * Nat.factorial (A + B - j) * (A ^ A * B ^ B) := by
  obtain ⟨d, rfl⟩ := Nat.exists_eq_add_of_le h
  have hAB : j + d + B - j = B + d := by omega
  have hfac : Nat.factorial (j + d) ≤ Nat.factorial j * (j + d) ^ d := by
    have := factorial_le_factorial_mul_pow (j := j) (A := j + d) (by omega)
    simpa using this
  have hBd : Nat.factorial B * B ^ d ≤ Nat.factorial (B + d) :=
    le_trans (Nat.mul_le_mul_left _ (Nat.pow_le_pow_left (by omega) d))
      Nat.factorial_mul_pow_le_factorial
  have hsplit : (j + d) ^ (j + d) = (j + d) ^ d * (j + d) ^ j := by
    rw [← pow_add]; ring_nf
  rw [hAB]
  calc Nat.factorial (j + d) * Nat.factorial B * ((j + d) ^ j * B ^ (B + d))
      = (Nat.factorial (j + d) * (Nat.factorial B * B ^ d)) * ((j + d) ^ j * B ^ B) := by rw [pow_add]; ring
    _ ≤ (Nat.factorial j * (j + d) ^ d * Nat.factorial (B + d)) * ((j + d) ^ j * B ^ B) :=
        Nat.mul_le_mul_right _ (Nat.mul_le_mul hfac hBd)
    _ = Nat.factorial j * Nat.factorial (B + d) * ((j + d) ^ d * (j + d) ^ j * B ^ B) := by ring
    _ = Nat.factorial j * Nat.factorial (B + d) * ((j + d) ^ (j + d) * B ^ B) := by rw [hsplit]

/-- **`k` is a mode of `Bin(L, k/L)`**, with denominators cleared. -/
theorem binPMF_mode_nat {L k j : ℕ} (hk : k ≤ L) (hj : j ≤ L) :
    L.choose j * (k ^ j * (L - k) ^ (L - j))
      ≤ L.choose k * (k ^ k * (L - k) ^ (L - k)) := by
  have hkL : k + (L - k) = L := by omega
  -- the factorial-level statement, from `mode_half` in whichever orientation
  have hcore : Nat.factorial k * Nat.factorial (L - k) * (k ^ j * (L - k) ^ (L - j))
      ≤ Nat.factorial j * Nat.factorial (L - j) * (k ^ k * (L - k) ^ (L - k)) := by
    rcases le_total j k with hjk | hjk
    · have := mode_half (A := k) (B := L - k) (j := j) hjk
      rwa [hkL] at this
    · have := mode_half (A := L - k) (B := k) (j := L - j) (by omega)
      rw [show L - k + k = L by omega, show L - (L - j) = j by omega] at this
      calc Nat.factorial k * Nat.factorial (L - k) * (k ^ j * (L - k) ^ (L - j))
          = Nat.factorial (L - k) * Nat.factorial k * ((L - k) ^ (L - j) * k ^ j) := by ring
        _ ≤ Nat.factorial (L - j) * Nat.factorial j * ((L - k) ^ (L - k) * k ^ k) := this
        _ = Nat.factorial j * Nat.factorial (L - j) * (k ^ k * (L - k) ^ (L - k)) := by ring
  -- multiply up by the two factorial products and cancel
  have hc1 : L.choose j * (Nat.factorial j * Nat.factorial (L - j)) = Nat.factorial L := by
    rw [← Nat.choose_mul_factorial_mul_factorial hj]; ring
  have hc2 : L.choose k * (Nat.factorial k * Nat.factorial (L - k)) = Nat.factorial L := by
    rw [← Nat.choose_mul_factorial_mul_factorial hk]; ring
  have hpos : 0 < (Nat.factorial j * Nat.factorial (L - j)) * (Nat.factorial k * Nat.factorial (L - k)) := by positivity
  refine Nat.le_of_mul_le_mul_left ?_ hpos
  calc (Nat.factorial j * Nat.factorial (L - j)) * (Nat.factorial k * Nat.factorial (L - k))
        * (L.choose j * (k ^ j * (L - k) ^ (L - j)))
      = (L.choose j * (Nat.factorial j * Nat.factorial (L - j)))
        * (Nat.factorial k * Nat.factorial (L - k) * (k ^ j * (L - k) ^ (L - j))) := by ring
    _ = Nat.factorial L * (Nat.factorial k * Nat.factorial (L - k) * (k ^ j * (L - k) ^ (L - j))) := by rw [hc1]
    _ ≤ Nat.factorial L * (Nat.factorial j * Nat.factorial (L - j) * (k ^ k * (L - k) ^ (L - k))) :=
        Nat.mul_le_mul_left _ hcore
    _ = (L.choose k * (Nat.factorial k * Nat.factorial (L - k)))
        * (Nat.factorial j * Nat.factorial (L - j) * (k ^ k * (L - k) ^ (L - k))) := by rw [hc2]
    _ = (Nat.factorial j * Nat.factorial (L - j)) * (Nat.factorial k * Nat.factorial (L - k))
        * (L.choose k * (k ^ k * (L - k) ^ (L - k))) := by ring


/-! ## From the integer statement to the mass function

Dividing `binPMF_mode_nat` by `L^L` is the whole content; the care is only in
the `ℕ`-subtraction casts. -/

private lemma pow_div_split {L j : ℕ} (hj : j ≤ L) (a b c : ℝ) :
    (a / c) ^ j * (b / c) ^ (L - j) = (a ^ j * b ^ (L - j)) / c ^ L := by
  rw [div_pow, div_pow, div_mul_div_comm, ← pow_add]
  congr 2
  omega

/-- **`k` is a mode of `Bin(L, k/L)`.** -/
theorem binPMF_mode {L k : ℕ} (hL : 0 < L) (hk : k ≤ L)
    (h0 : (0 : ℝ) ≤ (k : ℝ) / (L : ℝ)) (h1 : (k : ℝ) / (L : ℝ) ≤ 1)
    (j : Fin (L + 1)) :
    (binPMF L ((k : ℝ) / (L : ℝ)) h0 h1).p j
      ≤ (binPMF L ((k : ℝ) / (L : ℝ)) h0 h1).p ⟨k, by omega⟩ := by
  have hLR : (0 : ℝ) < (L : ℝ) := by exact_mod_cast hL
  have hLpow : (0 : ℝ) < (L : ℝ) ^ L := by positivity
  have hcompl : 1 - (k : ℝ) / (L : ℝ) = ((L - k : ℕ) : ℝ) / (L : ℝ) := by
    rw [Nat.cast_sub hk]
    field_simp
  have hjL : (j : ℕ) ≤ L := by omega
  -- the integer statement, cast to `ℝ`
  have hnat := binPMF_mode_nat (L := L) (k := k) (j := (j : ℕ)) hk hjL
  have hcast : ((L.choose (j : ℕ) : ℕ) : ℝ)
        * (((k : ℕ) : ℝ) ^ (j : ℕ) * (((L - k : ℕ) : ℝ)) ^ (L - (j : ℕ)))
      ≤ ((L.choose k : ℕ) : ℝ)
        * (((k : ℕ) : ℝ) ^ k * (((L - k : ℕ) : ℝ)) ^ (L - k)) := by
    exact_mod_cast hnat
  simp only [binPMF_apply, hcompl]
  rw [mul_assoc, mul_assoc, pow_div_split hjL, pow_div_split hk,
    ← mul_div_assoc, ← mul_div_assoc, div_le_div_iff_of_pos_right hLpow]
  exact hcast

/-- The method-of-types lower bound with its hypothesis discharged. -/
theorem types_lower_bound_at {L k : ℕ} (hL : 0 < L) (hk : k ≤ L)
    (h0 : (0 : ℝ) ≤ (k : ℝ) / (L : ℝ)) (h1 : (k : ℝ) / (L : ℝ) ≤ 1) :
    1 / ((L : ℝ) + 1)
      ≤ (L.choose k : ℝ) * ((k : ℝ) / (L : ℝ)) ^ k
        * (1 - (k : ℝ) / (L : ℝ)) ^ (L - k) :=
  types_lower_bound h0 h1 (k := ⟨k, by omega⟩) (binPMF_mode hL hk h0 h1)


/-! ## Poisson-binomial domination

The moment generating function of a sum of independent Bernoulli variables is
dominated by that of the binomial with the same mean:

    ∏ (1 + pᵢ·c)  ≤  (1 + q·c)^L,      q = (∑ pᵢ)/L

for `c ≥ -1`, which at `c = e^t - 1` is the statement that replacing the `pᵢ`
by their average can only increase the MGF.  This is the step that lets the
Chernoff bound for a Poisson binomial be read off the binomial one, and it is
the piece Mathlib does not have.

It is AM-GM with uniform weights, `Real.geom_mean_le_arith_mean_weighted`. -/

/-- **AM-GM with uniform weights**: a product is at most the `L`-th power of
the arithmetic mean. -/
theorem prod_le_pow_avg {L : ℕ} (hL : 0 < L) (z : Fin L → ℝ) (hz : ∀ i, 0 ≤ z i) :
    ∏ i, z i ≤ ((∑ i, z i) / (L : ℝ)) ^ L := by
  have hLR : (0 : ℝ) < (L : ℝ) := by exact_mod_cast hL
  set P : ℝ := ∏ i, z i with hP
  have hPnn : 0 ≤ P := Finset.prod_nonneg fun i _ => hz i
  have hw : ∀ i ∈ (Finset.univ : Finset (Fin L)), (0 : ℝ) ≤ (L : ℝ)⁻¹ := by
    intro i _; positivity
  have hw' : ∑ _i : Fin L, ((L : ℝ)⁻¹) = 1 := by
    rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
    field_simp
  have hgm := Real.geom_mean_le_arith_mean_weighted Finset.univ
    (fun _ => (L : ℝ)⁻¹) z hw hw' (fun i _ => hz i)
  have hprod : ∏ i : Fin L, z i ^ ((L : ℝ)⁻¹) = P ^ ((L : ℝ)⁻¹) := by
    rw [hP, ← Real.finsetProd_rpow _ _ (fun i _ => hz i)]
  have harith : ∑ i : Fin L, (L : ℝ)⁻¹ * z i = (∑ i, z i) / (L : ℝ) := by
    rw [← Finset.mul_sum]; ring
  rw [hprod, harith] at hgm
  have hpow := pow_le_pow_left₀ (Real.rpow_nonneg hPnn _) hgm L
  rwa [← Real.rpow_natCast (P ^ ((L : ℝ)⁻¹)) L, ← Real.rpow_mul hPnn,
    inv_mul_cancel₀ (ne_of_gt hLR), Real.rpow_one] at hpow

/-- The Poisson-binomial MGF domination, the specialisation of
`prod_le_pow_avg` at `zᵢ = 1 + pᵢc`.  Both `pᵢ ≤ 1` and `c ≥ -1` are needed,
and only to make every factor nonnegative. -/
theorem prod_one_add_mul_le {L : ℕ} (hL : 0 < L) (p : Fin L → ℝ) (c : ℝ)
    (hp0 : ∀ i, 0 ≤ p i) (hp1 : ∀ i, p i ≤ 1) (hc : -1 ≤ c) :
    ∏ i, (1 + p i * c) ≤ (1 + ((∑ i, p i) / (L : ℝ)) * c) ^ L := by
  have hLR : (0 : ℝ) < (L : ℝ) := by exact_mod_cast hL
  have hz : ∀ i : Fin L, 0 ≤ 1 + p i * c := by
    intro i
    rcases le_total 0 c with hcs | hcs
    · have : 0 ≤ p i * c := mul_nonneg (hp0 i) hcs
      linarith
    · have : c ≤ p i * c := by nlinarith [hp1 i, hp0 i]
      linarith
  refine (prod_le_pow_avg hL _ hz).trans_eq ?_
  congr 1
  rw [Finset.sum_add_distrib, Finset.sum_const, Finset.card_univ, Fintype.card_fin,
    nsmul_eq_mul, mul_one, ← Finset.sum_mul]
  field_simp

/-! ## The Poisson-binomial law

Indexed by the *set* of successes rather than by `Fin L → Bool`, so that a
single lemma — `Finset.prod_add` — supplies both totality and the generating
identity

    ∑_T P(T) z^|T|  =  ∏ (1 - pᵢ + pᵢ z).

The generating identity is the whole content of the Chernoff argument: the
left side has `z^k · P[S = k]` among its (nonnegative) terms, and the right
side is what `prod_one_add_mul_le` dominates. -/

/-- The law of `L` independent Bernoulli trials, indexed by the success set. -/
noncomputable def poissonBinom {L : ℕ} (p : Fin L → ℝ)
    (hp0 : ∀ i, 0 ≤ p i) (hp1 : ∀ i, p i ≤ 1) : FinPMF (Finset (Fin L)) where
  p := fun T => (∏ i ∈ T, p i) * ∏ i ∈ Finset.univ \ T, (1 - p i)
  nonneg := fun T => by
    have h1 : ∀ i ∈ Finset.univ \ T, (0 : ℝ) ≤ 1 - p i := by
      intro i _; linarith [hp1 i]
    exact mul_nonneg (Finset.prod_nonneg fun i _ => hp0 i) (Finset.prod_nonneg h1)
  total := by
    have h := Finset.prod_add (fun i => p i) (fun i => 1 - p i)
      (Finset.univ : Finset (Fin L))
    simp only [Finset.powerset_univ] at h
    calc ∑ T : Finset (Fin L), (∏ i ∈ T, p i) * ∏ i ∈ Finset.univ \ T, (1 - p i)
        = ∏ i : Fin L, (p i + (1 - p i)) := h.symm
      _ = 1 := by simp

@[simp] lemma poissonBinom_apply {L : ℕ} (p : Fin L → ℝ)
    (hp0 : ∀ i, 0 ≤ p i) (hp1 : ∀ i, p i ≤ 1) (T : Finset (Fin L)) :
    (poissonBinom p hp0 hp1).p T
      = (∏ i ∈ T, p i) * ∏ i ∈ Finset.univ \ T, (1 - p i) := rfl

/-- **The generating identity.**  `Finset.prod_add` again, now with `pᵢ z` in
place of `pᵢ`. -/
theorem poissonBinom_mgf {L : ℕ} (p : Fin L → ℝ)
    (hp0 : ∀ i, 0 ≤ p i) (hp1 : ∀ i, p i ≤ 1) (z : ℝ) :
    ∑ T : Finset (Fin L), (poissonBinom p hp0 hp1).p T * z ^ T.card
      = ∏ i, (1 - p i + p i * z) := by
  have h := Finset.prod_add (fun i => p i * z) (fun i => 1 - p i)
    (Finset.univ : Finset (Fin L))
  simp only [Finset.powerset_univ] at h
  calc ∑ T : Finset (Fin L), (poissonBinom p hp0 hp1).p T * z ^ T.card
      = ∑ T : Finset (Fin L),
          (∏ i ∈ T, p i * z) * ∏ i ∈ Finset.univ \ T, (1 - p i) := by
        refine Finset.sum_congr rfl fun T _ => ?_
        rw [poissonBinom_apply, Finset.prod_mul_distrib, Finset.prod_const]
        ring
    _ = ∏ i : Fin L, (p i * z + (1 - p i)) := h.symm
    _ = ∏ i, (1 - p i + p i * z) := Finset.prod_congr rfl fun i _ => by ring
/-! ## The Chernoff step

`z^k · P[S = k]` is a sub-sum of the generating function — the terms of the
full sum are nonnegative for `z ≥ 0` — so the generating identity and
`prod_one_add_mul_le` bound it in one line each. -/

theorem poissonBinom_chernoff {L : ℕ} (hL : 0 < L) (p : Fin L → ℝ)
    (hp0 : ∀ i, 0 ≤ p i) (hp1 : ∀ i, p i ≤ 1) (k : ℕ) {z : ℝ} (hz : 0 ≤ z) :
    z ^ k * (poissonBinom p hp0 hp1).prob (fun T => T.card = k)
      ≤ (1 + ((∑ i, p i) / (L : ℝ)) * (z - 1)) ^ L := by
  classical
  have hPnn := (poissonBinom p hp0 hp1).nonneg
  have step1 : z ^ k * (poissonBinom p hp0 hp1).prob (fun T => T.card = k)
      ≤ ∑ T : Finset (Fin L), (poissonBinom p hp0 hp1).p T * z ^ T.card := by
    have hcongr : ∑ T ∈ Finset.univ.filter (fun T : Finset (Fin L) => T.card = k),
            z ^ k * (poissonBinom p hp0 hp1).p T
        = ∑ T ∈ Finset.univ.filter (fun T : Finset (Fin L) => T.card = k),
            (poissonBinom p hp0 hp1).p T * z ^ T.card := by
      refine Finset.sum_congr rfl fun T hT => ?_
      rw [Finset.mem_filter] at hT
      rw [hT.2]; ring
    rw [FinPMF.prob, Finset.mul_sum, hcongr]
    exact Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _)
      (fun T _ _ => mul_nonneg (hPnn T) (pow_nonneg hz _))
  refine step1.trans ?_
  rw [poissonBinom_mgf]
  have heq : ∏ i, (1 - p i + p i * z) = ∏ i, (1 + p i * (z - 1)) :=
    Finset.prod_congr rfl fun i _ => by ring
  rw [heq]
  exact prod_one_add_mul_le hL p (z - 1) hp0 hp1 (by linarith)
/-! ## Optimising the free parameter

At `z = k(1-q) / (q(L-k))` the Chernoff bound collapses to

    q^k (1-q)^(L-k) L^L / (k^k (L-k)^(L-k)),

which is `exp(-L · D(k/L ‖ q))` written without an exponential.  Keeping it in
that form is what lets `types_lower_bound_at` cancel against it exactly. -/

private lemma opt_algebra {L k : ℕ} (hk0 : 0 < k) (hkL : k < L) {q z : ℝ}
    (hq0 : 0 < q) (hq1 : q < 1)
    (hzval : z = (k : ℝ) * (1 - q) / (q * ((L : ℝ) - k))) :
    (1 + q * (z - 1)) ^ L
      = q ^ k * (1 - q) ^ (L - k) * (L : ℝ) ^ L
          / ((k : ℝ) ^ k * ((L : ℝ) - k) ^ (L - k)) * z ^ k := by
  subst hzval
  have hkR : (0 : ℝ) < (k : ℝ) := by exact_mod_cast hk0
  have hkLR : (k : ℝ) < (L : ℝ) := by exact_mod_cast hkL
  have hLk : (0 : ℝ) < (L : ℝ) - k := by linarith
  have hqc : (0 : ℝ) < 1 - q := by linarith
  have hsum : (L - k) + k = L := by omega
  have hbase : 1 + q * ((k : ℝ) * (1 - q) / (q * ((L : ℝ) - k)) - 1)
      = (1 - q) * (L : ℝ) / ((L : ℝ) - k) := by
    field_simp
    ring
  rw [hbase]
  simp only [div_pow, mul_pow]
  rw [show ((1 : ℝ) - q) ^ L = (1 - q) ^ (L - k) * (1 - q) ^ k by
        rw [← pow_add, hsum],
      show ((L : ℝ) - k) ^ L = ((L : ℝ) - k) ^ (L - k) * ((L : ℝ) - k) ^ k by
        rw [← pow_add, hsum]]
  have h1 : ((1 : ℝ) - q) ^ k ≠ 0 := pow_ne_zero _ (ne_of_gt hqc)
  have h2 : (q : ℝ) ^ k ≠ 0 := pow_ne_zero _ (ne_of_gt hq0)
  have h3 : ((k : ℝ)) ^ k ≠ 0 := pow_ne_zero _ (ne_of_gt hkR)
  have h4 : ((L : ℝ) - k) ^ k ≠ 0 := pow_ne_zero _ (ne_of_gt hLk)
  have h5 : ((L : ℝ) - k) ^ (L - k) ≠ 0 := pow_ne_zero _ (ne_of_gt hLk)
  field_simp

/-- **The Chernoff bound for a Poisson binomial**, in the exponential-free form
`q^k (1-q)^(L-k) L^L / (k^k (L-k)^(L-k))`. -/
theorem poissonBinom_prob_le {L k : ℕ} (hL : 0 < L) (hk0 : 0 < k) (hkL : k < L)
    (p : Fin L → ℝ) (hp0 : ∀ i, 0 ≤ p i) (hp1 : ∀ i, p i ≤ 1)
    {q : ℝ} (hqval : (∑ i, p i) / (L : ℝ) = q) (hq0 : 0 < q) (hq1 : q < 1) :
    (poissonBinom p hp0 hp1).prob (fun T => T.card = k)
      ≤ q ^ k * (1 - q) ^ (L - k) * (L : ℝ) ^ L
          / ((k : ℝ) ^ k * ((L : ℝ) - k) ^ (L - k)) := by
  have hkR : (0 : ℝ) < (k : ℝ) := by exact_mod_cast hk0
  have hkLR : (k : ℝ) < (L : ℝ) := by exact_mod_cast hkL
  have hLk : (0 : ℝ) < (L : ℝ) - k := by linarith
  have hqc : (0 : ℝ) < 1 - q := by linarith
  have hz : (0 : ℝ) < (k : ℝ) * (1 - q) / (q * ((L : ℝ) - k)) :=
    div_pos (mul_pos hkR hqc) (mul_pos hq0 hLk)
  have hch := poissonBinom_chernoff hL p hp0 hp1 k (le_of_lt hz)
  rw [hqval, opt_algebra hk0 hkL hq0 hq1 rfl] at hch
  refine le_of_mul_le_mul_right ?_ (pow_pos hz k)
  linarith [hch]
/-! ## Cancelling the Chernoff bound against the method of types

`types_lower_bound_at` says `C(L,k) k^k (L-k)^(L-k) / L^L ≥ 1/(L+1)`, which is
exactly the reciprocal of the factor the Chernoff bound carries.  The two
therefore combine with no exponential and no entropy anywhere: the Poisson
binomial is pointwise at most `L+1` times the binomial of the same mean. -/

private lemma types_ratio {L k : ℕ} (hL : 0 < L) (hk0 : 0 < k) (hkL : k < L) :
    (L : ℝ) ^ L / ((k : ℝ) ^ k * ((L : ℝ) - k) ^ (L - k))
      ≤ ((L : ℝ) + 1) * (L.choose k : ℝ) := by
  have hkR : (0 : ℝ) < (k : ℝ) := by exact_mod_cast hk0
  have hkLR : (k : ℝ) < (L : ℝ) := by exact_mod_cast hkL
  have hLk : (0 : ℝ) < (L : ℝ) - k := by linarith
  have hLR : (0 : ℝ) < (L : ℝ) := by linarith
  have hL1 : (0 : ℝ) < (L : ℝ) + 1 := by linarith
  have hk : k ≤ L := le_of_lt hkL
  have h0 : (0 : ℝ) ≤ (k : ℝ) / (L : ℝ) := by positivity
  have h1 : (k : ℝ) / (L : ℝ) ≤ 1 := by rw [div_le_one hLR]; linarith
  have hD : (0 : ℝ) < (k : ℝ) ^ k * ((L : ℝ) - k) ^ (L - k) := by positivity
  have hLp : (0 : ℝ) < (L : ℝ) ^ L := by positivity
  have ht := types_lower_bound_at hL hk h0 h1
  have hcompl : 1 - (k : ℝ) / (L : ℝ) = ((L : ℝ) - k) / (L : ℝ) := by field_simp
  rw [hcompl, mul_assoc, pow_div_split hk] at ht
  rw [div_le_iff₀ hD]
  have h := mul_le_mul_of_nonneg_left ht (le_of_lt hL1)
  rw [mul_one_div, div_self (ne_of_gt hL1)] at h
  have h2 := mul_le_mul_of_nonneg_right h (le_of_lt hLp)
  rw [one_mul] at h2
  refine h2.trans_eq ?_
  field_simp

/-- **`lem:structured-route-domination`, pointwise form.**  A sum of `L`
independent Bernoulli variables puts at most `L+1` times as much mass on any
single value as the binomial with the same mean does. -/
theorem poissonBinom_le_binom {L k : ℕ} (hL : 0 < L) (hk0 : 0 < k) (hkL : k < L)
    (p : Fin L → ℝ) (hp0 : ∀ i, 0 ≤ p i) (hp1 : ∀ i, p i ≤ 1)
    {q : ℝ} (hqval : (∑ i, p i) / (L : ℝ) = q) (hq0 : 0 < q) (hq1 : q < 1) :
    (poissonBinom p hp0 hp1).prob (fun T => T.card = k)
      ≤ ((L : ℝ) + 1) * ((L.choose k : ℝ) * q ^ k * (1 - q) ^ (L - k)) := by
  have hb := poissonBinom_prob_le hL hk0 hkL p hp0 hp1 hqval hq0 hq1
  have hr := types_ratio hL hk0 hkL
  have hqk : (0 : ℝ) ≤ q ^ k * (1 - q) ^ (L - k) :=
    mul_nonneg (pow_nonneg (le_of_lt hq0) k) (pow_nonneg (by linarith) _)
  refine hb.trans ?_
  have hsplit : q ^ k * (1 - q) ^ (L - k) * (L : ℝ) ^ L
        / ((k : ℝ) ^ k * ((L : ℝ) - k) ^ (L - k))
      = (q ^ k * (1 - q) ^ (L - k))
        * ((L : ℝ) ^ L / ((k : ℝ) ^ k * ((L : ℝ) - k) ^ (L - k))) := by
    ring
  rw [hsplit]
  calc (q ^ k * (1 - q) ^ (L - k))
        * ((L : ℝ) ^ L / ((k : ℝ) ^ k * ((L : ℝ) - k) ^ (L - k)))
      ≤ (q ^ k * (1 - q) ^ (L - k)) * (((L : ℝ) + 1) * (L.choose k : ℝ)) :=
        mul_le_mul_of_nonneg_left hr hqk
    _ = ((L : ℝ) + 1) * ((L.choose k : ℝ) * q ^ k * (1 - q) ^ (L - k)) := by ring

end Spin
