import SpinCodes.UniformMixer.RegionalContinuumError

/-!
# Rounding routed bit sites to 64-bit encoder steps

These arithmetic lemmas quantify the displacement introduced when a bit
site is assigned to its enclosing 64-bit step. They will be summed over
the ordered support to compare discrete and continuum gap vectors.
-/

noncomputable section

namespace Spin.UniformMixer

theorem firstGap_rounding (R : ℕ) (hR : 0 < R) (x : ℕ) :
    |((x / 64 : ℕ) : ℝ) / R - (x : ℝ) / (64 * R)| ≤
      1 / (R : ℝ) := by
  have hRp : (0 : ℝ) < R := by exact_mod_cast hR
  have hRn : (R : ℝ) ≠ 0 := ne_of_gt hRp
  have hmodNat := Nat.mod_lt x (by omega : 0 < 64)
  have hmod : ((x % 64 : ℕ) : ℝ) ≤ 64 := by
    exact_mod_cast (Nat.le_of_lt hmodNat)
  have hdecomp : (x : ℝ) = 64 * ((x / 64 : ℕ) : ℝ) +
      ((x % 64 : ℕ) : ℝ) := by
    exact_mod_cast (Nat.div_add_mod x 64).symm
  have heq : ((x / 64 : ℕ) : ℝ) / R - (x : ℝ) / (64 * R) =
      -(((x % 64 : ℕ) : ℝ) / (64 * R)) := by
    rw [hdecomp]
    field_simp
    ring
  rw [heq, abs_neg, abs_of_nonneg (by positivity)]
  apply (div_le_iff₀ (by positivity : (0 : ℝ) < 64 * R)).mpr
  have hmul : (1 : ℝ) / R * (64 * R) = 64 := by
    field_simp
  rw [hmul]
  exact hmod

theorem middleGap_rounding (R : ℕ) (hR : 0 < R)
    (x y : ℕ) (hstep : x / 64 + 1 ≤ y / 64) :
    |((y / 64 - (x / 64 + 1) : ℕ) : ℝ) / R -
      ((y : ℝ) - (x : ℝ)) / (64 * R)| ≤
      3 / (R : ℝ) := by
  have hRp : (0 : ℝ) < R := by exact_mod_cast hR
  have hRn : (R : ℝ) ≠ 0 := ne_of_gt hRp
  have hcast : ((y / 64 - (x / 64 + 1) : ℕ) : ℝ) =
      ((y / 64 : ℕ) : ℝ) - ((x / 64 : ℕ) : ℝ) - 1 := by
    rw [Nat.cast_sub hstep]
    push_cast
    ring
  have heq : ((y / 64 - (x / 64 + 1) : ℕ) : ℝ) / R -
      ((y : ℝ) - (x : ℝ)) / (64 * R) =
      (((y / 64 : ℕ) : ℝ) / R - (y : ℝ) / (64 * R)) -
      (((x / 64 : ℕ) : ℝ) / R - (x : ℝ) / (64 * R)) -
      1 / R := by
    rw [hcast]
    ring
  rw [heq]
  have ht := firstGap_rounding R hR y
  have hx := firstGap_rounding R hR x
  have htri (a b c : ℝ) : |a - b - c| ≤ |a| + |b| + |c| := by
    calc
      _ ≤ |a - b| + |c| := by
        simpa only [sub_zero, zero_sub, abs_neg] using
          (abs_sub_le (a - b) 0 c)
      _ ≤ |a| + |b| + |c| := by
        have hab := abs_sub_le a 0 b
        simp only [sub_zero, zero_sub, abs_neg] at hab
        linarith
  have hh := htri
    (((y / 64 : ℕ) : ℝ) / R - (y : ℝ) / (64 * R))
    (((x / 64 : ℕ) : ℝ) / R - (x : ℝ) / (64 * R))
    (1 / R)
  rw [abs_of_nonneg (by positivity : (0 : ℝ) ≤ 1 / R)] at hh
  have hthree : 3 / (R : ℝ) = 1 / R + 1 / R + 1 / R := by ring
  rw [hthree]
  linarith

theorem lastGap_rounding (R : ℕ) (hR : 0 < R)
    (x : ℕ) (hstep : x / 64 + 1 ≤ R) :
    |((R - (x / 64 + 1) : ℕ) : ℝ) / R -
      (1 - (x : ℝ) / (64 * R))| ≤
      2 / (R : ℝ) := by
  have hRp : (0 : ℝ) < R := by exact_mod_cast hR
  have hRn : (R : ℝ) ≠ 0 := ne_of_gt hRp
  have hcast : ((R - (x / 64 + 1) : ℕ) : ℝ) =
      (R : ℝ) - ((x / 64 : ℕ) : ℝ) - 1 := by
    rw [Nat.cast_sub hstep]
    push_cast
    ring
  have heq : ((R - (x / 64 + 1) : ℕ) : ℝ) / R -
      (1 - (x : ℝ) / (64 * R)) =
      -(((x / 64 : ℕ) : ℝ) / R - (x : ℝ) / (64 * R)) -
      1 / R := by
    rw [hcast]
    field_simp
    ring
  rw [heq]
  have hx := firstGap_rounding R hR x
  have htri := abs_sub_le
    (-(((x / 64 : ℕ) : ℝ) / R - (x : ℝ) / (64 * R)))
    0 (1 / (R : ℝ))
  simp only [sub_zero, zero_sub, abs_neg] at htri
  rw [abs_of_nonneg (by positivity : (0 : ℝ) ≤ 1 / R)] at htri
  have htwo : 2 / (R : ℝ) = 1 / R + 1 / R := by ring
  rw [htwo]
  exact htri.trans (add_le_add_left hx _)

/-- Sum of absolute gap errors in an ordered support. The second argument
is the next encoder step; the third is the preceding bit position. -/
def supportGapError (R : ℕ) :
    List (Fin (64 * R)) → ℕ → ℕ → ℝ
  | [], prev, prevBit =>
      |(((R - prev : ℕ) : ℝ) / R) -
        (1 - (prevBit : ℝ) / (64 * R))|
  | x :: xs, prev, prevBit =>
      |(((x.val / 64 - prev : ℕ) : ℝ) / R) -
        ((x.val : ℝ) - (prevBit : ℝ)) / (64 * R)| +
        supportGapError R xs (x.val / 64 + 1) x.val

theorem supportGapError_bound {R prev : ℕ}
    {xs : List (Fin (64 * R))} (h : ValidSupport R prev xs)
    (hR : 0 < R) :
    ∀ prevBit : ℕ,
      ((prev = 0 ∧ prevBit = 0) ∨ prev = prevBit / 64 + 1) →
      supportGapError R xs prev prevBit ≤
        3 * (xs.length + 1 : ℕ) / (R : ℝ) := by
  induction h with
  | nil hprev =>
      intro prevBit hrelation
      rcases hrelation with ⟨hp, hb⟩ | hp
      · subst prevBit
        rw [hp]
        simp [supportGapError, div_self (by exact_mod_cast
          (Nat.ne_of_gt hR) : (R : ℝ) ≠ 0)]
        positivity
      · rw [hp]
        have hh := lastGap_rounding R hR prevBit (by omega)
        have hRp : (0 : ℝ) < R := by exact_mod_cast hR
        have h23 : 2 / (R : ℝ) ≤ 3 / R := by
          have hbase : 0 ≤ 1 / (R : ℝ) := by positivity
          have heq : 3 / (R : ℝ) = 2 / R + 1 / R := by ring
          rw [heq]
          linarith
        simpa [supportGapError] using hh.trans h23
  | @cons prev x xs hgap hinterior htail ih =>
      intro prevBit hrelation
      have hhead :
          |(((x.val / 64 - prev : ℕ) : ℝ) / R) -
            ((x.val : ℝ) - (prevBit : ℝ)) / (64 * R)| ≤
          3 / (R : ℝ) := by
        rcases hrelation with ⟨hp, hb⟩ | hp
        · subst prev
          subst prevBit
          have hh := firstGap_rounding R hR x.val
          have h13 : 1 / (R : ℝ) ≤ 3 / R := by
            have hbase : 0 ≤ 1 / (R : ℝ) := by positivity
            have heq : 3 / (R : ℝ) = 1 / R + 2 * (1 / R) := by ring
            rw [heq]
            linarith
          simpa using hh.trans h13
        · rw [hp]
          exact middleGap_rounding R hR prevBit x.val (by omega)
      have htailBound := ih x.val (Or.inr rfl)
      change
        |(((x.val / 64 - prev : ℕ) : ℝ) / R) -
          ((x.val : ℝ) - (prevBit : ℝ)) / (64 * R)| +
          supportGapError R xs (x.val / 64 + 1) x.val ≤
          3 * (xs.length + 1 + 1 : ℕ) / (R : ℝ)
      calc
        _ ≤ 3 / (R : ℝ) +
            3 * (xs.length + 1 : ℕ) / (R : ℝ) :=
          add_le_add hhead htailBound
        _ = 3 * (xs.length + 1 + 1 : ℕ) / (R : ℝ) := by
          push_cast
          ring

theorem good_support_gap_error {R : ℕ} (hR : 0 < R)
    (T : Finset (Fin (64 * R))) (hgood : ¬ badRegionalSites T) :
    supportGapError R (T.sort (· ≤ ·)) 0 0 ≤
      3 * (T.card + 1 : ℕ) / (R : ℝ) := by
  have h := supportGapError_bound
    (goodRegionalSites_validSupport hR T hgood) hR 0
    (Or.inl ⟨rfl, rfl⟩)
  simpa only [Finset.length_sort] using h

end Spin.UniformMixer
