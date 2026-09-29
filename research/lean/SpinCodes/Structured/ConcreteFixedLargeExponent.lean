import SpinCodes.Structured.ConcreteFixedLargeLaplace

/-! Closed form of the concave tilted exponent maximum. -/
noncomputable section
namespace Spin.Structured.Placement
open Set

def largeExponent (sig tau k : ℝ) : ℝ :=
  -k + tau*sig/(1-sig) + k*Real.log (k*(1-sig)/(tau*sig))

theorem log_affine_envelope {c tau k y : ℝ} (hc : 0 < c) (htau : 0 < tau)
    (hk : 0 < k) (hy : 0 ≤ y) :
    -tau*y + k*Real.log (1+c*y) ≤ -k+tau/c+k*Real.log (k*c/tau) := by
  have ht : 0 < 1+c*y := by positivity
  have ht0 : 0 < k*c/tau := by positivity
  have h := Real.log_le_sub_one_of_pos (div_pos ht ht0)
  rw [Real.log_div (ne_of_gt ht) (ne_of_gt ht0)] at h
  have hh := mul_le_mul_of_nonneg_left h hk.le
  have he : k*((1+c*y)/(k*c/tau)-1) = tau/c+tau*y-k := by
    field_simp
  rw [he] at hh
  linarith

/-- The unrestricted tangent maximum also bounds the unit interval. -/
theorem largeExponent_bounds {sig tau k y : ℝ} (hsig0 : 0 < sig) (hsig1 : sig < 1)
    (htau : 0 < tau) (hk : 0 < k) (hy : 0 ≤ y) :
    -tau*y + k*Real.log (1-y+y/sig) ≤ largeExponent sig tau k := by
  have hc : 0 < (1-sig)/sig := div_pos (sub_pos.mpr hsig1) hsig0
  have h := log_affine_envelope hc htau hk hy
  have ht : 1+(1-sig)/sig*y = 1-y+y/sig := by field_simp; ring
  have h1 : tau/((1-sig)/sig) = tau*sig/(1-sig) := by field_simp
  have h2 : k*((1-sig)/sig)/tau = k*(1-sig)/(tau*sig) := by ring
  simpa only [ht,h1,h2,largeExponent] using h

/-- Paper parameters, retaining the exact logarithm for interval certification. -/
def paperLargeExponent : ℝ := largeExponent (127/250) (133/125) (4/3)

theorem paperLargeExponent_bound {Q : ℕ} (hQ : 3 ≤ Q) (y : ℝ) (hy : y ∈ Icc (0:ℝ) 1) :
    -(133/125:ℝ)*y + (1+1/(Q:ℝ))*Real.log (1-y+y/(127/250)) ≤ paperLargeExponent := by
  apply exponent_le_three hQ (by norm_num) (by norm_num) ?_ y hy
  intro z hz
  have h := largeExponent_bounds (sig := (127/250:ℝ)) (tau := (133/125:ℝ))
    (k := (4/3:ℝ)) (by norm_num) (by norm_num) (by norm_num) (by norm_num) hz.1
  convert h using 1 <;> norm_num [paperLargeExponent]

#print axioms log_affine_envelope
#print axioms largeExponent_bounds
#print axioms paperLargeExponent_bound
end Spin.Structured.Placement
