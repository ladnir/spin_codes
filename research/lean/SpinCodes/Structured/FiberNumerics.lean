import SpinCodes.Structured.FiberNumericsDefs
import SpinCodes.Structured.PolyIdentity
import SpinCodes.Structured.ConcreteFiberBounds

/-! The executable Fourier sums compute the mathematical Krawtchouk sums.
This file does not assume or certify a numerical spectrum. -/

namespace Spin.Structured.FiberNumerics
open Finset

theorem list_sum_range (f : ℕ → ℤ) (n : ℕ) :
    ((List.range n).map f).sum = ∑ i ∈ Finset.range n, f i := by
  induction n with
  | zero => simp
  | succ n ih => simp [List.range_succ, ih, Finset.sum_range_succ]

theorem kraw_eq (j w : ℕ) : kraw j w = Spin.krawtchouk 128 j w := by
  simp only [kraw, Spin.krawtchouk, polyChoose_eq, list_sum_range]

theorem transform_eq (spectrum : List ℕ) {j : ℕ} {values : List ℤ}
    (hv : (List.range 129).map (kraw j) = values) (f : ℤ → ℤ) :
    transform spectrum values f =
      ∑ w : Fin 129, (spectrum.getD w 0 : ℤ) * f (Spin.krawtchouk 128 j w) := by
  rw [transform, list_sum_range]
  rw [Fin.sum_univ_eq_sum_range (fun w : ℕ =>
    (spectrum.getD w 0 : ℤ) * f (Spin.krawtchouk 128 j w)) 129]
  apply Finset.sum_congr rfl
  intro w hw
  have hw' : w < 129 := Finset.mem_range.mp hw
  rw [← hv, List.getD_eq_getElem ((List.range 129).map (kraw j)) _ (by simpa using hw')]
  simp only [List.getElem_map, List.getElem_range, kraw_eq]

theorem signedSum_eq (spectrum : List ℕ) {j : ℕ} {values : List ℤ}
    (hv : (List.range 129).map (kraw j) = values) :
    signedSum spectrum values =
      ∑ w : Fin 129, (spectrum.getD w 0 : ℤ) * Spin.krawtchouk 128 j w :=
  transform_eq spectrum hv id

theorem squareSum_eq (spectrum : List ℕ) {j : ℕ} {values : List ℤ}
    (hv : (List.range 129).map (kraw j) = values) :
    squareSum spectrum values =
      ∑ w : Fin 129, (spectrum.getD w 0 : ℤ) * (Spin.krawtchouk 128 j w) ^ 2 := by
  simpa only [squareSum, pow_two] using transform_eq spectrum hv (fun k => k * k)

theorem absSum_eq (spectrum : List ℕ) {j : ℕ} {values : List ℤ}
    (hv : (List.range 129).map (kraw j) = values) :
    absSum spectrum values =
      ∑ w : Fin 129, (spectrum.getD w 0 : ℤ) * |Spin.krawtchouk 128 j w| := by
  simpa only [absSum, Int.natCast_natAbs] using transform_eq spectrum hv (fun k => (k.natAbs : ℤ))

end Spin.Structured.FiberNumerics
