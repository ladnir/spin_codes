import SpinCodes.Numeric.BAEvalDefs

namespace Spin.Numeric.Fix

/-- The complete residual at a fixed positive generating-function witness. -/
def fResidual (lg : LogFn) (u A B W slope intercept : Fix) : Option Fix :=
  (fgObjL lg u A).bind fun g =>
  (fpiBestClampL lg A B).bind fun p =>
  (fpiBestClampL lg B W).bind fun q =>
  some (sub (add (add g p) q) (add (mul slope W) intercept))

/-- Center the complete residual before multiplying by coordinate widths.
This retains cancellation between derivatives of adjacent accumulator terms. -/
def fResidualCentered (lg : LogFn) (u A B W slope intercept : Fix) : Option Fix :=
  if feasOK A B && feasOK B W then
    (fgObjL lg u (midPt A)).bind fun g =>
    (fpiEvalL lg (midPt A) (midPt B)).bind fun p =>
    (fpiEvalL lg (midPt B) (midPt W)).bind fun q =>
    (flogIWL lg u).bind fun lu =>
    (fDpaL lg A B).bind fun pa =>
    (fDpcL lg A B).bind fun pb =>
    (fDpaL lg B W).bind fun qb =>
    (fDpcL lg B W).bind fun qw =>
    let center := sub (add (add g p) q) (add (mul slope (midPt W)) intercept)
    some (add (add (add center
      (mul (sub pa lu) (sub A (midPt A))))
      (mul (add pb qb) (sub B (midPt B))))
      (mul (sub qw slope) (sub W (midPt W))))
  else none

/-- The witness also chooses the bound: nonnegative selects the centered
bound, negative selects the component bound. Its absolute value encodes u. -/
def majorantLeaf (lg : LogFn) (slope intercept : Fix)
    (witness : Int) (A B W : Fix) : Bool :=
  infeasible A B W ||
    (decide (witness ≠ 0) &&
      match (if 0 < witness then fResidualCentered else fResidual)
          lg (sc (Int.ofNat witness.natAbs)) A B W slope intercept with
      | some v => decide (v.hi < 0)
      | none => false)

end Spin.Numeric.Fix

