import SpinCodes.Turbo.SparseFixedDefs22

/-!
# Core-only regional binary-powering certificate

Each trace step stores its squared stage and, when the bit is set, the stage
multiplied by the base. Every stage is a materialized array, so the kernel
checks one product per stored stage and never recomputes earlier stages.
The exponent is read MSB-first with an implicit leading one (the base).
-/
namespace Spin.Turbo

@[inline] def sparseRegionalStep22 (e : Nat) (step : Bool × Array Int × Array Int) : Nat :=
  2 * e + if step.1 then 1 else 0

def sparseRegionalLast22 (prev : Array Int) : List (Bool × Array Int × Array Int) → Array Int
  | [] => prev
  | (bit, sq, next) :: tail => sparseRegionalLast22 (if bit then next else sq) tail

def sparseRegionalPowerCheck22 (scale : Int) (base prev : Array Int) :
    List (Bool × Array Int × Array Int) → Bool
  | [] => true
  | (bit, sq, next) :: tail =>
    sparseIntPolyMulUpperCheck22 scale prev prev sq &&
      (if bit then sparseIntPolyMulUpperCheck22 scale sq base next else true) &&
      sparseRegionalPowerCheck22 scale base (if bit then next else sq) tail

/-- Full regional certificate: positive scale, exact exponent R, the chain,
and normalization of the final stage by choose(16R,d). -/
def sparseRegionalCertificateCheck22 (R : Nat) (scale : Int) (base : Array Int)
    (trace : List (Bool × Array Int × Array Int)) (result : Array Int) : Bool :=
  decide (0 < scale ∧ 0 < R ∧ trace.foldl sparseRegionalStep22 1 = R) &&
    sparseRegionalPowerCheck22 scale base base trace &&
    sparseIntRegionNormalizeCheck22 R (sparseRegionalLast22 base trace) result

end Spin.Turbo
