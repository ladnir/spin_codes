/-!
# Core-only integer sparse arithmetic checks

Audit: Numeric.FixedDefs and the old split certificate blocks keep kernel
checks free of Mathlib; soundness belongs in the companion module. The
fixed11-coordinate inner sum and row tests are explicitly unrolled. Flat
arrays give direct access without nested function-literal lookup. The
endpoint scale is explicit because Chernoff/outer-count amplification
requires finer absolute precision than the logarithm grid.
-/
namespace Spin.Turbo

@[inline] def sparseAll11Bool22 (f : Nat → Bool) : Bool :=
  f 0 && f 1 && f 2 && f 3 && f 4 && f 5 && f 6 && f 7 && f 8 && f 9 && f 10

@[inline] def sparseAll3Bool22 (f : Nat → Bool) : Bool := f 0 && f 1 && f 2

@[inline] def sparseSum11Int22 (f : Nat → Int) : Int :=
  f 0 + f 1 + f 2 + f 3 + f 4 + f 5 + f 6 + f 7 + f 8 + f 9 + f 10

@[inline] def sparseIntMatrixGet22 (a : Array Int) (i j : Nat) : Int :=
  a.getD (11 * i + j) 0

@[inline] def sparseIntPolyMatrixGet22 (a : Array Int) (d i j : Nat) : Int :=
  a.getD (121 * d + 11 * i + j) 0

@[inline] def sparseIntRowsGet22 (a : Array Int) (v j : Nat) : Int :=
  a.getD (11 * v + j) 0

def sparseIntMatrix22 (a : Array Int) (i j : Fin 11) : Int :=
  sparseIntMatrixGet22 a i.val j.val

def sparseIntPolyMatrix22 (a : Array Int) (d : Nat) (i j : Fin 11) : Int :=
  sparseIntPolyMatrixGet22 a d i.val j.val

def sparseIntRows22 (a : Array Int) (v : Nat) (j : Fin 11) : Int :=
  sparseIntRowsGet22 a v j.val

def sparseIntMatrixNonnegCheck22 (a : Array Int) : Bool :=
  sparseAll11Bool22 fun i => sparseAll11Bool22 fun j => decide (0 ≤ sparseIntMatrixGet22 a i j)

def sparseIntMatrixMulUpperCheck22 (scale : Int) (a b u : Array Int) : Bool :=
  sparseAll11Bool22 fun i => sparseAll11Bool22 fun j => decide
    (sparseSum11Int22 (fun k => sparseIntMatrixGet22 a i k * sparseIntMatrixGet22 b k j) ≤
      scale * sparseIntMatrixGet22 u i j)

@[inline] def sparseIntPolyDot22 (a b : Array Int) (r s i j : Nat) : Int :=
  sparseSum11Int22 fun k => sparseIntPolyMatrixGet22 a r i k * sparseIntPolyMatrixGet22 b s k j

def sparseIntPolyProduct22 (a b : Array Int) (d i j : Nat) : Int :=
  match d with
  | 0 => sparseIntPolyDot22 a b 0 0 i j
  | 1 => sparseIntPolyDot22 a b 0 1 i j + sparseIntPolyDot22 a b 1 0 i j
  | 2 => sparseIntPolyDot22 a b 0 2 i j + sparseIntPolyDot22 a b 1 1 i j + sparseIntPolyDot22 a b 2 0 i j
  | _ => 0

/-- All degree0..2 coefficients of an ordered polynomial product. -/
def sparseIntPolyMulUpperCheck22 (scale : Int) (a b u : Array Int) : Bool :=
  sparseAll3Bool22 fun d => sparseAll11Bool22 fun i => sparseAll11Bool22 fun j =>
    decide (sparseIntPolyProduct22 a b d i j ≤ scale * sparseIntPolyMatrixGet22 u d i j)

/-- choose(16R,d) for only the three sparse occupancies. -/
def sparseRegionDenominator22 (R d : Nat) : Nat :=
  match d with
  | 0 => 1
  | 1 => 16 * R
  | 2 => (16 * R) * (16 * R - 1) / 2
  | _ => 0

def sparseIntRegionNormalizeCheck22 (R : Nat) (u result : Array Int) : Bool :=
  sparseAll3Bool22 fun d => sparseAll11Bool22 fun i => sparseAll11Bool22 fun j =>
    decide (sparseIntPolyMatrixGet22 u d i j ≤
      (sparseRegionDenominator22 R d : Int) * sparseIntPolyMatrixGet22 result d i j)

/-- Every one of the65 shell rows is checked, including zero/default rows. -/
def sparseIntQ1StepUpperCheck22 (scale : Int) (r0 r1 prev next : Array Int) : Bool :=
  (List.range 65).all fun v => sparseAll11Bool22 fun j => decide
    (sparseSum11Int22 (fun i => sparseIntRowsGet22 prev v i * sparseIntMatrixGet22 r0 i j) +
      (if v = 0 then 0 else
        sparseSum11Int22 (fun i => sparseIntRowsGet22 prev (v - 1) i * sparseIntMatrixGet22 r1 i j)) ≤
      scale * sparseIntRowsGet22 next v j)

def sparseIntQ1InitialCheck22 (scale : Int) (v : Array Int) : Bool :=
  (List.range 65).all fun degree => sparseAll11Bool22 fun j => decide
    ((if degree = 0 ∧ j = 0 then scale else 0) ≤ sparseIntRowsGet22 v degree j)

/-- Positive Collatz vector for one whole joint regional marker matrix. -/
def sparseIntCollatzCheck22 (scale rhoNum rhoDen : Int) (m vector : Array Int) : Bool :=
  decide (0 < scale ∧ 0 < rhoDen ∧ 0 ≤ rhoNum) &&
    ((sparseAll11Bool22 fun i => decide (0 < vector.getD i 0)) &&
      sparseAll11Bool22 fun i => decide
        (rhoDen * sparseSum11Int22 (fun j => sparseIntMatrixGet22 m i j * vector.getD j 0) ≤
          scale * rhoNum * vector.getD i 0))

/-- MSB-first exponent encoded by a short square/multiply trace. -/
def sparsePowerTraceExponent22 (trace : List (Bool × Int)) : Nat :=
  trace.foldl (fun e step => 2 * e + if step.1 then 1 else 0) 0

def sparsePowerTraceLast22 (initial : Int) : List (Bool × Int) → Int
  | [] => initial
  | (_, next) :: tail => sparsePowerTraceLast22 next tail

def sparseIntPowerTraceCheck22 (scale base prev : Int) : List (Bool × Int) → Bool
  | [] => true
  | (bit, next) :: tail =>
    decide (0 ≤ prev ∧ 0 ≤ next ∧
      (if bit then scale * scale * next ≤ prev * prev * base
        else scale * next ≤ prev * prev)) && sparseIntPowerTraceCheck22 scale base next tail

/-- Endpoint use is restricted to the low-tilt q1/q2 family. The audited
retained tilts have lambda*d < 168 (less than 243 binary exponent bits);
the sufficient first 3/first 2 families need less than 122 bits. A 512-bit
absolute trace must finish strictly positive, as this check requires.
High-tilt tail witnesses use a separate explicit exponent/rescaling path. -/
def sparseIntPowerCertificateCheck22 (scale base : Int) (exponent : Nat)
    (trace : List (Bool × Int)) : Bool :=
  decide (0 < scale ∧ 0 ≤ base ∧ sparsePowerTraceExponent22 trace = exponent ∧
    0 < sparsePowerTraceLast22 scale trace) && sparseIntPowerTraceCheck22 scale base scale trace

end Spin.Turbo


