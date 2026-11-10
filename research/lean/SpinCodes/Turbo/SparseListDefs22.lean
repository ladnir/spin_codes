import SpinCodes.Turbo.SparseRowBaseDefs22

/-!
# Core-only list-structured integer checks

Indexed access into large literals is not constant time in the kernel; the
flat-array checks exhaust kernel memory at endpoint size. These checks
traverse lists in lockstep instead. Matrices are lists of 11 rows of 11
entries; degree-≤2 polynomial matrices are lists of three matrices; q1 rows
are lists of 65 row vectors.
-/
namespace Spin.Turbo

def sparseZero11L22 : List Int := [0,0,0,0,0,0,0,0,0,0,0]

/-- Row vector times matrix: Σ_i x_i · (row i of M). -/
def sparseVecMatL22 (x : List Int) (M : List (List Int)) : List Int :=
  (x.zip M).foldl (fun acc p => List.zipWith (· + ·) acc (p.2.map (p.1 * ·))) sparseZero11L22

def sparseVecAddL22 (a b : List Int) : List Int := List.zipWith (· + ·) a b

def sparseLeScaledL22 (s : Int) (lhs rhs : List Int) : Bool :=
  (List.zipWith (fun a b => decide (a ≤ s * b)) lhs rhs).all id

def sparseLen11L22 (x : List Int) : Bool := decide (x.length = 11)

def sparseMatShapeL22 (M : List (List Int)) : Bool :=
  decide (M.length = 11) && M.all sparseLen11L22

def sparsePolyShapeL22 (P : List (List (List Int))) : Bool :=
  decide (P.length = 3) && P.all sparseMatShapeL22

def sparseMatMulL22 (A B : List (List Int)) : List (List Int) :=
  A.map fun row => sparseVecMatL22 row B

def sparseMatAddL22 (A B : List (List Int)) : List (List Int) :=
  List.zipWith sparseVecAddL22 A B

def sparseMatLeScaledL22 (s : Int) (L R : List (List Int)) : Bool :=
  (List.zipWith (sparseLeScaledL22 s) L R).all id

def sparsePolyProdL22 (a b : List (List (List Int))) (d : Nat) : List (List Int) :=
  match d with
  | 0 => sparseMatMulL22 (a.getD 0 []) (b.getD 0 [])
  | 1 => sparseMatAddL22 (sparseMatMulL22 (a.getD 0 []) (b.getD 1 []))
      (sparseMatMulL22 (a.getD 1 []) (b.getD 0 []))
  | 2 => sparseMatAddL22 (sparseMatAddL22 (sparseMatMulL22 (a.getD 0 []) (b.getD 2 []))
      (sparseMatMulL22 (a.getD 1 []) (b.getD 1 [])))
      (sparseMatMulL22 (a.getD 2 []) (b.getD 0 []))
  | _ => []

/-- `a`, `b` are assumed shaped (checked where they are produced); `u` is
shape-checked here. -/
def sparsePolyMulCheckL22 (s : Int) (a b u : List (List (List Int))) : Bool :=
  sparsePolyShapeL22 u &&
    sparseMatLeScaledL22 s (sparsePolyProdL22 a b 0) (u.getD 0 []) &&
    sparseMatLeScaledL22 s (sparsePolyProdL22 a b 1) (u.getD 1 []) &&
    sparseMatLeScaledL22 s (sparsePolyProdL22 a b 2) (u.getD 2 [])

/-- `u_d ≤ choose(16R,d) · result_d`, entrywise. -/
def sparseNormalizeCheckL22 (R : Nat) (u result : List (List (List Int))) : Bool :=
  sparsePolyShapeL22 result &&
    sparseMatLeScaledL22 (sparseRegionDenominator22' R 0) (u.getD 0 []) (result.getD 0 []) &&
    sparseMatLeScaledL22 (sparseRegionDenominator22' R 1) (u.getD 1 []) (result.getD 1 []) &&
    sparseMatLeScaledL22 (sparseRegionDenominator22' R 2) (u.getD 2 []) (result.getD 2 [])
where
  sparseRegionDenominator22' (R d : Nat) : Int :=
    match d with
    | 0 => 1
    | 1 => (16 * R : Nat)
    | 2 => ((16 * R) * (16 * R - 1) / 2 : Nat)
    | _ => 0

def sparseRegionalLastL22 (prev : List (List (List Int))) :
    List (Bool × List (List (List Int)) × List (List (List Int))) → List (List (List Int))
  | [] => prev
  | (bit, sq, next) :: tail => sparseRegionalLastL22 (if bit then next else sq) tail

def sparseRegionalStepL22 (e : Nat)
    (step : Bool × List (List (List Int)) × List (List (List Int))) : Nat :=
  2 * e + if step.1 then 1 else 0

def sparseRegionalPowerCheckL22 (s : Int) (base prev : List (List (List Int))) :
    List (Bool × List (List (List Int)) × List (List (List Int))) → Bool
  | [] => true
  | (bit, sq, next) :: tail =>
    (sparsePolyMulCheckL22 s prev prev sq &&
      (if bit then sparsePolyMulCheckL22 s sq base next else true)) &&
      sparseRegionalPowerCheckL22 s base (if bit then next else sq) tail

/-- Base stage dominates choose(16,d)·row_d·mult, with an exact shape. -/
def sparseRegionalBaseCheckL22 (row : Nat → Nat → Nat → Int) (mult : Int)
    (base : List (List (List Int))) : Bool :=
  sparsePolyShapeL22 base &&
    (List.range 3).all fun d => (List.range 11).all fun i => (List.range 11).all fun k =>
      decide (sparseBinom16Small22 d * row d i k * mult ≤
        ((base.getD d []).getD i []).getD k 0)

def sparseRegionalCertificateCheckL22 (R : Nat) (s : Int) (row : Nat → Nat → Nat → Int)
    (mult : Int) (base : List (List (List Int)))
    (trace : List (Bool × List (List (List Int)) × List (List (List Int))))
    (result : List (List (List Int))) : Bool :=
  decide (0 < s ∧ 0 < mult ∧ 0 < R ∧ trace.foldl sparseRegionalStepL22 1 = R) &&
    sparseRegionalBaseCheckL22 row mult base &&
    sparseRegionalPowerCheckL22 s base base trace &&
    sparseNormalizeCheckL22 R (sparseRegionalLastL22 base trace) result

/-! q1 rows: 65 shells of 11 coordinates. -/

def sparseQ1InitialRowsL22 (s : Int) : List (List Int) :=
  (s :: List.replicate 10 0) :: List.replicate 64 sparseZero11L22

def sparseQ1StepCheckL22 (s : Int) (r0 r1 prev next : List (List Int)) : Bool :=
  decide (next.length = 65) && next.all sparseLen11L22 &&
    (next.zip (prev.zip (sparseZero11L22 :: prev))).all fun t =>
      sparseLeScaledL22 s (sparseVecAddL22 (sparseVecMatL22 t.2.1 r0)
        (sparseVecMatL22 t.2.2 r1)) t.1

/-! Scaled matrix squaring with explicit binary exponents (tail powers). -/

def sparseMatScaleL22 (c : Int) (A : List (List Int)) : List (List Int) :=
  A.map fun row => row.map (c * ·)

/-- One stage `(down, up, U')`: 2^down · (U·U) ≤ s · 2^up · U' entrywise. -/
def sparseSquareStageCheckL22 (s : Int) (U : List (List Int))
    (stage : Nat × Nat × List (List Int)) : Bool :=
  sparseMatShapeL22 stage.2.2 &&
    sparseMatLeScaledL22 (s * 2 ^ stage.2.1)
      (sparseMatScaleL22 (2 ^ stage.1) (sparseMatMulL22 U U)) stage.2.2

def sparseSquareChainCheckL22 (s : Int) (U : List (List Int)) :
    List (Nat × Nat × List (List Int)) → Bool
  | [] => true
  | st :: rest => sparseSquareStageCheckL22 s U st && sparseSquareChainCheckL22 s st.2.2 rest

def sparseSquareChainLastL22 (U : List (List Int)) :
    List (Nat × Nat × List (List Int)) → List (List Int)
  | [] => U
  | st :: rest => sparseSquareChainLastL22 st.2.2 rest

end Spin.Turbo
