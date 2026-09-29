/-
Fixed-point arithmetic: the *definitions*.

`DECISIONS.md` D2 records the measurement that forced this file: the Lean
kernel cannot reduce `ℚ` arithmetic through `decide`, but it reduces `Int`
arithmetic without trouble.  The box covers of `eq:structured-ba-dense-tail`
and `eq:ba-spectrum-majorant` are thousands of checks each, so they have to be
discharged by kernel reduction, which means integers.

A value is represented by a pair of integer numerators over the fixed implicit
denominator `scale`: `Fix.mk lo hi` stands for the set of reals `x` with
`lo/scale ≤ x ≤ hi/scale`.  Every operation rounds *outward*, so the enclosure
is always conservative; `Fixed.lean` proves that against `ℝ`.

This file imports nothing (beyond the automatic `Init`), for the same reason
`Structured/PolyCertDefs.lean` does: a data module that pulls in Mathlib makes
the certificate checks an order of magnitude slower and, at the sizes involved
here, exhausts memory.  Keep it that way — all reasoning lives next door.
-/

namespace Spin.Numeric

/-- Decimal digits of fixed-point precision. -/
def prec : Nat := 30

/-- The fixed-point denominator.  Written as a literal rather than `10 ^ prec`
so that kernel reduction never has to run `Nat.pow`. -/
def scale : Int := 1000000000000000000000000000000

/-! ## Directed integer division

`Int`'s `/` is Euclidean, which for a positive divisor is floor division.  That
gives rounding down directly, and rounding up by negating twice. -/

/-- Round `n / d` down (`d > 0`). -/
def fdiv (n d : Int) : Int := n / d

/-- Round `n / d` up (`d > 0`). -/
def cdiv (n d : Int) : Int := -((-n) / d)

/-! ## Order helpers

Defined here rather than taken from an instance, so that the kernel reduces
them by `Int.decLe` alone and no `LinearOrder` structure is unfolded. -/

def imin (a b : Int) : Int := if a ≤ b then a else b
def imax (a b : Int) : Int := if a ≤ b then b else a

/-! ## The enclosure type -/

/-- An enclosure of a real number: `lo/scale ≤ x ≤ hi/scale`. -/
structure Fix where
  lo : Int
  hi : Int
  deriving DecidableEq, Repr

namespace Fix

/-- The exact value `n / scale`. -/
def sc (n : Int) : Fix := ⟨n, n⟩

/-- The exact integer `n`. -/
def ofInt (n : Int) : Fix := ⟨n * scale, n * scale⟩

/-- An outward-rounded enclosure of `p / q`, for `q > 0`. -/
def ofFrac (p q : Int) : Fix := ⟨fdiv (p * scale) q, cdiv (p * scale) q⟩

def add (a b : Fix) : Fix := ⟨a.lo + b.lo, a.hi + b.hi⟩

def neg (a : Fix) : Fix := ⟨-a.hi, -a.lo⟩

def sub (a b : Fix) : Fix := ⟨a.lo - b.hi, a.hi - b.lo⟩

/-- Interval multiplication: outward-rounded hull of the four corner products.
The products live at scale `scale^2`, so both ends are divided back down. -/
def mul (a b : Fix) : Fix :=
  let p1 := a.lo * b.lo
  let p2 := a.lo * b.hi
  let p3 := a.hi * b.lo
  let p4 := a.hi * b.hi
  ⟨fdiv (imin (imin p1 p2) (imin p3 p4)) scale,
   cdiv (imax (imax p1 p2) (imax p3 p4)) scale⟩

/-- Division by a positive integer. -/
def divInt (a : Fix) (k : Int) : Fix := ⟨fdiv a.lo k, cdiv a.hi k⟩

/-- Reciprocal of an enclosure that stays strictly positive (`0 < lo`).

`1/y` for `y ∈ [lo/scale, hi/scale]` runs over `[scale/hi, scale/lo]`, whose
numerators at scale `scale` are `scale^2 / hi` and `scale^2 / lo`. -/
def inv (a : Fix) : Fix := ⟨fdiv (scale * scale) a.hi, cdiv (scale * scale) a.lo⟩

/-- Division, as multiplication by the reciprocal.  Sound when the divisor's
lower end is strictly positive; reusing `mul` means the corner analysis is not
repeated. -/
def div (a b : Fix) : Fix := mul a (inv b)

def pow (a : Fix) : Nat → Fix
  | 0 => ofInt 1
  | n + 1 => mul (pow a n) a

/-- The exact zero. -/
def zero : Fix := ⟨0, 0⟩

/-- The symmetric enclosure `[-r.hi, r.hi]`: everything of absolute value at
most `r.hi / scale`.  Used to absorb a series remainder. -/
def pm (r : Fix) : Fix := ⟨-r.hi, r.hi⟩

/-! ### The logarithm series

Evaluated by Horner rather than term by term.  The naive form recomputes
`y ^ i` at every term and so costs `O(n^2)` multiplications; at the fifty-odd
terms needed here that is the difference between a tenth of a second and
several seconds *per argument*, which the box covers cannot absorb.

    hornerAux y m i  =  y * (1/i + y * (1/(i+1) + ...))   -- m factors
-/
def hornerAux (y : Fix) : Nat → Nat → Fix
  | 0, _ => zero
  | m + 1, i => mul y (add (ofFrac 1 (i : Int)) (hornerAux y m (i + 1)))

/-- `∑_{i=1}^{n} y^i / i`, the partial sum whose negation approximates
`log (1 - y)`. -/
def series (n : Nat) (y : Fix) : Fix := hornerAux y n 1

/-- `2 * rho^(n+1)`.  This bounds the series remainder `|y|^(n+1)/(1-|y|)`
whenever `|y| ≤ rho ≤ 1/2`, since then `1/(1-|y|) ≤ 2`. -/
def remB (n : Nat) (rho : Fix) : Fix := mul (ofInt 2) (pow rho (n + 1))

/-- One half, exactly: the reduction window `[3/4, 3/2]` is the widest a
power-of-two shift can guarantee, and on it `|1 - z| ≤ 1/2`. -/
def half : Fix := ofFrac 1 2

/-! ### The odd series

`log x = 2·atanh z` at `z = (P-Q)/(P+Q)`, which is the same series with the
even terms cancelled.  The reduction window widens to `[2/3, 3/2]` — still
reachable by a power-of-two shift, since its ratio `9/4` exceeds `2` — and
there `|z| ≤ 1/5`, so every term buys a factor of `25` rather than `2`.  Eight
multiplications do what forty-eight did.

Horner in `z^2`, with the denominators `1, 3, 5, …` running upward. -/
def oddGo (z2 : Fix) : Nat → Nat → Fix
  | 0, _ => zero
  | m + 1, i => add (ofFrac 1 (i : Int)) (mul z2 (oddGo z2 m (i + 2)))

/-- `∑_{k<m} z^(2k+1)/(2k+1)`. -/
def oddSeries (m : Nat) (z : Fix) : Fix := mul z (oddGo (mul z z) m 1)

/-- Remainder bound for the odd series when `|z| ≤ 1/c` and `3 ≤ c`:
`2·|z|^(2m+1)/(1-|z|) ≤ 3·(1/c)^(2m+1)`. -/
def remOdd (c : Int) (m : Nat) : Fix := mul (ofInt 3) (pow (ofFrac 1 c) (2 * m + 1))

/-- `log 2`, from the odd series at `z = 1/3`, since `(1+1/3)/(1-1/3) = 2`.

This cannot go through `flogA`: `2/1` is outside the reduction window, and
reducing it needs `log 2`. -/
def log2 : Fix :=
  add (mul (ofInt 2) (oddSeries 21 (ofFrac 1 3))) (pm (remOdd 3 21))

/-- `Real.log (P/Q)` by the odd series, for `P/Q` inside `[2/3, 3/2]`. -/
def flogA (P Q : Int) (m : Nat) : Fix :=
  add (mul (ofInt 2) (oddSeries m (ofFrac (P - Q) (P + Q)))) (pm (remOdd 5 m))

/-- `2 ^ k`, spelled out so the kernel reduces it by `Int` multiplication
alone. -/
def ipow2 : Nat → Int
  | 0 => 1
  | n + 1 => 2 * ipow2 n

/-- Enclosure of `Real.log (p / q)`, reducing the argument by `2 ^ k` and
truncating the series at `n` terms.

Sound when `0 < p`, `0 < q` and the shifted argument lands in the reduction
window, `3 * q ≤ 4 * p * 2^k` and `2 * p * 2^k ≤ 3 * q`.  Those are `Int`
comparisons, so a certificate discharges them by `decide` alongside the bound
it is actually checking. -/
def flog (p q : Int) (k n : Nat) : Fix :=
  let y := ofFrac (q - p * ipow2 k) q
  sub (add (neg (series n y)) (pm (remB n half)))
      (mul (ofInt (k : Int)) log2)

/-- Logarithm of an enclosure that stays strictly positive.

`Real.log` is monotone, so the enclosure of `log` on `[lo, hi]` is bounded by
the logarithms of the two endpoints — and those are logarithms of *rationals*,
which `flog` already handles.  The two shifts are supplied separately because
the endpoints can straddle a power of two. -/
def flogI (a : Fix) (klo khi n : Nat) : Fix :=
  ⟨(flog a.lo scale klo n).lo, (flog a.hi scale khi n).hi⟩

/-- `true` only when every real in the enclosure is negative. -/
def isNeg (a : Fix) : Bool := decide (a.hi < 0)

/-- `true` only when every real in the enclosure is nonnegative. -/
def isNonneg (a : Fix) : Bool := decide (0 ≤ a.lo)

end Fix

end Spin.Numeric
