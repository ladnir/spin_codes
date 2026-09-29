/-
Fixed-point evaluators for the BA outer exponent: the *definitions*.

Import-free apart from `FixedDefs`, for the usual reason — this is the module
the box covers reduce through.

**How the reduction shift is chosen.**  `flog` needs a power-of-two shift
landing in the window `[3/4, 3/2]`, and threading that through `h`, `π` and `g`
as a parameter would mean a dozen shift arguments per box.  Instead
`shiftFor` searches for one, and the result is wrapped in `Option`: the window
conditions are *checked* by the definition, not assumed.

Nothing here proves the search succeeds, and nothing needs to.  If it finds a
good shift the guard passes and the enclosure is sound; if it does not, the
answer is `none` and the box check fails visibly.  That is the right failure
mode — a heuristic that silently returned a wrong enclosure would not be.
-/
import SpinCodes.Numeric.FixedDefs

namespace Spin.Numeric

namespace Fix

/-- Bounded search for the least `k` with `2·den ≤ 3·num·2^k`. -/
def shiftGo (num den : Int) : Nat → Nat → Nat
  | 0, k => k
  | m + 1, k => if 2 * den ≤ 3 * (num * ipow2 k) then k else shiftGo num den m (k + 1)

/-- Bounded search for the least `j` with `2·p ≤ 3·q·2^j`. -/
def shiftDownGo (p q : Int) : Nat → Nat → Nat
  | 0, j => j
  | m + 1, j => if 2 * p ≤ 3 * (q * ipow2 j) then j else shiftDownGo p q m (j + 1)

/-- The reduction window is `[2/3, 3/2]`, so an argument may need scaling in
either direction — `G(u) ≥ 1` always, so the downward case is not hypothetical.
Scale the numerator up when the argument is below the window, the denominator
up when it is above.  The window's ratio `9/4` exceeds `2`, so a power-of-two
shift always reaches it.  128 doublings is far more than
`scale = 10^30 < 2^100` can need. -/
def shiftPair (p q : Int) : Nat × Nat :=
  if 2 * q ≤ 3 * p then (0, shiftDownGo p q 128 0) else (shiftGo p q 128 0, 0)

/-- `Real.log (p/q)` for arbitrary positive `p, q`, with the two-sided shift
searched for; `none` if the search fails.

`log (p/q) = log ((p·2^k)/(q·2^j)) - k·log 2 + j·log 2`.

Nothing here proves the search succeeds, and nothing needs to.  If it finds a
good pair the guard passes and the enclosure is sound; if it does not, the
answer is `none` and the box check fails visibly.  That is the right failure
mode — a heuristic that silently returned a wrong enclosure would not be. -/
def flogKJ (p q : Int) (k j n : Nat) : Option Fix :=
  if 0 < p * ipow2 k ∧ 0 < q * ipow2 j ∧
      2 * (q * ipow2 j) ≤ 3 * (p * ipow2 k) ∧ 2 * (p * ipow2 k) ≤ 3 * (q * ipow2 j) then
    some (add (sub (flogA (p * ipow2 k) (q * ipow2 j) n) (mul (ofInt (k : Int)) log2))
              (mul (ofInt (j : Int)) log2))
  else none

/-- The same with the shift pair supplied by the search. -/
def flogQ (p q : Int) (n : Nat) : Option Fix :=
  flogKJ p q (shiftPair p q).1 (shiftPair p q).2 n

/-! ### The log oracle

Every evaluator below takes its logarithms from a function `Int → Option Fix`
rather than computing them.  Two instances matter: `directLog`, which computes
each one, and `treeLog`, which looks them up in a table verified once.

Across the dense-tail cover the same logarithm is asked for about eighteen
times, so the table is what makes the cover affordable.  Abstracting over the
source rather than duplicating the evaluators keeps one soundness proof. -/

/-- A source of logarithm enclosures: `lg p` encloses `log (p / scale)`. -/
abbrev LogFn := Int → Option Fix

/-- Compute each logarithm on demand. -/
def directLog (n : Nat) : LogFn := fun p => flogQ p scale n

/-- A lookup table keyed by numerator.

Soundness needs no search-tree invariant: `find` returns a value only at a node
whose key it has *compared equal*, so a verified node is all that is required.
The ordering is there for speed alone, and a mis-ordered tree can only make
lookups fail, never lie. -/
inductive LogTree where
  | leaf : LogTree
  | node : Int → Fix → LogTree → LogTree → LogTree

def LogTree.find : LogTree → Int → Option Fix
  | .leaf, _ => none
  | .node key val l r, k =>
      if k = key then some val
      else if k < key then LogTree.find l k else LogTree.find r k

/-- Every node carries a logarithm `flogQ` agrees with. -/
def LogTree.check (n : Nat) : LogTree → Bool
  | .leaf => true
  | .node key val l r =>
      decide (flogQ key scale n = some val) && LogTree.check n l && LogTree.check n r

def treeLog (t : LogTree) : LogFn := fun p => t.find p

/-- `Real.log` of an enclosure: the two endpoints, by monotonicity. -/
def flogIWL (lg : LogFn) (a : Fix) : Option Fix :=
  (lg a.lo).bind fun L =>
  (lg a.hi).bind fun H =>
  some ⟨L.lo, H.hi⟩

/-- The exact one. -/
def one : Fix := ofInt 1

/-! ### `x log x`

Evaluated as a unit, so that the feasibility boundary `c = a/2` — where the
factor is zero and the logarithm is `log 0` — is an ordinary point rather than
a failure.

*Upper* end: the maximum of a convex function on a box is at an endpoint, and
at `0` the endpoint value is exactly `0`.

*Lower* end: the tangent at the box midpoint, written `x·(log m + 1) - m`.
This is exact at `x = m` and needs no constant for the minimum at `1/e`.  The
factored form matters: evaluating the equivalent `x·log m + x - m` would let
`x` appear twice as independent intervals and lose an order of magnitude on a
wide box. -/

/-- `1/e`, to thirty places.  This is a *choice of expansion point*, not a
bound: `xlogx_tangent` holds at every positive `m`, so nothing about this
constant needs proving.  But it is the choice that makes the tangent bound
exact rather than merely sound — the tangent at `1/e` is horizontal, and `1/e`
is where `t log t` bottoms out. -/
def einv : Int := 367879441171442321595523770161

/-- Where to take the tangent: `1/e`, clamped into the box.  Right of `1/e`
the tangent at `lo` bottoms out at `f lo`; left of it the tangent at `hi`
bottoms out at `f hi`; and if `1/e` is inside, the tangent there is flat at
the true minimum.  So this recovers the exact lower bound in all three cases,
without a monotonicity lemma. -/
def tangentPt (a : Fix) : Int := imax (imax a.lo (imin einv a.hi)) 1

/-- `t log t` at the single rational point `p / scale`. -/
def xlPointL (lg : LogFn) (p : Int) : Option Fix :=
  if p = 0 then some zero
  else (lg p).bind fun l => some (mul (sc p) l)

/-- `t log t` over an enclosure. -/
def fxlogxL (lg : LogFn) (a : Fix) : Option Fix :=
  (xlPointL lg a.lo).bind fun L =>
  (xlPointL lg a.hi).bind fun H =>
  if a.hi = 0 then
    -- then `0 ≤ lo ≤ hi = 0`, the enclosure is the single point `0`, and
    -- `0 log 0 = 0`.  This is the feasibility boundary `c = a/2`.
    some zero
  else
    (lg (tangentPt a)).bind fun lm =>
    some ⟨(sub (mul a (add lm one)) (sc (tangentPt a))).lo, imax L.hi H.hi⟩

/-- Binary entropy `-(x log x + (1-x) log (1-x))`. -/
def fhEntL (lg : LogFn) (x : Fix) : Option Fix :=
  (fxlogxL lg x).bind fun t1 =>
  (fxlogxL lg (sub one x)).bind fun t2 =>
  some (neg (add t1 t2))

/-- The accumulator exponent in the verifier's division-free form. -/
def fpiEvalL (lg : LogFn) (a c : Fix) : Option Fix :=
  (fxlogxL lg c).bind fun t1 =>
  (fxlogxL lg (sub one c)).bind fun t2 =>
  (fxlogxL lg (sub c (divInt a 2))).bind fun t3 =>
  (fxlogxL lg (sub (sub one c) (divInt a 2))).bind fun t4 =>
  (fxlogxL lg (sub one a)).bind fun t5 =>
  some (add (sub (sub (add (add (mul a log2) t1) t2) t3) t4) t5)

/-- The Golay weight enumerator `G(u)`.

In `v = u^4` this is `1 + 759v² + 2576v³ + 759v⁴ + v⁶`, and Horner in `v`
evaluates it in six multiplications where `u^8, u^12, u^16, u^24` cost sixty.
`u` is always a point here, so repeating `v` costs no width. -/
def fGolayG (u : Fix) : Fix :=
  let u2 := mul u u
  let v := mul u2 u2
  let v2 := mul v v
  add one (mul v2 (add (ofInt 759)
    (mul v (add (ofInt 2576) (mul v (add (ofInt 759) v2))))))

/-- The objective `(1/24) log G(u) - a log u` whose infimum is `g(a)`. -/
def fgObjL (lg : LogFn) (u a : Fix) : Option Fix :=
  (flogIWL lg (fGolayG u)).bind fun lgG =>
  (flogIWL lg u).bind fun lu =>
  some (sub (divInt lgG 24) (mul a lu))

/-! ### The centered bound

`π` on a box, expanded about the box midpoint, with the two remainders carried
by enclosures of the partial derivatives over the whole box.  This is the
verifier's `p_taylor_upper`, and it is the one bound the cover cannot do
without — the plain interval extension of `π` does not converge.

The derivative expressions are written so that each interval variable appears
once; the dependency lesson from `fxlogx` applies here too. -/

/-- The midpoint of an enclosure, as an exact point. -/
def midPt (A : Fix) : Fix := sc (fdiv (A.lo + A.hi) 2)

/-- `∂π/∂a = ln 2 + (ln(c - a/2) + ln(1 - c - a/2))/2 - ln(1 - a)`. -/
def fDpaL (lg : LogFn) (A C : Fix) : Option Fix :=
  (flogIWL lg (sub C (divInt A 2))).bind fun l1 =>
  (flogIWL lg (sub (sub one C) (divInt A 2))).bind fun l2 =>
  (flogIWL lg (sub one A)).bind fun l3 =>
  some (sub (add log2 (divInt (add l1 l2) 2)) l3)

/-- `∂π/∂c = ln c - ln(1-c) - ln(c - a/2) + ln(1 - c - a/2)`. -/
def fDpcL (lg : LogFn) (A C : Fix) : Option Fix :=
  (flogIWL lg C).bind fun l1 =>
  (flogIWL lg (sub one C)).bind fun l2 =>
  (flogIWL lg (sub C (divInt A 2))).bind fun l3 =>
  (flogIWL lg (sub (sub one C) (divInt A 2))).bind fun l4 =>
  some (add (sub (sub l1 l2) l3) l4)

/-- The centered bound for `π` on the box `A × C`. -/
def fpiCenteredL (lg : LogFn) (A C : Fix) : Option Fix :=
  (fpiEvalL lg (midPt A) (midPt C)).bind fun p0 =>
  (fDpaL lg A C).bind fun da =>
  (fDpcL lg A C).bind fun dc =>
  some (add (add p0 (mul da (sub A (midPt A)))) (mul dc (sub C (midPt C))))

/-! ### Intersecting two enclosures

Two sound enclosures of the same real can be intersected.  The cover uses this
to take whichever of the plain and centered bounds is better on a given box,
without having to decide in advance which that is. -/

def meet (a b : Fix) : Fix := ⟨imax a.lo b.lo, imin a.hi b.hi⟩

/-- The feasibility guard of the centered bound, as a decision. -/
def feasOK (A C : Fix) : Bool :=
  decide (0 < 2 * C.lo - A.hi) && decide (0 < 2 * scale - 2 * C.hi - A.hi)
    && decide (A.hi < scale) && decide (0 < C.lo) && decide (C.hi < scale)

/-- The centered bound with its guard checked rather than assumed. -/
def fpiCenteredG (lg : LogFn) (A C : Fix) : Option Fix :=
  if feasOK A C then fpiCenteredL lg A C else none

/-- The better of the plain and centered bounds, whichever are available.

Boxes touching the feasibility boundary have no centered bound — the
derivatives are not defined there — but they do have the plain one, which
`fxlogx` handles because `0 log 0` is an ordinary value for it.  Boxes well
inside get both, and their intersection. -/
def fpiBestL (lg : LogFn) (A C : Fix) : Option Fix :=
  match fpiEvalL lg A C, fpiCenteredG lg A C with
  | some p, some q => some (meet p q)
  | some p, none => some p
  | none, some q => some q
  | none, none => none

/-! ### Tables split across modules

A single table literal does not scale: 23,541 entries elaborate to a 246 MB
`.olean` and the `decide` over them reached 16 GB resident before being
stopped.  Splitting the entries across several `LogTree`s, each in its own
module, keeps every literal small enough to elaborate and check in bounded
memory.

The trees are not merged — the *oracles* are composed.  A lookup walks the
chunk list until one answers, which costs a constant factor and nothing in
soundness: `treeLogN_find` still ends at a node whose key was compared equal
in some verified chunk. -/

def treeLogN : List LogTree → LogFn
  | [], _ => none
  | t :: ts, p => (t.find p).orElse (fun _ => treeLogN ts p)

/-- Every chunk's entries are the logarithms they claim to be. -/
def checkAll (n : Nat) : List LogTree → Bool
  | [] => true
  | t :: ts => t.check n && checkAll n ts

/-! ### The split tree

The certificate is the branch-and-bound's own tree, not a flat list of leaves.
That matters for the proof rather than for the data: "these boxes cover the
original" is then structural induction — at each node the two children differ
in one coordinate and every point falls into one of them — instead of a
tiling argument over thousands of boxes. -/

inductive BoxTree where
  | leaf : Int → BoxTree
  | split : Nat → Int → BoxTree → BoxTree → BoxTree

/-- Narrow a coordinate to the part below the cut. -/
def setHi (X : Fix) (m : Int) : Fix := ⟨X.lo, m⟩

/-- Narrow a coordinate to the part above the cut. -/
def setLo (X : Fix) (m : Int) : Fix := ⟨m, X.hi⟩

/-- Walk the tree, narrowing the box at each split and testing every leaf. -/
def checkTree (leafOK : Int → Fix → Fix → Fix → Bool) :
    BoxTree → Fix → Fix → Fix → Bool
  | .leaf u, A, B, W => leafOK u A B W
  | .split 0 m l r, A, B, W =>
      checkTree leafOK l (setHi A m) B W && checkTree leafOK r (setLo A m) B W
  | .split 1 m l r, A, B, W =>
      checkTree leafOK l A (setHi B m) W && checkTree leafOK r A (setLo B m) W
  | .split _ m l r, A, B, W =>
      checkTree leafOK l A B (setHi W m) && checkTree leafOK r A B (setLo W m)

/-! ### Clamping to the feasible part

A box may straddle the feasibility boundary: `c - a/2` then runs negative over
part of it, and `xlogx` of a negative enclosure is simply unavailable — the
logarithm does not exist there.  Without clamping, every such box is
unbounded, and since the root box `β ∈ [0,1]` straddles at almost every `α`,
the search never closes.

The feasible points all have `c - a/2 > 0`, so the lower end may be raised to
zero.  That is exactly what the original verifier's
`max(ZERO, b_lo - a_hi/2)` does.  Soundness then needs the feasibility
hypothesis, which is precisely what `DenseClaim` supplies. -/

/-- Raise the lower end to zero. -/
def clampLo (X : Fix) : Fix := ⟨imax 0 X.lo, X.hi⟩

/-- `π` on a box, with the two derived quantities clamped to their feasible
part.  Sound under the feasibility hypothesis only. -/
def fpiEvalClampL (lg : LogFn) (A C : Fix) : Option Fix :=
  (fxlogxL lg C).bind fun t1 =>
  (fxlogxL lg (sub one C)).bind fun t2 =>
  (fxlogxL lg (clampLo (sub C (divInt A 2)))).bind fun t3 =>
  (fxlogxL lg (clampLo (sub (sub one C) (divInt A 2)))).bind fun t4 =>
  (fxlogxL lg (sub one A)).bind fun t5 =>
  some (add (sub (sub (add (add (mul A log2) t1) t2) t3) t4) t5)

/-- The better of the clamped plain bound and the centered bound. -/
def fpiBestClampL (lg : LogFn) (A C : Fix) : Option Fix :=
  match fpiEvalClampL lg A C, fpiCenteredG lg A C with
  | some p, some q => some (meet p q)
  | some p, none => some p
  | none, some q => some q
  | none, none => none

/-! ### The dense-tail leaf test

A leaf of the cover is accepted for one of two reasons, and the certificate
has to be able to give either.

*The box misses the feasibility region.*  `π` is `-∞` off the region
`a/2 ≤ c ≤ 1-a/2`, so the claim is vacuous there.  `infeasible` is the
negation of the verifier's `intersects_feasible_region`, and it is four `Int`
comparisons.

*The bound is below the threshold.*  `g` at the leaf's own witness `u`, plus
the two `π` bounds, with `fpiBestL` choosing per box. -/

/-- No point of the box satisfies both accumulator constraints. -/
def infeasible (A B W : Fix) : Bool :=
  decide (2 * B.hi < A.lo) || decide (2 * scale - A.lo < 2 * B.lo)
    || decide (2 * W.hi < B.lo) || decide (2 * scale - B.lo < 2 * W.lo)

/-- The leaf test.  `thr` is the threshold numerator: the paper's
`-7.68·10⁻⁸` is `thr = -76800000000000000000000` at `scale = 10^30`. -/
def leafOK (lg : LogFn) (thr : Int) (u : Int) (A B W : Fix) : Bool :=
  infeasible A B W ||
    (decide (0 < u) &&
      (match fgObjL lg (sc u) A, fpiBestClampL lg A B, fpiBestClampL lg B W with
       | some g, some p1, some p2 => decide ((add (add g p1) p2).hi < thr)
       | _, _, _ => false))

/-! ### The direct instances

These are the shapes the soundness pins name; each is the corresponding
oracle-taking evaluator at `directLog n`. -/

def flogIW (a : Fix) (n : Nat) : Option Fix := flogIWL (directLog n) a
def xlPoint (p : Int) (n : Nat) : Option Fix := xlPointL (directLog n) p
def fxlogx (a : Fix) (n : Nat) : Option Fix := fxlogxL (directLog n) a
def fhEnt (x : Fix) (n : Nat) : Option Fix := fhEntL (directLog n) x
def fpiEval (a c : Fix) (n : Nat) : Option Fix := fpiEvalL (directLog n) a c
def fgObj (u a : Fix) (n : Nat) : Option Fix := fgObjL (directLog n) u a
def fpiCentered (A C : Fix) (n : Nat) : Option Fix := fpiCenteredL (directLog n) A C

end Fix

end Spin.Numeric
