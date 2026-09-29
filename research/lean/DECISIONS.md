# Decision log

## D1 — Interval arithmetic: build a purpose-built `ln` enclosure, do not depend on `girving/interval`

**Date:** 2026-09-23 (task T5)

**Question.** The BA tail cover (`eq:structured-ba-dense-tail`, 6,747 boxes), the
concave majorant (`eq:ba-spectrum-majorant`, 10,721 box checks) and the
positive-occupation cover all need rigorous enclosures of expressions built from
`h(x) = -x ln x - (1-x) ln(1-x)` on rational boxes. Do we import an interval
arithmetic library or build what we need?

### What was evaluated

**`girving/interval`** — conservative interval arithmetic in Lean, Apache-2.0,
actively maintained (last push 2026-08-10), 64+64-bit *software* floating point
(Lean's `Float` is untrusted, so they implemented their own), with `Approx`
soundness typeclasses and conservative `exp`, `log`, `pow`, `sqrt`, `sincos`,
plus complex `Box`. It is the right library for this class of problem and it is
built by someone who used it for verified Mandelbrot bounds.

**Mathlib** has no outward-rounded interval arithmetic. But it does have the
piece that is actually hard:

```
Real.abs_log_sub_add_sum_range_le (h : |x| < 1) (n : ℕ) :
  |(∑ i ∈ range n, x^(i+1)/(i+1)) + log (1 - x)| ≤ |x|^(n+1) / (1 - |x|)
```

— an explicit, proved remainder bound for the log series. It also has
`log_two_gt_d9` / `log_two_lt_d9` for argument reduction.

### Decision

**Build a purpose-built rational `ln` enclosure on Mathlib's remainder bound.**

Reasons, in order of weight:

1. **Version gap.** `girving/interval` is pinned to `leanprover/lean4:v4.27.0-rc1`;
   this project is on `v4.34.0`. Adopting it means either downgrading the whole
   project onto a release candidate seven versions back, or porting the library
   across seven versions of Mathlib API churn. Neither is proportionate.
2. **Scope mismatch.** This project needs `ln` on rationals in a bounded range.
   It does not need `exp`, `sin`, `cos`, `sqrt`, complex boxes, or a
   floating-point layer. A general library is mostly surface we would not use
   but would have to maintain.
3. **Mathlib supplies the hard part.** With the remainder bound above plus
   argument reduction (`ln q = k ln 2 + ln r`, `r ∈ [2/3, 4/3]`, so `|1-r| ≤ 1/3`
   and the series converges like `3^{-n}`), the enclosure is arithmetic on `ℚ`.
   `ln 2` itself comes from the same series at `x = 1/2`.
4. **Trust surface.** No third-party dependency; the axiom base stays exactly
   Mathlib's, which `scripts/check.sh` already enforces.

### What would change this

If the project later needs `exp`, powers, or trigonometric enclosures at scale —
the IMT transfer bounds (T9) involve `z^k` for real `z`, and the continuum limit
(T11) may need more — then porting `girving/interval` becomes the better trade,
because reimplementing those is no longer a small job. Revisit at T9.

### Evidence

`SpinCodes/Numeric/LogBounds.lean` implements the two-sided form and
demonstrates the precision actually reached; see the worked bounds there.
Precision required by the certificates is around `1e-9`
(the majorant's largest residual is `-1.70e-9`), which the series reaches with a
few tens of terms at `|x| ≤ 1/3`.


## D2 — Box covers must use integer fixed-point, not `ℚ`: `decide` cannot reduce rational arithmetic

**Date:** 2026-09-23 (task T7)

**Finding.** The Lean kernel cannot evaluate `ℚ` arithmetic through `decide`.
Measured directly:

| expression | `decide` |
|---|---|
| `Nat.gcd 123456 789012 = 12` | reduces |
| `(123456789 : Int) * 987654321 - 1 < …` | reduces |
| `Rat.blt 1 2 = true` | reduces |
| `(1/3 + 1/6 : ℚ) = 1/2` | **stuck** |
| `(2 : ℚ) + 3 = 5` | **stuck** |

So it is not `Nat.gcd` that blocks (that reduces fine) — it is `ℚ`'s
arithmetic instances. Integer arithmetic reduces without trouble, which is why
the sparse Collatz certificate (T3) worked: it was deliberately scaled to `ℤ`.

**Consequence.** The plan for T7 — computable `ℚ` interval arithmetic, whole
cover discharged by one `decide` — does not work. The options are:

1. **Integer fixed-point with directed rounding.** Represent values as `ℤ`
   numerators over a fixed implicit denominator, with outward rounding on
   multiplication and division. Kernel-reducible, and the same technique that
   already succeeded in T3.
2. `norm_num` per box. Rejected: 6,747 + 10,721 boxes, each needing several
   log enclosures, is hours of elaboration and makes the project unbuildable
   in practice.
3. Compiled evaluation. Rejected: enlarges the axiom closure, which
   `scripts/check.sh` exists to prevent.

**Decision: option 1.**

**This partially re-opens D1.** D1 judged `girving/interval` not worth its
version cost because "Mathlib already has the hard part" — the log series
remainder bound. That was true of the *mathematical* layer, and `LogBounds`
plus `RatLog` deliver it. But the box covers need a *computable* layer, and
building one means fixed-point arithmetic with directed rounding — which is
precisely the bulk of what `girving/interval` is, and precisely why it ships
its own software floating point instead of using `ℚ`.

D1's reasoning stands for what it covered (the one-off constants of T6h and
T6j, which `norm_num` handled comfortably). Its cost estimate for the covers
was wrong. Revisit the import-vs-build question before T9, as D1 already
flagged, now with this extra evidence on the build side of the ledger.

### Outcome (T7a, 2026-09-23)

Option 1 works. `Numeric/FixedDefs.lean` + `Numeric/Fixed.lean` implement it,
and the kernel reduces the whole thing through `decide` — verified on
`ofFrac`, `pow`, and an 8,000-element list check.

Measured throughput: **~21 ms per box** at 8 operations (8,000 boxes, five
multiplications and three subtractions each, 2m49s). Enough for the covers
*if* the log enclosures are tabulated rather than recomputed per box; a log is
roughly 75 operations, which at this rate would put a single cover over half
an hour. See the note on T7a2.

## D3 — If a certificate table is too big for one module, chunk it and compose the oracles

**Date:** 2026-09-23 (task T7b3h)

**Observation.** The dense-tail logarithm table has 23,541 entries. As a
single `LogTree` literal it is a 2.78 MB source file, and checking it has so
far taken 80 CPU-minutes — against a 17-minute projection from the per-lookup
and per-computation costs measured in isolation. The small cover's 3,530-entry
table produced a **36.9 MB** `.olean`, about 10 KB per entry, which is the
warning sign: cost per entry is not constant, it grows with the size of the
enclosing term.

**Fallback, if the single-module table proves impractical.** Split the table
across several modules and compose the oracles rather than the trees:

```
def treeLog2 (t u : LogTree) : LogFn := fun p => (t.find p).orElse (fun _ => u.find p)
```

with `Oracle (treeLog2 t u)` following from `t.check n` and `u.check n`
separately. A `k`-way split costs `O(k)` extra per lookup, which at `k = 4` is
nothing against a depth-15 descent, and it buys three things:

1. per-module elaboration cost bounded, so the build reports progress instead
   of disappearing for hours;
2. a failure localises to one chunk;
3. `.olean` sizes stay in a range the toolchain handles comfortably.

This is the same move that `Structured/PolyCertDefs.lean` and the seven
`SparseData/Row*.lean` files made for T3, for the same reason, and it should
probably be the default shape for any table beyond a few thousand entries —
including T7c's.

**Not applied yet.** The current build is progressing steadily (CPU accruing,
568 MB resident, no thrash), and editing any Lean source now would invalidate
the `.olean` it is producing and discard the work already done. Decide once it
finishes or fails.

### D3, confirmed

`DenseTailLogs.olean` came out at **246,553,752 bytes** — 10.5 KB per entry,
exactly the scaling the small cover's 36.9 MB predicted. And the 93 minutes it
took bought no checking at all: that module holds only `def logTable`, while
`logs_ok` lives in `DenseTail.lean`. Ninety-three minutes to elaborate and
serialise a literal.

Separating the two costs clarifies what is actually unavoidable:

* **Checking** the table is 23,541 logarithm evaluations at ~42 ms, once in
  the elaborator and once in the kernel — about 33 minutes. That is a floor,
  and it is acceptable.
* **Carrying** the table as one 2.78 MB literal is the part that is not, and
  it is what D3's chunking removes.

So the chunked design should be adopted whether or not the current run
finishes: it does not reduce the floor, but it removes an overhead several
times larger than the floor, and it bounds memory — the module now running is
at 2.3 GB and climbing, against 568 MB for the table alone.

### D3, applied — and the tree needs splitting too

The monolithic run was stopped after both data modules had built
(`11,766 s` each) and `DenseTail.lean` reached **16.3 GB resident**. That is
not a slow build, it is a build that was going to take the machine down, so it
was killed rather than waited out.

The memory tells us something the timings did not. It is not only the table
literal: `tree_ok` is a single `decide` over 13,513 leaves, and that is its own
unbounded object. The small cover's 1,416-leaf tree checked in 831 s at modest
memory, so the working unit is on the order of a thousand leaves, not fourteen
thousand.

So the chunking applies twice:

* **The table** splits into several `LogTree`s, composed by `treeLogN`, which
  walks the chunk list until one answers. Soundness is unchanged —
  `treeLogN_find` still ends at a node whose key was compared equal inside a
  chunk that `checkAll` verified. Implemented.
* **The tree** splits at a fixed depth, each subtree proved in its own module
  against its own narrowed box, and the root assembled from them. `checkTree`
  on a `.split` is definitionally the conjunction of its children, so the
  assembly is a rewrite rather than a new argument. Not yet implemented.

The second is the one that matters for memory, and it also makes the build
report progress instead of disappearing: sixteen modules of about 850 leaves
each, cached independently, rather than one four-hour all-or-nothing step.

## D4 — The majorant cover is blocked on cost, and what would unblock it

**Date:** 2026-09-24 (task T7c)

**Measured.** Driving the paper's own verification — a 3-D interval cover of
`g(a) + π(a,c) + π(c,x) - (s·x + intercept) < 0` per support, over the
support's active `x`-interval — with the same `Fix` bound Lean checks:

    segment 14 (slope 0.05)    closed: 64,107 nodes, 32,054 leaves
    segment  0 (slope 1.666)   >212,000 nodes, did not close

Segment 14 alone is 2.4 times the entire dense-tail cover. Fifteen left
segments plus the central one would put the total in the high hundreds of
thousands of leaves if every segment cost what segment 14 did. See the
correction below: that extrapolation is not trustworthy, and the figures it
produced (three days of compute, 18 GB of `.olean`) should not be relied on.

**Two things this is not.** It is not a defect in the paper — the certificate
reports `status: proved` over the same region. And it is not a looseness in the
`Fix` bound: on the dense tail that bound reproduced the mpmath partition leaf
for leaf (13,513 both ways).

**Why it is expensive.** Superseded — see the correction below. The
near-tangency explanation recorded here was wrong.

**Correction (2026-09-24, same day, after re-measurement).** Three claims above
were wrong, and they compounded in the same direction — all three inflated the
apparent cost, which is what made an external server look necessary.

1. *"Nineteen left supports."* The certificate has **fifteen** `left_segments`.
   The leaf-count extrapolation was scaled by the wrong factor.

2. *"The certificate's own numbers are consistent with what I measure."* They
   are not, and the comparison was unit-confused. The certificate's **101,940
   boxes cover the entire claim** — every segment plus the central constant.
   My 32,054 leaves covered **one segment**. The published search is therefore
   roughly thirty times more efficient than the one here, which inverts the
   conclusion: the cost is an artifact of my search, not a property of the
   problem.

3. *"The supports are near-tangent, so the supremum is attained in the
   interior."* Re-measuring the true margin on every segment gives roughly
   `2e-5` throughout, tightest at segment 4 (`-1.695e-5`). Segment 0, the one
   that would not close, has the **most** slack (`-9.94e-5`). Cost does not
   track tightness at all, so near-tangency does not explain it. The
   certificate's `largest_accepted_upper_natural = -2.2e-9` is an
   accepted-box artifact, not evidence of a near-tangency.

Also measured, and relevant: the `Fix` bound has **no precision floor** here —
shrinking a box around the argmax drives the enclosure width to ~`8e-13`, far
below the margin. So the blockage is not precision either.

**What the certificate does not contain.** The file is 3,702 bytes and records
only summary statistics: the fifteen segments' slopes and intercepts, the box
counts, `status: proved`. **The boxes themselves are not stored, and neither is
the search procedure.** So the published certificate cannot be re-checked from
its own contents, and there is no cover available to transcribe into a
`BoxTree`. This is precisely the gap the Lean replay is meant to close, and it
means the cover must be regenerated here rather than imported.

**Revised diagnosis.** T7c is a search-engineering problem, not a compute
problem and not a mathematical one. A working cover of about 100k boxes
demonstrably exists.

Two leads were then tested on segment 0, instrumented with *resolved volume*
(the summed volume of closed and pruned boxes — exact progress for a DFS, where
node counts say nothing):

* **Feasibility-tightened start** — begin from `c ≤ 2·max(x)`, `a ≤ 2·max(c)`
  rather than `[0,1]^2`. Large throughput gain (infeasible boxes prune
  cheaply), **no** convergence gain.
* **Raised depth cap** — `MAJ_MAXDEPTH=120`. **No effect whatsoever.** The
  observed depth stays pinned at 56 under both caps, so `MAXDEPTH` was never
  binding and its proximity to the certificate's `deepest_box = 55` was
  coincidence.

What the volume trace actually shows:

    nodes= 40000  resolved=0.222412  depth<=56  eta=0.18h
    nodes= 60000  resolved=0.320357  depth<=56  eta=0.16h
    nodes= 80000  resolved=0.320640  depth<=56  eta=0.25h
    nodes=100000  resolved=0.321095  depth<=56  eta=0.29h
    nodes=120000  resolved=0.321887  depth<=56  eta=0.37h

Progress **stalls at 32%** while the ETA climbs monotonically. Extrapolating
the remaining 68% at the observed rate gives ~1.7e7 nodes, about ten days at
~20 nodes/s — and because the rate is decelerating that is a floor, not an
estimate. The true shape is boxes that fail to certify until they are ~`2^-56`
across, i.e. on the order of `1e16` of them. **The search does not terminate in
practice, at any budget.**

**The gap is in the leaf bound, and it is enormous.** ~`1e16` boxes here versus
101,940 for the published cover: about eleven orders of magnitude on the same
inequality. No split heuristic closes a gap of that size — three were tried
(`width`, `trial`, `collapse`) and all failed. The per-box bound in
`majorant.py` is the natural interval extension of `g(a) + π(a,c) + π(c,x)`,
which loses badly to dependency: `c` occurs in two `π` terms and `a` in two
places, so each variable's range is inflated independently in every occurrence.
A mean-value / centred form, or exploiting monotonicity in `x` to evaluate at
an endpoint rather than over the interval, is where the missing orders of
magnitude must come from. Untested.

**Not yet ruled out:** that the published cover was produced from a different
formulation of the inequality (for instance the two-dimensional exchanged-suprema
form below), in which case 101,940 is not a target this search can reach at all
and the reduction becomes necessary rather than optional.

**What would unblock it.** Exchanging the suprema,

    sup_x [a_BA(x) - s·x] = sup_{a,c} { g(a) + π(a,c) + sup_x [π(c,x) - s·x] }

makes the inner supremum one-dimensional in `x`, with a closed form reachable
through the generating function `∑_j C(j,m) z^j = z^m/(1-z)^{m+1}`. That drops
the cover to two dimensions and replaces the tangency in `x` with an exact
maximisation.

It is deliberately **not** being done here. The loop's job is to formalize the
paper's argument, and this would substitute a different one. Recording it so
the option is on the record rather than lost.

**Consequence.** T7c is marked BLOCKED and the queue moves to T8. Everything
`eq:ba-spectrum-majorant` feeds into stays unformalized; that is stated rather
than papered over.


## D5 — T7c search diagnosis corrected (2026-09-25)

The original three-dimensional method is practical. The prior nontermination
claim and the `1e16`-box extrapolation in D4 were not justified. Search depth
counts splits across coordinates; it does not say that every coordinate has
width `2^-depth`. Resolved volume also does not give a reliable linear ETA
when verification cost varies strongly across the region.

The original generator was present at
`research/workstreams/paper_architecture/certificates/single_sampled_ba_rm2sub/certify_golay_ba_concave_majorant.py`.
Its `Box.taylor_upper` centers the **whole objective**. The stalled Lean search
omitted this bound. Adding its fixed-point counterpart retains cancellation
in `-log(u)+Dpa`, `Dpc+Dpa`, and `Dpc-slope`. The new search also uses the
original absolute-width split rule and memoizes exact logarithm evaluations.

The paper uses `imt_asymptotic/d11/OUTER_REFINED.json`, with 39 active supports;
the old Lean search loaded the earlier 31-support file. The refined artifact
stores partition witnesses for its five new left supports. The inherited
supports still need their covers regenerated. The new generator regenerates
all 19 active nonconstant left supports, rounds segment endpoints outward,
and records every split and generating-function witness.

Result: 108,045 tree nodes, 54,032 leaves, maximum depth 56, approximately
22 seconds of search. All 387 independent numerical modules passed kernel
checking. The combined bound, closed-boundary identities, analytic central
bound, exact active-line comparisons, and reflection have also been proved
in Lean. No exchanged-suprema reduction or change to the paper was needed.

The old claims that there was no source generator, that this search could
not terminate at any budget, and that a different mathematical argument might
be necessary must not be used as current status. See `STATUS.md` for the
final assembly check and `CLOSURE_AUDIT.md` for the remaining theorem inputs.
