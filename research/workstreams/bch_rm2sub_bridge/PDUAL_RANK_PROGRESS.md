# P-dual rank refinement: case 15 closed, case 7 remains

Follow-up: `PDUAL_CASE7_SEARCH.md` records bounded case-7 searches and the
exact fixed-degree-five reduction. Case 7 remains open.

The new exact witnesses close Fourier case 15. Case 7 remains at rank 30,
so the distance-32 claim for the complete P dual is still unproved. The
unconditional K=2^28 BCH/RM2Sub certificate remains at 42.510028 bits.
The 43.585033-bit result in `LOW_SHELL_ROI.md` remains conditional.

There is a new narrower theorem: an explicitly defined dimension-117
subcode of the extended P dual has minimum distance at least 32. Consequently,
every possible weight-30 word of the full P dual has nonzero seventh Fourier
value. This does not exclude those words or give their number.

## Construction and new exact statement

Let P0 be the primitive binary BCH code of length 255 and designed distance
37 used in the existing certificates. Let P be its even extension, of
dimension 131. Coordinates of length-256 words are indexed by F256, with
the extra coordinate at zero. The field uses modulus 0x14D and primitive
element alpha=2, as in the frozen construction.

For a binary vector c=(c_x) indexed by F256, define

\[
F_7(c):=\sum_{x\in\mathbb F_{256}}c_x x^7,
\qquad
E:=\{c\in P^\perp:F_7(c)=0\}.
\]

**Checked statement.** The binary code E has length 256, dimension 117,
and minimum distance at least 32. The independently audited complete cover
still proves only minimum distance at least 30 for P dual itself.

The map F7 on P dual has binary rank eight. The checker verifies this by
eliminating its values on a basis of P dual, which has dimension 125.
Thus its kernel E has dimension 117. This is a specific subcode statement,
not a replacement BCH outer or a change to the encoder.

## What changed in the rank search

The old witness search used one fixed nonzero Fourier value as its append
pivot. The new search can alternate among several values that a branch
has already required to be nonzero. It also starts from prefixes of checked
old witnesses, instead of reconstructing every prefix from rank one.

The successful case-15 partition is exhaustive:

| Case-15 branch | Checked rank | Result |
|---|---:|---|
| F(23)=0 | 31 | Closed |
| F(23) is nonzero | 31 | Closed; append pivots 15 and 23 both used |

Here F(e)=sum_(j=0)^254 c_j alpha^(ej) for a word in P0 dual. In case 15,
F(15) is nonzero and every earlier allowed cyclotomic coset has value zero.
The nonzero branch reuses a valid prefix from the earlier Q-dual refinement.
The checker verifies the new path against the P-case assumptions directly;
it does not infer a P-specific bound merely from the source file's name.

The retained case-7 partition has three leaves:

| Case-7 branch | Checked rank | Status toward rank 31 |
|---|---:|---|
| F(15)=0 and F(23)=0 | 31 | Closed |
| F(15)=0 and F(23) is nonzero | 30 | Open |
| F(15) is nonzero | 30 | Open |

All fourteen original cases after case 15 already have sufficient checked
witnesses. The complete refined partition therefore has 19 leaves.
Its minimum rank is 30; after excluding case 7, its minimum rank is 31.
The independent audit reconstructs the first-nonzero partition from the
generator polynomial and checks both children of every split.

## Why a multi-pivot witness proves a rank lower bound

Fix a binary word c in a branch and let S be its support. For each exponent e,
write v_e=(alpha^(ej)) indexed by j in S. These vectors lie in F256^|S|.
Their coordinate sums are F(e).

Start with v_b for a value F(b) known to be nonzero. A stored step (k,s,p)
first transforms every exponent e to 2^k e+s modulo 255. Applying the field
automorphism to every coordinate preserves linear independence. Multiplying
each coordinate j by alpha^(sj) also preserves linear independence.

The checker requires every transformed exponent to lie in the branch's
zero set. Every vector in their span then has coordinate sum zero.
The next pivot p must have F(p) nonzero, so v_p is outside this span.
Appending v_p therefore increases the independent set size by one.
The pivot may change between steps; its nonzero premise must hold in that
same branch. A checked rank-31 path implies |S|>=31, and evenness implies
|S|>=32.

Beam width, tie-breaking, and time limits affect discovery only. The verifier
uses sets, finite-field operations, and exact integer comparisons. It neither
trusts the search score nor assumes that the beam exhausts all valid paths.

## Transfer of the kernel result to length 256

For the punctured dual, F(7)=0 excludes the first Fourier case. Every remaining
case now has rank at least 31. Its words are even, so the punctured kernel
has minimum distance at least 32.

The checker then constructs a basis of P dual from the 124 punctured dual
generator shifts, each with zero extra coordinate, and the length-256 all-one
word. It checks their independence, evenness, and orthogonality to P.

The already checked affine coordinate generators preserve P and hence P dual.
On the new dual basis, exact field calculations verify

\[
F_7(T_{a,b}c)=a^7F_7(c)
\quad\text{for }(a,b)=(1,1),(2,0),
\]

where T_(a,b) moves coordinate x to ax+b. Linearity extends these identities
to every word of P dual. Both generators therefore preserve E. Since 2 is
primitive, they generate all translations and nonzero scalings.

Suppose E had a nonzero word of weight below 32. It has a zero coordinate.
Translate that coordinate to the extra position and puncture it. The result
is a nonzero word of P0 dual with the same weight and F(7)=0. This contradicts
the punctured kernel bound. Hence the extended kernel E has distance at least
32, as claimed.

## Bounded effort and negative findings

This turn ran 21 serial rank searches, with about 249 seconds of total
reported search time. Three saved witnesses reached rank 31. Increasing the
beam to 5000 while starting at rank 20 spent its budget before reaching the
known rank-30 endpoint. Those timed-out outputs do not replace stronger
retained witnesses.

Further splits on cosets 23 and 29 did not close the remaining case-7
branches. The search stopping message `no further state` refers only to the
retained beam. It is not an exhaustive impossibility result.

A separate small counterexample check exhausted combinations of up to three
of the 124 punctured dual generator shifts. The smallest weight found was
52. This is a restricted search, not evidence of the absence of weight 30
among general linear combinations. No weight-30 counterexample was found.

## Next work

Focus on the residual case F(7) nonzero. Preserve the completed case-15 tree
and do not rerun all sixteen original cases. A useful next bounded attempt
would combine a different rank argument or case refinement with a targeted
search for an actual weight-30 word. Any counterexample must pass direct
code-membership and weight checks.

If rank-only refinements continue to stall, the alternative is a counting
bound for the residual weight-30 shell, or an LP extension using the newly
proved kernel E. Neither alternative has a measured payoff yet. The current
result does not justify applying the zero-weight-30 LP premise, and it does
not validate a reference-spectrum heuristic.

## Artifacts and replay

`generated/pdual_partial_cover_v1.json` embeds the complete 19-leaf witness
forest, the audit, and source hashes in about 105 KB. Replay needs neither
the old search directory nor an optimizer:

```text
python -B pdual_case_audit.py --output generated/pdual_partial_cover_v1.json --verify
python -B -m unittest test_pdual_refinement test_low_shell_roi test_curve_spectrum
```

The tests check recorded paths, exhaustive split assumptions, rejection of
tampered ranks and missing branches, and valid multi-pivot transitions on
all 31 nonzero length-five binary words. Discovery receipts remain under
`generated/pdual_refinement_v1/`. All changes are in this workstream; frozen
certificate sources and the other agent's worktree were not modified.
