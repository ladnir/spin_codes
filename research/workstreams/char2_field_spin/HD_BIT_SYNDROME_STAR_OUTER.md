# Rank-eight star-quadratic outer: exact local checks

This note records two new local mixer candidates. Neither has a regional
distance provider, a native implementation, or a whole-code certificate.
The existing H* implementation and all its receipts remain unchanged.
The coupled binary coefficient variant is the preferable candidate from
these local checks; its runtime has not been measured.

The separate [nonbinary syndrome proposal](HD_NONBINARY_SYNDROME_PROPOSAL.md)
is exploratory. It replaces one star check by a field-valued quadratic and
proves quiet minimum weight at least six, at a larger untimed circuit cost.
Its nonquiet spreading has not been checked; it does not supersede either
of the tested local variants recorded here.

## Construction and quiet subspace

Use the AES field F of order 256 and the existing sixteen labels
`u_i = 06 + a0 + 0c*a1 + 50*a2 + b0*a3`, where `a` is the four-bit
expansion of i. Let V have rows `(1,a0,a1,a2,a3)`, and let B be the
existing invertible five-by-five coefficient matrix of H*:

`H* = I + V B V^T`.

Add the three columns `f=(a0*a1,a0*a2,a0*a3)` and put W=[V,f].
For a symmetric invertible binary matrix C, define

`H_C = I + W diag(B,C) W^T = H* + f C f^T`.

The two checked choices are

```
identity:  1 0 0       coupled:  0 1 0
           0 1 0                 1 0 1
           0 0 1                 0 1 1
```

Both C and I+C are invertible for the coupled choice. All products of two
columns of W have binary degree at most three: the quadratic columns share
a0. Summing any such monomial over all sixteen labels gives zero. Thus
`W^T W=0` and `H_C^2=I`. The checked coefficient determinant is 0x19.

W spans arbitrary affine functions separately on the two a0 halves.
Consequently it has rank eight, and the quiet subspace is

`ker(W^T) = ker(affine_3^T) on the even labels`
`            direct-sum ker(affine_3^T) on the odd labels`.

Each half is the GF256 extension of the binary [8,4,4] extended Hamming
code. Both variants therefore have the same eight-dimensional quiet
subspace and quiet minimum weight four. Rank eight is maximal for a
sixteen-dimensional square-zero correction: its image must lie in its
kernel.

For both candidates, use independent uniform nonzero D1 and D2 entries,
sampled once at setup for each group and RS column. The outer column map is
`D2 H_C D1`. This is a new ensemble, not a change in how an old receipt is
named. Conditional on a complete output mask, D2 still supplies independent
nonzero amplitudes.

## Exact rank and quiet counts

The helper computes ranks over GF2, equivalently over its GF256 extension.
Its half decomposition is checked against the actual eight-column matrix.
The following entries count labelled active row sets, not affine orbits.

| Arity | Rank distribution for W |
|---|---|
|7|4:16, 5:448, 6:4704, 7:6272|
|8|4:2, 5:128, 6:1764, 7:7840, 8:3136|
|9|5:16, 6:448, 7:4704, 8:6272|

For comparison, H* has distributions `4:240,5:11200` at arity 7,
`4:30,5:12840` at arity 8, and `5:11440` at arity 9.
The two full a0 halves retain quiet dimension four; the stronger syndrome
does not remove those labelled classes.

Exact full-support counts follow by subset Möbius inversion in each
eight-point half, then multiplication. The full quiet weight enumerators
are saved in both local summary receipts. Some comparisons are:

|Exact column weight|H* quiet words|Either star variant|
|---|---:|---:|
|4|35,700|7,140|
|6|54,521,040|3,627,120|
|8|303,930,208,350|8,337,892,590|
|9|46,544,944,857,600|0|
|All weights|256^11|256^8|

Quiet counts never increase, since W includes the original five
syndromes. Weight nine here refers to one column's actual input support;
it does not exclude quiet columns in a message with nine active RS rows.

## Sparse-output check and its limitation

Every entry of both new matrices is nonzero. The check is especially
simple for the current labels: H* has neither a zero entry nor an entry
equal to 1. Adding the binary matrix f C f^T can therefore create neither
zero entries nor a new singleton hole.

The helper also exhausts all 120 source pairs and all 255 nonzero amplitude
ratios. Counts are per fixed labelled pair. That pair remains fixed
across RS columns; this calculation does not resample it.

The identity-star choice has maximum eight holes. Twenty-four projective
pair rays attain that maximum. For example, source pair (3,5), at ratio
`u_3/u_5`, has output mask 5437 and output weight eight. Involution then
gives a nonquiet arity-eight input whose output has weight two. This is
not a distance counterexample: the exact amplitude ray has to be counted
through the RS word. It is an adverse local feature that the proof cannot
replace by a quiet-kernel improvement.

The coupled choice restores maximum six holes. The complete projective
input-arity counts for outputs of weight two are

`{10:17, 12:71, 14:37, 15:1460, 16:29015}`.

These sum to 30,600. There are no such rays at arity eight or nine.
Together with the singleton check, this proves output weight at least
three for every nonzero input supported on at most nine coordinates.
It does not establish that weight three is attained. No triple-amplitude
enumeration was run.

Sixteen labelled pairs nevertheless have a worse maximum hole count than
the same pair under H*. Thus the old labelled pair laws do not transfer,
even though the global maximum remains six. Both summary receipts retain
all pair-law classes and every such worsening pair.

## Circuit cost

The shared circuit preserves the original pair and quarter sums. The
three quadratic syndromes need five additional XORs: two reuse odd-input
pair sums already needed for the a0 syndrome. Incorporating them into
the odd-half output deltas needs three more XORs.

The even and odd correction halves use 7 and 8 XORs respectively, the same
fifteen correction XORs as the old sixteen-point evaluation. The coupled
choice additionally forms `C*b=(b2,b1+b3,b2+b3)`, costing two XORs.

|Mixer|Fixed field products|Binary XORs|
|---|---:|---:|
|H*|8|66|
|Identity star|8|74|
|Coupled star|8|76|

These are exact scalar-circuit primitive counts, not measured native
instruction counts or a runtime prediction. Neither candidate changes
the diagonal keys, group size, RS transforms, output layout, or number of
passes. Register lifetimes and code generation remain native obligations.

## Validation, pins, and next gate

Eight focused tests pass in 1.745 seconds. They cover the factorization,
Gram-zero identity, involution, all matrix entries, 128 binary basis inputs
and 32 deterministic random inputs for each variant, rank decomposition,
an exhaustive binary quiet enumerator, GF256 rank/count checks, complete
pair regressions, and the exact primitive counts. An independent agent
reviewed the algebra, circuit, and count interfaces without rerunning a
provider. The two exact summaries were regenerated together in 1.76 seconds.

Files and SHA256 pins:

- `theory/hd_bit_syndrome_star_outer.py`:
  `bb7c767017eb8e9a62d7843c8eca2f27be0c313fdfe9f8adcb79007456d19345`
- `theory/hd_bit_syndrome_star_outer_tests.py`:
  `26f255eb500e1924504e6f39dad79ecb7d625f66951bf3639f22804709af519f`
- `theory/hd_bit_syndrome_star_identity_exact.json`:
  `da9681b67c7cd611123a4f4393d428485600e982ba288dd727d86a328e45a8ed`
- `theory/hd_bit_syndrome_star_coupled_exact.json`:
  `49e8df203482b84c492f58bf6883dc3b6e600a459d9094984439b34a0075d58c`

The receipts carry distinct construction identities ending in
`-identity` and `-coupled`, as well as the exact C matrix. Each records
`whole_code_certificate=false` and `old_receipts_transfer=false`.

Reproduce from the theory directory:

```
C:/Python314/python.exe -m unittest hd_bit_syndrome_star_outer_tests -v
C:/Python314/python.exe -c "import json,hd_bit_syndrome_star_outer as s; print(json.dumps(s.exact_summary('coupled'),indent=2))"
```

The next construction-specific proof gate should use the coupled variant
and retain its actual nonquiet fiber counts. The correction is affine on
each eight-point half, so its zero set is a pair of affine flats, including
empty flats. There are 52 choices per half, giving at most 2704 patterns.
This offers a finite counting interface; it is not yet an implemented
full kernel or a valid replacement for H*'s old 308-flat table. Actual cross
ranks, inside fibers, labelled chronology, and all RS constraints must
be recomputed. The current exact-empty-refresh improvement is independently
useful and should be tested on the unchanged encoder before a native pivot.
