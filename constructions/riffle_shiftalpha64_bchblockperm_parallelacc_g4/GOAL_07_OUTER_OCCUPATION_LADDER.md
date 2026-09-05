# Goal 07: outer occupation ladder

## Question

Can the fixed-profile kernel from Goal 06 absorb the number of outer words as
the number of occupied field symbols increases?

Let (n=16{,}386) and (q=2^{64}). The field code is an
([n,n-2,3]_q) MDS code. Define (A_s^{\mathrm{MDS}}) as the number of its
codewords with exactly (s) nonzero symbols. The exact formula is

\[
A_s^{\mathrm{MDS}}
=
\binom ns
\sum_{j=0}^{s-3}(-1)^j\binom sj
\left(q^{s-2-j}-1\right).
\tag{1}
\]

The first comparison charges (A_s^{\mathrm{MDS}}) against the compressed
inner bound for (s) equal BCH weights. This comparison is diagnostic. An
equal-weight profile does not cover the other BCH-weight profiles.

## Equal-weight diagnostics

At bad-output threshold (D=188{,}766), the dyadic ladder gives:

| (s) | \(\log_2 A_s^{\mathrm{MDS}}\) | inner exponent, weight 22 | deficit | inner exponent, weight 64 | deficit |
|---:|---:|---:|---:|---:|---:|
| 4 | 179.415 | 39.198 | 140.217 | 131.851 | 47.564 |
| 8 | 480.700 | 82.128 | 398.572 | 267.667 | 213.033 |
| 16 | 1075.742 | 188.935 | 886.807 | 540.065 | 535.677 |
| 32 | 2250.299 | 458.298 | 1792.001 | 1086.426 | 1163.873 |
| 64 | 4567.838 | 1073.490 | 3494.349 | 2186.029 | 2381.809 |
| 128 | 9139.143 | 2389.294 | 6749.849 | 4440.987 | 4698.156 |

A positive deficit means that the shell multiplicity exceeds the certified
inner exponent. The deficit increases along both tested profiles. The
weight-64 profile is much better than the minimum-weight substitution, but it
does not close these rungs.

These values do not refute the construction. They show that one fixed BCH
weight cannot serve as the outer-to-inner interface.

## Spectrum compression with the MDS equations

For a nonzero field value (z), let (h(z)) be the weight of its BCH
encoding. For a coefficient tilt (x>0), define

\[
a_x(z)
:=
\frac{x^{-h(z)}}{\binom{128}{h(z)}},
\qquad
a_x(0):=0.
\tag{2}
\]

Fix an outer support (S) of size (s). The two MDS equations parametrize
the codewords supported inside (S) by a vector

\[
\mathbf z\in\mathbb F^{s-2}.
\]

Each coordinate is a surjective linear form (L_i(\mathbf z)). Hölder's
inequality gives

\[
\sum_{\mathbf z\in\mathbb F^{s-2}}
\prod_{i=1}^{s}a_x(L_i(\mathbf z))
\le
q^{s-3}\sum_{z\in\mathbb F}a_x(z)^s.
\tag{3}
\]

The exact BCH spectrum evaluates the moment on the right side:

\[
\sum_{z\in\mathbb F}a_x(z)^s
=
\sum_{h>0} A_h
\left(\frac{x^{-h}}{\binom{128}{h}}\right)^s.
\tag{4}
\]

Equations (3) and (4) collapse every BCH-weight profile in one occupation
shell to one scalar moment. They retain both outer parity equations.

The unweighted version of (3) costs zero visible bits at (s=3) and (s=4)
at the displayed precision. Thus the obstruction is not the ordinary MDS
multiplicity count.

## Result of the scalar spectrum envelope

After optimizing the inner and BCH tilts, the complete shell bounds remain
trivial:

| (s) | independent-value envelope | MDS-Hölder envelope |
|---:|---:|---:|
| 3 | (2^{205.867}) | (2^{183.223}) |
| 4 | (2^{272.163}) | (2^{302.586}) |

The smaller envelope depends on (s); neither closes the first shell.

At (s=4), the optimized Hölder moment is dominated almost equally by BCH
weights 22 and 128:

| BCH weight | log contribution to the moment |
|---:|---:|
| 22 | -371.135 bits |
| 128 | -371.193 bits |
| 24 | -389.456 bits |
| 26 | -407.449 bits |

The unique all-one BCH word competes with the complete minimum shell. One
scalar coefficient tilt must cover both endpoints. It therefore discards the
central shape of the BCH spectrum.

Hölder also discards the alignment among the coordinate forms (L_i). For
example, it charges an (s)-fold endpoint product whenever one marginal
endpoint is large. The outer equations may forbid most or all such aligned
tuples for a given support.

## Consequence

Block occupation is too little state. Exact BCH weight is too much state. The
next interface should retain a coarse BCH-weight band.

A first banding should separate:

1. the low shell, beginning at weights 22 and 24;
2. the central spectrum around weight 64;
3. the high endpoint, including the unique weight-128 word.

Each band can use its own coefficient tilt. The outer calculation should
retain only the number of symbols in each band. This adds two independent
count coordinates, not one coordinate per BCH weight.

The high endpoint should also be checked directly against the first MDS
equation before it enters a moment bound. In characteristic two, an odd
number of identical all-one preimages cannot sum to zero. This simple check
already removes some products charged by the scalar Hölder envelope.

## Reproduction

Run the two optimized low-occupation shells sequentially:

```powershell
python scripts/analyze_riffle_shiftalpha64_outer_occupation.py --occupation 3 --scaled-cost 66.4 --zero-scale 23.7 --weight-tilt 1.507 --outer-envelope holder --optimize --optimizer-maxiter 100
python scripts/analyze_riffle_shiftalpha64_outer_occupation.py --occupation 4 --scaled-cost 88.6 --zero-scale 31.5 --weight-tilt 1.637 --outer-envelope holder --optimize --optimizer-maxiter 100
```

The analysis script has SHA-256
`DAE9D463BA02EF82CC9D5D3104072BB7FB05B8955B08AC3FF0289375A7FE3E70`.
The summary receipt is
`receipts/goal07_outer_occupation_ladder.json`.

## Next goal

Implement a three-band outer envelope. Start with occupation sizes 3 and 4,
where the outer solution spaces have dimensions one and two. Compare the band
envelope with exact low-shell multiplier intersections. Then climb through
8, 16, and 32 occupied symbols. If packet support becomes the computational
bottleneck, apply the support-integral compression proposed in Goal 06.

