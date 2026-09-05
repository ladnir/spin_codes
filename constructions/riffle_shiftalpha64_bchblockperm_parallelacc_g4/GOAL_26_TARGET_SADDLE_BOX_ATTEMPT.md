# Goal 26: determine the target saddle box

## Result

A rigorous target saddle box cannot yet be selected. The target-length outer
bulk sum has not been bounded, so there is no proved set of packet histograms
whose complement has negligible first-moment contribution.

Two numerical calibrations nevertheless change the proof plan.

1. The small-instance saddle range from Goal 25 does not contain the natural
   target centers. At output weight (D=76000), the center has fourth saddle
   coordinate approximately (-6.408), whereas the earlier calibration ended
   near (-3).
2. Covering every positive histogram is also the wrong replacement. A packet
   count can equal one, which drives its saddle coordinate toward minus
   infinity as the block length grows.

The proof should first separate rare packet coordinates from the interior.
The within-BCH permutation supplies a strong uniform lower-tail gate for the
number (h_4) of weight-four packets. This gate makes such a separation
plausible, but it does not by itself perform the outer-code sum.

## What a required box means

Fix a boundary input weight (H=2D). Let (mathcal R_D) be a set of packet
histograms. A box is required by the proof only after establishing


\[
\sum_{\mathbf h\notin\mathcal R_D}
C(\mathbf h)
\Pr[\operatorname{wt}(y)\le D\mid\mathbf h]
\le \varepsilon_{\rm out}.
\tag{1}
\]

The Fourier certificate then needs to cover the saddles of the histograms in
(mathcal R_D). Neither Goal 17 nor Goal 18 establishes (1) for the target
double-parity outer code. Goal 17 treats one sparse zero-parity shell. Goal 18
states the bulk interface but does not sum it.

Therefore a numerical range observed on selected histograms is a calibration
box, not the box required by the proof.

## Perron-limit calibration

Set (x_0=1) and (x_k=e^{\ell_k}). The probe optimizes the Perron limit
of the boundary matrix from Goal 21. Its center is the packet histogram of a
uniform binary word conditioned to have weight (H). This center is not the
proved outer histogram law; it is a scale diagnostic.

At (N=524352) packets, the centers are:

| (D) | center counts ((h_0,h_1,h_2,h_3,h_4)) | saddle ((\ell_1,\ell_2,\ell_3,\ell_4)) |
|---:|---|---|
| 76000 | ((388091.3,121290.7,14215.1,740.4,14.5)) | ((-1.602,-3.204,-4.806,-6.408)) |
| 77000 | ((386497.9,122508.0,14561.6,769.3,15.2)) | ((-1.595,-3.190,-4.784,-6.379)) |
| 188743 | ((237097.0,208151.2,68527.0,10026.8,550.2)) | ((-1.053,-2.106,-3.158,-4.211)) |

The optimizer matches every requested packet fraction to less than
(2\mathbin\cdot10^{-8}). These are floating-point Perron-limit results.
They omit the finite bridge endpoint correction.

Fixing each nonzero packet count to one, while maximum-entropy tilting the
other counts under the constraints on (N) and (H), produces coordinates
from approximately (-12.863) through (-0.664) over the three displayed
distances. This is not an extremal search. It suffices to show why a theorem
over all positive histograms cannot use the Goal 25 calibration box.

More generally, if (h_k) stays fixed while (N) grows, then the tilted mean
(h_k/N) tends to zero. The corresponding coordinate leaves every fixed
compact set. A compact-family local limit theorem therefore applies only
after a rare-coordinate split.

## A uniform BCH gate for (h_4)

Fix one permuted BCH block of weight (r), and let (J_r) count its
weight-four packets. Its expectation is

\[
\mathbb E[J_r]
=32\frac{\binom{124}{r-4}}{\binom{128}{r}}
=32\frac{(r)_4}{(128)_4}.
\tag{2}
\]

Every nonzero extended-BCH word has weight at least 22. The ratio
(\mathbb E[J_r]/r) is increasing for (r\ge22). Consequently, conditional
on any fixed outer word of total binary weight (H), independence of the
within-block permutations gives

\[
\mathbb E[h_4]
\ge
H\frac{32\mathbin\cdot21\mathbin\cdot20\mathbin\cdot19}
{128\mathbin\cdot127\mathbin\cdot126\mathbin\cdot125}
=0.0009973753281H.
\tag{3}
\]

At (D=76000), equation (3) gives (mathbb E[h_4]\ge151.601). The
uniform-word center above gives only 14.5 because it ignores the BCH block
minimum distance.

The complete one-block law also has an exact coefficient formula. For
(0<s<1), define

\[
G_r(s)
:=
\frac1{\binom{128}{r}}
\sum_{j=0}^{32}
\binom{32}{j}
[x^{r-4j}](1+4x+6x^2+4x^3)^{32-j}s^j.
\tag{4}
\]

Let

\[
c(s):=\max_{r:A_r>0,\ r>0}\frac{\log G_r(s)}r.
\]

For every fixed outer word of weight (H), Chernoff's inequality gives

\[
\Pr[h_4\le J]
\le \exp\bigl(Hc(s)-J\log s\bigr).
\tag{5}
\]

The exact one-block counts in (4) reproduce (\binom{128}{r}) for every BCH
weight. Numerical optimization of (s) gives the following diagnostics at
(H=152000):

| Threshold | log-base-two upper bound from (5) |
|---:|---:|
| (h_4\le0) | -220.02 |
| (h_4\le16) | -144.77 |
| (h_4\le32) | -101.55 |
| (h_4\le64) | -47.19 |
| (h_4\le96) | -17.11 |
| (h_4\le128) | -2.83 |

Weight 22 is the maximizing BCH weight at every displayed optimum. The
counts are exact, but most displayed optimized logarithms are not
outward-rounded proof certificates.

The (h_4 <= 64) row has a separate exact rational witness. Set (s=27/64),
and let (u) be a 64-bit dyadic upper envelope satisfying (G_r(s) <= u^r)
for every nonzero BCH weight. Exact integer arithmetic verifies every local
inequality and

\[
(27/64)^{-64}u^{152000}\le2^{-47}.
\tag{6}
\]

Thus the per-word (h_4 <= 64) gate at (D=76000) is a certified 47-bit
statement, not merely the numerical -47.19-bit optimum.

Equation (5) is a conditional packetization bound for one fixed outer word.
It cannot be used alone to discard (h_4\le64) from the first moment. The
number of outer words and their accumulator probabilities still occur in
(1). The gate instead identifies a useful coordinate split for the joint
outer-inner sum.

## Target face audit

The (H=16) proof gym exposed only three open face graphs after Goals 23 and
24. Across target-length boundary weights, every nonempty subset of packet
weights one through four can occur; a fixed value of (D) may remove some
faces by divisibility. There are fifteen faces containing packet weight zero.

The exact audit gives:

- six faces already covered by the exact formulas in Goals 23 and 24;
- three open faces audited in Goal 25;
- six additional target faces;
- all fifteen face graphs are primitive.

The six additional faces are

\[
\{0,1,2\},\ \{0,1,4\},\ \{0,2,3\},\ \{0,3,4\},
\ \{0,1,2,4\},\ \{0,2,3,4\}.
\]

Their presence does not recreate a histogram-by-histogram case ladder. A
single finite-state Fourier argument can be parameterized by the active face.
The audit shows that primitivity is available in every case. It also records
the exact lattice index at lengths six through twenty; the index is not
uniformly two on the lower-dimensional faces.

## Revised proof shape

The target proof should not certify one large four-dimensional box first.
It should use the following order.

1. Bound the target outer boundary sum in bands that retain at least the
   active face and (h_4).
2. Choose a threshold (J) from that joint contribution bound.
3. For (h_4<J), extract the (x_4) coefficient explicitly and apply a
   lower-dimensional local bound to the remaining variables.
4. For (h_4\ge J), certify the ordinary four-dimensional bridge on the
   compact saddle region selected by the outer sum.
5. Apply the same rare-coordinate rule recursively if another packet count
   approaches a face.

The fourth boundary component is especially simple: its nonzero entries lie
on (a+b=4) and all have multiplicity one. Tracking at most (J) uses a
block lift of the five-state transfer matrix rather than a table over all
histograms. Thus the split preserves the compression gained in Goals 20--25.

## Conclusion

The attempted target box did not close Goal 25. It found that Goal 25 was
premature and that its small-shell face list was incomplete at target length.
The obstruction is now narrower: construct a target outer boundary envelope
that retains one rare-coordinate count. Once that envelope determines (J),
both the compact Fourier box and the lower-dimensional remainder become
well-defined proof obligations.

## Reproduction

Run the probes sequentially:

```powershell
python scripts/probe_riffle_target_boundary_saddle_box.py
python scripts/probe_riffle_boundary_h4_gate.py
python scripts/audit_riffle_target_boundary_faces.py
```

They write:

- `receipts/goal26_target_boundary_saddle_calibration.json`;
- `receipts/goal26_boundary_h4_gate.json`;
- `receipts/goal26_target_boundary_faces.json`.

The next goal should derive a target-length outer boundary envelope indexed
by an (h_4) band. It should compose the exact BCH packetization polynomial
with an outer occupation bound before optimizing the accumulator saddle.
