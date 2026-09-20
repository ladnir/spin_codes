# Riffle DP-2Lap g=4: global block-moment route

## Purpose

The 24-node minimum-distance route leaves 57 component supports and scales
poorly. This checkpoint replaces a minimum-distance condition on every local
word with one cumulative condition over the complete 32,772-node interval.
The construction is unchanged and remains deterministic.

## Exact block-moment lemma

Let (N=32{,}772), (K=N/12=2{,}731), and (M=16N=192K). Fix a
nonzero character (chi\in\mathbb F_2^{64}) and a nonzero packet value
(v\in\mathbb F_2^4). At node (t), define the signed slot sum

\[
S_t(\chi,v)
=\sum_{s=0}^{15}
(-1)^{\langle\chi,T^{t+1}(v\mathbin{\ll}4s)\rangle}.
\]

For block (j), define

\[
R_j(\chi,v)=\sum_{t=12j}^{12j+11}S_t(\chi,v),
\qquad 0\le j<K.
\]

The normalized full-interval character sum is

\[
b_v(\chi)=\frac1{192K}\sum_{j=0}^{K-1}R_j(\chi,v).
\]

The following two bounds would prove the coarse character caps required by
the support-33 terminal-zero gate:

\[
\sum_{j=0}^{K-1}R_j(\chi,v)^2\le K\cdot96^2
\quad(v\ne15),
\tag{BM-1}
\]

and

\[
\sum_{j=0}^{K-1}R_j(\chi,15)^2\le K\cdot120^2.
\tag{BM-15}
\]

Indeed, Cauchy--Schwarz gives

\[
|b_v(\chi)|
\le
\frac{\sqrt{K\sum_j R_j(\chi,v)^2}}{192K}.
\]

Thus (BM-1) gives (|b_v(\chi)|\le1/2), and (BM-15) gives
(|b_{15}(\chi)|\le5/8). There is no endpoint term because 12 divides
32,772 exactly.

## Equivalent pair-difference character bound

Each block contains 192 packet cells. Let \(\mathcal D_v\) be the multiset
of contribution differences

\[
z_v(p)\mathbin\oplus z_v(q)
\]

over every unordered pair of distinct cells (p,q) in the same 12-node
block, and over all (K) blocks. Its size is

\[
|\mathcal D_v|=K\binom{192}{2}=18{,}336K.
\]

Expanding the square gives the exact identity

\[
\sum_j R_j(\chi,v)^2
=192K+2\sum_{d\in\mathcal D_v}(-1)^{\langle\chi,d\rangle}.
\]

Consequently, (BM-1) is equivalent to the one-sided Fourier cap

\[
\frac1{|\mathcal D_v|}
\sum_{d\in\mathcal D_v}(-1)^{\langle\chi,d\rangle}
\le\frac{47}{191},
\]

and (BM-15) is equivalent to the cap (74/191). This formulation exposes
the global permutation structure: \(\mathcal D_v\) is a union of 18,336
orbits of length 2,731 under the 12-step transpose recurrence.

## Refuted node-moment precursor

A first attempt applied Cauchy--Schwarz to the individual node sums. It
proposed

\[
\frac1N\sum_t S_t(\chi,v)^2\le64
\quad(v\ne15)
\]

and the analogous cap 100 for value 15. Exact Walsh enumeration proves these
bounds on every pure irreducible component. The largest pure-component values
are 36 for non-15 values and 100 for value 15; both occur in the degree-one
component. The degree-18 and degree-20 maxima are below 16.51.

The mixed-component claim is false. For

\[
v=7,
\qquad
\chi=\mathtt{0xe6f85de5120c7879},
\]

the character has support (C_1\oplus C_2), and

\[
S_t=(-10,10,2)
\]

with period three. Its node second moment is exactly 68. Its normalized full
character bias is only (1/24), and every 12-node block sum equals 8. The
counterexample shows that node energy discards essential cancellation at the
small component periods. It does not refute the desired character cap or the
block-moment lemma.

## Block-moment evidence

The block-moment probe evaluates every reported character exactly. Its search
set contains:

- the prior full-character-bias witnesses;
- the exact maximizer from each pure irreducible component;
- all XOR combinations of those seven component maximizers, followed by
  coordinate ascent from the best combination;
- the exact (C_{18}\oplus C_{20}) 12-node counterexample;
- the dimension-61 and full-dimension witnesses from the 10,000-restart
  24-node distance search; and
- pseudorandom coordinate-ascent restarts.

No block-moment violation was found. The largest reported non-15 block mean
square is 5,184, attained by the degree-one character for values 1 and 9.
The proposed cap is 9,216. For value 15, the degree-one character attains
14,400 exactly, so (BM-15) is tight.

The two known obstructions have small cumulative block energy. The
(C_1\oplus C_2) node-moment counterexample has block mean square 64. The
(C_{18}\oplus C_{20}) local-distance counterexample has block mean square
approximately 197.666789, even though one block has magnitude 100.

These values are refutation evidence only. Coordinate ascent can miss narrow
linear subspaces, as the low-rank short-window checks demonstrated. The probe
does not prove (BM-1) or (BM-15).

## Structural split suggested by the evidence

Let

\[
L=C_1\oplus C_2.
\]

The component root orders are one and three, so the 12-step recurrence acts
as the identity on (L). All eight characters in (L) can therefore be
handled by exact enumeration. The remaining problem has a nonzero component
in

\[
C_4\oplus C_9\oplus C_{10}\oplus C_{18}\oplus C_{20}.
\]

This split isolates the precise mechanism that refuted the node-moment
lemma. It also prevents a constrained refutation search from falling back to
the easy low-period extremizers.

Marginal component bounds alone do not prove a mixed-component Fourier cap.
A joint distribution can have small marginal coefficients and a large mixed
coefficient. Any proof must therefore control the pair-difference distribution
jointly or establish a conditional bound after fixing the three low-mode
bits.

## Next bounded checkpoint

The next proof task is to compute exact Walsh spectra of the 12-node
pair-difference distribution on (L) and on each pure high component. In
parallel, the refutation task is to search under the explicit constraint that
the high-component projection is nonzero. The checkpoint succeeds if it
either finds a block-moment counterexample or produces exact component
margins and a credible conditional inequality for mixed components.

The route should be abandoned if high-component constrained search approaches
the cap or if the exact pair-distribution calculation recreates the same
support-by-support explosion. Otherwise, (BM-1) and (BM-15) are the current
best global proof targets.

## Artifacts

- Exact pure-component node moments:
  `explorations/riffle_dp_g4_component_second_moment_certificate.json`
- Node-moment counterexample and replay:
  `explorations/riffle_dp_g4_global_second_moment_probe.json`
- Block-moment diagnostic:
  `explorations/riffle_dp_g4_global_block_moment_probe.json`
- Sources:
  `scripts/certify_riffle_dp_g4_component_second_moment.py`,
  `scripts/probe_riffle_dp_g4_global_second_moment.py`, and
  `scripts/probe_riffle_dp_g4_global_block_moment.py`

