# Proof status

The candidate is promising at relative distance \(0.09\), but it does not yet
have a certificate. The exact calculations cover both proposed extremal
packet-rank profiles. The remaining issue is to prove that these extremes
bound every intermediate profile.

## Exact inner reduction

Two state classes suffice in the ideal field-checksum model. Class \(Z\)
contains the deterministic zero state. Class \(U\) contains a fresh uniform
word in \(\mathbb F_2^{64}\), including zero.

For a fixed nonzero epoch output, the four independent field coefficients
make the next state uniform in \(\mathbb F_2^{64}\). Thus every nonzero output
enters class \(U\).

For a fresh incoming state and fixed epoch input, write
\(B_\ell=A(X_\ell)\). The epoch output is

\[
(B_0+v,B_1+v,B_2+v,B_3+v),
\]

where \(v\) is uniform in \(\mathbb F_2^{64}\). For \(0<z<1\), its weight
moment is at most

\[
\left(\frac{1+z^4}{2}\right)^{64}.
\]

The exact diagnostics retain the outer candidate bits instead of applying
this input-independent bound. A 16-state dynamic program follows the four
accumulator bits at each packet position. It records the output-weight moment
and whether the complete epoch output is zero. Exhaustive enumeration for
three packet positions agrees with this recurrence to
\(1.35\mathbin\cdot10^{-15}\).

## Current numerical evidence

The calculation uses the modeled regular spectrum of a
\([256,128,38]\) outer code, \(2^{20}\) message bits, and target distance
\(0.09\cdot2^{21}\). The following rows are floating-point diagnostics, not
certificates.

| Active outer blocks | Packet-rank geometry | Margin (bits) |
|---:|:---|---:|
| 1 | one rank-1 group | 55.96 |
| 2 | one rank-2 group | 126.02 |
| 3 | one rank-3 group | 189.85 |
| 4 | one rank-4 group | 209.41 |
| 4 | four rank-1 groups in one lane | 210.05 |
| 8 | two rank-4 groups | 512.70 |
| 8 | eight rank-1 groups in one lane | 514.03 |
| 12 | three rank-4 groups | 716.40 |
| 12 | four rank-3 groups | 717.10 |
| 12 | six rank-2 groups | 717.84 |
| 12 | twelve rank-1 groups in one lane | 718.59 |

The exact packed calculation covers every profile \(4q+r\), not merely the
rows in this table. Its pointwise and summed regular-word margins are both
55.96 bits because the one-active-block row dominates all other rows.

The second exact extremum fills the 2,048 packet groups one lane at a time.
It therefore covers maximally split profiles from zero through 8,192 active
blocks. Its minimum margin is also 55.96 bits. Among uniform-rank profiles
with the same total rank, fuller packets are consistently worse after tilt
optimization.

The original support-only relaxation fails by about 859,290 bits because it
replaces all candidate bits in a nonempty zero-state epoch by one output bit.
That failure diagnoses the relaxation, not the construction.

## Failed one-extreme lemmas

Packing is not entrywise monotone at the epoch level. An exhaustive audit of
2,092 toy profiles with four packet positions refutes the claim that one
packed group always dominates every split profile of the same total rank.
For example, four rank-1 packets in one lane have a larger nonzero-output
moment than one rank-4 packet at \(z=0.8\).

This local counterexample is compatible with the full-region evidence.
Splitting a group creates more occupied epochs after the packet permutation.
Each additional nonzero epoch applies another independent 64-bit checksum.
Over a complete 2,048-packet region, this refresh benefit slightly exceeds
the same-lane cancellation loss. The packed rank-4 profile is then worse.

Packing alone also fails as a pointwise scalar bound at every tilt. In a toy
four-epoch model, four separated rank-1 groups exceed the packed scalar
moment by 1.19 bits at \(z=0.2\). Thus the proof cannot discard the split
extremum before optimizing the Chernoff tilt.

The corrected diagnostic takes the larger of the packed and maximally split
moments at each tilt and minimizes only afterward. For active-block counts
1 through 12, this calculation still gives 55.96 summed bits. The packed
profile is active at every optimizing tilt in this sparse range.

## Proof target

For a rank profile \(\rho\), let \(R_\rho(z)\) be its exact two-state region
matrix after averaging the random placement of its occupied packet groups.
Let \(P_a\) and \(S_a\) be the maximally packed and maximally split profiles
of total rank \(a\). The next lemma should prove

\[
e_Z^\mathsf{T}R_\rho(z)^{256}\mathbf 1
\;\leq\;
\max\left\{
e_Z^\mathsf{T}R_{P_a}(z)^{256}\mathbf 1,
e_Z^\mathsf{T}R_{S_a}(z)^{256}\mathbf 1
\right\}
\]

for every \(0<z\leq1\) and every profile \(\rho\) of total rank \(a\).
This is a scalar statement after a complete region and all 256 repetitions;
it does not assert the false entrywise matrix order.

An exhaustive toy audit checked all 2,092 mask/rank profiles of total rank at
most four. The two-extreme inequality held at all four tested tilts, with a
maximum apparent violation of \(1.17\mathbin\cdot10^{-12}\) bits at \(z=1\),
consistent with floating-point noise. A useful proof may expose the number
of occupied epochs as an intermediate variable and condition on it.

After that comparison, the remaining certificate work is to handle the
exceptional all-one word, replace the conjectured outer spectrum if needed,
and add outward-rounded arithmetic.
