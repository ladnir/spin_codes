# Finite (k=2^{20}) independent-row setup

## Frozen theorem interface

This note fixes the probability space for the first finite certificate. It
does not claim that the complete 40-bit distance proof is closed.

Set

\[
  B=240,\qquad L=8832,\qquad N_+=BL=2{,}119{,}680,
  \qquad K_+=LB/2=1{,}059{,}840.
\]

The message space is the subspace of \(\mathbb F_2^{K_+}\) in which the final
11,264 coordinates are zero. The first \(2^{20}\) coordinates carry the
message. This convention matches the implemented transposed interface.

For each row \(i\in[L]\), setup independently performs the following steps.

1. Apply ten copies of the extended binary Golay \([24,12,8]\) encoder.
2. Sample \(\rho_{i,1}\gets S_{240}\), apply it, and apply the terminated
   rate-one accumulator.
3. Sample \(\rho_{i,2}\gets S_{240}\), apply it, and apply a second terminated
   rate-one accumulator.
4. Repeat steps 2 and 3 with fresh independent permutations until the
   resulting Golay--BA-3 code has no nonzero word of weight outside
   \([23,217]\). Denote this event by \(\mathcal G_{240}^{23}\).
5. Sample an independent route permutation \(\sigma_i\gets S_{240}\).

The rejection step defines the exact conditional law used by the proof. The
mathematical sampler has no failure return and terminates almost surely. The
proved lower bound

\[
  \Pr[\mathcal G_{240}^{23}]\ge
  \mathtt{0x1.fcc11f1547a88p-1}
\]

gives at most 1.00638 expected BA trials per row. This is not yet an efficient
setup algorithm: no efficient decision procedure for \(\mathcal G_{240}^{23}\) is
known. A capped sampler, a seeded pseudorandom sampler, or an unconditioned
sampler has a different probability space and is not covered by this note.

After the row maps, arrange their outputs as an \(L\)-by-\(B\) matrix, apply
every \(\sigma_i\), and transpose. For each region \(j\in[B]\), sample an
independent \(\pi_j\gets S_L\). All BA, route, and region permutations are
mutually independent except for the explicit per-row conditioning on
\(\mathcal G_{240}^{23}\).

Serialize the 240 permuted regions into 16,560 consecutive epochs
\(X_r\in\mathbb F_2^{128}\). Let

\[
  A:\mathbb F_{2^{19}}\to\mathbb F_2^{128},\qquad
  C:\mathbb F_2^{128}\to\mathbb F_{2^{19}}
\]

be the fixed audited RM2Sub-S19 constituents. Independently sample
\(\alpha_r\gets\mathbb F_{2^{19}}^*\) for every epoch, set \(Q_0=0\), and
compute

\[
  Y_r=X_r+A(Q_r),\qquad
  Q_{r+1}=\alpha_rQ_r+C(X_r).
\]

All 16,560 words \(Y_r\) form the output. There is no terminal-state
constraint and no terminal-state symbol. The implementation evaluates the
transpose in reverse storage order; binding that circuit to this recurrence
remains obligation O6.

## Finite target

For the resulting encoder \(E\), define

\[
  Z_d(E)=\left|\{x\ne0:\operatorname{wt}(E(x))\le d\}\right|,
  \qquad d=233{,}164.
\]

The target is

\[
  \Pr[Z_d(E)>0]\le2^{-40}.
\]

The probability is over the conditional BA draws, route permutations, region
permutations, and nonzero field multipliers just specified. Shortening only
restricts the parent message space, so a first-moment bound for the full
\(K_+\)-dimensional parent is conservative for this target.

## Established finite component

`certify_golay_ba_rm2sub_finite_one_active.py` proves, with one-sided
binary64 arithmetic, that messages occupying exactly one outer row contribute
at most

\[
  \mathtt{0x1.6f3c66666f370p-47}<2^{-46}
\]

to \(\mathbb E[Z_d]\). The receipt covers every permitted outer weight from
23 through 217 and includes the factor of 8,832 possible active rows. The
revised receipt has suffix `w23_217`.

`certify_golay_ba_rm2sub_finite_q2_64.py` proves that occupations 2 through
64 together contribute at most

\[
  \mathtt{0x1.4dc4b8b530941p-1}\cdot2^{-84}<2^{-84}.
\]

The union of occupations 1 through 64 therefore remains below \(2^{-46}\).
These are proved components, not a complete distance certificate.

## Performance evidence and limits

The implementation task measured the independent-row transposed online map at
10.407942 ms on one pinned Zen 4 core. The B=256 Structured SPIN baseline in
the same binary took 10.696591 ms, so the independent-row map was 2.699%
faster. The independent setup occupied 18,547,208 bytes and the workspace
occupied 38,154,240 bytes.

The external performance receipt has SHA-256
`4b612e86c267a914b6cb5f6b231579f2ceb8dcfd7aa5981ee0b5807ce4fa01ee`.
It measures neither ordinary encoding nor conditional setup. It uses seeded
pseudorandom schedules rather than the exact uniform conditional law above.
Consequently it supports the online performance comparison, but not the
finite probability claim or a complete linear-time setup claim.

## Open obligations

- Give an efficient exact or authenticated implementation of the
  \(\mathcal G_{240}^{23}\) conditioning step.
- Outward-certify occupations 65 through 8,832; this dense range is the main
  mathematical gap.
- Complete the RM2Sub implementation-equivalence and source-receipt audit.
- Measure ordinary encoding and conditional setup if the final theorem claims
  concrete costs for them.
