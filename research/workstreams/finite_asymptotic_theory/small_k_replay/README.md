# Exact-spectrum replay at smaller message lengths

This folder compares fixed repeated constituents with authenticated weight
spectra. The construction uses one fixed constituent, repeated in every outer
row, followed by uniform routing and RandomStepConv. The experiment varies the
total message length and inner memory. Its comparison set is restricted to
exact-spectrum BCH, exact-spectrum RM, and a random rate-half outer model.

The constituents are:

- extended BCH ([32,16,8]);
- extended BCH ([128,64,22]); and
- RM(4,9) ([512,256,32]).

For total message length \(k=2^e\), a rate-half \([B,K]\) constituent is
repeated \(L=k/K\) times. The output length is \(N=BL=2k\). The distance
cutoff is

\[
D=\lceil0.109N\rceil.
\]

One constituent is fixed and reused.  Each row-coordinate permutation, each
region permutation, and each RandomStepConv map is sampled independently
once and fixed for the resulting code.

`evaluate_exact_spectra_q1_phase.py` uses every exact spectrum coefficient.
It does not replace BCH or RM by a random-code spectrum. Its numerical
results are binary64 diagnostics. Cases that clear 40 bits at occupation one
are candidates for an all-occupation replay, not distance certificates.

The best deterministic all-proof-regime candidate found so far is RM(4,9) at
\(k=2^{13}\) and memory \(M=22\). Its diagnostics give 40.717 bits for
occupation one and 84.277 bits for occupation two. The corresponding random
rate-half first-moment benchmark has 59.967 bits. At smaller \(k\), some
low-occupation RM bounds exceed 40 bits, but the random benchmark does not;
those cases are therefore poor candidates for the intended random-like
all-occupation proof.

The corrected banded Bernoulli run splits the all-zero candidate-input atom
before applying the Chernoff factor. It remains vacuous: live compositions
made entirely from the lowest RM weight band dominate. The earliest-pivot
reduction is also vacuous because reducing every active row to one pulse
leaves an inner-annihilation floor of order \(2^{-M}\). These are failures of
the attempted upper bounds, not evidence that the code lacks the target
distance.

The random comparison now also includes the requested constituent model: one
uniform full-rank \([B,B/2]\) map, sampled once and reused in every row. It is
used only at matched block sizes. At \(k=2^{13}\) and \(M=12\), neither the
random \([32,16]\) nor the random \([128,64]\) constituent clears Q1; their
margins are -10.529 and 6.180 bits. The corresponding exact BCH margins are
-10.134 and 7.181 bits. Random \([512,256]\) has Q1 and Q2 margins of 91.019
and 180.187 bits, but this is an apples-to-apples control only for RM(4,9),
whose matched Q1 margin is 17.048 bits. Occupations \(Q\ge3\) remain open for
the random 512-bit control because the rank and linear-relation type of the
active local messages must be retained. The older global random rate-half
calculation is kept only as an ideal reference.

See `SMALL_K_EXACT_SPECTRUM_STATUS.md` for the current boundary and the
remaining proof obligations.

`MATCHED_CONSTITUENT_CURVES.md` gives the controlled comparison requested for
10% distance. It fixes \(M(k)=\lceil\log_2 k\rceil+2\) for every curve and
compares each authenticated constituent only with a same-size reused random
constituent. The reported values are Q1 margins, not full certificate margins.

`RATE_HALF_FAMILY_CURVES.md` extends that comparison across every accepted
nontrivial rate-half BCH-derived and RM constituent below block length 1024,
and across random full-rank blocks of lengths 8, 16, 32, 64, 128, 256, 512,
and 1024. `SPECTRUM_SOURCES.md` authenticates the added tables and gives the
exact stopping rules. The BCH-derived ladder stops at length 128 because no
complete rate-half spectrum is available at the next sizes. The RM ladder
stops at length 512 because there is no rate-half RM code at length 1024.
The random ladder first clears the 40-bit Q1 screen at block length 256. No
structured curve clears it; RM(4,9) remains closest at 38.925 bits.

`RM2SUB_CALIBRATION_AND_FAMILY_REPLAY.md` repeats the experiment after
calibrating RM2Sub against RandomStepConv. The selected epoch length is 64.
The neutral schedule is (s(e)=\max(7,e-4)), which matches the earlier
RandomStepConv persistence exponent (M(e)=e+2) for (e\ge11). Under this
schedule, RM(4,9) has at least 40 binary64 Q1 bits for (14\le e\le18), with
a 41.462-bit maximum at (e=16). Extended BCH ([128,64,22]) reaches only
33.966 bits. These remain occupation-one diagnostics, not full distance
certificates. `audit_rm2sub_calibration_replay.py` exactly rechecks the
selected A/B maps and the receipt structure.

`TEN_PERCENT_PARAMETER_STUDY.md` defines the next optimization study. Its
primary target is 10% distance and 40 bits of total failure margin. Distances
near 10.9% or 11% are bonuses and do not drive parameter selection. The study
separates a fully certified frontier from a conditional frontier under one
explicit realized outer-spectrum cap.

The first tranche evaluates RM2Sub epochs 128 and 256 at `k=2^14` and
`k=2^16`. The exact RM(4,9) constituent is the only 40-bit Q1 survivor. Every
tested `t128` and `t256` state clears the corrected Q2 screen. The initial
`t128_s13` elimination was invalid because it combined aggregate results from
disjoint tilt grids instead of minimizing pointwise before aggregation. These
are binary64 occupation screens, not distance certificates.
`rm2sub_primary_frontier_summary_d100.json` is the corrected consolidated
receipt. `rm2sub_primary_tranche_audit.json` records a passing exact audit of
the selected inner maps and the summary structure.

`RM2SUB_PERSISTENCE_UNCONFOUNDED.md` corrects the initial interpretation of
the larger-epoch screen. A complete pointwise Q2 grid gives 85.947 bits for
`t128_s13`; the earlier negative value resulted from combining aggregate
outputs from disjoint witness grids. A literal nested `t=128` family shows
smooth behavior from `s=13` to `s=14`. At fixed `k=2^16` and persistence 20,
both Q1 and Q2 margins decrease from `t=64` to 128 to 256.

`RM2SUB_Q_LADDER.md` extends the matched-persistence comparison through
occupation eight and the strongest `t64_s14` sparse calculation through
occupation 32. A positive two-colour recurrence handles every mixture of
regular and all-one RM rows. Its direct-enumeration audit passes. The
`t64_s14` diagnostic remains positive through Q=29, with 30.683 bits there,
and becomes vacuous at Q=30. This locates the next exploratory task without
imposing a 40-bit acceptance gate.

The first bridge probes reduce that gap. Jointly optimized dense spectrum
bands close every pure face of weight at least 64 at Q=29. The minimum-weight
dense face becomes positive at Q=53, while its exact sparse counterpart stays
positive through Q=54.

`RM2SUB_TWO_BAND_BRIDGE.md` now covers every composition supported on the
minimum-weight band and one other spectrum band for Q=30 through 52. The
union of the displayed mixed families has 1172.882 diagnostic bits.

`RM2SUB_REFINED_BAND_BRIDGE.md` completes the diagnostic multiband bridge.
The refined groups are weights 1--95, 96--416, and 417--512. Their complete
union over Q=30 through 256 checks 2,857,249 compositions and has 1228.638
bits. The consolidation audit passes. Combining the strongest receipts over
Q=1 through 256 gives 42.577982 diagnostic bits, with Q=1 weakest. Directed
outward rounding remains open.
