# BCH [256,128] with RM2Sub: occupation-one integration

**Current result:** `T64_S20_FULL_CLOSURE.md` closes the full distance
target for the selected t64_s20 map: K=2^20, rate1/2, relative distance
greater than0.1, with setup-failure bound below2^-50 under fresh independent
nonzero field multipliers. All8192 occupancies passed certification and
replay; `generated/larger_coverage_full.json` checks the exact union.
The t128_s15 discussion below is retained as the earlier, separate case.

**Earlier t128_s15 status:** `ZERO_STATE_FIRST_MOMENT_OBSTRUCTION.md`
proves that the actual bad-message expectation at Q2620 exceeds 2^19600.
Thus the unconditional first-moment route cannot close t128_s15. This does
not by itself refute the requested setup-failure probability. Verified
partial coverage through Q1655 remains valid; the full goal remains open.
`PAIR_FOURIER_AND_PARITY.md` records the latest shared-type certificates,
global pair-parity mixing, and the exact-weight-80 second-moment family.
These are dependency bounds, not a full setup-failure probability result.
`WEIGHTED_PARITY_AND_JENSEN.md` adds a polynomial weighted-parity theorem
and a sharper actual first-moment denominator for the exact-weight80 family.
`PAIR_ALIAS_TABLE_AND_CENTRAL_WITNESS.md` certifies73600 shared-region types
and a central-count witness that restricts only the counted bad messages,
not the encoder. Its full second moment remains open.

This directory joins the existing deterministic BCH spectrum bounds to an
RM2Sub inner. It changes neither the completed random-matrix proof nor the
other agent's worktree. Inputs from that worktree are local snapshots.

The first milestone covers messages with exactly one nonzero BCH row at
message length `2^20`. The follow-up in `OCCUPATION_TWO.md` also covers two
nonzero rows for t128_s15, with an upper bound below 2^-91. `OCCUPATION_THREE.md`
adds a bound below 2^-126 for three rows. `GENERAL_OCCUPANCIES.md` gives a
kernel-aware formula for all occupancies and certifies every occupancy from
4 through 16. `ADAPTIVE_OCCUPANCY_RANGE.md` removes composition enumeration
and certifies every occupancy from 17 through 64.
`OA29_SHELL_CAPS_AND_TIGHTENING.md` extends verified coverage through 128
and derives stronger deterministic shell caps. `SHARED_RANGE_1024.md` uses
those caps and separate group probabilities to certify every occupancy
through 1024. `CONSTANT_ROW_SPLIT.md` extends complete coverage through
Q1655 using stronger exact shell caps and shared witnesses. It also certifies a broad class
with many all-one rows and the all-ordinary composition at Q8192. The
certified union has contribution below 2^-49; full coverage above Q1655
remains open. The derivation below concerns Q1.

## Construction and probability space

Use the fixed code C from
`../../bch_spectrum_work/bch_spectrum_codex_bundle/RANDOM_INNER_M22_CLOSURE.md`.
Its parameters are [256,128,d >= 38]; its spectrum is not assumed known.
The code is the previously audited five-dimensional syndrome restriction of
the extended designed-distance-37 BCH code over GF(256), with modulus 0x14D.
The retained outer certificate supplies the spectrum constraints used below.

Split the message into L = 8192 rows of 128 bits. Encode each row with C.
Independently permute the 256 coordinates of each encoded row. Transpose
into 256 regions of L bits, and independently permute each region. Serialize
in region order, obtaining N = 2^21 bits.

Fix one of the three selected map pairs in `inputs/`. In column-vector
notation, A maps F_2^s into F_2^t and B = A^T maps F_2^t into F_2^s.
The audited maps have full rank, BA = 0, and distinct nonzero B columns.
Identify F_2^s with a fixed basis of F_{2^s}. For each successive t-bit epoch,
emit Y and update the state q by

\[
q_0:=0,\qquad Y_i:=X_i+Aq_i,\qquad
q_{i+1}:=\alpha_iq_i+BX_i.
\]

Each multiplier alpha_i is sampled independently from F_{2^s} minus zero.
The state continues across region boundaries; there is no final flush.
The setup consists of all permutations and multipliers, sampled once and
shared by every message. Probability statements below concern this setup.
They do not replace the multipliers with a seeded generator or reuse them
across epochs.

Let Z_1 count messages with exactly one nonzero row whose output weight is
at most H = 209716. This cutoff is retained from the random-matrix comparison;
it is one larger than floor(0.1N). We bound Pr[Z_1 > 0] by E[Z_1].

## Why retain an activation state?

When q = 0 and X = e_J, the update gives q' = B(e_J). The fresh multiplier
does not randomize that state. The selected maps have only t possible
activation states, whereas there are 2^s - 1 nonzero field elements.

The imported finite two-state calculation labels this transition live, then
uses a near-uniform live-state average. Its stated distribution invariant
does not cover activation. A separate domination argument might still justify
its numerical bound; this integration does not assume one.

The other worktree's dense-occupation analysis already distinguishes a
deterministic state. We make that distinction explicit in this finite Q1
calculation too. Historical `*_screen.json` files are unverified two-state
diagnostics. Their certification entry point is disabled. Only filenames
containing `_activation_` are current results.

## An activation-aware epoch bound

Fix z in (0,1), and write M = 2^s - 1 and kappa = M/(M-1). Let a_w be the
number of nonzero words of weight w in im(A), and d_A their minimum weight.
Define

\[
m_0(z):=\frac1M\sum_w a_wz^w,\qquad
m_1(z):=\frac1{tM}\sum_w a_w
  \bigl(wz^{w-1}+(t-w)z^{w+1}\bigr).
\]

The second expression averages the effect of flipping one uniform coordinate.
Use three classes of normalized state laws:

- Z: the point mass at zero;
- D: any distribution on nonzero states;
- L: a distribution on nonzero states with point probabilities at most kappa/M.

These classes describe weighted measures, not an exact three-state Markov
chain. At every step, decompose the measure weighted by z to the accumulated
output weight into nonnegative masses times laws in these classes. A transfer
entry bounds the mass assigned to its destination class.

For a zero input and for X = e_J with J sampled from [t], respectively, use
the row-to-column transfers

\[
T_0=\begin{pmatrix}
1&0&0\\
0&0&z^{d_A}\\
0&0&\kappa m_0
\end{pmatrix},\qquad
T_1=\begin{pmatrix}
0&z&0\\
z^{d_A-1}/M&0&z^{d_A-1}\\
\kappa m_1/M&0&\kappa m_1
\end{pmatrix}.
\]

To justify the transfers, first consider a nonzero entering state q.
Conditioned on q and X, alpha q is uniform nonzero and independent of the
emitted weight. If X = 0, the next state is uniform nonzero. If X = e_J,
its probability of becoming zero is 1/M. Conditioned on survival, its law is
uniform outside zero and B(e_J), hence belongs to L. Mixtures over q and J
preserve the L bound, including after output weighting.

A D state emits at least d_A bits on zero input, and at least d_A - 1 bits
on weight-one input. An L law bounds the corresponding moments by kappa m_0
and kappa m_1. For T_1, we bound the surviving mass by the full moment;
retaining the additional termination entry only enlarges the bound.
The Z row follows directly from initialization and nonzero B columns.
These facts establish the weighted-measure invariant by induction.

## Region averaging and the BCH interface

A region has E = L/t epochs. For one nonzero BCH row, a selected region has
exactly one input bit, uniform among its L positions. Other regions are zero.
Thus its marked epoch and within-epoch coordinate are independent and uniform.
The entering state is independent of this region's permutation. Define

\[
R_0:=T_0^E,\qquad
R_1:=\frac1E\sum_{j=0}^{E-1}T_0^jT_1T_0^{E-1-j}.
\]

Let e_Z = (1,0,0), and let 1 denote the three-entry column of ones. For an
outer word of weight w, its region support is uniform among w-subsets. Hence

\[
p_w(z):=\min\left\{1,\frac{z^{-H}}{\binom{256}{w}}
[x^w]e_Z(R_0+xR_1)^{256}\mathbf1\right\}
\]

bounds its bad-output probability. Every coefficient is nonnegative, so
interval arithmetic can evaluate the expression without cancellation.
Choose a separate tilt for each w and set c_w := L p_w(z_w).
If A_w(C) is the outer spectrum, then

\[
\mathbb E[Z_1]\le\sum_w A_w(C)c_w.
\]

This is where the existing BCH work enters. Complement symmetry gives
A_w(C) = A_{256-w}(C). The existing exact LP certificate bounds a positive
weighted sum of the three paired shells w = 38,40,42. Let its positive paired
coefficients be b_w and its upper bound be J. Set

\[
\rho:=\max_{w\in\{38,40,42\}}
\frac{c_w+c_{256-w}}{b_w}.
\]

The six-shell contribution is at most rho J. All other nonzero shells use
the existing deterministic caps U_w. Their sum with rho J gives the bound
U_1 on E[Z_1]. Only the old *spectrum objective* is reused: no random-matrix
transfer probability or higher-occupation bound is substituted for RM2Sub.

## Checked results

The implementation searches z = exp(-exp(j/10)) for integer j from -120 to 0.
Binary64 selects tilts; 256-bit Arb recomputes their positive transfers.
Exact rational upper endpoints feed the BCH aggregation. An independent
512-bit Arb implementation checks all 92 coefficients for each map pair.

| t | s | d_A | Diagnostic -log2(U_1) | Exact checked inequality |
|---:|---:|---:|---:|---:|
| 64 | 16 | 16 | 49.8589352990 | U_1 < 2^-49 |
| 128 | 15 | 48 | 49.4813480343 | U_1 < 2^-49 |
| 256 | 14 | 112 | 48.8336361505 | U_1 < 2^-48 |

These are upper bounds on failure for occupation one, not estimates of the
actual failure probability. All three leave more than eight bits below the
2^-40 target for this contribution. For t128_s15, the follow-up also checks
the sum of U_Q for Q=1 through 64 is below 2^-49. The full distance theorem
still requires a bound for occupations 65 through 8192 for that configuration.

## Reproduction and merge boundaries

All new code and notes belong to this directory. No production encoder,
shared script, existing BCH certificate, or file in worktree `3061` was edited.
The map manifest records source and snapshot hashes. Runtime imports use our
local BCH bundle only; the other worktree is unnecessary after snapshotting.

Dependencies are Python, NumPy, and python-flint. From this worktree root:

```powershell
python -B workstreams/bch_rm2sub_bridge/test_activation_bridge.py
python -B workstreams/bch_rm2sub_bridge/verify_bridge.py
```

Also replay the retained outer certificate from its bundle directory:

```powershell
python -B code/verify_bch_m22_closure.py
```

That verifier additionally checks the old random-matrix theorem. Its success
does not transfer the old higher-occupation bounds to this construction.

For a fresh local snapshot, run `bridge.py snapshot --source PATH`, where
PATH is the other worktree root. Then run `activation_bridge.py screen` and
`activation_bridge.py certify`, each with `--configuration t64_s16`, or either
other configuration. Outputs are write-once; existing outputs cause failure.

The input and generated directories are ignored by Git. The existing BCH
repository snapshot is also notes-only. Therefore these small scripts and
notes are not a self-contained computational certificate without the retained
local BCH evidence and selected map inputs. Do not commit the experimental
directories merely to make replay paths resolve.

Next: extend the adaptive bound above Q=64 and share witnesses across longer
occupancy intervals; see `ADAPTIVE_OCCUPANCY_RANGE.md`. Retain t64_s16 as the
strongest Q1 bound among the three configurations. No performance comparison
has been run here.
