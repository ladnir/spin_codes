# Complete finite parameter grid

The finite grid includes Q1, adaptive Q2..64, composition-preserving
Q2/Q3/Q4, and a typed cover of Q65..L on every native tuple. The producers
retain signed bounds, counting fallbacks, and authenticated witnesses.
Exact-region composition boxes provide stronger selected RM results.
See [`GRID_FINDINGS.md`](GRID_FINDINGS.md) for the final coverage audit,
parameter comparisons, and remaining numerical slack.

The goal is to evaluate the RM2Sub parameter landscape for complete
exact outer spectra and a modest random-outer reference range. BCH-256 is
excluded from this grid. Its separately completed bounded-spectrum result
remains a comparison point; it does not determine the fitted BCH curve.

## Finite scope

| Outer family | Constituent lengths | Evidence used |
|---|---|---|
| BCH | 8, 32, 64, 128 | Complete exact spectra |
| Reed--Muller | 8, 32, 128, 512 | Complete exact spectra, RM(1,3) through RM(4,9) |
| Random | 8, 16, 32, 64, 128, 256, 512, 1024 | Uniform full-rank rate-half subspace, sampled once and reused |

Every constituent has dimension half its length. RM(5,11) has only partial
low-weight information in `../rm511_outer/`, so it is outside this grid.
Only four complete BCH spectra are available in the curated inventory.

For each constituent, use message exponents 16, 18, 20, 22, 24 and epoch
sizes t=64,128,256. For each t, retain every state dimension
log2(t)+1 <= s <= 20. These are all prefixes of the same frozen s=20
generator chain used by the pilot. The smallest prefix is the complete
affine RM subcode. This gives 39 map configurations and 3,120 candidate
tuples. Exactly 3,108 tuples are native: length-1024 random outers at
message exponent 16 and t=256 have only 128 outer rows, less than an epoch,
so those twelve tuples are recorded as inapplicable.

The cutoff is floor(N/10), and every signed margin is retained. A 40-bit
margin is a comparison threshold, not a filter that deletes failed screens.
One outer constituent is reused across rows. Row-coordinate permutations
and region permutations are independent; state carries across regions.
Output precedes the fresh nonzero-multiplier state update, without a flush.

## Q1 production and evidence

The grid reuses 1,012 authenticated preferred Q1 observations and fills
2,096 missing tuples. Each map is a restart boundary. Completed batches
have a CSV and a final manifest containing source and CSV hashes. Pending
files never count as results. Existing map identities are preserved so
that occupation comparisons join exactly.

The first new batch uses the retained NumPy coefficient implementation.
`run_complete_q1_grid_native.py` resumes with an equivalent C++ recurrence.
The native implementation uses fixed three-state contractions, two reusable
buffers and log arithmetic. It does not discard small positive terms by
converting whole region matrices to ordinary floating-point probabilities.
Independent dense-matrix and extreme-log tests compare it to NumPy.
The source, build script and binary identity are recorded in the receipts.

The witness grid is log-surprisal -16,-15.9,...,0. Boundary witnesses remain
visible. In particular, very small states can give weak bounds with a
boundary optimum; that is numerical search slack, not evidence of code
failure. Additional witness refinement can improve an upper bound without
changing the code or the parameter tuple.

From this directory, run numerical producers sequentially:

```text
powershell -ExecutionPolicy Bypass -File build_activation_q1.ps1
python run_complete_q1_grid_native.py
python register_complete_grid.py --register --require-complete
python build_landscape_db.py
python query_landscape.py export landscape_export.csv
```

The runner skips verified complete batches. It refuses a concurrent run
through an exclusive lock. After a hard interruption, remove that lock only
after confirming the producer has stopped. A final CSV without a final
manifest requires inspection; it is not silently accepted or overwritten.
`register_complete_grid.py` without flags reports current partial coverage.

## Higher occupations

`activation_occupation.py` implements an activation-aware epoch envelope
using the complete A and kernel spectra. Its three state classes are zero,
arbitrary nonzero, and density-bounded nonzero. For epoch input weight j,
the kernel mass k_j/choose(t,j) determines whether zero stays zero. The
nonzero emission moment is averaged over uniform weight-j inputs; taking
the maximum over actual A weights bounds an arbitrary nonzero state.
Both valid termination-mass bounds are intersected before transfer.

Truncated positive matrix-polynomial powering gives the exact
without-replacement region coefficient formula. Its cost depends on the
requested occupation ceiling; truncation preserves lower coefficients.
No independence approximation to positions in the original code is used.

For deterministic outer multiplicities or simultaneous caps, partition the
nonzero weights into bands. For each band g choose a Bernoulli parameter
p_g and a root rho_g satisfying

    rho_g^B >= max_(w in g) A_w / (choose(B,w) p_g^w (1-p_g)^(B-w)).

The all-one band uses its exact endpoint law. Starting with K_j^(0)=R_j,
the entrywise recurrence

    K_j^(q+1) = max_g rho_g ((1-p_g) K_j^(q) + p_g K_(j+1)^(q))

bounds every fixed assignment of q bands. With G bands, multiplying the
terminal moment by choose(L,q) G^q z^(-H) covers all such assignments.
Every requested integer occupation is evaluated; isolated passing samples
do not cover the integers between them. The implementation is nearest
binary64, and a valid formula through q=L is distinct from a completed run.

The transfer follows the activation-aware general-occupation and adaptive
counting derivations in the ba80 task's `bch_rm2sub_bridge` workstream. This
implementation has no runtime dependency on that worktree. Exact GF(16)
tests cover every epoch weight, including kernel weights four and eight,
and every pair of epoch weights. Separate tests check region coefficients,
truncation, adaptive domination and band densities.

For a random full-rank [B,D] outer, Q1 uses its exact expected spectrum.
For higher occupations, `random_spectrum_caps` provides simultaneous
integer shell caps using Markov's inequality and an equal allocation of a
setup-failure budget 2^-60 across B shells. A shell cap is the minimum of
the deterministic support/mass bounds and floor(B 2^60 E[A_w]). On the
good setup event these caps apply to the one constituent reused everywhere.
The bad setup term is charged once in the final union. It must not be
charged once per occupation or hidden inside an ordinary diagnostic margin.
Products of expected spectra are not used for this purpose.

## Completion criteria and next work

The first production batch, `activation_occupation_grid_sparse_v1`, has
450 rows for Q2..16 at message exponent 16, t=64,128,256, s=12,18, and
constituent lengths 128 and 512. A second batch with one shell per band,
`activation_occupation_grid_singleton_v1`, adds 60 rows at t=64 and
length 512. Both include fixed exact spectra and conditional random caps.
Neither batch closes these parameters. For example, the singleton-band
RM(4,9) Q2 bounds at t=64 are about -125.8 and -123.3 bits for s=12 and
s=18, while the exact-support Q2 calculation gives 83.7 and 89.1 bits.
Thus narrower bands alone do not remove the adaptive bound's substantial
slack. Sparse refinement should retain outer-weight compositions more
closely before these bounds are used to select state sizes.

The occupation receipts are
imported with `register_occupation_grid.py DIRECTORY`; this validates their
source hashes and appends explicit catalog entries. Database schema 4
stores each conditional random setup event and its failure budget.

`complete_grid_coverage.csv` lists every candidate, its native status and
the actually evaluated occupation coverage. It separates Q1 completion,
Q2 completion and full-occupation completion. It cannot promote a partial
screen or an exact outer spectrum to an outward distance certificate.

The remaining work is to finish the whole-grid dense cover and refine
promising state/epoch choices with the exact-region composition method.
The coefficient-free typed cover handles large L, while the composition
method avoids its conditioning loss for selected finite occupations.
All occupation ranges must be accounted for
before reporting a full union. Numerical failures and loose bounds remain
part of the landscape rather than being interpreted as code impossibility.
`parameter_cost_frontiers.py` now compares state/epoch cost at explicit
occupation coverage. `EXACT_FAMILY_PROJECTIONS.md` gives exact-spectrum
Q1 fits, largest-size holdouts, and estimates beyond the known lengths.
These reports must be refreshed against the completed grid. Outward
replay is a separate step for selected implementation candidates.
