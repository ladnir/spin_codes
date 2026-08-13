# Race assignment: close the frozen `g=4` Riffle proof

You are competing independently against Codex to produce the first complete
end-to-end certificate for the frozen `g=4` Riffle construction. Both
contestants start from the repository state recorded on 2026-08-13.

This is not a request for another high-level audit. Work autonomously on the
mathematics, implementation, diagnostics, and final verification. You may
change proof machinery and proof-search code. You may not change the frozen
construction and still claim victory in this race.

Write your running conclusions and final submission to

```text
CLAUDE_G4_RACE_REPORT.md
```

Your isolated Windows worktree is

```text
C:\Users\peter\repo\permute_conv-claude-g4
branch: codex/g4-race-claude
```

Do not edit `C:\Users\peter\repo\permute_conv`; that worktree belongs to the
other contestant. Both race branches must initially resolve to the baseline
commit supplied with this prompt.

Label intermediate claims as **CERTIFIED**, **PROVED**, **DIAGNOSTIC**,
**CANDIDATE**, or **REFUTED**. Do not promote a binary64 result to a proof.

## Finish line

The frozen parameters are

```text
N = 2^21
K = 2^20
d = floor(0.09 N) = 188743
packet width g = 4
packet count M = N/g = 524288
```

Let `Z_d` count nonzero messages whose encoded word has Hamming weight at
most `d`. A winning submission must give a complete, reproducible,
outward-rounded certificate for

```text
E[Z_d] <= 2^-40.
```

The expectation is over the setup randomness declared by the construction:
independent lane bijections, the global packet permutation, and recursive
state permutations. Do not insert an unstated independence assumption.

The certificate must cover every feasible integer packet profile

```text
a = (a0,a1,a2,a3,a4),
sum_j a_j = M,
sum_j j a_j >= 21.
```

There are exactly

```text
binom(M+4,4) - 717
```

such profiles. A uniform profile allocation would require an upper bound of
`-111.41506501642425` bits per profile. You may instead use the sharper
cell-local aggregation already implemented.

Closing the current hard profile is necessary but not sufficient. Closing a
finite source-cell queue is also insufficient. Victory requires the complete
integer profile domain and a final outward upper endpoint no greater than
`-40`.

The other contestant must be able to rerun and adversarially audit the final
certificate. The race clock stops at submission, but victory occurs only
after that audit passes.

## Frozen construction

Both local-code layers remain the binary systematic extended BCH code
`[128,64,22]`. The parameter `g=4` changes only the permutation atom between
the two layers. It does not replace either BCH code.

The outer layout uses the existing three-band construction and its exact
24-dimensional graph word in 128 holes. The inner transform is the existing
systematic BCH recursive chain with the 64-bit accumulator state. Preserve
the existing setup-randomness model.

The packet-profile orbit size is

```text
Q_4(a) = M! / product_j a_j! * product_j binom(4,j)^a_j.
```

This orbit formula depends on independent within-packet lane randomization.
The global packet permutation alone does not justify it. Read the
transitivity lemma and implementation invariant in
`explorations/riffle_group_chain_proof.md` before changing any profile bound.

Construction changes may be recorded separately as useful follow-up ideas.
They do not close this race target.

## Starting proof architecture

The current proof framework has the following sound components:

1. The domain is partitioned into exact positive-support strata.
2. Each stratum has an exact integer hull and rational simplicial mesh.
3. An exact rational BSP may refine each source cell independently.
4. Each leaf uses one fixed witness or one fixed rational convex mixture.
5. Convexity reduces each leaf proof to outward checks at its vertices.
6. Dual column generation prices every eligible finite-atlas witness.
7. Cell-local profile-count bounds feed an outward global log-sum-exp.
8. The final verifier independently reloads and hardens selected witnesses.

Do not weaken these invariants for convenience. In particular:

- Different witnesses at different vertices do not certify a convex cell.
- A pointwise minimum of convex witness bounds need not be convex.
- A sparse zero-fugacity witness cannot leave its exact support face.
- Volume equality alone does not prove that a mesh has no holes or overlaps.
- A rounded tuning profile is discovery data, not proof of a rational vertex.
- A finite candidate cap in the mixture LP is not a complete atlas price scan.

The complete `g=2` proof is the model for hybrid convex-cell and discrete
singleton accounting. Its optimized outward result is
`-62.7194713852` bits. The relevant artifacts are listed in
`PROOF_STATUS.md`.

## Current `g=4` frontier

The best available discovery ledger processes the leading 1024 source cells.
It selects the best of four exact BSP trees per cell and reattaches leaves over
5587 frozen witnesses.

```text
largest updated source-cell term       +16514.50512918122
largest unprocessed source-cell term   +12078.982162066526
leading cell                           s1fr033931
leading rounded profile                [429359,24531,38395,28039,3964]
```

The leading profile has physical weight `201294`. Its current decomposition
is approximately

```text
converged inner probability            -267186.9845
optimized exact-graph linear outer     +268485.2337
combined                                 +1298.2492
uniform profile target                    -111.4151
deficit                                   1409.6643
```

The converged inner value came from a 20,000-step power trajectory on Peach.
It has not yet been packaged as a canonical local artifact. Reproduce or copy
it before relying on it.

The following avenues have already been tested at this profile:

| Test | Observed gain or result |
| --- | --- |
| Extend power iteration from 320 to 20,000 steps | About `116.9` bits gained; inner MGF converges near `30.8161` bits. |
| Replace the adversarial graph-hole tax by exact graph conditioning | About `491.4` bits gained. |
| Reoptimize the exact-graph linear-BL outer tilt | About `0.005` additional bits. |
| Share systematic BCH columns across inner rows | About `0.035` bits beyond the row-specific bound. |
| Retain every four-packet drive-pattern multiset at the inherited tilt | No material gain after the convergence correction. |
| Optimize only the positive 65-state Collatz test vector | The component has no plausible 1400-bit gain unless the operator or tilt changes; its complete converged MGF is about `30.8` bits. |

These measurements trigger a proof-family redesign. Do not spend the first
phase on more Powell depth, more power iterations, or another arbitrary state
vector unless you identify a new mathematical coupling.

## Candidate routes already visible

You are free to choose another sound route. These candidates are supplied so
you start from the same information as Codex.

### A. Exact band-pair-spectrum outer

Exact split spectra are available at

```text
out/ebch85_band01_split_spectrum.csv
out/ebch86_band12_split_spectrum.csv
```

A candidate proof can combine the exact `(0,1)` and `(1,2)` membership
probabilities with rank-one domination and a degree-4 Finner inequality. The
goal is a packet-profile outer bound that retains more three-band structure
than the current linear-BL outer.

An unfinished prototype exists at
`scripts/probe_packet_group_g4_pair_spectrum_outer.py`. Treat it as untrusted.
Before using any output, check:

1. the direction of every rank-one constraint;
2. the degree-4 Finner exponent;
3. the band and tile multiplicities;
4. the packet coefficient extraction;
5. the `t=1` total-mass invariant, which must upper-bound `2^K`;
6. graph-hole and puncture treatment; and
7. outward evaluation after freezing optimized parameters.

The prototype currently omits graph holes and has never been run on the hard
profile. Repairing it is one possible opening, not an assigned route.

### B. Exact drive-pattern inner with joint retuning

`scripts/probe_packet_group_pattern_exact_transfer.py` retains all 4845
multisets of four packet-drive patterns. Its inherited-tilt result did not
improve the leading profile. Jointly retuning the outer and exact-pattern
inner may expose a different saddle point.

The present Python path is expensive. If you choose this route, cache
fugacity-independent structure or implement a fixed-layout native kernel
before a broad search. Keep the independent outward implementation separate.

### C. Outer-profile realizability

The current outer bounds may count packet profiles with far greater
multiplicity than the frozen outer code permits. A sound support-conditioned
enumerator, rank argument, or graph-conditioned shaped spectrum could remove
the leading profile without strengthening the inner bound.

A heuristic failure to find the profile is not an exclusion proof. State the
exact lemma that converts any computation into an upper bound on profile
multiplicity.

### D. Hybrid regional proof

A stronger branch need not dominate globally. It may certify only the
remaining BSP leaves or support strata. Preserve exact ownership or use a safe
overcount when combining the new branch with the existing ledger.

## Required work phases

### Phase 1: reproduce and challenge the frontier

1. Read `PROOF_STATUS.md` and the mathematical notes named below.
2. Verify the hashes of the canonical local artifacts.
3. Reproduce the leading-profile decomposition.
4. Identify the exact inequality responsible for the remaining deficit.
5. Test at least one qualitatively stronger bound on that profile.

Do not launch a broad campaign until the stronger bound demonstrates a
plausible gain on the leading profile.

### Phase 2: convert the gain into reusable witnesses

1. Freeze every optimized parameter and source digest.
2. Separate discovery arithmetic from outward certificate arithmetic.
3. Reattach or replay the affected BSP cells with the stronger witnesses.
4. Measure whether a new profile family becomes dominant.
5. Automate residual extraction and witness generation when the pattern is
   repeatable.

### Phase 3: complete the global cover

1. Cover all source cells and every exact-support stratum.
2. Enumerate terminal lattice profiles exactly when convex coverage is not
   economical.
3. Bind all source artifacts by SHA-256.
4. Run the independent outward verifier.
5. Compute the complete cell-local union with directed rounding.

### Phase 4: adversarial self-audit

Before submission, try to invalidate the result. Check profile counts,
support eligibility, mesh or BSP coverage, mixture weights, graph
conditioning, terminal ownership, logarithmic rounding, and source hashes.
Run the final verifier from a clean process.

## Required reading

Read in this order:

1. `PROOF_STATUS.md`
2. `explorations/riffle_group_chain_proof.md`
3. `explorations/g2_triangular_cover_math.md`
4. `explorations/g4_delaunay_cover_math.md`
5. `CLAUDE_SECOND_OPINION_REPORT.md`
6. `CLAUDE_G4_CLOSURE_SESSION_REPORT.md`

Then inspect these proof-search and verification paths:

```text
scripts/probe_packet_group_g4_cell_bsp.py
scripts/probe_packet_group_g4_cell_bsp_batch.py
scripts/reattach_packet_group_g4_cell_bsp_batch.py
scripts/certify_packet_group_g4_cell_bsp.py
scripts/build_packet_group_g4_support_regular_mesh.py
scripts/certify_packet_group_g4_anchor_mesh.py
scripts/probe_packet_group_joint_inner_opt.py
scripts/probe_packet_group_pattern_exact_transfer.py
scripts/probe_packet_group_exact_graph_linear_bl.py
scripts/probe_packet_group_exact_graph_puncture_outer.py
scripts/upgrade_packet_group_exact_graph_outer.py
scripts/probe_packet_group_g4_pair_spectrum_outer.py
```

The native kernels are diagnostic accelerators:

```text
scripts/packet_group_shared_drive_native.cpp
scripts/packet_group_native.py
scripts/packet_group_drive_stratified.py
```

The outward verifier must not trust their floating-point output.

## Baseline artifacts

Verify at least these files before modifying their descendants:

| Artifact | SHA-256 |
| --- | --- |
| `out/g4_cell_bsp_top1024_best_collatz_merged_round4_reattached.json` | `69ff23973c9d38a1a93ac0378b4b38e834b83b0b913004340bb8870ad027fd17` |
| `out/g4_cell_bsp_top1024_best_collatz_merged_round4_reattached_failures.json` | `de996e74ce8622926e20deacf916a897d5f27f7bf1f6cd112d9f2261f6891930` |
| `out/g4_leader_429359_joint_fullcollatz_multistart.json` | `f7951ff6e116de9d81411b9abd214642f2fed9d6490735ec7ddc098a8fc31675` |
| `out/g4_leader_429359_joint_fullcollatz_exact_graph.json` | `1d7b5271ca477784f9e40b4796604030f76a4777798086387187690be0b90f1e` |
| `out/g4_leader_429359_exact_graph_linear_optimized.json` | `b9141dad0e74e170822a9c3981b240ee4e1c91ad4663a5a216f21366c31d5ee2` |

The merged BSP artifact references additional sources under
`/tmp/permute-conv-g2/out/` on Peach. Their expected hashes are embedded in
the artifact. If a source is missing locally, copy it from Peach or regenerate
it. Do not substitute a similarly named file.

## Compute and implementation rules

- Never run two benchmarks or long proof-search jobs simultaneously.
- Use all available cores for one coordinated batch when parallelism helps.
- Benchmark a native replacement against the reference before relying on its
  performance.
- Preserve a Python or exact-arithmetic reference for every accelerated
  diagnostic kernel.
- Use deterministic seeds and record complete commands.
- Checkpoint long jobs and write partial artifacts atomically.
- A failed diagnostic is information, not a proof that the construction
  fails.

## Peach access and shared-compute protocol

Peach is the Linux compute host used for the larger witness searches. The
configured connection is

```text
host: peach48.devcore4.com
port: 9022
user: prindal
Claude workspace: /tmp/permute-conv-g4-claude
Codex workspace: /tmp/permute-conv-g4-codex
shared historical artifacts: /tmp/permute-conv-g2/out
```

The supported local connection helper is

```text
C:\Users\peter\.codex\skills\connect-to-peach\scripts\connect-peach.ps1
```

Open an interactive shell from PowerShell with

```powershell
powershell -ExecutionPolicy Bypass -File "$env:USERPROFILE\.codex\skills\connect-to-peach\scripts\connect-peach.ps1"
```

Run one remote command with

```powershell
powershell -ExecutionPolicy Bypass -File "$env:USERPROFILE\.codex\skills\connect-to-peach\scripts\connect-peach.ps1" -Command "cd /tmp/permute-conv-g4-claude && git status --short"
```

The helper binds the configured PuTTY key and the expected host-key
fingerprint. Do not bypass its host-key check. Exact raw connection and
WinSCP transfer commands are documented in

```text
C:\Users\peter\.codex\skills\connect-to-peach\references\connection-details.md
```

For reference, transfers use
`C:\Users\peter\OneDrive\tools\WinSCP.com`, the configured `.ppk`, and the
host-key fingerprint from that reference file. Prefer the documented command
instead of inventing `scp` or converting the key.

The established Python environment on Peach uses the repository scripts plus
the local dependency directory:

```bash
cd /tmp/permute-conv-g4-claude
PYTHONPATH=scripts:/tmp/permute-conv-g2/pydeps OMP_NUM_THREADS=1 python3 <script> <arguments>
```

Set the script's explicit worker count for parallel batch work. Previous
campaigns used 16 workers effectively. Inspect `lscpu` before changing that
policy. `OMP_NUM_THREADS=1` prevents each worker from creating a nested BLAS
thread pool.

All contestants share Peach. Before a long run, acquire the common lock
directory atomically. A suitable one-shot command has this shape:

```bash
set -eu
lock=/tmp/.g4-proof-race-long-job-lock
if ! mkdir "$lock" 2>/dev/null; then
  echo "another long g=4 job owns the lock"
  cat "$lock/owner" 2>/dev/null || true
  exit 75
fi
printf 'owner=claude\nstarted_utc=%s\ncommand=%s\n' "$(date -u +%FT%TZ)" '<command>' > "$lock/owner"
cleanup() {
  rm -f "$lock/owner"
  rmdir "$lock"
}
trap cleanup EXIT INT TERM
cd /tmp/permute-conv-g4-claude
PYTHONPATH=scripts:/tmp/permute-conv-g2/pydeps OMP_NUM_THREADS=1 python3 <script> <arguments>
```

Do not remove a lock merely because it exists. First inspect its owner and
confirm that the recorded job is no longer running. Short read-only probes,
hash checks, compilation checks, and syntax checks do not require the long-job
lock. A benchmark, large atlas pass, tuning batch, or full verifier
regeneration does require it.

Copy every winning artifact back to the Windows workspace. Verify its SHA-256
on both hosts before citing it. Preserve remote checkpoints until the opposing
audit finishes.

Do not edit `/tmp/permute-conv-g2`. That older checkout is shared historical
state. Read or copy missing source artifacts from its `out/` directory into
your isolated workspace, then record their original and copied hashes.

## Winning submission

Your final report must contain:

1. the new lemma or inequality, with its probability space and assumptions;
2. a proof or a precise reduction to independently checked finite arithmetic;
3. the complete certificate architecture;
4. exact commands needed to regenerate the final artifact;
5. wall-clock and hardware information for long runs;
6. final artifact paths and SHA-256 hashes;
7. the outward interval for `log2 E[Z_d]`;
8. a list of every remaining assumption or authenticated external input;
9. an adversarial soundness audit; and
10. a concise explanation of why the result applies to the frozen
    construction.

If you find a genuine construction failure, submit the explicit witness and a
verifier instead. A low-weight codeword of weight at most `188743` is decisive.

Do not stop after improving the leading profile. Continue until the complete
outward certificate passes, a decisive counterexample is verified, or a
specific mathematical obstruction prevents further progress.
