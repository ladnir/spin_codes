# Source migration map

## Migration rule

The new paper is organized around the three SPIN families. Existing text is
reused only when its mathematical object matches the new construction or when
it provides a proof technique that can be restated honestly.

Use four actions:

- **reuse:** the statement and object already match;
- **adapt:** the mathematics is useful but notation, scope, or construction
  must change;
- **source only:** consult the proof, but do not migrate its theorem claim;
- **exclude:** historical material that does not support the new spine.

## Compiled paper

| Source | New destination | Action | Reason or constraint |
| --- | --- | --- | --- |
| `intro.tex` | Introduction | source only | Its motivation is useful, but its named construction results are obsolete for SPIN. |
| `prelim.tex` | Preliminaries | adapt | Retain only notation used by the new interface. |
| `framework.tex`, bad-codeword count and Markov step | SPIN framework | reuse | This is the exact first-moment core. Generalize the index from weight `w` to class `tau` before specializing back to `w`. |
| `framework.tex`, uniform-slice lemma | SPIN framework | reuse | Applies to Accumulator SPIN and Random SPIN. Do not apply it to the structured interleaver. |
| `framework.tex`, coarse outer and inner assumptions | Construction sections or appendix | adapt | Introduce envelopes only when an instantiation needs them. |
| `innerAcc.tex`, accumulator definition | Accumulator SPIN | reuse | Normalize notation and boundary conventions. |
| `innerAcc.tex`, exact enumerator | Accumulator SPIN | reuse | Central transparent example of the framework. |
| `innerAcc.tex`, tail and integration theorem | Accumulator SPIN | adapt/source only | The tail is reusable. The combined theorem uses the old sliding outer and does not prove the random-block result. |
| `randomDenseConstruction.tex` | Random SPIN | source only | Its three-part organization is useful, but its outer is not the SPIN random-block outer. |
| `outerDense.tex` | Random-B outer analysis | source only | Techniques may help, but the sliding convolutional outer is a different ensemble. |
| `innerDenseScalar.tex` | Random SPIN inner | adapt | Reuse the recurrence analysis after fixing the Random SPIN inner definition. |
| `integrationDense.tex` | Random SPIN proof | source only | Reuse exponent and range-splitting techniques. Do not migrate the `0.109` theorem unchanged. |
| `localCodeStructured.tex` | Structured SPIN | source only | Replace the old RM/EBCH construction story with the current Structured SPIN family and reference member. |
| `localCodeOuter.tex` | Known-spectrum outer and certificate format | adapt | Reuse block-spectrum bookkeeping and exact finite-sum techniques. |
| `localCodeInner.tex` | Structured inner | source only | The current frozen inner is RM2Sub-S19 and needs its own forward definition. |
| `localCodeCertificate.tex` | Finite certificates | source only | Reuse certificate organization and outward arithmetic patterns, not the old construction theorem. |
| `localCodeProjections.tex` | Related designs or appendix | exclude from main | Projection heuristics are not part of the core SPIN claim. |

## Block-Accumulate source paper

| Source | New destination | Action | Reason or constraint |
| --- | --- | --- | --- |
| `BA_paper/blockAcc.tex` | Accumulator SPIN attribution and construction relation | adapt and cite | Describe Accumulator SPIN as the one-stage truncation of the BA architecture. |
| `BA_paper/block.tex`, random block enumerator | Accumulator SPIN outer analysis | adapt and cite | Import the local enumerator and the independent-block global formula. The theorem conditions each block map to be injective and therefore uses its own local generating function. |
| `BA_paper/turbo.tex`, serial-concatenation enumerator | Accumulator SPIN proof integration | source only | Reuse only after the final independent-versus-shared block ensemble is fixed. |

## Finite and asymptotic theory workstream

These sources live in
`C:/Users/peter/.codex/worktrees/3061/permute_conv/workstreams/finite_asymptotic_theory/`.

| Source | New destination | Action | Reason or constraint |
| --- | --- | --- | --- |
| `ACCUMULATOR_SPIN_ASYMPTOTIC.md` | Accumulator SPIN theorem and proof | reuse with normalized notation | Proves the canonical independent-injection ensemble, including the sharp accumulator contraction and the rate-half point `B>=3 log_2 N`, `delta=0.0037`. |
| `certify_accumulator_spin_asymptotic.py` | Reproducibility support | reuse as certificate script | Gives an outward interval check of the stated strict parameter margin. |
| `RANDOM_SPIN_PROOF_AUDIT.md` | Random SPIN theorem and proof appendix | reuse with normalized notation | Proves the all-random-tap, independent-injection ensemble by an exact transfer matrix plus sparse/linear/top input-weight bounds. The archived fixed-tap candidate is explicitly excluded. |
| `certify_random_spin_linear_exponent.py` | Reproducibility support | reuse as certificate script | Certifies the rate-half point `B>=17 log_2 N`, `m=ceil((51/50)log_2 N)`, `delta=0.11002`. |
| `SINGLE_SAMPLED_BA_RM2SUB_D11.md` | Canonical Structured SPIN instantiation and theorem | reuse with normalized notation | Certifies the one-sampled Golay--BA-3/RM2Sub-S19 family at rate `1/2`, relative distance `0.11`, and linear ordinary and transposed work. The current schedule uses `b_m>=(39/4) log_2 N_m`. |
| `SINGLE_SAMPLED_BA_RM2SUB_CERTIFICATE_MANIFEST.json` and named verifier scripts | Structured SPIN certificate appendix and artifact bundle | imported and hash verified | The paper-owned snapshot binds the outer good event, occupation split, weight-coupled fixed-occupation check, outward endpoints, arithmetic semantics, and source hashes. |

The two theorem sources define the random SPIN baselines because they have the same
outer--permutation--inner composition. Their outer setup samples a fresh
uniform injection in each diagonal position.

## Adjacent prior constructions

| Source | New destination | Action | Reason or constraint |
| --- | --- | --- | --- |
| `C:/Users/peter/repo/gen-BAA/main.tex` | Accumulator SPIN related work | adapt and cite | Chosen-Block BAA is the direct precedent for repeating one fixed constituent code. Its deterministic Golay, Reed--Muller, and BCH-subcode instances are a related option, not the random SPIN ensemble. |
| Couteau--Rindal--Raghuraman, *Silver* | Random SPIN attribution | cite and compare | Earlier structured recursive encoder obtained by solving a banded lower-triangular system; convolutional in form, but not the iid time-varying Random SPIN ensemble. |
| Raghuraman--Rindal--Tanguy, *Expand-Convolute Codes* | Random SPIN attribution | cite and compare | Closest prior architecture for a sparse outer followed by independently sampled time-varying convolution taps. Its expander outer differs from Random SPIN. |

## Planning documents

| Source | New destination | Action | Reason or constraint |
| --- | --- | --- | --- |
| `SPIN_NAMING.md` | Names throughout | reuse | Controls the paper-facing family names and single-interleaver terminology. |
| `RIFFLE_NEXT_WORK_ROADMAP.md`, Track B | Global paper spine | reuse | Supplies the three-family progression. |
| Same roadmap, Tracks C--E | Finite lengths, asymptotics, and complexity | adapt | These are research targets, not achieved claims. |
| `PAPER_RESTRUCTURE_PLAN.md` | Proof-source index | source only | It organizes the previous manuscript, not the restarted SPIN paper. |
| `PAPER_INVENTORY.md` | Proof-source index | source only | Use its status notes to avoid accidental claim migration. |

## Frozen Structured SPIN record

All paths in this table are relative to
`constructions/riffle_parityfanout31x33_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s19/`.

| Source | New destination | Action | Reason or constraint |
| --- | --- | --- | --- |
| `MAIN_CODE_FREEZE.md` | Frozen reference member | reuse | Authoritative scope and current proof gates for one member. It is not the family definition. |
| `README.md`, construction formula | Structured SPIN reference definition | adapt | Restate the forward outer, permutation, and inner maps in mathematical order. |
| `README.md`, ParityFanout map | Known-spectrum outer | adapt | Derive the exact transition from base weight to expected transformed weight. |
| `PROOF_STATUS.md` | Finite certificate status | reuse with labels | Preserve the modeled-source and nearest-binary64 qualifications. |
| `PERFORMANCE.md` | Implementation and performance | reuse | Preserve machine, trial, and measurement conditions. |
| `SOURCE_MANIFEST.json` or the manifest named by the freeze | Reproducibility appendix | reuse | Bind the reference-member proposition to exact sources and hashes. |
| Frozen source header defining the outer and routing | Mathematical definition and implementation correspondence | adapt | Extract the forward map. Keep packing and transposed evaluation out of the construction definition. |
| Frozen RM2Sub-S19 header | Structured-inner definition | adapt | State the forward recurrence before mapping it to optimized code. |
| Correctness and benchmark receipts | Evaluation | reuse as measured evidence | These receipts do not prove minimum distance. |

## Adjacent ParityFanout artifacts

| Source | New destination | Action | Reason or constraint |
| --- | --- | --- | --- |
| `constructions/riffle_parityshear12_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s19/receipts/parityfanout31x33_outer256_d38_expected_spectrum.json` | Known-spectrum outer | adapt/source only | Records the expected transformed spectrum and Bernoulli envelope under independent fanout sampling. The input spectrum is modeled. |
| Adjacent first-moment ledgers | Structured finite evaluation | reuse with `DIAGNOSTIC` label | Preserve their range partition and reported margins. Do not call them certificates until outward evaluation and the actual source spectrum are supplied. |

## Material intentionally left out of the new spine

The following material may remain in historical appendices or exploration
notes, but it should not shape the main narrative:

- the random sliding dense outer as if it were the SPIN random-block outer;
- the RM/EBCH full-split construction as if it were Structured SPIN;
- projection-only BCH spectra;
- sampled-grid or nearest-binary64 results stated as theorems;
- implementation details introduced before the mathematical forward map;
- the legacy `Riffle` names as paper-facing family names.

## Completed migration passes

The first drafting pass migrated these units:

1. the exact first-moment proof from `framework.tex`;
2. the uniform-slice lemma from `framework.tex`;
3. the accumulator definition and enumerator from `innerAcc.tex`;
4. the random-inner recurrence and local lemmas from `innerDenseScalar.tex`;
5. the frozen reference-member facts from its freeze and README;
6. the ParityFanout transition rule;
7. the evidence-status and certificate conventions;
8. the frozen implementation and performance facts.

The next proof pass migrated the two independent-block asymptotic theorems
from the finite/asymptotic theory workstream. The third pass migrated the
one-sampled Golay--BA-3/RM2Sub-S19 theorem as the canonical scalable Structured
SPIN instantiation. The frozen ParityFanout ledger remains diagnostic and is
not a source for that theorem.
