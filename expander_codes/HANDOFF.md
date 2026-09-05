# Expander-code project handoff

Date: 2026-09-03

This document records the state of the expander-code paper, certificates, and
libOTe implementation. It is intended to let a new chat continue without
reconstructing the earlier discussion.

## Start here

The active paper project is:

```text
C:\Users\peter\.codex\worktrees\30ff\permute_conv\expander_codes
```

The active manuscript entry point is:

```text
C:\Users\peter\.codex\worktrees\30ff\permute_conv\expander_codes\paper\main.tex
```

The latest compiled paper is:

```text
C:\Users\peter\.codex\worktrees\30ff\permute_conv\expander_codes\output\exact-enumerator-bounds-for-expander-codes.pdf
```

The dedicated libOTe worktree is:

```text
C:\Users\peter\.codex\worktrees\expander-codes\libOTe
```

It is on branch `codex/expander-codes` at base commit
`d0e6161526ea2db85b474d0485fc2f3cb267cd88`. Its expander-code changes are
uncommitted.

The paper worktree is detached at
`89406e04ffce2f43873885ce90b919906a178e94`. Almost all active manuscript,
certificate, and script work after that commit is uncommitted. Do not reset or
clean either worktree.

## Project objective

The project replaces the loose state-process analysis in the published
Expand--Accumulate (EA) and Expand--Convolute (EC) proofs with exact
input--output enumerators. The intended audience is cryptographers rather than
coding theorists. Define coding terminology at first use and review new paper
prose under Controlled Writing for Cryptography.

The paper separates three contributions:

1. Tighter analysis of the original binary EA and EC ensembles.
2. New regular-expander variants with better finite parameters.
3. Labeled prime-field constructions with certified rate-one-half parameters.

The final first-moment argument is still a union bound. The improvement comes
from retaining the exact recursive-map enumerator and separating support
ranges, rather than assigning one loose estimate to every message.

The current work focuses on minimum distance. It does not provide a decoding
algorithm or a complete application-security reduction.

## Canonical terminology and constructions

Use **left-regular expander** when only the left degree is fixed. Use
**two-sided regular expander** when both left and right degrees are fixed.
Use **regular expander** only when the distinction is already clear.

The two-sided construction divides the `n` expander outputs into `d_L`
regions. Each message coordinate has one edge in every region. Within each
region, every output coordinate has exactly `d_R` incident edges. At rate one
half, `d_L = 2 d_R`.

For binary EC, `d_R` must be odd. If `d_R` is even, the all-one message maps
to zero before the invertible recursive map.

The optimized binary encoder uses the recurrence

```text
y[t] = u[t] + y[t-m] + sum_{j=1}^{m-1} b[t,j] y[t-j].
```

The oldest tap is fixed to one. Earlier discussions and the implementation
call this **wrapped convolution**, following the paper terminology. The time
index itself is not cyclic; out-of-range terms are omitted.

The prime-field EC construction is different. Every expander edge receives an
independent label from `F_p^*`. Every time-and-lag feedback coefficient is
sampled independently from the full field `F_p`. Sampling feedback taps only
from `F_p^*` does not give the reduced-state law used by the proof.

## Manuscript state

The manuscript has 34 compiled pages. Its title is *Exact Enumerator Bounds
for Expand--Accumulate and Expand--Convolute Codes*. The source is split into:

```text
paper/01-introduction.tex
paper/02-constructions.tex
paper/03-enumerator-method.tex
paper/04-binary-ea.tex
paper/05-binary-ec.tex
paper/06-regular-expanders.tex
paper/07-prime-fields.tex
paper/08-certificates.tex
paper/09-discussion.tex
paper/appendix-*.tex
```

`paper/PLAN.md` contains the contribution and section plan.
`paper/CWC_REVIEW.md` records the completed section-by-section writing review.
`TODO.md` and the red `\TODO{...}` markers in the TeX sources record work that
the manuscript deliberately does not claim.

Build the paper from `expander_codes/paper` with:

```powershell
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The historical permute-convolute paper was restored at:

```text
C:\Users\peter\.codex\worktrees\30ff\permute_conv\enumerator_paper
```

Do not turn that directory into the new manuscript. An archival copy also
exists under `expander_codes/notes/previous_draft`.

The current manuscript does not yet contain the new Goldilocks-specific
certificate or implementation measurements.

## Certified binary results

All probabilities below are over the sampled code. The quantity
`-log_2(U)` is the **failure exponent**, not a computational security level.

At rate one half, the current two-sided regular EC frontier is:

| left/right degree | memory | length `n` | cutoff `L` | failure exponent |
|---:|---:|---:|---:|---:|
| 18/9 | 4 | 2,097,270 | 230,742 | 21.6746 bits |
| 14/7 | 6 | 2,097,270 | 230,742 | 26.0875 bits |
| 10/5 | 15 | 2,097,150 | 230,729 | 32.4990 bits |
| 6/3 | 79 | 2,097,150 | 230,729 | 21.4222 bits |

The `10/5`, memory-15 point is the current binary implementation target. Its
cutoff is `0.9999309371` times the binary rate-one-half GV distance. The
`6/3`, memory-79 certificate shows that much smaller degree is possible, but
at an unattractive convolution cost. The simple proxy `d_L+m` is minimized by
the `14/7`, memory-6 point, with value 20.

The corresponding certificate files are under `expander_codes/results`:

```text
binary_biregular_ec_rate_half_d18_m4_gv.json
binary_biregular_ec_rate_half_d14_m6_gv.json
binary_biregular_ec_rate_half_d10_m15_gv.json
binary_biregular_ec_rate_half_d6_m79_gv.json
```

The structural checker and interval verifier for the `10/5` profile run as:

```powershell
python scripts/check_binary_biregular_ec_certificate.py `
  results/binary_biregular_ec_rate_half_d10_m15_gv.json
python scripts/binary_biregular_ec_certificate.py verify `
  results/binary_biregular_ec_rate_half_d10_m15_gv.json
```

Run certificate verifiers serially. The degree-`6/3` verifier is especially
expensive and must not run concurrently with a benchmark.

Two small `10/5`, memory-15 certificates support the current padding floor:

| message size `k` | code length `n` | distance cutoff | failure exponent |
|---:|---:|---:|---:|
| 525 | 1,050 | 10% of `n` | about 20.2526 bits |
| 585 | 1,170 | 10% of `n` | about 24.0318 bits |

The Silent OT configuration pads any positive requested output count. The
mathematical floor is `k=525`. PPRF partition alignment makes `k=585` the
smallest physical protocol instance observed in the integration tests. Do not
claim that every padded length has its own finite certificate. The
implementation comment explicitly identifies this limitation.

For comparison, the manuscript also includes:

- left-regular binary EA with degree 64 and failure exponent 20.0864 bits;
- left-regular binary EC with degree 28, memory 9, and failure exponent
  22.2273 bits;
- Bernoulli EA and EC improvements for the unchanged original ensembles.

The two-sided regular `62/31` EA candidate fails at message support two.
Two-sided regularity should not be claimed to improve every construction.

## Certified prime-field results

For `p = 2^127-1`, rate one half, and the exact floored `p`-ary GV cutoff, the
current frontier is:

| map | memory | left/right degree | length `n` | cutoff `L` | failure exponent |
|---|---:|---:|---:|---:|---:|
| EA | 1 | 30/15 | 2,097,150 | 1,032,064 | 29.9715 bits |
| EC | 3 | 28/14 | 2,097,144 | 1,032,062 | 27.6486 bits |
| EC | 4 | 26/13 | 2,097,134 | 1,032,057 | 71.4724 bits |
| EC | 6 | 24/12 | 2,097,144 | 1,032,062 | 39.5003 bits |
| EC | 12 | 22/11 | 2,097,150 | 1,032,064 | 27.2951 bits |

The manuscript currently treats `22/11`, memory 12 as the flagship because it
minimizes expander degree in this certified table. The practical implementation
scan instead favored `26/13`, memory 4.

The proof uses projective message counting, support grouping, constraint
traces, and a singleton-free-region refinement. Floating-point optimization
chooses decimal markers. The verifier fixes those markers and recomputes every
bound with outward-rounded Arb arithmetic through `python-flint`.

Do not restore the obsolete degree-18 or degree-24, memory-1 certificates.
They applied a projective cap after an averaged transition probability, which
is not valid when different projective messages impose different equations.

## Goldilocks result and scalar extension

The Goldilocks prime is

```text
p = 2^64 - 2^32 + 1 = 18446744069414584321.
```

The same `26/13`, memory-4 profile has a new certificate at the exact floored
Goldilocks `p`-ary GV cutoff:

| parameter | value |
|---|---:|
| `k` | 1,048,567 |
| `n` | 2,097,134 |
| cutoff `L` | 1,015,822 |
| relative cutoff `L/n` | 0.4843858332371703 |
| certified bound `U` | about `1.89001286438e-26` |
| failure exponent | about 85.4517344129 bits |

The certificate is:

```text
results/prime_field_biregular_ec_goldilocks_rate_half_d26_m4_singleton_gv.json
```

Verify it with:

```powershell
python scripts/prime_field_biregular_ec_singleton_certificate.py verify `
  results/prime_field_biregular_ec_goldilocks_rate_half_d26_m4_singleton_gv.json
```

`prime_field_biregular_ec_singleton_certificate.py generate` now accepts
`--prime`. It recomputes the field-specific GV cutoff, clamps field-marker
domains, and reoptimizes the stored markers. A marker-only `rebind` command
also exists, but reusing the large-field markers gave a useless Goldilocks
bound. The checked-in Goldilocks certificate was regenerated and optimized;
it is not the failed rebind result.

Suppose `E/F_p` is a field extension and the generator matrix has entries in
`F_p`. Scalar extension preserves minimum distance exactly. Write an extension
message in an `F_p` basis as `x = a + u b`. Then

```text
xG = aG + u(bG).
```

The support of `xG` is the union of the supports of `aG` and `bG`. It is at
least the base-field minimum distance. Embedding a minimum-weight base-field
codeword into `E` gives the reverse inequality. Therefore the Goldilocks
certificate also proves the same distance for a quadratic scalar extension,
provided all generator coefficients remain in the Goldilocks base field.

This observation does not remove the need for a genuine extension-field type
in the rest of a VOLE implementation.

## libOTe implementation

The binary transpose encoder is:

```text
libOTe/Tools/ExConvCode/RegularEcCode.h
```

It implements two-sided regular regional incidence and the binary wrapped
convolution. The hot convolution loop uses compile-time tap recursion and
packed coefficient masks. Do not replace it with runtime polymorphism or a
generic dynamic loop without measuring the generated code.

The generic field transpose encoder is:

```text
libOTe/Tools/ExConvCode/RegularEcFieldCode.h
```

Its template is:

```cpp
RegularEcFieldCode<LeftDegree, Memory, Scalar>
```

`Scalar` is the type of edge labels and feedback coefficients. The context
supplies field arithmetic, sampling, vectors, and copies. Each edge label is
sampled from `Scalar^*`; each feedback coefficient is sampled from `Scalar`.
The encoder supports values from a larger module when the context can multiply
those values by `Scalar`.

The binary `10/5`, memory-15 encoder is integrated into Silent OT through the
new enum value:

```cpp
MultType::RegularEc10x5x15
```

Relevant integration files are:

```text
libOTe/TwoChooseOne/ConfigureCode.h
libOTe/TwoChooseOne/Silent/SilentOtExtSender.cpp
libOTe/TwoChooseOne/Silent/SilentOtExtReceiver.cpp
```

The adapter distinguishes the physical PPRF prefix from the padded EC code
length. It expands the physical prefix and fills the public padding with
zeros. Both separate and packed choice encodings use the new encoder.

Silent VOLE currently rejects this binary `MultType`, just as it rejects the
binary quasi-cyclic encoder. Those checks are in:

```text
libOTe/Vole/Silent/SilentVoleSender.h
libOTe/Vole/Silent/SilentVoleReceiver.h
```

The encoder tests are registered through the existing test frontend:

```text
libOTe_Tests/RegularEcCode_Tests.cpp
libOTe_Tests/RegularEcCode_Tests.h
libOTe_Tests/UnitTests.cpp
```

They compare transpose encoding with materialized generator matrices over
`F_2`, `Fp31`, and `GF(2^128)`. They also test a two-coordinate Goldilocks
module and verify exact left and right degrees.

The Silent OT tests are part of
`OtExt_Silent_provedCodes_Test` in `libOTe_Tests/SilentOT_Tests.cpp`. That test
covers ExAcc and regular EC, EC dimension rounding, and packed choices.

There must be no expander-code-specific `main` function or standalone test
executable in `libOTe_Tests`. The earlier standalone test and benchmark entry
points were removed after review. libOTe uses `frontend/main.cpp` as the single
frontend entry point.

The EC benchmarks now live in:

```text
frontend/RegularEcBench.cpp
frontend/RegularEcBench.h
frontend/benchmark.h
```

The intended commands are:

```text
frontend_libOTe -bench -regularEc
frontend_libOTe -bench -regularEc -all
frontend_libOTe -bench -regularEc -fp31 -profile 26/13/m4
frontend_libOTe -bench -regularEc -goldilocks
```

## Measurements already obtained

At approximately `n=2^21`, the binary `10/5`, memory-15 block encoder took
about 32 ms.

The Fp31 profile scan produced:

| profile | setup | tables | encode | convolution | expander |
|---|---:|---:|---:|---:|---:|
| 22/11, m=12 | 360 ms | 268 MiB | 276 ms | 129 ms | 146 ms |
| 24/12, m=6 | 317 ms | 236 MiB | 226 ms | 62 ms | 163 ms |
| 26/13, m=4 | 275 ms | 236 MiB | 214 ms | 42 ms | 171 ms |
| 28/14, m=3 | 296 ms | 244 MiB | 222 ms | 31 ms | 188 ms |

The `26/13`, memory-4 profile was the practical winner in this scan.

For Goldilocks at `k=1,048,567` and `n=2,097,134`, the raw encoder benchmark
reported:

| workload | time | approximate output rate |
|---|---:|---:|
| one Goldilocks coordinate | 74 ms | 14.17 million/s |
| two Goldilocks coordinates | 190 ms | 5.52 million/s |
| fused API for extension value plus base value | 258 ms | 4.06 million/s |
| schedule construction | 331 ms | n/a |

The schedule tables occupied about 372 MiB. These are encoder-only numbers.
They exclude PPRF generation, base VOLE, communication, and other protocol
costs.

`dualEncode2` currently invokes the two convolution and expansion kernels in
sequence. It shares the configured schedule but does not yet fuse arithmetic
inside the hot loops.

Never run two benchmarks at the same time.

## Goldilocks extension implementation status

The benchmark represents an extension value as:

```cpp
FVec<Goldilocks, 2>
```

This is a product module, not a quadratic field. It gives the correct encoder
cost when every generator coefficient lies in the Goldilocks base field,
because scalar multiplication acts independently on the two coordinates. Do
not describe `FVec<Goldilocks,2>` as `F_{p^2}`.

There is no genuine Goldilocks quadratic-extension type in this libOTe
worktree or in `C:\Users\peter\repo\gen-BAA`.

A natural implementation is:

```text
Goldilocks2 = Goldilocks[u] / (u^2 - 7).
```

An exact Legendre-symbol check found that 7 is a quadratic nonresidue modulo
the Goldilocks prime. This check has not yet been added as a source test.

`CoeffCtxGoldilocks` is sufficient for the current encoder core. It is not a
complete context for an end-to-end Silent VOLE with
`F=FVec<Goldilocks,2>, G=Goldilocks`: some generic context methods assume that
the value type itself is `Goldilocks`. `CoeffCtxFVec<Goldilocks,2>` does not
solve this problem because it models a product and correctly reports that the
two-coordinate value type is not a field.

The next implementation should either:

1. add a genuine `Goldilocks2` type and a context supporting
   `F=Goldilocks2, G=Goldilocks`; or
2. add a carefully scoped module context if the surrounding VOLE needs only
   base-field scalar multiplication.

The first option is conceptually cleaner for a VOLE over a genuine extension
field. The EC encoder itself does not need extension-by-extension
multiplication when its matrix entries remain in the base field. Other parts
of the protocol may need it.

## Verification and build status

Before the frontend cleanup, the following checks passed:

- materialized-generator tests for the binary and field encoders;
- normal and packed-choice Silent OT tests;
- the Goldilocks raw encoder benchmark;
- rigorous verification of the Goldilocks certificate.

After the frontend cleanup, these translation units compiled successfully:

```text
frontend/RegularEcBench.cpp
libOTe_Tests/RegularEcCode_Tests.cpp
libOTe_Tests/SilentOT_Tests.cpp
libOTe_Tests/UnitTests.cpp
```

`git diff --check` passes apart from line-ending conversion warnings.

A complete `frontend_libOTe` link was attempted in
`out/build/regular-ec-silent`. The build initially required regular and sparse
DPF feature flags because the existing monolithic frontend includes those
types. It then stopped in the unrelated `RingLpn_Tests.cpp`, where
`convertToOle` is unavailable when RingLPN is disabled. No EC compile or link
error was observed. The registered post-cleanup tests have therefore compiled,
but they have not yet been executed through a newly linked frontend binary.

The ignored build cache in `out/build/regular-ec-silent` was reconfigured with
`ENABLE_REGULAR_DPF=ON` and `ENABLE_SPARSE_DPF=ON`. No source files were
changed to work around the unrelated frontend configuration failures.

## Dirty files in the libOTe worktree

At handoff, the intended changes are:

```text
M  frontend/benchmark.h
M  libOTe/TwoChooseOne/ConfigureCode.h
M  libOTe/TwoChooseOne/Silent/SilentOtExtReceiver.cpp
M  libOTe/TwoChooseOne/Silent/SilentOtExtSender.cpp
M  libOTe/Vole/Silent/SilentVoleReceiver.h
M  libOTe/Vole/Silent/SilentVoleSender.h
M  libOTe_Tests/CMakeLists.txt
M  libOTe_Tests/SilentOT_Tests.cpp
M  libOTe_Tests/SilentOT_Tests.h
M  libOTe_Tests/UnitTests.cpp
?? frontend/RegularEcBench.cpp
?? frontend/RegularEcBench.h
?? libOTe/Tools/ExConvCode/RegularEcCode.h
?? libOTe/Tools/ExConvCode/RegularEcFieldCode.h
?? libOTe_Tests/RegularEcCode_Tests.cpp
?? libOTe_Tests/RegularEcCode_Tests.h
```

Review the actual `git status` before editing because the list may change
after this handoff.

## Recommended continuation order

1. Review and commit the current paper and libOTe work in separate commits.
   Preserve the detached paper worktree state by creating a branch first if
   desired.
2. Build libOTe with its standard all-feature frontend configuration. Run
   `RegularEcCode_encode_test`, `RegularEcCode_config_test`, and
   `OtExt_Silent_provedCodes_Test` through `frontend_libOTe`.
3. Run each registered EC benchmark serially and confirm that moving it into
   the frontend did not change measurements.
4. Implement `Goldilocks2` and its coefficient context, then add algebraic
   tests for `u^2=7` and field inversion.
5. Integrate `RegularEcFieldCode<26,4,Goldilocks>` into Silent VOLE for
   `F=Goldilocks2` and `G=Goldilocks`.
6. Benchmark the complete VOLE. Separate schedule construction, PPRF, EC
   encoding, base VOLE, and communication costs.
7. Decide whether the Goldilocks theorem and measurements belong in the main
   paper or an implementation section. If included, preserve the distinction
   between the distance theorem, scalar extension, and protocol performance.
8. Return to the open paper TODOs: recent repeat--accumulate related work,
   certificate-schema unification, remaining degree--memory points, wrapping
   versus nonwrapping comparison, and an application-level reduction.

## Suggested prompt for the next chat

```text
Continue the expander-code project using
C:\Users\peter\.codex\worktrees\30ff\permute_conv\expander_codes\HANDOFF.md.
Read the complete handoff before acting. The paper worktree and the dedicated
libOTe worktree are dirty, so preserve all existing changes. First inspect both
git statuses. Then continue from the recommended continuation order. Do not
add standalone main functions to libOTe_Tests; tests and benchmarks must run
through frontend/main.cpp. Never run two benchmarks concurrently.
```
