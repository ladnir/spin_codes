# Weight-40/42 sampler validation

Updated: 2026-09-04

The fast endpoint samplers and independent replay implementations passed their
structural tests and two fixed-size preflights. This completes the implementation
milestone. It establishes neither proposed shell cap.

The mathematical experiment is defined in HIGHER_SHELL_ENDPOINT_REDUCTION.md.
The implementation samples a syndrome in S={0,...,31}, then a subset of 20 or
22 distinct nonzero field elements. It solves the locator equations and counts
agreements. A singular system counts as a non-hit; the sampler never replaces
it with a new trial.

## Independent implementations

`code/test_bch_higher_endpoints.cpp` provides both C++ paths, specialized at
compile time for each weight. It leaves all weight-38 certificate sources intact.

- The primary path uses Newton divided differences for the first remainder,
  polynomial shifts for subsequent remainders, and closed one- or two-variable
  solves. Eight explicit AVX2 lanes count agreements through aligned term tables.
- Replay separately interpolates every remainder with an incremental algorithm.
  It solves the small system by Gaussian elimination and counts agreements with
  eight explicit scalar Horner chains.
- Python independently solves the full 20-by-20 or 22-by-22 system, counts roots
  across the field, and checks the decoding of retained random bytes.

The C++ paths share field tables. An independent carryless-product implementation
checks all 65,536 multiplication entries. The hot paths use fixed-size arrays
and compile-time specialization, with no per-trial heap allocation.

The validation driver is `code/run_higher_endpoint_preflight.py`. Before drawing
random bytes, it hashes the sources, executable, fixtures, and successful kernel
check. It retains all CNG tape blocks. Sampling and replay run sequentially under
the existing BCH experiment mutex. The driver refuses an existing output
directory and cannot accept a statistical shell claim.

## Completed checks

Both C++ paths passed all 4,096 Python fixtures:

- 2,048 known endpoints, covering weights 40 and 42 and both zero and nonzero
  syndromes;
- 2,039 nonsingular non-endpoints;
- nine singular systems.

The known endpoints came from 32 retained witness records and multiple subsets
of their roots. They are positive controls, not independent random samples of
the code's spectrum.

The preflight used 65,536 fresh trials at each weight:

| Check | Weight 40 | Weight 42 |
|---|---:|---:|
| Endpoint hits | 0 | 0 |
| Singular systems, counted as non-hits | 314 | 227 |
| Python-verified audit records | 288 | 288 |
| Singular records among those audits | 32 | 32 |
| Retained tape bytes | 2,097,152 | 2,097,152 |

Full replay matched the agreement histograms, singular counts, byte accounting,
and audit records. The independent Python checks passed all 576 audit records.
The executable also rejected attempts to overwrite an existing result, select
an unsupported weight, use a negative trial budget, or specify a wrong protocol.

These sample counts are far below the prospective statistical budgets.
The zero-hit observations therefore remain implementation diagnostics.

## Reproduction and retained state

The receipt is
`generated/higher_endpoint_preflight_20260904/preflight.json`.
The pre-sampling manifest has SHA-256
`204611cf359e56adeca52350535015f0c7d827e49a0bcf503b18eccf8b5253c8`.
The same directory retains the fixtures, check output, both tapes, full results,
and audit records. Recheck it with:

    python -B code/run_higher_endpoint_preflight.py verify generated/higher_endpoint_preflight_20260904

This verification was rerun successfully. The existing weight-38 certificate
also passed its verifier, including input hashes and all 407 audit records.
Neither verification reruns the long weight-38 experiment.

The files listed in the new manifest are now frozen inputs. Preserve them if
extending the workflow; put the full statistical driver in a new source file.

## Subsequent full experiment (completed)

The planned milestone below has now completed with fresh tapes, zero hits,
matching full replay, and 1,150 independent Python audit checks. See
RANDOM_MODEL_NARROW_CERTIFICATE.md and
generated/higher_endpoint_certificate_20260904/certificate.json. The following
paragraphs preserve the preflight's original plan, not current pending work.

Implement a separate fixed-budget statistical driver, freeze its acceptance
rules, and execute the full experiments sequentially with complete replay.
The current prospective budgets are 436,628,033 trials for weight 40 and
139,431,904 for weight 42, at ideal-IID false-accept error at most 2^-41 each.
Use fresh tapes rather than treating the preflight as an adaptive prefix.

Any statistical acceptance retains the CNG randomness qualification and differs
from a deterministic spectrum theorem. The existing weight-38 test plus those
two error allocations would give familywise false-accept error at most 2^-39.
That error is separate from the target SPIN setup failure probability.
At the preflight milestone, Q2 outward certification and Q>=3 remained open.
Q2 and Q3 have since been certified; the current remaining occupation tail is Q>=4.
