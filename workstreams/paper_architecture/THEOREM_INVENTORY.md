# SPIN theorem inventory

## Status labels

- **AVAILABLE:** a suitable statement and proof already exist.
- **ADAPT:** useful proof material exists, but the statement or construction
  must change for SPIN.
- **TARGET:** the paper needs this result; no complete proof is asserted here.
- **DIAGNOSTIC:** numerical or modeled evidence only.

These labels describe source readiness. They are not theorem environments and
do not promote a result beyond its current evidence.

## Framework

| ID | Object | Required content | Status | Placement |
| --- | --- | --- | --- | --- |
| F1 | Definition: SPIN encoder | `E_N=I_N o Pi_N o O_N`; inputs, output, setup randomness, admissible lengths, and rate | AVAILABLE; written in main | main |
| F2 | Definition: interface class | A class `tau` retaining the outer information needed by routing and the inner transfer | AVAILABLE; written in main | main |
| F3 | Theorem: class-indexed first moment | `E[Z_d]=sum_tau A_tau^out Q_tau(d)` under an explicit construction distribution | AVAILABLE; written in main | main |
| F4 | Lemma: uniform-slice law | A uniform permutation maps every fixed weight-`w` word uniformly over its Hamming slice | AVAILABLE | main, short proof |
| F5 | Corollary: weight-indexed first moment | `E[Z_d]=sum_w A_w^out p_w^in(d)` for the uniform-interleaver setting | AVAILABLE | main |
| F6 | Corollary: envelope bound | Replace exact terms by rigorous upper bounds and add any outer failure event | ADAPT | main |
| F7 | Lemma: admissible-length wrapper | Padding, shortening, or neighboring block sizes extend a theorem beyond its native divisibility conditions | TARGET | appendix |

The class-indexed theorem should be proved by linearity of expectation. It
must state the probability space over outer setup, interleaver setup, and inner
setup. The weight-indexed theorem is a specialization, not the universal
framework.

## Accumulator SPIN

| ID | Object | Required content | Status | Placement |
| --- | --- | --- | --- | --- |
| A1 | Definition: random-block outer | Sample an independent uniform linear injection in each of the `M=N/B` positions | AVAILABLE; written in main | main |
| A1a | Lemma: uniform random-block enumerator | Give the exact expected local input--output enumerator for an unconditioned uniform `K_B`-by-`B` matrix | AVAILABLE from BA codes; written in main | main |
| A1b | Formula: BA independent-block enumerator | Count total input and output weights by the number of active input blocks when the diagonal matrices are sampled independently | AVAILABLE from BA codes; written in main | main |
| A3 | Definition: accumulator | State the prefix-sum recurrence, boundary convention, and output length | AVAILABLE | main |
| A4 | Lemma: exact accumulator enumerator | Count inputs of weight `w` producing output weight `h` | AVAILABLE | main |
| A5 | Lemma: sharp accumulator contraction | For output fraction `delta`, prove `p_{N,h}(floor(delta N)) <= rho(delta)^h`, where `rho(delta)=2 sqrt(delta(1-delta))` | AVAILABLE; written in main | main |
| A6 | Theorem: Accumulator SPIN distance | For independent uniform injections, prove positive relative distance under an explicit logarithmic block schedule | AVAILABLE; written in main; rate-half point `B>=3 log_2 N`, `delta=0.0037` is interval certified | main |
| A7 | Proposition: Accumulator SPIN cost | State setup storage, encoding work, and memory in the selected cost model | AVAILABLE for the direct encoder | scaling section |

Accumulator SPIN is the one-interleaver--accumulator truncation of the
Block-Accumulate architecture. A1a and A1b are credited to that work. The
independent-block injection proof from the finite/asymptotic theory workstream
now proves A6. The repeated-constituent option is related work, not a theorem
target in the present SPIN paper.

Chosen-Block BAA is the direct source for repeating one fixed constituent in
all outer positions. It treats certified deterministic blocks. The SPIN paper
cites that option but analyzes independent random-block injections.

## Random SPIN

| ID | Object | Required content | Status | Placement |
| --- | --- | --- | --- | --- |
| R1 | Definition: all-random-tap recursive inner | Sample an independent `m`-bit tap vector at each time; use zero initial state and non-wrapping termination | AVAILABLE; written in main | main |
| R2 | Lemmas: exact transfer matrix and random-inner exponents | Give the exact bivariate transform, Perron-root identity, linear-range exponent, and sparse-range bound | AVAILABLE; statements in main, proofs in appendix |
| R3 | Theorem: Random SPIN distance | Combine independent uniform injections with R2 | AVAILABLE; written in main | main |
| R4 | Concrete asymptotic parameter point | Rate `1/2`, `B>=17 log_2 N`, `m=ceil((51/50)log_2 N)`, `delta=0.11002` | AVAILABLE; outward interval margins certified | main |
| R5 | Proposition: Random SPIN cost | Give encoding work and state size under the theorem schedule | AVAILABLE for the direct encoder | scaling section |

The proved recurrence uses fully sampled time-varying taps. It replaces the
archived fixed-tap candidate, whose proof audit found a reversed monotonicity
step and an invalid independence decomposition. Silver is the earlier source
for a structured recursive encoder of convolutional form. Expand-Convolute is
the closer prior architecture for independently sampled time-varying
convolution taps, but it uses a sparse expander outer.

## Structured SPIN

| ID | Object | Required content | Status | Placement |
| --- | --- | --- | --- | --- |
| S1 | Definition: Structured SPIN family | Define admissible native lengths, the selector `theta_m`, and the outer, interleaver, and inner rules | AVAILABLE; written in main | main |
| S2 | Proposition: reference-member recovery | Show that `m=21` and the recorded parameters recover the frozen encoder bound to the manifest | TARGET | main statement; verification appendix |
| S3 | Lemma: ParityFanout transition | Derive the exact output-weight law from input weight under uniform disjoint `S,T` | AVAILABLE; written for the finite specialization | appendix |
| S4 | Corollary: expected transformed spectrum | Apply S3 to a certified base spectrum or envelope | TARGET | main |
| S5 | Certificate: base outer spectrum | Authenticate the actual base constituent used by the frozen member | TARGET | artifact and appendix |
| S6 | Definition: structured interleaver | Define the factored permutation as one map and state all setup randomness | AVAILABLE; written in main | main |
| S7 | Lemma: structured route law | Reduce the routed profile to active-position fraction and mean constituent weight at subexponential cost | AVAILABLE; focused proof written | appendix; exhaustive witnesses in artifact |
| S8 | Definition: structured inner | State the mathematical forward recurrence and termination independently of its transposed implementation | AVAILABLE; written in main | main |
| S9 | Theorem: structured-inner transfer | Bound the low-output probability for every class produced by S7 | AVAILABLE; three-state, fixed-occupation, and growing-sparse proofs written | appendix; exhaustive witnesses in artifact |
| S10 | Theorem: conditional finite Structured SPIN | Combine S4, S7, and S9 through F3 for explicit finite parameters | TARGET | main |
| S11 | Certificate: frozen reference distance | Outward-rounded complete ledger for `N=2^21` using the actual outer spectrum | DIAGNOSTIC at present | main result after completion |
| S12 | Theorem: asymptotic Structured SPIN | One-sampled Golay--BA-3 outer, structured route, and RM2Sub-S19 inner have rate `1/2` and relative distance `0.11` on the native schedule | AVAILABLE; written in main; current certified block constant `39/4` | main |
| S13 | Proposition: Structured SPIN complexity | Evaluate the outer, routing, and inner costs under `theta_m` | AVAILABLE for any fixed base constituent with a linear-time tiled encoder | main |
| S14 | Theorem: linear-time Structured SPIN | Exhibit a schedule satisfying S12 with `O(N)` ordinary and transposed encoding | AVAILABLE; included in S12 | main |
| S15 | Lemma: requested-length wrapper | Zero extend the largest preceding native member to any sufficiently large requested output length with `o(1)` rate and distance loss | AVAILABLE; written in main | main |

S11 must not be stated as complete while its source spectrum is modeled or its
arithmetic is nearest binary64. This finite-member restriction does not limit
S12, whose one-sampled Golay--BA-3 certificate is independent of ParityFanout.

## Certificates and implementation evidence

| ID | Object | Required content | Status | Placement |
| --- | --- | --- | --- | --- |
| C1 | Definition: certificate schema | Parameters, ensemble, class partition, exact or outward bounds, verifier, and artifact hashes | TARGET | finite-certificate section |
| C2 | Lemma: certificate acceptance | A verifier accepting a ledger with total below one implies existence of a code with the stated distance | TARGET | main, short proof |
| I1 | Implementation correspondence | Each optimized stage computes the mathematical frozen reference map | ADAPT from correctness records | appendix |
| I2 | Correctness report | Independent dense-inner and end-to-end tests, with scope stated | AVAILABLE as measured evidence | evaluation |
| I3 | Performance report | Operation ledger, machine conditions, trial count, and median runtime | AVAILABLE as measured evidence | evaluation |

Correctness tests and benchmarks are evidence about the implementation. They
do not discharge S5, S7, S9, or S11.

## Minimal theorem spine

The initial paper should organize its principal statements in this order:

1. F3: exact first-moment theorem;
2. A4, A5, and A6: accumulator mechanism and asymptotic family;
3. R2, R3, and R4: stronger random inner and asymptotic family;
4. S1: Structured SPIN family definition;
5. S7, S9, and S12--S15: scalable routing, transfer, distance, complexity,
   and requested-length wrapper;
6. S3--S5 and S10--S11: ParityFanout specialization and its still-open finite
   certificate.

The random construction lines use independent uniform injections throughout.
Repeated fixed constituents remain part of the BA-related-work discussion,
not part of the principal theorem spine.
