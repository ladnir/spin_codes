# SPIN terminology and notation

## Family names

| Canonical term | Meaning | Avoid |
| --- | --- | --- |
| SPIN codes | Single-Permutation INterleaved codes | lowercase `spin codes`; using SPIN for a construction without one global interleaver |
| Accumulator SPIN | independently sampled uniform-injection blocks, uniform interleaver, accumulator inner | accumulator warmup as the formal name |
| Block-Accumulate (BA) codes | prior architecture with a random block-diagonal outer and one or more interleaver--accumulator stages | presenting Accumulator SPIN without BA attribution |
| truncated BA code | Accumulator SPIN viewed as retaining the first interleaver--accumulator stage of BA | omitting BA attribution |
| Chosen-Block BAA | prior specialization that repeats one certified constituent code in all outer positions | generalized BA as a paper title; claiming that its deterministic constituent is sampled like the SPIN random block |
| Random SPIN | independently sampled uniform-injection blocks, uniform interleaver, and all-random time-varying recursive inner | dense+dense; random dense construction |
| Structured SPIN | block-diagonal constituent outer, factored structured interleaver, and batch-state structured inner | Fast SPIN; Riffle as a paper-facing synonym |
| frozen reference member | Structured SPIN at `N=2^21`, `B=256`, `K_B=128`, `t=128`, `s=19` | frozen family; final Structured SPIN |
| scalable Structured SPIN instantiation | one sampled BA-3 constituent reused diagonally, structured route, and RM2Sub-S19 inner; the base constituent and base length are parameters | treating the current Golay certificate as an architectural restriction |
| current certified scalable instantiation | the scalable interface with extended Golay base length `a=24` and certified schedule `b_m>=(39/4)log_2 N_m` | transferring its theorem to another base constituent without a new spectrum certificate |
| finite-length specialization | a member that changes constituent or routing parameters to optimize one target length, such as the frozen ParityFanout member | using it as the definition of the scalable family |

The paper may introduce another family name later only if a new component
combination supports a distinct theorem or design objective.

## Components

| Canonical term | Meaning | Avoid |
| --- | --- | --- |
| outer encoder | first map `O_N` in the SPIN composition | outer code when the map and image must be distinguished |
| random-block outer | one independently sampled uniform linear injection `F_2^{K_B} -> F_2^B` in each of the `M=N/B` diagonal positions | one shared sampled block; unconditioned random matrices |
| BA random-block outer | one independently sampled unconditioned random matrix in each diagonal position | conditioned uniform injections |
| repeated constituent | optional related architecture that reuses one fixed local code in every position, as in Chosen-Block BAA | the Accumulator SPIN or Random SPIN outer |
| chosen-block outer | one fixed constituent code reused in all outer positions | implying that the constituent must be random or deterministic without stating which ensemble is active |
| shared sampled constituent | one sampled local code reused in every outer position; its setup expectation remains outside powers of its realized enumerator | treating different diagonal positions as independent samples |
| local map | injective map `G_{m,j}:F_2^{k_m} -> F_2^{b_m}` in one Structured SPIN diagonal position | constituent without stating whether setup is shared or local |
| base constituent | rate-half map `H_a:F_2^(a/2) -> F_2^a` tiled inside one scalable BA-3 local map | fixing `a=24` in the family definition; calling every base constituent Golay |
| known-spectrum outer | structured local rule with a rigorous fixed spectrum, expected spectrum, transition law, or sufficient envelope | inferring a spectrum from minimum distance alone |
| interleaver | the single global permutation `Pi_N` | calling its implementation factors separate interleavers |
| uniform interleaver | a uniform permutation of `[N]` | random permutation when the distribution is relevant |
| structured interleaver | the declared factored SPIN permutation distribution | uniform permutation; uniform Hamming-slice routing |
| accumulator inner | rate-one prefix-sum recursion | accumulator code without its boundary convention |
| all-random-tap recursive inner | non-wrapping recurrence whose full `m`-bit tap vector is sampled independently at each time; inner ensemble of the proved Random SPIN theorem | fixed-tap candidate; time-invariant convolution |
| fixed-tap candidate | archived recurrence that fixes one tap coefficient; not used in the current theorem because its proof audit found two gaps | Random SPIN theorem recurrence |
| Silver | prior code whose encoder solves a banded lower-triangular system recursively | claiming that Silver samples iid time-varying convolution taps |
| Expand-Convolute | prior code using a sparse expander outer and an iid time-varying randomized recursive convolution | treating its expander outer as the Random SPIN block outer |
| structured inner | the batch-state recurrence `Y_i=X_i+A(Q_i)`, `Q_{i+1}=alpha_i Q_i+C(X_i)` with declared field multipliers | Riffle inner as the paper-facing term; causal recursive inner |

## Outer-spectrum terms

| Canonical term | Meaning | Avoid |
| --- | --- | --- |
| base outer constituent | local code or map before ParityFanout | complete Structured SPIN outer |
| ParityFanout transform | invertible sparse output map defined by sampled disjoint sets `S,T` | treating it as an unstructured random outer |
| ParityFanout transition law | conditional distribution of output weight given input weight | exact enumerator of every sampled fanout realization |
| expected transformed spectrum | base spectrum propagated through the transition law | fixed-realization spectrum |
| modeled base spectrum | even-floor source model used by the current receipt | BCH spectrum; certified spectrum |
| certified outer envelope | rigorous upper envelope derived from the actual base constituent and declared transforms | nearest-binary64 modeled receipt |

ParityFanout is part of the known-spectrum architecture because its transition
law is analytic. The current finite distance ledger remains diagnostic because
its base spectrum is modeled and its final arithmetic is not outward rounded.

## Probability and evidence terms

| Canonical term | Meaning | Avoid |
| --- | --- | --- |
| construction distribution | joint distribution of all randomized setup choices | randomness without naming its source |
| realized encoder | one fixed outcome of the construction distribution | ensemble |
| expected spectrum | expectation over the declared construction randomness | spectrum of every realization |
| finite certificate | exact or outward-rounded complete first-moment bound below one | numerical evidence; sampled grid |
| diagnostic | modeled, sampled, or nearest-rounding evidence | theorem; certificate |
| implementation evidence | correctness tests, checksums, operation counts, and benchmarks | mathematical proof |

Whenever the paper takes an expectation, it should identify whether the
randomness comes from the outer, interleaver, inner, or several independent
components.

## Core notation

| Symbol | Type and meaning |
| --- | --- |
| `N` | complete outer-output and inner-output length |
| `K` | message dimension of the complete encoder |
| `B` | outer block length |
| `K_B` | message dimension of one outer block |
| `M=N/B` | number of outer-block positions |
| `O_N` | outer encoder |
| `Pi_N` | one global interleaver |
| `I_N` | recursive inner encoder |
| `E_N=I_N o Pi_N o O_N` | complete SPIN encoder |
| `m` | positive index of a native Structured SPIN member; for the current certified selector, `L_m=128m` |
| `theta_m` | selector for the Structured SPIN parameters at native index `m` |
| `b_m,k_m,L_m` | Structured SPIN outer length, local dimension, and number of diagonal positions, with `L_m=N_m/b_m` |
| `G_{m,j}` | realized injective local map in diagonal position `j` |
| `A_m,C_m` | fixed Structured SPIN state-output and input-state maps, with `C_m o A_m=0` |
| `Q_i` | Structured SPIN state before inner step `i` |
| `alpha_i` | independent nonzero multiplier in `F_{2^{s_m}}` for inner step `i` |
| `tau` | interface class passed from the outer analysis to routing and inner analysis |
| `A_tau^out` | expected number of nonzero messages whose outer output has class `tau` |
| `Q_tau(d)` | conditional probability of final output weight at most `d` from class `tau` |
| `A_w^out` | expected outer weight spectrum in the weight-indexed specialization |
| `p_w^in(d)` | uniform-slice inner tail in the uniform-interleaver specialization |
| `Z_d` | number of nonzero messages producing final weight at most `d` |

Use `:=` for definitions and deterministic assignments. Use sampling notation
only after the sample space or distribution has been defined.

## Single-interleaver convention

For Structured SPIN, write the implementation factorization as

\[
\Pi_N=\Pi_{\rm reg}\circ T\circ\Pi_{\rm blk}.
\]

The three factors are stages of one interleaver. Use “block shuffle,”
“transpose,” and “region shuffle” for the factors. Do not describe Structured
SPIN as a multiple-interleaver code.

## Claim language

Use these verbs consistently:

- **define** for a new object or notation;
- **sample** for a randomized setup operation;
- **fix** for an arbitrary realization held constant;
- **prove** only for a complete mathematical argument;
- **certify** only for an independently checkable exact or outward ledger;
- **measure** for performance and empirical correctness evidence;
- **model** for an unproved analytic approximation or surrogate.

It is valid to state that the current certified scalable Structured SPIN
ensemble has relative distance `0.11` with high probability. Attribute that
statement to the one-sampled Golay--BA-3/RM2Sub-S19 theorem. Do not transfer the
claim to another base constituent or to the frozen ParityFanout member without
the corresponding certificate.
