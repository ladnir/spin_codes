# SPIN paper outline

## Scope of this restart

This outline starts from the intended SPIN families rather than the current
compiled manuscript. Existing proofs remain source material. They do not
determine the paper's construction spine.

The paper studies families indexed by admissible lengths. Accumulator SPIN and
Random SPIN are asymptotic families. Structured SPIN is also intended to be a
family, naturally indexed by powers of two. The frozen `N=2^21` construction is
one reference member, not the definition of Structured SPIN.

## One construction interface

For each admissible length `N`, a SPIN encoder has the form

\[
E_N = I_N\circ \Pi_N\circ O_N.
\]

Here `O_N` is the outer encoder, `Pi_N` is one permutation of the complete
outer output, and `I_N` is a recursive inner encoder. An implementation may
factor `Pi_N` into several routing stages. Those stages still define one
interleaver.

Let `B` be the outer block length, `K_B` its message dimension, and
`M=N/B`. A random-block outer independently samples a `K_B`-to-`B` injection
for each diagonal position.

The paper uses two outer choices:

1. **random-block outer:** independently sample one uniform `K_B`-to-`B`
   injection in each of the `M` positions;
2. **known-spectrum outer:** use a structured local rule whose fixed or
   ensemble-expected spectrum is analytically controlled.

The paper uses two interleaver choices:

1. **uniform interleaver:** sample one uniform permutation of `[N]`;
2. **structured interleaver:** sample the declared factored SPIN routing rule.

The paper uses three inner choices:

1. **accumulator:** the rate-one prefix-sum recursion;
2. **random recursive inner:** a declared random convolutional recursion;
3. **structured inner:** the recursive local transform used by Structured
   SPIN.

Only three combinations receive paper-facing names:

| Family | Outer | Interleaver | Inner | Primary result |
| --- | --- | --- | --- | --- |
| Accumulator SPIN | random blocks | uniform | accumulator | asymptotic linear distance |
| Random SPIN | random blocks | uniform | random recursive | stronger asymptotic distance |
| Structured SPIN | known-spectrum | structured | structured | rate `1/2`, relative distance `0.11`, linear time |

Other combinations may appear as proof bridges or experiments. They are not
named families unless they later support a distinct theorem.

## Common first-moment framework

The framework should begin with the bad-codeword count, not with separate
outer and inner assumptions.

Fix a target distance `d`. For one sampled encoder `E_N`, define

\[
Z_d:=\bigl|\{x\ne 0:\operatorname{wt}(E_N(x))\le d\}\bigr|.
\]

Then

\[
\Pr[Z_d>0]\le \mathbb E[Z_d].
\]

Choose an interface class `tau` that retains exactly the information needed
after the outer stage. Let

\[
A_\tau^{\rm out}
:=\mathbb E[\#\{x\ne0:O_N(x)\text{ has class }\tau\}],
\]

and let `Q_tau(d)` be the conditional probability that the interleaver and
inner encoder map a fixed outer object of class `tau` to output weight at most
`d`. When the construction distribution makes `Q_tau(d)` depend only on
`tau`, linearity of expectation gives

\[
\mathbb E[Z_d]=\sum_\tau A_\tau^{\rm out}Q_\tau(d).
\]

For a uniform interleaver, the class can often be the outer Hamming weight
`w`. The existing slice-to-tail framework is then recovered with

\[
\mathbb E[Z_d]=\sum_w A_w^{\rm out}p_w^{\rm in}(d).
\]

For the structured interleaver, `tau` may need to record block occupation or a
small routing profile. The paper must derive that route law. It must not assume
that a structured interleaver is uniform on a Hamming slice.

This exact finite identity is the common proof interface. Asymptotic envelopes
and computer-assisted bounds are later ways to bound the same sum.

## Proposed paper spine

### 1. Introduction

State the design question: how much randomness and structure can be removed
from a serially concatenated code while retaining linear distance and fast
encoding?

Introduce the three families as a progression:

1. Accumulator SPIN isolates the first-moment mechanism.
2. Random SPIN improves the distance through a stronger recursive inner.
3. Structured SPIN replaces the random components by explicit fast
   components.

State separately what is already established and what remains a target. Do
not place the frozen `N=2^21` parameters in the definition of the family.

### 2. SPIN framework

Define the encoder interface, setup randomness, rate, and admissible lengths.
Then state the class-indexed first-moment identity.

Specialize the identity to a uniform interleaver and outer weight. Retain the
existing uniform-slice lemma and the exact weight-sum formula. Introduce only
the envelopes required by the next construction.

Main text:

- encoder and setup interfaces;
- bad-codeword count;
- exact class-indexed first moment;
- uniform-slice specialization.

Appendix:

- optional generating-function forms;
- alternative envelopes;
- routine rounding and truncation lemmas.

### 3. Accumulator SPIN

Define the random-block outer ensemble by sampling the local injections
independently. State the conditioning on injectivity explicitly.

Define the accumulator and give its exact input-output enumerator. Use the
uniform interleaver to obtain the slice-to-tail probability. Combine the
independent outer generating function with the accumulator tail.

The main result should be an asymptotic positive relative-distance theorem.
Finite numerical examples are optional and secondary.

Material available now:

- accumulator definition;
- exact accumulator enumerator;
- accumulator tail calculations;
- the general first-moment reduction.

Completed result:

- rate-half relative distance `0.0037` with `B>=3 log_2 N`;
- an outward interval check of the strict parameter margin;
- direct encoding cost `Theta(N log N)`.

### 4. Random SPIN

Keep the random-block outer and uniform interleaver from Accumulator SPIN. Replace
only the accumulator by the random recursive inner.

Define the inner setup distribution and recurrence before stating its tail
bound. Then prove the stronger asymptotic distance theorem and state the
encoding cost under the selected memory schedule.

Material available now:

- a detailed analysis of a random dense recursive inner;
- episode, termination, and large-deviation bounds;
- an existing asymptotic integration argument using a different outer.

Completed result:

- an all-random-tap recursive-inner definition and exact transfer matrix;
- rate-half relative distance `0.11002` with `B>=17 log_2 N` and
  `m=ceil((51/50)log_2 N)`;
- outward interval checks of the sparse- and linear-range margins.

### 5. Structured SPIN

Define direct Structured SPIN members at admissible native lengths. A selector
`theta_m` chooses the local maps, routing parameters, fixed inner maps, and
state-field basis. The scalable outer also selects an even base length `a` and
a certified rate-half base constituent `H_a`. It tiles `H_a` to an outer-block
length `b_m` divisible by `a`, then applies two permutation--accumulator stages.
The current certificate takes `a=24`, `L_m=128m`, and
`b_m>=(39/4) log_2(L_m b_m)`. Zero extension from the largest preceding native
member supplies every sufficiently large requested length. The `N=2^21`
member is introduced only after the scalable theorem and is not the family
definition.

#### 5.1 Known-spectrum outer

Define the block-diagonal map from realized local injections `G_{m,j}`. State
whether the local maps are fixed, share one sampled constituent, or use
independent setup. The manuscript records the distinct enumerator formulas
for repeated, independently sampled, and shared-sampled constituents.

ParityFanout belongs here. For an input word of weight `w`, let

\[
a=|S\cap\operatorname{supp}(x)|,
\qquad
b=|T\cap\operatorname{supp}(x)|.
\]

If `a` is even, ParityFanout preserves the weight. If `a` is odd, it changes
the weight to `w+|T|-2b`. Uniform disjoint-set sampling gives an exact
transition law from the base spectrum to the expected transformed spectrum.

The current `31x33` receipt uses this exact transformation with a modeled base
spectrum. The transform is analytic; the present instantiated spectrum is not
yet a certificate for the actual base constituent.

For the scalable member, parameterize the base constituent as a rate-half map
`H_a:F_2^(a/2)->F_2^a` with a certified spectrum. Tile it to any outer-block
length divisible by `a`, then apply two independently sampled
permutation--accumulator stages. Sample this BA-3 constituent once and reuse it
in all `L_m` diagonal positions. Keep the expectation over its two BA
permutations outside the `q`th power of its realized enumerator. The current
certificate instantiates this interface with extended Golay `[24,12,8]`.

#### 5.2 Structured interleaver

Index coordinates by `(outer block, local coordinate)`. Define the independent
block shuffles, deterministic transpose, and independent region shuffles as
factors of one permutation `Pi_m`. The construction now states the exact
region sizes and conditional support symmetries. For the current certified
member, the certificate reduces the routed profile to the active-position
fraction and the mean relative constituent weight, at subexponential cost.

#### 5.3 Structured inner

The manuscript now states the forward recurrence independently of the
optimized transposed implementation:

\[
Y_i=X_i+A(Q_i),\qquad Q_{i+1}=\alpha_iQ_i+C(X_i).
\]

It defines the field-valued state, independent nonzero multipliers, fixed
maps, zero initial state, discarded final state, and proves inner
bijectivity. The scalable certificate supplies the three-state RM2Sub-S19
transfer bound consumed by the first-moment sum.

#### 5.4 Complete encoder

Define `E_m=I_m o Pi_m o O_m`. Deduce injectivity and rate `k_m/b_m` from
the injective local maps, permutation route, and bijective inner.

#### 5.5 Scalable constituent family

Define the parameterized BA-3 constituent interface, the
block-shuffle--transpose--region-shuffle route, and fixed RM2Sub-S19 parameters
`(t,s)=(128,19)`. Then instantiate the current certificate with a one-sampled
Golay--BA-3 constituent and state its theorem at the native lengths:

\[
\Pr[d_{\min}(E_m)\le \lfloor0.11N_m\rfloor]=o(1),
\]

with rate `1/2` and `O(N_m)` ordinary and transposed encoding work. The block
constant `39/4` is the current safe certificate value; changing the base
constituent or tightening its analysis changes the schedule, not the route or
inner architecture. Then state the requested-length wrapper with `o(1)` rate
and distance loss.

#### 5.6 Finite-length specialization

At the requested target `N=2^21`, retain the BA-3 outer and RM2Sub-S19 inner
while leaving the base constituent and its length under selection. Bind the
selected construction to a new manifest. Report correctness, distance, and
performance as separate forms of evidence.

### 6. Finite-length certificates

Give one certificate format shared by all SPIN variants:

\[
\sum_\tau \widehat A_\tau^{\rm out}\widehat Q_\tau(d)<1,
\]

where every hatted quantity is an exact value or a rigorous outward bound.
The certificate records parameters, distributions, class partitions,
arithmetic semantics, and artifact hashes.

The main text should show the partition and final bound. Detailed ledgers and
verification procedures belong in the appendix or artifact bundle.

### 7. Asymptotic scaling and complexity

Treat the three families separately.

- Accumulator SPIN: prove a small positive relative distance.
- Random SPIN: prove a larger relative distance.
- Structured SPIN: define the parameterized BA-3 constituent interface, state
  the currently certified Golay--BA-3/RM2Sub-S19 theorem on its native
  schedule, then apply the requested-length wrapper.

Define the cost model before claiming linear time. For Structured SPIN, write

\[
T(N)=\frac{N}{B}C_{\rm out}(B)
    +\frac{N}{t}C_{\rm in}(t,s)
    +O(N).
\]

Every scalable member with a fixed linear-time base constituent achieves
linear ordinary and transposed work.
The modular interface remains important because later constituents can improve
finite-length behavior or the logarithmic block constant without changing the
route or inner definitions.

### 8. Implementation and performance

Describe the optimized transposed evaluation after the mathematical
construction. Explain how each implementation stage realizes the forward
map. Report correctness checks, operation counts, benchmark conditions, and
the frozen checksum.

Performance evidence does not support a distance theorem. Distance receipts
do not establish implementation correctness or runtime.

## Drafting status and next order

The first manuscript wave filled these parts:

1. the SPIN interface and three-family overview;
2. the exact finite first-moment identity;
3. the accumulator definition and exact enumerator;
4. the random recursive inner definition and reusable local lemmas;
5. the Structured SPIN modular definition and current certified scalable theorem;
6. the ParityFanout transition law;
7. the certificate format and evidence-status convention;
8. the implementation and benchmark description.

The second wave added complete random-block theorems for Accumulator SPIN and
Random SPIN. The paper now adopts independent uniform injections as the
canonical random-block outer. Repeated fixed constituents remain related work.
The third wave integrated the one-sampled Golay--BA-3/RM2Sub-S19 certificate.
The fourth wave expanded that result into a focused proof appendix covering
the exact outer spectrum, one-sample selection, route domination, finite and
continuum transfer bounds, all occupation regimes, and theorem completion.
The exhaustive rational witnesses remain in a hash-bound artifact bundle.
The frozen ParityFanout member remains a separate finite-length optimization;
its modeled spectrum and floating ledger remain diagnostic.
