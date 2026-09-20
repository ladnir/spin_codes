# Reconciliation with the preliminary first-moment framework

**Status:** superseded by the clean restart recorded in `OUTLINE.md` and
`MERGE_SUMMARY.md`. Retain this file only as a record of the earlier
reconciliation pass; do not merge it into the new paper spine.

## Role of the top-level framework

The top-level `framework.tex` is a good starting point. It already isolates the
two quantities used by a first-moment proof:

- the number of outer words in each weight class;
- the probability that routing and inner encoding map such a word below the
  target weight.

The final paper should preserve this mechanism. It need not preserve the
current notation, theorem scope, or weight-only interface.

The exact identity is finite because it counts words at a selected length.
Its main role is to prove asymptotic theorems for Accumulator SPIN and Random
SPIN after the paper supplies parameter schedules and uniform envelopes.
Finite instantiations of those two families are optional consequences.

## Content to preserve

The following parts can remain structurally unchanged:

1. The complete encoder is a serial composition of an outer code, one
   interleaver, and a length-preserving inner map.
2. A uniform interleaver sends a fixed weight-`h` word to the uniform Hamming
   slice of weight `h`.
3. A nonnegative random variable counts low-weight images of nonzero outer
   words.
4. Markov's inequality bounds the probability that this count is nonzero.
5. Outer-spectrum and inner-tail envelopes can be substituted after the exact
   reduction.

The accumulator and random-convolution sections can continue to consume the
weight-indexed specialization.

## Semantic issues to correct

### 1. Separate realized and expected spectra

The current source uses one spectrum symbol for two objects:

- the number of weight-`h` words in a realized outer code;
- the expectation of that number over an outer-code ensemble.

The exact fixed-code theorem needs the first object. The ensemble corollary
needs the second. The final notation should distinguish them.

For a realized outer code `C`, define

\[
  N_h(C)
  :=\bigl|\{c\in C:\operatorname{wt}(c)=h\}\bigr|.
\]

For a random outer code `C`, define

\[
  A_h^{\mathrm{out}}
  :=\mathbb E_C[N_h(C)].
\]

This distinction removes the apparent double expectation in the current
envelope corollary.

### 2. State the exact theorem before its envelopes

The most general exact identity is wordwise. It does not require an
interleaver symmetry or a class partition.

Fix a realized outer code `C`. For a sampled interleaver `Pi` and inner map
`In`, define

\[
  X_d(C,\Pi,\mathsf{In})
  :=
  \bigl|\{c\in C\setminus\{0\}:
  \operatorname{wt}(\mathsf{In}(\Pi(c)))\le d\}\bigr|.
\]

Define the bad event directly by

\[
  \mathsf{Bad}_d
  :=\{X_d\ge1\}.
\]

For `d >= 0`, this event includes every nonzero kernel word. If the restricted
map is injective, the event is equivalent to final minimum distance at most
`d`.

Let

\[
  q_C(c;d)
  :=
  \Pr_{\Pi,\mathsf{In}\mid C}
  [\operatorname{wt}(\mathsf{In}(\Pi(c)))\le d].
\]

Then

\[
  \Pr[\mathsf{Bad}_d]
  \le \mathbb E[X_d]
  =
  \mathbb E_C
  \left[
    \sum_{c\in C\setminus\{0\}}q_C(c;d)
  \right].
\]

This statement remains valid for a joint distribution of the outer code,
interleaver, and inner map. Dependence is represented by the conditional law
inside `q_C`.

### 3. Make the class-indexed form an envelope theorem

A class partition is useful only after the exact wordwise reduction. For each
realized code `C`, let

\[
  \tau_C:C\setminus\{0\}\to\mathcal T
\]

be a finite class map, and define the realized class count

\[
  N_\tau(C)
  :=\bigl|\{c\in C\setminus\{0\}:\tau_C(c)=\tau\}\bigr|.
\]

Suppose a deterministic number `P_tau(d)` satisfies

\[
  q_C(c;d)\le P_{\tau_C(c)}(d)
\]

for every supported outer realization `C` and every nonzero `c in C`. Define

\[
  A_\tau^{\mathrm{out}}
  :=\mathbb E_C[N_\tau(C)].
\]

The exact theorem then gives

\[
  \Pr[\mathsf{Bad}_d]
  \le
  \sum_{\tau\in\mathcal T}
  A_\tau^{\mathrm{out}}P_\tau(d).
\]

Use equality only when the conditional tail is exactly constant within each
class. In certificate statements, `P_tau(d)` is normally an upper bound, so
the class-indexed expression is an inequality.

### 4. Recover the preliminary weight formula as a corollary

Assume the following conditions:

- `Pi` is independent of `C` and `In`;
- `Pi` is uniform in `S_N`;
- `In` has the declared inner-ensemble distribution.

For `U_h` uniform on the weight-`h` Hamming slice, define

\[
  p_h^{\mathrm{in}}(d)
  :=
  \Pr_{U_h,\mathsf{In}}
  [\operatorname{wt}(\mathsf{In}(U_h))\le d].
\]

The uniform-slice lemma gives

\[
  q_C(c;d)=p_h^{\mathrm{in}}(d)
  \qquad
  \text{whenever }\operatorname{wt}(c)=h.
\]

Therefore

\[
  \mathbb E[X_d]
  =
  \sum_{h=1}^{N}
  A_h^{\mathrm{out}}p_h^{\mathrm{in}}(d).
\]

This is the main interface for Accumulator SPIN and Random SPIN.

Structured SPIN instead uses, at each native admissible length
`N_m=L_m b_m`, a class map that records the information preserved by its
factored routing distribution. For the canonical scalable selector, the
certificate uses the active-position fraction and mean relative constituent
weight; the class set and envelopes vary with `m` through `L_m` and `b_m`.

### 5. Keep outer injectivity separate

The framework begins with an outer code `C`, not an arbitrary message-to-outer
map. A construction theorem must therefore do one of the following:

- define an injective outer encoder whose image is `C`;
- condition the outer ensemble on full rank;
- add failure of outer injectivity to the construction's bad event.

This choice belongs in the random-block outer definition. It must not be left
implicit in the generic first-moment theorem.

### 6. Preserve block independence

The random-block outer independently samples one uniform injection in every
diagonal position. Let `M=N/B`, and let `W_i(z)` denote the local weight
enumerator in position `i`. Conditional on the sampled injections, the outer
has enumerator

\[
  \prod_{i=1}^M W_i(z).
\]

Averaging factors across independent positions. When the local distributions
are identical, the expected enumerator is

\[
  \bigl(\mathbb E[W_1(z)]\bigr)^M.
\]

Repeating one sampled constituent would instead produce
`E_G[W_G(z)^M]`. Chosen-Block BAA motivates that related option, but it is not
the Accumulator SPIN or Random SPIN ensemble analyzed in this paper.

### 7. Define the inner interface before its property

The current framework introduces the inner as a matrix and later describes its
randomness. The revised framework should first specify:

- domain and codomain;
- setup randomness;
- persistent state, when applicable;
- initialization and terminal-state convention;
- whether the sampled map is independent of the interleaver.

Only then should the paper define `q_C(c;d)` or `p_h^in(d)`.

## Recommended theorem sequence

The framework section should use the following order.

### Definition 1: SPIN experiment

Define the outer code or outer ensemble, interleaver distribution, inner
ensemble, and sampled final map.

### Definition 2: bad event and low-weight counter

Define `X_d` and `Bad_d`. Interpret zero outputs as kernel witnesses.

### Theorem 1: exact wordwise first moment

State the conditional-probability identity and the Markov upper bound. This is
the semantic core inherited from `framework.tex`.

### Corollary 2: class-indexed envelope

Introduce classes only when they compress the wordwise sum. This corollary is
the common interface for all SPIN variants.

### Lemma 3: uniform-slice law

State the permutation counting lemma.

### Corollary 4: weight-indexed uniform-interleaver formula

Recover the preliminary outer-spectrum and inner-tail formula.

### Corollary 5: cutoff and piecewise envelopes

State the outer cutoff event and allow different envelopes on different
weight or profile ranges.

### Definition 6 and Lemma 7: certificate interface and aggregation

Define exact/outward certificate rows and prove that a complete nonnegative
cover upper-bounds the class sum.

## Mapping from current labels

| Current item | Final role | Recommendation |
| --- | --- | --- |
| `sec:sc-framework` | Common SPIN framework | Keep or alias the section label. |
| `lem:interleaver-uniform-slice` | Uniform-slice lemma | Preserve the label if the statement is unchanged. |
| `eq:pw-def` | Uniform-slice tail | Rename notation to `p_h^in(d)` and place it after the inner interface. |
| `ass:outer-interface` | Optional asymptotic envelope | Move after the exact finite theorem. Do not make it the framework's primary object. |
| `ass:inner-interface` | Optional asymptotic envelope | Replace by construction-specific transfer theorems. |
| `thm:sc-framework` | Weight-indexed specialization | Preserve as an alias or corollary; assign a new label to the wordwise theorem. |
| `cor:framework-envelope` | Cutoff and piecewise envelope | Retain after the realized/expected spectrum notation is corrected. |

## Relationship to the three SPIN variants

| Variant | Class map | Routing law | Framework output |
| --- | --- | --- | --- |
| Accumulator SPIN family | independent random-block expected global weight `h` | uniform global interleaver | proved asymptotic theorem from the factored injection generating function and accumulator contraction |
| Random SPIN family | independent random-block expected global weight `h` | uniform global interleaver | proved asymptotic theorem from the factored injection generating function and random-recursive transfer |
| Structured SPIN family | active-position fraction and mean constituent weight on native lengths `N_m=L_m b_m` | one factored structured interleaver `Pi_m` per member | certified rate-`1/2`, distance-`0.11`, linear-time theorem for the one-sampled Golay--BA-3/RM2Sub-S19 selector |

The scalable certificate uses the active-position fraction and mean relative
constituent weight as its structured interface. Its proof is independent of
the frozen ParityFanout member, whose modeled ledger remains finite diagnostic
evidence rather than a source of asymptotic uniformity.

## Provisional status

This reconciliation is complete with respect to the current top-level
framework. The finite/asymptotic theory workstream supplies both random-SPIN
theorems and the one-sampled Golay--BA-3/RM2Sub-S19 Structured SPIN theorem.
The latter uses native admissible lengths and a neighboring-length wrapper for
requested powers of two.
