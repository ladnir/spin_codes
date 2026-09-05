# SPIN codes: next-work roadmap

The frozen Structured SPIN implementation is the main code. Its legacy artifact
name is ParityFanout-31x33 S19. New design work must use a distinct SPIN variant
name.

## Track A: package and clean the frozen implementation

1. Start from the verified source snapshot in the frozen construction folder.
2. Preserve the snapshot and create a separate cleanup working copy.
3. Reproduce the independent correctness checks and frozen checksum locally
   and on Peach.
4. Separate construction parameters, generated schedules, implementation,
   correctness tests, and benchmarks.
5. Replace the patch scripts as the canonical build path, but preserve them as
   provenance.
6. Refactor in small behavior-preserving steps. Re-run correctness after every
   step and the serialized benchmark after every hot-path change.
7. Finish with a source manifest, build command, test command, benchmark
   command, operation ledger, and performance receipt.

No proof or algorithm change belongs in this track.

## Track B: establish the paper spine

The paper should define one serial-concatenation interface and instantiate it
three times.

1. **Framework.** Define the outer spectrum, interleaver, inner transfer law,
   and the exact finite-length first-moment theorem.
2. **Accumulator SPIN.** Use a random block outer, a global random
   permutation, and one accumulator. Prove linear distance with a small
   constant. The exact accumulator enumerator makes the proof transparent.
3. **Random SPIN.** Keep the random block outer and random
   permutation, but replace the accumulator with a random convolution. Prove a
   substantially larger distance and explain its computational cost.
4. **Structured SPIN.** Replace the outer, permutations, and inner in
   that order. State precisely which spectrum and transfer properties each
   replacement preserves. End with the frozen code, finite certificate, and
   benchmark.
5. **Scaling and complexity.** Separate finite-length validity, asymptotic
   distance, failure-margin growth, and encoder complexity.

The current `framework.tex`, `innerAcc.tex`, and random dense-construction
files already provide the first three sections' proof skeletons. The final
section should be rewritten around the frozen construction rather than the old
RM/EBCH full-split candidate.

## Track C: support every finite length

First prove the framework for arbitrary finite parameters. State construction
theorems for admissible lengths that satisfy the block and step divisibility
conditions. Then add one rounding lemma that handles all other lengths by
padding, shortening, or two adjacent block sizes. Do not duplicate boundary
logic in every construction proof.

The finite theorem should return an explicit bound

\[
  \sum_h A_h^{\mathrm{out}}p_h^{\mathrm{in}}(d),
\]

not promise one fixed distance or margin for every small length. A parameter
selector and certificate generator can evaluate that bound for any requested
finite length.

## Track D: prove asymptotic families

Treat the three constructions separately.

- The accumulator construction should yield a small positive relative
  distance when the outer block size is proportional to `log n`.
- The random-convolution construction should yield a larger relative distance.
- The structured construction should first be stated conditionally on uniform
  outer-spectrum and inner-transfer hypotheses. Concrete constituent families
  can then discharge those hypotheses.

With logarithmic outer blocks, the worst messages activate one outer block and
contain only `Theta(log n)` shuffled positions. Their late-placement event has
only `2^{-Theta(log n)}` suppression. Thus the natural global margin is often
`Theta(log n)`, even when linear-weight message classes have an
`Theta(n)` margin. This distinction should be explicit.

## Track E: investigate a linear-time family

Fix the cost model before changing constituents. For outer block size `B`,
inner step size `t`, and state size `s`, write the cost as

\[
  \frac{n}{B}C_{\mathrm{out}}(B)
  + \frac{n}{t}C_{\mathrm{in}}(t,s)
  + O(n)
\]

for routing and permutations.

A dense random `B/2`-to-`B` outer costs `Theta(B^2)` per block. With
`B=Theta(log n)`, it gives `Theta(n log n)` work. A block family with a
linear-size encoder, `C_out(B)=O(B)`, removes that factor.

This is necessary but may not be sufficient. If the proof requires
`s,t=Theta(log n)`, the inner also needs `C_in(t,s)=O(t)`. Otherwise a
quasilinear local transform produces `O(n log log n)` or worse. The first
analysis task is therefore to determine whether the frozen inner can retain
constant `s,t` as `B` grows. If it can, only the outer must change. If it
cannot, both constituents need linear-size local circuits.

The first candidate outer for **Linear SPIN** should be a linear-time BA, accumulate, or
systematic expander construction with a provable spectrum envelope. Any such
replacement is a new variant; it does not alter the frozen main code.

## Dependency order

Tracks A and B can start independently. Track C should stabilize the theorem
interfaces before the final paper proofs are written. Track D uses those
interfaces. Track E can begin with a cost audit, but construction changes
should wait until Track D identifies the weakest sufficient constituent
properties.
