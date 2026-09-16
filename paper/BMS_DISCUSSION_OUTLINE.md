# BMS discussion: outline and remaining formalization

Working subsection: `scaling_complexity.tex`, `sec:bms-question`, under
“Finite lengths, asymptotics, and complexity.” Do not promote a resolution
claim to the introduction until the model correspondence below is checked.

## Narrative

1. State the obstruction: constant outer memory and sublinear inner memory
   cannot give an asymptotically good depth-two serial concatenation.
2. Paraphrase the question in Bazzi–Mahdian–Spielman, Section IV: is
   logarithmic outer memory enough? Cite the existing `bms09` entry.
3. Use Accumulator SPIN for the basic positive result: logarithmic outer
   blocks, one uniform global permutation, a one-bit accumulator, rate 1/2,
   relative distance greater than 0.0037 with probability 1-o(1).
4. Use Random SPIN to explain how stronger recursive mixing gives distance
   above 0.11002, with logarithmic inner state and O(N log N) direct work.
5. Use Structured SPIN for the stronger result: logarithmic outer blocks,
   one structured global permutation, fixed 128-bit steps and 19-bit state,
   relative distance above 0.11, and O(N) work. This is near the rate-1/2
   GV benchmark; do not silently turn the numerical theorem into an
   arbitrary-epsilon or exact-GV theorem.
6. Explain that constant-size field arithmetic is constant work per step.
   Constant recursive state does not mean constant setup storage, a single
   time-invariant transition rule, or constant memory for the full encoder.

## Formal bridge to check

- Match BMS's bit-input, constant-output-per-transition convention.
  Our direct description uses block encoders; a buffering implementation
  is not automatically a literal instance of their stated interface.
- Give a bit-serial outer realization with O(log N) state bits. Track
  buffering and any block-position counter. Treat local map descriptions
  as fixed code data, consistently with the automata model.
- Specify startup delay and termination. If O(log N) extra coordinates or
  fixed termination inputs are used, prove the resulting rate and distance
  statements instead of assuming that arbitrary shortening preserves them.
- The accumulator is already a one-bit, time-invariant inner. It is the
  cleanest construction on which to establish the correspondence first.
- For Structured SPIN, explain constant-size serialization of each 128-bit
  step and its fixed buffers. Check where BMS permits time-varying maps:
  their explicit discussion of this extension is in the repeat-convolute
  analysis, so do not silently treat it as the definition for every result.
  If a position counter is internalized, distinguish its O(log N) bits
  from the constant recursive data state.
- Keep the distinction between local permutations inside the outer and
  the single global permutation between the two constituents.
- State any resulting corollary with its precise model and asymptotic
  parameters. Do not attach a new efficient-decoding claim.

## Attribution and priority

Accumulator SPIN uses the truncated Block-Accumulate architecture already
credited in the paper. Check that work and later literature on the BMS
question before saying “first,” “resolves a long-standing open problem,”
or claiming novelty for the architecture itself. Providing a positive
construction in the requested regime and establishing historical priority
are separate claims.

## Sources checked

- Bazzi, Mahdian, Spielman, *The Minimum Distance of Turbo-Like Codes*,
  Sections I-A, III-A, and IV:
  https://www.cs.yale.edu/homes/spielman/Research/tc.pdf
- The paper's Accumulator SPIN, Random SPIN, and Scalable Structured SPIN
  theorems, and the structured-inner recurrence.

Next step: prove the automata correspondence for Accumulator SPIN, then
extend the discussion to the structured inner and audit priority. Keep the
introduction and contribution bullets unchanged during this step.
