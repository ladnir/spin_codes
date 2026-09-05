# Initial g=4 result

The first evaluated point fixes

\[
n=2^{20},\qquad B=1024,\qquad g=4,\qquad \sigma=15,
\]

and asks for relative distance \(0.09\). The packed-group calculation gives

\[
\lambda=195.28\text{ bits}.
\]

The comparable Riffle S-Stripe RandomStepConv point gave \(110.03\) bits.
TransposePacketShuffle therefore gains about \(85.25\) bits in this model.

The dominant term has one active outer block. It contributes about \(2^{523}\)
candidate messages, while the inner low-weight event has upper bound
\(2^{-718.28}\). Their product is about \(2^{-195.28}\). This occupation is
modeled exactly and does not invoke the active-group packing lemma.

For this canonical case, every one of the 1024 transposed rows contains one
distinguished packet among 512 packet positions. That packet is nonzero with
probability one half. Its position is independently uniform in each row.
Thus the active block cannot place all its impulses near the end of the full
inner stream. The former bulk late-start mechanism is broken into many
independently located row-level impulses.

The two-active-block term is the second largest displayed term. Its
pointwise exponent is \(-437.35\), about 242 bits below the dominant term.
All remaining displayed terms are smaller.

The packing-domination lemma is now proved in
`proof/PACKING_DOMINATION.md`. The proof compares the nonzero-packet count
before and after one packing move. It then uses the uniform row permutation
to couple the packed support inside the original support. The averaged inner
process is monotone under this inclusion.

The exploratory output Chernoff parameters were selected from a finite grid.
The certificate checker treats them as fixed inputs and performs no
optimization. It evaluates all 2048 occupations with outward rounding and
proves

\[
\log_2\mathbb E[\#\text{ bad nonzero messages}]
\le -195.283900426038041150079401119.
\]

Thus a sampled setup has a nonzero codeword of weight at most 188743 with
probability at most \(2^{-195.283900426038}\). Except with that probability,
the minimum distance is at least 188744, which exceeds 9% of the output
length.

An independent interval calculation evaluates the one- and two-block terms
at fixed Chernoff parameters. It proves the respective pointwise margins

\[
195.28390042685\quad\text{and}\quad 437.35106328497
\]

to the displayed precision. The interval calculation does not use numerical
optimization.

Artifact:

- `receipts/g4_b1024_sigma15_initial.json`.
- `receipts/g4_b1024_sigma15_sparse_interval.json`.
- `receipts/g4_b1024_sigma15_full_interval.json`.

Recommended next goal: audit the ensemble assumptions against the intended
encoder cost model, then locate the best certified distance for this same
parameter point before changing the outer constituent.
