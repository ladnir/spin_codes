# Remove one outer symbol randomizer without losing the bound

Exactly one outer symbol randomizer per group can be replaced by a fixed
invertible map. The complete expected outer output measure is unchanged,
so the existing first-moment distance bound transfers without loss.
The statement covers both certified constructions:

| Outer group | Parallel rows r | Symbol alphabet Q | Fixed native forward map |
|---|---:|---:|---|
| GF16 RS[16,8] | 4 | 2^16 | Q16^T at symbol 0 |
| GF16 RS[8,4] | 8 | 2^32 | Identity at symbol 0 |

Here Q16 is the fixed four-by-four coordinate transpose, with index map
c↦4(c mod4)+floor(c/4). The old native reverse kernel retains Q16 when it
omits the field multiplier. The new reverse kernel uses the identity.
Any one fixed symbol position works; it may differ deterministically by
group. Keeping an unused setup draw does not change the conclusion.

## Equality of expected measures

Fix a linear row code C0 over GF16. Form C from r parallel copies of C0,
viewed as vectors of symbols in GF16^r. Let Q=16^r. For a nonempty exact
symbol support S, write N_S for the number of words of C with that support.
The zero message is excluded throughout.

A matrix H in GL(r,GF16) mixes the rows of C. It preserves C and each exact
symbol support. For any nonzero symbols a,b in GF16^r, some such H maps a
to b. Consequently, for each j in S and every nonzero a, exactly
N_S/(Q-1) codewords have support S and value a at coordinate j.

Fix one distinguished coordinate j. At every other coordinate, sample an
independent invertible randomizer with uniform nonzero fixed-input images.
At j, use any fixed invertible binary map F. Fix an output x with exact
support S, and put h=|S|.

If j is not in S, its zero symbol is unchanged; the expected number of
messages producing x is N_S/(Q-1)^h. If j is in S, the input symbol must be
F^-1*x_j. The preceding count gives N_S/(Q-1) compatible codewords. Each
reaches the other h-1 output labels with probability (Q-1)^-(h-1). Thus in
both cases the expected multiplicity at x is

\[
\mu(x)=\frac{N_S}{(Q-1)^h}.
\]

This is exactly the expected multiplicity when all symbol randomizers are
independent and transitive. The argument requires neither an MDS row code
nor alignment between a randomizer's field basis and the GF16 row basis.
It uses row-code symmetry to count complete output vectors, not to infer
joint independence from marginals.

Group setups remain independent, so their expected measures multiply.
Apply the unchanged independent packet shuffles, regional routing, and
inner maps to either measure. The expected number of low-weight outputs
is identical. Markov's inequality therefore gives the same complete
first-moment setup-failure bound. Individual messages need not have the
same output distributions under the two ensembles.

For the old RS16/s15 construction the bound remains 62.04117326862128 bits.
For the new RS8x8/s15 construction it remains 61.8840972765688264 bits.
These are transfers of the respective complete receipts; no deterministic
seed is certified by either statement.

## Why a second omission is different

Two unrandomized labels need not be independent in the codeword count.
On a minimum MDS support, the shortened code is one-dimensional over the
extension alphabet. Fixing one nonzero symbol determines the second by a
fixed GF16 scalar. The second label is therefore not uniform independently
of the first.

A general but weak pointwise bound is still available. Fixing a second
label cannot increase the count after the first label is fixed. If both
omitted positions are active, this gives

\[
\mu_{\text{omit two}}(x)\le(Q-1)\mu_{\text{all mixed}}(x).
\]

When at most one omitted position is active, equality with the all-mixed
measure still holds. The worst-case product bound across q active groups
loses q*log2(Q-1) bits. For old 16-bit symbols, the q=1 bound would retain
about 46 bits, but the critical q=41 contribution would lose about 656 bits.
For the new 32-bit symbols even the q=1 bound would drop to about 30 bits.
These are limitations of this coarse domination bound, not low-distance
attacks. No additional omission is covered by the lossless lemma.

`test_omit_one.py` exhaustively checks a small analogue: two parallel
GF4 RS[4,2] rows, all 256 messages, and each possible omitted coordinate.
It verifies exact equality of expected pointwise measures. On a minimum
support, omitting two mixers gives a fifteen-fold pointwise inflation,
matching the coarse factor Q-1 for that analogue.

## Sharing randomness is not the same modification

A shared GF16-linear row map at every symbol commutes with the row code
and merely permutes its messages. It does not randomize the outer code set.
The same is true for shared extension-field scalars when their basis is
aligned with that GF16-linear action. This invariance is not automatic for
an arbitrary binary tower basis.

More generally, sharing an arbitrary binary map can preserve each
one-symbol marginal while introducing joint correlations. Neither the
lossless omission lemma nor the independent-symbol certificate covers
that replacement. Its joint expected measure must be derived separately.
