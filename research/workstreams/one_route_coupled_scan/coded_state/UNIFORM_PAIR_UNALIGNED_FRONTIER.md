# Clean episodes without an aligned first activation

The fixed three-shear model now closes the complete h=6,a=0,c=6
no-full-reset branch at 46% distance, with failure at most 2^-78.
The probability is over the fixed setup, and the event quantifies over
all nonzero messages in this branch. This is not a whole-code certificate.
The construction and the certified uniform-feedback control are unchanged.

The [aligned branch](UNIFORM_PAIR_ALIGNED_FRONTIER.md) remains separately
proved. This continuation reuses its local factors and regional dependence
correction without changing its verifier.

## What changes in the count

Here all six activations occur after the first state. They consume the
six defective batches, leaving no internal defective batch. The five
closed early episodes are therefore clean, with lengths l_i from 1 to 8.
None is forced to begin at a region boundary.

Let d=max_i l_i and E=sum_i l_i. Then d+4<=E<=5d. If p episodes cross
a regional boundary, there are m=5+p early fragments, with 0<=p<=5.
The terminal interval contributes r=ceil((T-E+1)/16) regions. Treating
these terminal regions as completely available is conservative.

Use the previously proved single-region factor for the 5-p noncrossing
episodes and the two-region factor for the p crossing episodes. Take
coefficientwise length envelopes over 1,...,d, but keep E in the route
cap. The actual longest episode still forces at least 2^d active RS rows.

The early fragment fractions sum to (E+5)/16. For each endpoint tuple,
the eligible-row bound uses m fragments, mean (E+5)/(16m), and threshold
129-r. Its factorial-moment bound and shared-region correction are the
same as in the aligned proof. Six intervals have fewer than 2^156
endpoint tuples. Requesting 256 bits in the eligible-row cap leaves
a common route exception of at most 2^-100 after the endpoint union.

For g duplicate region copies, set k_v=m+r-128 and k=k_v-g. The count
retains actual RS zeros Z, message dimension ks-Z, and virtual deficit
Z+gs. Sum Z from 0 to (k-1)s, then active-row counts s from 2^d to the
eligible-row cap. Pay the regional sharing correction and one global
density factor two. The endpoint tuple fixes the count parameters, so
their maxima require no second count-class union.

## Complete branch integration

Recomputing the retained verifier leaves exactly 40 live lengths,
T=1893,...,1932. Its worst exponent is +3003 at T=1931. Every other
length has good-route exponent at most -6990.

The new count covers all 40 residual lengths and 5711 count cells.
Its worst good-route exponent is -139, first at T=1925. Hence

    2^-6990 + 2^-80 + 2^-139 + 2^-100 <= 2^-78.

The old verifier supplies the 46% output-weight restriction. The new
clean-episode exclusion does not need an output-weight restriction.

## The first branch with an internal defect

The neighboring h=6,a=0,c=5 branch already passes the old verifier:
there are no residual lengths, and its worst good-route exponent is
-76454. It does not require another episode refinement.

For h=6,a=1,c=6, the old verifier leaves 34 lengths, T=1908,...,1941.
Its worst exponent is +4091 at T=1940. Five later activations consume
five defects, leaving at most one defective batch inside a live interval.

First restrict this branch to clean closed early episodes. The remaining
internal defect can occur in the terminal live interval, which the new
count leaves unrestricted. There are again five clean early episodes.
Ignoring the first episode's forced alignment only enlarges their class,
so the same five-factor count applies without modification.

All 34 residual lengths in this subclass close, across 9822 count cells.
The worst new good-route exponent is -74 at T=1941. The old verifier
covers other lengths with exponent at most -6244. Thus this subclass has

    2^-6244 + 2^-80 + 2^-74 + 2^-100 <= 2^-72.

This does not close the complete h=6,a=1,c=6 branch. Its remaining case
places the internal defective batch inside one of the five closed early
episodes. The other four closed episodes remain clean. The defective
episode can have length up to 17, so its input interval can touch three
regions. The older selected length-three, two-region charge-band test
does not cover this case.

The next discriminating test should join the defect-envelope identity
to this actual geometry. Retain that births and one-sided pairs can occur
in only one incoming batch. Allowing them independently at every step
would discard information supplied by the branch's defect budget.
Account for one-, two-, and three-region placement, then test T=1940
and T=1941 before a complete integration. Failure of a relaxed envelope
would be a proof limitation, not by itself a code counterexample.

## Reproduction and scope

`diagnose_pair_unaligned_frontier.py --legacy-branch` recomputes the
complete a=0,c=6 residual list. `--integrate` covers all its residuals.
`--clean-early-branch` checks only the stated a=1,c=6 subclass.
The legacy-only options `--initial` and `--episodes` inspect neighboring
branches without changing the new five-clean-episode count.

`--test` passes 576 comparisons of the five-factor length envelope against
explicit length tuples. Local factors and the multifragment correction
retain the tests in the aligned verifier. Runs were serial, used Python
with bytecode disabled, bounded caches, and 40-second internal limits.
The integrations took approximately 8 seconds each on the local machine;
these are proof-verifier times, not encoder benchmarks. Numerical work
began with more than 22 GiB free and wrote no result files.

Full-reset histories, general defect-rich histories, and extension-field
closure remain open. These results support continuing the proof; they do
not yet identify a structural requirement warranting an inner redesign.
