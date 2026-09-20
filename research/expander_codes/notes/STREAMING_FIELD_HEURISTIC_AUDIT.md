# Audit of the streaming field heuristic

## Scope

This note studies the compact field encoder in libOTe.  The encoder is a
heuristic implementation of regular Expand--Convolute.  It is not sampled from
the ensemble analyzed in the paper.

Let the code length be `n = Lq` and the message length be `k = Rq`.  Write a
message index as `(b,s)`, where `b` is in `[q]` and `s` is in `[R]`.  The
expander places one edge from each message index into each of the `L` regions.
The encoder then applies a nonzero label to each code coordinate and a
finite-memory convolution.  The current implementation regenerates the
convolution coefficients from the public seed at every position.

The labels and the convolution are invertible diagonal and triangular maps.
They therefore cannot repair a rank defect in the expander.

## Defect in the unsigned striped scatter

The first streaming implementation used one cyclic shift `a[r,s]` for each
region and slot.  Its edge from `(b,s)` in region `r` ended at

```text
b + a[r,s] mod q.
```

All edges that ended at one code coordinate shared the same nonzero label.
Fix slot constants `c[0],...,c[R-1]` whose sum is zero.  Assign `c[s]` to every
message coordinate `(b,s)`.  Each expander coordinate then receives the shared
label times the same zero sum.  Thus the message maps to the all-zero word.

This argument gives an `(R-1)`-dimensional kernel for every seed and every
region size.  The 26/13 construction consequently had nullity at least 12.
The periodic convolution does not affect this conclusion.

## Current repair

For fields of odd characteristic, the implementation retains the cyclic
shifts and adds one sign to each region-slot pair.  Region zero uses only `+1`.
For `1 <= r < R`, region `r` uses `-1` in slot `r` and `+1` in every other
slot.  The remaining regions use sampled signs.

For a base-constant message, subtracting the region-zero equation from region
`r` gives `-2 c[r] = 0`.  The first `R` region equations therefore have full
rank in odd characteristic.  This calculation removes the explicit
constant-slot kernel.  It does not prove that the complete expander has full
rank for every seed.

At sufficiently large region sizes, configuration also rejects repeated slot
differences.  This check guarantees that two distinct message rows share at
most one expander coordinate.  The deployed region size 80,659 enables the
check.  It changes only configuration and does not add work to encoding.

In characteristic two, signs cannot distinguish edges.  The implementation
instead generates one compact permutation of all `k` message indices for each
region.  Bucketing a permutation into groups of size `R` gives exact right
degree.  A region of size one is rejected because every such unweighted
expander has rank one.

## Evidence

The libOTe unit tests perform exact rank calculations on 128 reduced 26/13
instances over `Fp31`. They also test 128 characteristic-two seeds with the
full-index permutation path. A deployed-size test configures 1,024 signed
26/13 expanders with region size 80,659. Every accepted schedule has distinct
slot differences in all 78 slot pairs. The transpose, in-place, paired,
degree, and deterministic-regeneration tests also pass.

The script `scripts/streaming_ec_heuristic_audit.py` provides an independent
reduced-size model. Over `F_127`, it tested 128 seeds at region sizes 2, 3, 5,
7, 11, 13, 17, and 29. All 1,024 signed expanders had full row rank. All 1,024
permutation expanders also had full row rank. The unsigned striped expander
had nullity 12 at every tested size. These finite experiments are a rank
check, not a distance estimate.

A second rank stress test used the collision-free striped model with shared
column labels.  All 128 degree-6/3 instances over `F_3` had full row rank at
region size 13.  All 128 degree-26/13 instances had full row rank at region
size 29.  At degree 4/2 and region size ten, however, two of eight seeds had
nullity one.  Giving every edge an independently generated sign removed the
defect in 128 of 128 tests.  This separates removal of the deterministic
constant-slot kernel from a universal rank guarantee: the current signs do the
former, while the finite tests only provide evidence for the latter at the
intended degrees.

On the current benchmark machine, the following measurements use the 26/13/m4
Goldilocks profile with `k = 1,048,567` and `n = 2,097,134`.

| expander schedule | schedule | scalar encode | extension encode | receiver pair |
|---|---:|---:|---:|---:|
| unsigned stripes, defective | 9.4 KiB | 0.040 s | 0.080 s | 0.134 s |
| full-index permutations | 10.0 KiB | 0.085 s | 0.165 s | 0.209 s |
| signed stripes, 256-periodic taps | 11.7 KiB | 0.052 s | 0.103 s | 0.136 s |
| signed stripes, current on-the-fly taps | 3.7 KiB | 0.057 s | 0.101 s | 0.122 s |

The current implementation uses random coordinate signs instead of full-field
coordinate labels. It applies paired signs during convolution and avoids a
separate pass over the receiver's three field elements. These timings are
engineering measurements, not statistical bounds.

### September 2026 stress round

The pair-collision audit sampled 10,000 unfiltered offset tables at region
size 80,659. It found no projective key of multiplicity three or more. A total
of 1,637 tables contained a repeated key of multiplicity two. The C++ encoder
rejects these tables because it requires distinct slot differences. The
1,024-seed libOTe test above checks the property after this rejection step.

An exact ternary experiment enumerated every nonzero message for 1,024
`[18,9]` codes. The reference ensemble had mean distance 3.9629 and minimum
distance two. The signed streaming model had mean distance 3.9404 and minimum
distance one; eight streaming samples had distance one. Region size three
cannot satisfy the six-region difference condition. This experiment therefore
warns against using the structured construction at such small region sizes.

A meet-in-the-middle search tested the smallest collision-free 6/3 geometry:
`[42,21]` codes over `F_3` with region size seven and memory four. All eight
reference generators and all eight streaming generators had full row rank.
The mean best weights found were 9.0 for the reference and 8.625 for streaming.
One seed used 64 randomized joins; the other seven used 32 joins each. These
weights are upper bounds on the sampled distances.

A deterministic search tested eight `[754,377]` codes over `F_127`, with
degree 26/13, region size 29, and memory four. Every reference and streaming
generator had full row rank. The mean lightest row-pair weights were 684.625
for the reference and 688.5 for streaming. The mean interval-search weights
were 371.75 and 371.875. Thus this search found no systematic loss from the
streaming structure at the intended degree.

The libOTe frontend also audits the generator that the C++ encoders actually
materialize. For each pair of generator rows, the audit tries every nonzero
projective coefficient ratio and records the lightest resulting word. All
eight generators had full row rank at region sizes 29 and 79 over `F_127`.
At region size 29, the mean pair weights were 684.0 for the reference and
667.5 for streaming. At region size 79, they were 1859.0 and 1810.875.
Shared coefficients therefore cause a visible pair-weight loss in these small
instances. The found words remain much heavier than the corresponding
rate-half GV cutoff, and the search covers only messages of weight two.
Moreover, neither small region size enables the deployed difference check.

The deployed topology admits a stronger rank audit because it is
quasi-cyclic. Here `q = 80,659 = 79 * 1,021`. Over the Goldilocks field,
`X^q-1` has one linear factor, two degree-39 factors, one degree-1,020 factor,
and six degree-13,260 factors. The audit evaluated the first four factors
exactly for eight C++ schedules. Every resulting 26-by-13 matrix had rank 13.

The same eight schedules passed two exhaustive checks over all 80,659 cyclic
frequencies. A numerical complex audit found no singular value below `1e-4`;
the smallest observed value was at least 0.98. An exact audit over the
auxiliary prime 2,150,530,259 found rank 13 at every frequency. The auxiliary
prime contains a primitive 80,659-th root of unity, so this second audit can
detect characteristic-independent dependencies at every frequency.

The Python topology sampler replays the C++ sampler exactly. For the fixed
seed `(0, 2^64-1)`, both implementations produce the offset digest
`0x5374e401c8ebc7eb`. The libOTe unit test fixes this digest as a regression
value. This comparison ties the spectral experiments to the deployed
schedule generator rather than to an independent model.

The final 30-run Goldilocks benchmark measured 0.057 seconds for scalar
encoding, 0.101 seconds for degree-two extension encoding, and 0.122 seconds
for the paired receiver path. The streaming schedule occupied 3.7 KiB.

## Remaining questions

The current odd-field topology is quasi-cyclic. At a `q`-th root of unity `z`,
its expander becomes an `L` by `R` matrix whose entries have the form

```text
sign[r,s] * z^a[r,s].
```

A complete Goldilocks rank audit must test this matrix over every irreducible
factor of `X^q-1`. The deployed-size audit leaves the six degree-13,260
factors unchecked over Goldilocks. The numerical complex audit and the exact
auxiliary-prime audit are strong evidence against a structural kernel. They
do not prove full rank over those six Goldilocks extension fields.

The public coordinate signs come from one AES-round counter. They are not
secret or cryptographically random. The ideal fresh-tap calculation gives the
same zero-pattern distribution for every fixed nonzero label sequence. The
reduced experiments found no consistent distance penalty from random signs or
all-one labels. The implementation exposes neither ablation as another
encoder; its sole heuristic uses random signs. `STREAMING_LABEL_ABLATION.md`
gives the calculation, experiments, and scope limitation.

Finally, full rank is only a correctness gate.  The reduced-distance ablation
in `scripts/streaming_ec_distance_ablation.py` compares the proved regular
ensemble, signed stripes, shared labels, and periodic taps.  Rank-based
shortening is more informative than generator-row weights for this comparison.
