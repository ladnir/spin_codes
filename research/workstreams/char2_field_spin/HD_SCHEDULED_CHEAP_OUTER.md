# Cheap column mixing with fixed wiring

Status: checked native encoding at 3.442226 ms, selected sparse FLOAT
bounds, and a directed high-total-arity bound. There is no whole-code
certificate, library promotion, or paper change. The certified
3.798160-ms / 55.0927385-bit fallback and its pinned manifest are unchanged.

## Question and construction

The measured H* encoder takes 3.37 ms at K20. Adding one crossed mixing
stage takes 3.85 ms; adding two takes 4.35 ms. Those costs motivate
replacing arithmetic and changing local wiring before adding stages.
The measurements and independent native checks remain in
[the original cost record](HD_OUTER_PAIR_NATIVE.md).

The new candidate keeps sixteen RS[16,8] rows in each outer superblock.
They encode 128 input field symbols into a 16-by-16 array over GF(256).
It also keeps the global paired route, inner, and independent nonzero
diagonals before and after local mixing. It adds no crossed layer.

For column j, first apply its independent input diagonal D1. Send source
row i to P_j(i), apply H_lambda, then send output row r to Q_j(r). Apply
the independent output diagonal D2 at the resulting physical positions.
The fixed permutations are

    P_j(i) = inverse_GF16(i xor j), with inverse_GF16(0) = 0;
    Q_j(r) = r xor j.

GF16 labels use the polynomial x^4+x+1. These permutations act only on
the sixteen row labels inside a superblock. They do not change which
superblocks are active or add a global permutation.

Let b_k be the length-sixteen indicator of labels whose bit k is one.
The mixer is

    H_lambda = I + lambda * 1*1^T + sum_{k=0}^3 b_k*b_k^T.

Here lambda is AES-field byte 0xbc. It satisfies lambda^2=lambda+1, so
all coefficients lie in the embedded GF(4)={0,1,0xbc,0xbd}. Multiplication
and code symbols still use GF(256). The local circuit uses 57 XORs and
one fixed field multiplication per column, excluding both diagonals and
the wiring. The previous H* circuit used 66 XORs and eight field products.
These primitive counts are not full-encoder timing predictions.

All entries of H_lambda are lambda or lambda+1. Thus a singleton input
activates all sixteen outputs. With V=[1,b_0,b_1,b_2,b_3], V^T V=0 and
H_lambda-I=V diag(lambda,1,1,1,1) V^T. Consequently H_lambda is an
involution. Its quiet subspace is ker(V^T), whose minimum weight is four.
The varied column wiring changes which original input patterns meet this
subspace, and where their outputs occur in the inner's chronology.

## What changed in the proof

Every test below fixes the original set A of active RS rows across all
sixteen columns. It does not resample A per column or average schedules.
Final D2 supplies independent nonzero amplitudes conditional on support.
Each physical region still receives a subset of the same active group
set G. Thus the existing regional interface for Fmax(G)<=1 remains valid.

The selected source class has r=18 active superblocks and Fmax<=1. The
target remains K=2^20, rate 1/2, and forbidden output weight 734003.
The common regional measure is taken from the pinned mark-9927/10000
bank. Its small FLOAT scale lift is retained; it is not a directed bound.

The no-cross construction permits independent column factors. For fixed
boundary states, let J_b be the r-th moment of branch b under the common
positive matrix measure. Let n_b be its exact raw input count and w_b its
original input weight. Form the polynomial

    P(X) = sum_b n_b * J_b^(1/r) * X^(w_b)

entrywise. Raise each entry to the integer power r, then multiply the
sixteen polynomial matrices in physical chronological order. Keep only
degrees at least 9*a*r, where a=|A|. Coefficientwise Holder bounds the
expanded moments. This floor relaxes the separate weight-nine constraint
on each original RS row and each group replica.

Previously, two independent columns were multiplied before this Holder
step. Separating them is coefficientwise no worse: generalized Holder
bounds the sum over their intermediate state before the product of roots.
The new calculation therefore removes a proof relaxation without changing
the no-cross encoder. Independent small tests check this inequality and
noncommuting state chronology.

After the complete word tail, subtract 8*a*log2(255) from each r-th-root
logarithm. Sum the roots over all original A, then raise to r. The source
count and projective factor are included once through the common prefix.
Bounds for different arities may be combined at the root level only when
the mixer, schedule, diagonals, route, and regional measure are identical.

GF(4) relation kernels give exact column counts without enumerating all
GF(256) inputs. For x in GF(256)^a, the kernel of v -> sum_i v_i*x_i over
GF(4) determines every input and output zero. A kernel of codimension d
has exactly product_{h<d}(256-4^h) inputs. There are 44 kernel categories
at a=3 and 529 at a=4. Their ordering stays fixed, but a word may have a
different relation kernel in every column.

## Recorded comparisons

All entries are log2 upper bounds for the stated selected source class.
A positive entry means the bound is inadequate, not that a low-distance
codeword was found. None is a whole-code margin.

| Mixer and wiring | Added crossed layer | a=1 | a=2 | a=3 | a=4 |
|---|---|---:|---:|---:|---:|
| Hbin, input XOR | One symmetric layer | -2745.55 | -1390.73 | -2943.78 | Not run |
| H_lambda, input inverse | None, paired-column bound | -3956.52 | +70.12 | Not run | Not run |
| H_lambda, input inverse | None, single-column bound | -3982.74 | +69.27 | Not run | Not run |
| H_lambda, input inverse and output XOR | None, single-column bound | -3982.74 | -165.56 | -1934.58 | -574.84 |

The first row covers all 560 labelled triples, but is a different
construction from the final row. Its bound cannot be mixed with the
no-cross construction. The last row covers all 16 singleton sets, all
120 pairs, all 560 triples, and all 1,820 quadruples. The largest
four-row contribution is A={4,5,12,13}. These include the first quiet
supports; the complete four-row run took approximately 262 seconds.
Combining all four arities at the root level gives -165.5611964 bits for
this selected class, including heterogeneous groups with arities in1..4.

The subsequent five-row calculation covers all 4,368 original row sets.
Its log2 bound is -168.59973724857355. The marked-fiber bound suffices
except for 120 sets, where complete exact GF4 tails provide the refinement.
No requested refinement was skipped. Combining these roots with the
completed arities one through four gives -149.03605021531962.
This covers heterogeneous groups whose individual arities lie in1..5,
still only at r18/Fmax<=1 under the same common measure. It does not
cover every tuple with total arity at most90.

The all-five receipt is
`theory/hd_lambda_arity5_cover_complete_9927_float.json`, SHA256
`916f64d08573221db5e929155125ddd0e73a55522167cb06959f799b488d5ec7`.
Its coverage check finds all4,368 distinct masks, no invalid records,
and unchanged source and regional-artifact pins. Two360-second runs
and a16.50-second continuation completed the calculation. Intermediate
receipts are retained. The continuation reorders labelled sets for cache
reuse; it neither identifies translated sets nor multiplies an orbit
representative. All these arithmetic results remain FLOAT bounds.

Evidence, kept without overwriting earlier failures:

- `theory/hd_pair_scheduled_Hbin_9927_float.json`: all five prescribed
  low-arity Hbin schedules with one symmetric crossed layer.
- `theory/hd_pair_relation_a3_input_xor_9927_float.json`: all Hbin triples
  for that same input-XOR, one-crossed-layer construction.
- `theory/hd_pair_scheduled_lambda2_9927_float.json`: first H_lambda
  low-arity crossed-layer controls. Their lambda is AES byte 2.
- `theory/hd_scheduled_no_pair_9927_float.json`: no-cross paired-column
  controls. Lambda is still byte 2 in this receipt.
- `theory/hd_scheduled_column_9927_float.json`: no-cross column controls,
  explicitly using lambda=0xbc.
- `theory/hd_scheduled_column_extra_9927_float.json`: output-XOR control
  and its input-XOR alternative, also with lambda=0xbc.
- `theory/hd_lambda_column_a3_inverse_output_9927_float.json`: complete
  560-triple calculation for the final candidate, about 38.4 seconds.
- `theory/hd_lambda_column_a4_inverse_output_9927_float.json`: complete
  1,820-quadruple calculation for the same candidate.

The low-arity support counts for lambda=2 and lambda=0xbc are identical,
as checked by exact branch tests. This limited identity does not transfer
higher-arity results or identify the two full setup distributions.

The main new controllers are `hd_scheduled_column_control.py` and
`hd_lambda_column_control.py`. The latter can also keep weighted
consequences of the individual RS floors. For fixed positive row scores
c_i, replace branch degree by sum_i c_i*1[x_i!=0] and use floor
9*r*sum_i c_i. Scores must remain fixed across columns and replicas.
Take minima only after complete fixed-A bounds, never inside a state or
column sum. This option is implemented and unit-tested but was not needed
for the reported output-XOR results.

## Dense-input diagnostic

The same candidate has a useful dual description. Write the bare column
map as M_j=Q_j H_lambda P_j. It is orthogonal, so the dual of M_j F^A is
M_j F^(A-complement), with the same schedule and original complement.
For a=|A| and b=16-a, the exact raw column enumerator is

    K_A(t) = 256^(-b) * sum_{y in M_j F^(A-complement)}
             product_{y_i=0}(1+255*t_i) * product_{y_i!=0}(1-t_i).

This formula could handle a>=12 using the existing complement categories
of size at most four. First apply the regional majorant to the original
primal support polynomial, which has nonnegative coefficients. Only then
apply this identity under the common positive matrix measure. Both dual
factors are nonnegative at its atoms. Applying separate moment bounds to
the alternating expansion of (1-t)^r would be invalid.

The reusable dual kernel and a bounded controller were implemented, but
the proposed dense-source run was deliberately not launched. A positive
path in the existing common measure already rules out success for that
calculation. Choose the 0-to-1 transition at the first region, followed by
255 live-to-live transitions. Their atoms are t01 followed by 255 ones.
All entries of H_lambda are nonzero, so the first column enumerator on
this path is 256^(a-1)*(1+255*t01). Each other column contributes 256^a.

Let P be the saved common prefix and let e01,e11 be the logarithms of
these transition masses. This path alone gives the following lower floor
on the positive, RS-floor-dropped majorant, including all labelled A:

    P + r*log2(binomial(16,a))
      + r*((16*a-1)*8 + log2(1+255*t01) - 8*a*log2(255))
      + e01 + 255*e11.

For a=12,13,14,15,16, its approximate values are +332.20, +1454.41,
+2567.22, +3667.71, and +4748.52 bits. They are lower bounds on this
majorant, not lower bounds on code failure or evidence of bad codewords.
The saved arithmetic is `theory/hd_lambda_dual_common_measure_floor_float.json`.

Restoring the RS support floors cannot remove this obstruction from the
same ambient counting comparison. Restrict its inputs to all 16*a
original coordinates nonzero. Every row then has weight sixteen, and the
first output factor is at least t01. The resulting floor is

    P + r*log2(binomial(16,a))
      + r*(8*a*log2(255) + log2(t01)) + e01 + 255*e11.

It remains above +312 bits already at a=12. A different regional mark or
measure, or the existing uniform-amplitude comparison based on unique
preimages, is needed for these dense classes. Additional local mixing is
not a remedy for this specific comparison loss.

`theory/hd_lambda_dual_kernel.py` has seven tiny tests, including exact
GF(4) MacWilliams comparisons, noncommuting matrix measures, endpoint
atoms, and correct raw-mass normalization. Its controller drops original
input floors explicitly. Neither file supplies dense-class coverage here.

## Directed high-total-arity class

The invertibility comparison now covers a different selected class:
exactly 18 active superblocks with at least 247 active RS rows in total.
There are at most 288 such rows. The calculation includes every indexed
source support and all 99 even/odd pairing profiles, not only Fmax<=1.
It does not cover the complementary lower-total-arity class or other
numbers of active superblocks.

For this comparison, every fixed column map is invertible. Compare its
output with independent uniform field symbols, including zero. The
unique-preimage density factor is 256^(256*r)/255^(8*A0), with r=18 and
A0=247. Source-support counts remain outside the fractional power, and
the projective factor is applied once. No row-set or cohort multiplier
is added. The proof uses full independent nonzero D1 and D2.

`theory/hd_lambda_high_arity_certificate.py` runs fresh regional arithmetic
at mark 691/700 and alpha 1/5. It reuses reviewed arithmetic functions,
not another construction's saved certificate. Fresh 256- and 384-bit
runs agree on all 99 rational profile ceilings. A separate 512-bit sum
of those ceilings gives failure at most

    2^(-63942051987187/1000000000000),

or a 63.942051987187-bit margin for this class. The source-pinned receipts
and verification are in `theory/hd_lambda_high_arity_pinned/`. Earlier
receipts remain available but predate one added dependency pin. Use the
pinned directory for the current adapter, whose SHA256 is
`156ca16d6f18370eb2484e6a1bf74ef83cbede670b63b2ba03db27f465d2e036`.

This uses a different regional comparison from the sparse bounds. The
results may be combined only through a justified partition of source
classes, not by mixing their per-group root bounds.

### Testing the dense gap

A separate bounded control asks whether the uniform comparison can start
at total arity216, which includes eighteen groups of arity12. It rebuilds
the regional arithmetic at three prescribed marks, still with alpha1/5,
and sums all99 profiles for each mark. Fresh256- and384-bit runs agree
on every profile ceiling and every integer cutoff ceiling from216 to288.

| Mark | Log2 bound at cutoff216 | First cutoff reaching60 bits |
|---|---:|---:|
| 691/700 |+332.577878481012|247|
| 172/175 |+1084.25300234999|None through288|
| 171/175 |+2317.282456539453|None through288|

For a fixed mark, increasing the cutoff from216 to A0 improves the log2
bound by `(8/5)*(A0-216)*log2(255)`. These cutoffs describe the tested
witnesses, not a lower bound on the code's failure probability. At the
first mark, cutoff246 gives51.151086488213 bits, but still leaves the
arity12 case uncovered. Nearby tested marks do not repair that gap.
No union penalty is charged over marks: they are alternative bounds
on the same event, not disjoint events.

The controller is `hd_lambda_dense_bridge_control.py`, SHA256
`a9a844d33d3f7c74269e1dcc7206b3e02c76c0729a2e1ddc2a48153815b19dab`.
The fresh receipts are `theory/hd_lambda_dense_bridge_256.json` and
`theory/hd_lambda_dense_bridge_384.json`, with hashes
`6b93a7e7b194ca932118bdfca891efd1f6324149e99b18d580bebb33cb43597d`
and `062138e3c1e22a9c8a1c837ca1e44ab70deb5b11160028961a00f3d6dec0bba4`.
All297 profile ceilings and219 cutoff ceilings agree. The construction
binding, dependency pins, full profile inventory, and density adjustment
received independent review.

## Middle-arity fiber comparison

Fix the mapped input set A in one column, and let O be its complement.
Write C=H_lambda[O,A], c=rank(C), and a=|A|. The exact support enumerators
of im(C) on O and ker(C) on A have masses 256^c and 256^(a-c).
Call these two polynomials C_A and Z_A. Since H_lambda is an involution,
the inside outputs for each outside value form a translate of ker(C).
At coordinate markers in [0,1], its zero translate maximizes the weight
sum by nonnegative additive Fourier factors. Therefore the raw allowed
input enumerator satisfies K_A(t)<=C_A(t_O)Z_A(t_A).

`hd_lambda_fiber_outer.py` computes both labelled polynomials exactly by
zero-set ranks and integer Mobius inversion. It transports both masks
through Q_j and retains the original A across columns. The inequality
is pointwise, not coefficientwise, and initially drops RS input floors.

The first bounded control restored only nine nonempty columns per word.
It used the positive split K_A-1 <= (C_A-1)Z_A+(Z_A-1); zero coefficients
were omitted exactly, not subtracted from numerical upper bounds. The
empty branch retained its full matrix moment. This still failed:
fixed five-row examples gave log2 bounds +2857.08 and +2772.45, and
fixed eight-row examples gave +933.72 and +895.67. These are individual
fixed-A diagnostics, with no extrapolation to other row sets.
The receipt is `theory/hd_lambda_fiber_representatives_9927_float.json`.

The subsequent repair retains exact input activity on the zero-outside fiber
when H_lambda fixes that fiber pointwise. Every nonquiet branch may be
assigned all a input coordinates, since this only increases activity.
This preserves upper-tail bounds, not each coefficient separately.
Keeping original row floors this way is a proof refinement, not another
mixing stage or a change to the native encoder.

That repair now has bounded numerical tests. Reverse each column's
activity degree w to the deficit a-w, then retain total deficit at most
7*a*r during chronological composition. Every coefficient is positive,
so discarding larger deficits early is safe. This enforces the same
aggregate input-weight floor 9*a*r without computing the discarded tail.
It still relaxes the separate original-row floors.

| Fixed original A | First nonempty-column control | Repaired deficit tail | Exact GF4 deficit tail |
|---|---:|---:|---:|
| {0,1,2,3,4} |+2857.08|-604.40|-2751.62|
| {0,1,2,4,8} |+2772.45|-520.99|-2930.86|
| {0,1,2,3,4,5,6} |Not run|-1462.26|Not run|
| {0,1,2,3,4,5,6,7} |+933.72|-2155.55|Not run|
| {0,1,2,3,4,5,8,12} |+895.67|-1986.44|Not run|

These are fixed-A FLOAT diagnostics for the same selected r18/Fmax<=1
class, not all row sets. The repaired column bound also takes the minimum
of two valid bounds on the same nonquiet branch: the factored 2r Holder
bound and the direct sum over disjoint outside/inside output masks.

`hd_lambda_arity5_relations.py` supplies the exact five-row data. Its
12,278 GF4 kernel categories retain original input weights and all
determined output coordinates. One category has zero GF256 mass and is
skipped before taking logarithms. The exact raw mass is 256^5. Tiny tests
include exhaustive agreement over all 1,024 GF4 inputs for one column.

The measured diagnostic receipts are
`theory/hd_lambda_marked_middle_9927_float.json` and
`theory/hd_lambda_exact5_representatives_9927_float.json`.
The complete five-row result above uses the cheaper repaired bound and
exact counts where needed.

For later arities, full quiet identity can fail. There is a safe partial
replacement. Work in mapped coordinates before Q_j and define I_fix to contain
each i in A for which

    rank((H_lambda-I)[O union {i}, A]) = rank(C).

On ker(C), these coordinates satisfy x_i=y_i. For a nonzero quiet output
support T, assign (T intersect I_fix) union (A minus I_fix) as input support.
Its weight is a-|I_fix|+|T intersect I_fix|. This contains every possible original
input support. Keep the unique zero output at input weight zero. These
rules preserve upper-tail bounds even when the full quiet fiber is not
fixed pointwise. Physical masks must first be pulled back through Q_j;
weighted original-row scores also require P_j inverse.

The exception counts have a short affine-geometric derivation. At arity8,
only the14 hyperplanes with nonzero even-weight normal lose full identity.
At arities9,10,11, exceptions occur when the complement lies in an
eight-point affine hyperplane: there are240,840,1680 such mapped sets.
Their fixed-coordinate counts are1,2,3. At arity12 every mapped set is
exceptional:140 have an affine-plane complement and no fixed coordinates;
the other1680 retain four. These counts describe individual mapped
columns. They do not replace enumeration of the actual fixed original A
through its sixteen scheduled columns.

`hd_lambda_partial_marked_fiber_control.py` implements the partial rule
without changing the frozen full-identity helper. Exact GF4 tests check
exceptional quiet fibers, physical/original labels, and upper-tail
domination. The controller uses scalar Chernoff bounds, not the bounded
polynomial backend. Nine fixed-A tests completed in0.985seconds:

| Original A | Arity | Fixed-A log2 bound |
|---|---:|---:|
| {0,1,2,3,4,5} |6|-721.76|
| {0,1,2,4,8,15} |6|-677.42|
| {0,1,2,3,4,5,6} |7|-1395.12|
| {0,1,2,4,8,12,15} |7|-1259.44|
| {0,...,7} |8|-2150.44|
| {0,...,8} |9|-2485.26|
| {0,...,9} |10|-2112.67|
| {0,...,10} |11|-969.91|
| {0,...,11} |12|+180.60|

These are representative FLOAT controls, not all-set coverage. The
receipt is `theory/hd_lambda_partial_middle_9927_float.json`, SHA256
`c8c3beeb0fe509ce88ea57c13fee2bcc3339574a4fd76ae71d9fa4eecbaff22b`.
This coarse marker grid uses no refinement, so its arity-seven result
is weaker than the earlier optimized marker and exact-tail results.
At arity12 the best tested input mark is1; increasing it only weakens
the bound. This agrees with the dense-measure obstruction above.

Before enumerating every remaining middle-arity set, resolve that dense
gap. The current high-arity comparison starts at total arity247 and does
not cover eighteen groups of arity12, whose total arity is216. Nor does
its cutoff eliminate arities14..16 from the complementary calculation:
mixed tuples can contain those arities while their total remains below247.

## Checked native implementation and cost

`scheduled_lambda_native/` isolates the prescribed construction from the
old experiment and production code. The input and output schedules fold
into specialized loads and stores. The kernel retains explicit row
operations, final non-temporal stores, the same RS circuit, the same
global route, and the same inner. It adds no global data pass.

The normalized sampler preserves the full-independent-diagonal ensemble.
For each column, let a be the original row-zero D1 scalar. Replacing D1
by D1/a and D2 by a*D2 leaves the map unchanged. Conditional on a, every
remaining key is still independent uniform nonzero. The omitted D1 key
belongs to original row zero, not mixer row zero after P_j.

The baseline and candidate pass independent scalar-reference checks at
K=8192, 24576, and 2^20, including unaligned and aliasing cases. Additional
candidate checks cover 8192 binary basis directions, 16 random halves,
both RS cosets, and 256 diagonal-gauge cases. The natural-output build
also passes through K20. The candidate passes the small correctness suite
under ASan and UBSan, with leak detection enabled. The original baseline checksum is unchanged;
the candidate's full native-layout checksum is 17225622645250162225.

Matched serial measurements use Peach CPU2, GCC15, and precomputed setup.
Each entry is the median of 501 calls. Order is baseline/candidate,
candidate/baseline, baseline/candidate, under both shared benchmark locks.
Other users' host activity can still add noise.

| Variant | Round1 ms | Round2 ms | Round3 ms | Median of medians ms | Outer-only ms |
|---|---:|---:|---:|---:|---:|
| Unchanged H* |3.354983|3.348391|3.362177|3.354983|1.618532|
| Scheduled H_lambda, no cross |3.453258|3.442226|3.430756|3.442226|1.745299|

The candidate is about 2.6% slower than this matched baseline despite
fewer field products. Fixed wiring and generated instruction scheduling
therefore cannot be treated as free. It remains close enough to the
3.4-ms objective to prioritize the missing proof before more tuning.
Outer-only times are one subsequent 501-call run for each variant.

Native logs and scripts are retained in the isolated directory. Remote
scratch is `/tmp/hd-lambda-20261009-JkBHbA`. The configured WinSCP command
rejected its host-key argument; the transfer used the existing binary-safe
SSH method with the same pinned SHA256 host key. Native source hashes:

- `HdGroupedOuter.h`: `6f20a354c395e3bc008d64a2974d0c160dd1ff48918660d15572dccfc61efaa5`.
- `hd_grouped_forward.cpp`: `4b58e91aea91beb2ca6a9a5d6024b55d7a928d55443b0712704d498a6e9c9cc9`.
- `build.sh`: `72b2f35e79f43f5c89d67701c9ed4c2ea0683c2e9b41a159b1e1335634ca2f53`.

## Remaining work

Resolve the dense gap before a broad all-set enumeration. The next proof
should retain more information about the allowed column subspaces;
their dimension alone is not a pointwise probability reduction.
The existing positive dual formula offers a possible interface, but
the current sparse measure cannot close it at arity12.

For an arity12 column the dual dimension is four, so exact local dual
counts are available. The required correction is the weighted dual sum,
not merely a factor256^-4. At all-zero coordinate markers, both raw
enumerators equal1, which already rules out that bare factor. The
existing dual helper retains the full correction under the positive
nonzero-amplitude measure; its dense-path obstruction remains.
The successful uniform-amplitude provider has already averaged packet
values and applied fractional powers. Its final2x2 matrices do not
expose that same coordinate-product measure. A new comparison must insert
the support or dual constraint before those operations. Multiplying the
finished uniform matrices by a codimension factor is not justified.

One valid local interface is support-dependent density. Let C=M_j F^A,
let B be a parity-check matrix of C, and let n_C(S) count vectors of C
with exact support S. Independent nonzero D2 multipliers assign each
fixed target vector with support S raw counting mass

    n_C(S) / 255^|S| <= 255^(-rank(B[:,S])).

Choose pivot coordinates of B[:,S]. Each assignment of nonzero free
coordinates determines at most one assignment to the pivots, proving
the inequality. This is a raw local counting statement, not a normalized
whole-code failure bound. Its support weight must enter before the
regional provider's fractionalization and placement averaging. Fixed-A
bounds also need an explicit row-set union or norm sum; they do not
inherit the high-arity comparison's absence of a row-set factor.

The subsequent Renyi comparison below resolves the homogeneous arity12
case without changing the regional provider. The support constraint
enters through an exact likelihood moment before fractionalization.
This avoids the old positive-measure floor; it does not modify that
measure or claim its exact dual count was missing a codimension factor.

Then complete the low-total-arity comparison, extend pairing profiles,
and cover all active-group counts with directed arithmetic. At r18,
the current sparse gate covers four of99 ordered pairing profiles.
The previous certified construction's weakest saved weight was r226;
it is a useful later stress test, not a bound for this changed outer.
Exact GF4 enumeration grows rapidly beyond a=5, so do not launch an
unbounded brute-force extension. The high-total-arity comparison above
supplies only its stated class.

Retain matched serial measurements and independent scalar checks for any
further change. Lower field-product counts do not guarantee lower latency.

A read-only native review found a concrete implementation route in
`cross_pair_native/HdGroupedOuter.h`. Specialize each local column. Mixer
row r loads original row `inverse_GF16(r) xor j`; physical output row R
stores mixer row `R xor j`. This preserves the sixteen explicit physical
stores and avoids SIMD row-shuffle instructions or an extra global pass.
D1 keys must follow the original source row; D2 keys follow the physical
output row. In particular, the normalized input-row-zero gauge must not
be attached to mixer row zero after P_j. The expected tradeoff is lower
arithmetic versus larger specialized code and possible register pressure.
The checked implementation above tests this tradeoff directly.

The combined local regression suite passes 66 tests. Independent review
checked the fixed schedules, GF(4) counts, positive polynomial bounds,
normalizers, and the dense diagnostic. Subsequent native measurements and
the new directed class are recorded above; neither closes the full goal.

The subsequent focused regression passes37 tests, including resumable
coverage, exact five-row counts, partial-coordinate fibers, and the dense
bridge. Independent reviews cover the new wrappers and their scope.
No further native benchmark was needed for these proof-only changes.

## A directed twelve-row bridge through the outer density

The [Renyi derivation](HD_LAMBDA_RENYI.md) fixes a row set A across all
sixteen columns. It compares the allowed-input column distribution,
averaged over independent nonzero D2 multipliers, to uniform GF256^16.
Exact support counts give its likelihood's third moment I3. Holder's
inequality then converts an outer exponent3/10 to the existing uniform
provider's exponent1/5. All row-set choices are summed inside3/10;
indexed source sets remain outside, and projective normalization occurs
once with exponent3/10. No scalar correction is applied to the provider.

The serial exact inventory computes all1,820 four-coordinate complements
of H_lambda, using its full parity matrix including the diagonal I.
It took188.7542101seconds and found

    max I3 = 6209620547094272723822049200534401189384404922308296922792507977306624
             / 6161198947895343982399660568426761903988167460346780717372894287109375
           = approximately1.0078591195655953.

For r18, with every active superblock containing exactly twelve active
RS rows, the directed failure bound is2^-69.240061061278 at35% distance.
The row sets may differ across groups. Both fresh256/384-bit provider
replays agree on all99 profile ceilings and the union; summing rounded
profile ceilings at512bits gives the same final ceiling. This is not a
bound for every tuple of total arity216, since such a tuple may contain
a group of arity below twelve. Other r and mixed arities remain open.

Pinned sources and receipts in `theory/`:

- `hd_lambda_renyi_support.py`: c4662bb016f248faa7ffcc1972e01349bc66306c40762e12df9245fc04a3c03e.
- `hd_lambda_renyi_bridge.py`: 1076ab663e89835943a4f2fda982ac79ebfe63a9881a3a12afb56077cf607281.
- `hd_lambda_renyi_a12_inventory.json`: a5560281ebc57e0e1782aa5e40f653ac4547c9ec498d19edf16a5bbcb1ee56a9.
- `hd_lambda_renyi_a12_256.json`: d2862c36ea6cdf77699120ce0d6e29442b84b2094f5d879f64f99242000c18d4.
- `hd_lambda_renyi_a12_384.json`: 346bd05046fcf6f6ef15a56c2b762b129e7977267931cda4eb3737c791bb6d1b.

Independent review checked the probability comparison, normalizers,
full-parity rank calculation, exact Mobius inversion, and source binding.
Eleven focused tests pass. The native code and both certified fallbacks
are untouched. Next, apply a common moment order to a positive
mixed-arity coefficient sum; do not combine incompatible root bounds.

That mixed-arity step now passes. Exact third/fourth moments require134
support tables after a proved48-map local symmetry reduction, while all
2,517 labelled complements are retained. The arity12 third moments agree
exactly with the independent1,820-set sweep. The controller uses each
original row set's actual sixteen scheduled moment factors, then sums
positive polynomial coefficients over ordered arity assignments.
No whole-word symmetry or per-column row-set resampling is used.

The fourth-moment bound gives40.822846778659bits for r18, every group
of arity12..16, and total arity216..245. The earlier directed comparison
gives51.151086488213bits for all total arities246..288. Their512-bit
probability union gives40.82172502955bits. Both256/384-bit replays agree;
all99 source profiles are included. The remaining r18 class has total
arity at most245 and at least one group of arity at most11. Other r
remain open. The linked Renyi note records all receipts and exact scopes.
Sixty focused regression tests pass. This proof progress adds no native
work; candidate source hashes and the certified fallback manifest agree
with their previous values.
