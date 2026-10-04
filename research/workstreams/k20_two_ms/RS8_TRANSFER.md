# RS[8,4] with 32-bit column randomization at K20

Shortening the constituent halves its parity circuit, but doubles the
column randomizer width. The complete proof transfer closes at
**51.822466 bits of margin at 10% distance**. The isolated implementation
now measures about **2.86 ms**, versus roughly 2.96 ms for the retained
control; this is not yet a promoted library profile. Independent mathematical and
source review found no blocker. The retained RS[16,8]/s22 implementation
remains the control.

## Construction

There are 8192 groups. Each group maps 128 message bits to 256 code bits.
Interpret its input as eight rows, each containing four GF(16) symbols.
Encode each row with the same fixed systematic RS[8,4] generator. The
generator must be MDS. Collect the eight symbols at each code position
into a 32-bit column, using a fixed row/bit order.

Fix any binary basis of GF(2^32). Independently for every group and each
of its eight columns, sample a nonzero field coefficient. Multiply the
column by that coefficient. Binary adjoints of these multiplication maps
give the same expected-measure argument, provided the implementation
declares its forward and transposed orientations consistently.

Split each randomized column into eight four-bit packets. Independently
shuffle all 64 packets within each group. Route one packet from every
group to each of 64 regions, using an independent uniform group
permutation in every region.

Keep the exact certified22 inner: width 64, state dimension 22, and its
independent uniform GL(22,2) state updates. The state begins at zero,
persists across every physical step and region, and is not flushed.
Outer coefficients, group shuffles, regional permutations, and inner
updates are mutually independent. Setup is fixed across all messages.

This is the same group size, packet width, region count, and inner as the
retained K20 certificate. It changes only the outer distribution. It is
not the previously studied width-eight-packet route.

## Exact outer measure

Let Q=2^32. Grouping the eight RS rows gives the column-support enumerator
of an MDS[8,4] code over an alphabet of size Q. One can prove this without
choosing a field embedding. Shortening to a specified set of h columns
leaves max(h-4,0) GF(16) dimensions in each row. The eight rows therefore
have Q^(h-4) words when h>=4. Inclusion-exclusion gives the usual MDS
counts A_h, with minimum nonzero column support five.

For every fixed nonzero column, its independent randomizer gives a
uniform nonzero 32-bit vector. The packet-support polynomial is therefore

    P(z) = (1+15z)^8 - 1,
    C(z) = sum_(h=5)^8 A_h (P(z)/(Q-1))^h.

C_v is an expected number of local messages, not a probability or the
spectrum of every sampled code. Its sum is 2^128-1. For any fixed packet
support, all nonzero nibble labels have equal joint weight in this measure.
The independent local shuffle then makes support uniform given its size.
Thus each post-shuffle vector with v active packets has expected count

    C_v / (binomial(64,v) 15^v).

Independent setups between groups permit products of these measures.
No independence between messages under one setup is assumed.

For a specified h-column support, at most (Q-1)^(h-4) codewords have exactly
that support. Indeed, shortening produces an MDS code of dimension h-4;
its words are determined by any h-4 coordinates, each required nonzero.
Dividing by (Q-1)^h for the column randomizers bounds every complete
256-bit atom by (Q-1)^-4. Hence its uniform-input envelope satisfies

    beta_new = 2^256/(2^32-1)^4
             < 2^256/(2^16-1)^8 = beta_old.

## Complete certificate transfer

The retained certificate supplies exact outward endpoints for all active
group counts q=1,...,8192, including sparse support endpoints W_v.

For q=1, use the new exact fold 8192 sum_v C_v W_v. The cached values for
v<9 equal one, which is valid but deliberately crude.

The old outer has zero mass below packet support nine. Define F_bad as
the event that any group has any nonzero local word with support below
nine. Markov's inequality and a union over groups give

    Pr[F_bad] <= Delta = 8192 sum_(v=5)^8 C_v.

For v>=9, let r=max_v C_v/B_v, where B_v is the old RS[16,8] expected
support count. The maximum occurs at v=9. Conditional uniformity of
support and joint packet labels gives pointwise domination of the
restricted complete expected measure by r times the old measure.
Consequently, regular q=2 words cost at most r^2 times the retained q2
endpoint.

This does not condition setup on F_bad's complement. On that complement,
any q2 witness uses two local words of support at least nine. Its bad-word
indicator is bounded by the original indicator times those two support
indicators. Summing under the original independent setup distribution
then gives the product of truncated expected measures. Add Delta
separately. Counting the same rare words again in the unconditional q1
fold is conservative.

For q>=3, the improved uniform-input envelope preserves every retained
ordinary and routing-conditioned endpoint. Product domination holds for
each fixed active-group set and route before applying the inner. Clipping
and fractional moments therefore keep their original direction. The
active-group-set union remains outside each fractional power. No setup
conditioning or unproved rescaling of the retained endpoints is needed.

| Contribution | Margin in bits |
|---|---:|
| Rare local support five through eight | 51.890772646 |
| q1 exact shell fold | 56.257090200 |
| Regular q2 transfer | 116.743629984 |
| Retained q>=3 union | 81.503165538 |
| Complete union | **51.822466415** |

The resulting ideal setup has minimum distance at least 209716 except
with the stated probability. This does not certify the realized distance
of a seed or a heuristic implementation of the setup distribution.

## Reproduction and next step

[rs8_transfer.py](rs8_transfer.py) authenticates the retained complete
receipt and all 138 dependencies with the independent verifier. It also
pins the unchanged [field128 arithmetic helpers](field128_transfer.py).
The new calculation uses exact rational shell counts, ratios, and folds.
Integer conversion rounds every dyadic endpoint upward. Only display
logarithms use floating point; the 40-bit acceptance test is exact.

The ignored receipt is `tmp/k20-two-ms-rs8-transfer-v1.json`. Four tests
check an alternate MDS formula, minimum packet support and normalization,
the pointwise envelope at every column support, and scalar-map transitivity.
The frozen field128 producer and existing certificate sources were not
modified. No local operator census was run.

An independent review checked the shortened-code enumeration, the
information-set envelope, joint support/label averaging, and the rare-event
split. It also checked that multiplication adjoints give the required
transitive action in the implementation's fixed tower basis. That review
did not rerun numerical endpoints; the producer's exact result remains
separately reproducible.

    C:/Python314/python.exe -B -m unittest discover -s research/workstreams/k20_two_ms -p test_rs8_transfer.py -v
    C:/Python314/python.exe -B research/workstreams/k20_two_ms/rs8_transfer.py --output tmp/k20-two-ms-rs8-transfer-fresh.json

The root replay reproduced the receipt byte for byte: SHA256
`77f6e6c6a18f7856e647db6b7946d659f87ddc484313a1ea5bfc73c9fa6f6f4f`.
The prototype's scalar oracle, field arithmetic, binary adjoint, output
alignments, and in-place suffix checks pass at K2048 and K20. Its declared
forward orientation uses the adjoints of the implemented transpose
randomizers in the fixed GF32 tower basis.

Full timings include packing, randomization, routing, and conversion.
Two unseen setup seeds and ABBA ordering preserve the approximately 4%
gain. The [public-profile report](../../../spin/RS8.md)
records the subsequent promotion measurements. Next prioritize the remaining data movement;
the present result is far short of the 2 ms goal.
