# A larger outer for the byte-packet construction

The larger outer improves the sparse bound but does not close the middle range.
At K=65,536, the floating proposal still loses approximately 3,697 bits near
61 active outer groups. No implementation change is recommended from this gate.

## Construction and exact geometry

One group contains eight parallel GF(16) RS[16,8] rows. Its input has 256 bits
and its output has 512 bits. At each of the 16 symbol positions, concatenate
the eight row symbols into a 32-bit vector. These vectors form an MDS[16,8]
code over GF(2^32): the original generator remains MDS after extending scalars.

Independently at each symbol position, sample an invertible binary map whose
image of every fixed nonzero input is uniform among the nonzero 32-bit vectors.
A uniform GL(32,2) map or multiplication by a uniform nonzero GF(2^32) element
satisfies this fixed-input premise. These families need not have the same joint
law on several messages. Only the stated fixed-message law is used here.

Split each randomized symbol into four byte packets. Shuffle all 64 packets
within the group, then independently shuffle group slots in each of 64 regions.
Every group contributes one packet per region. The unchanged inner processes
eight packets per physical step, with a 16-bit state and fresh GL(2,GF(256)) updates.
State persists between regions; the initial state is zero and the final state is discarded.

| Quantity | Prior small outer | Tested larger outer |
|---|---:|---:|
| Group input/output bits | 128/256 | 256/512 |
| Parallel GF(16) rows | 4 | 8 |
| Randomized symbol bits | 16 | 32 |
| Groups at K=65,536 | 512 | 256 |
| Byte regions | 32 | 64 |
| Physical steps per region | 64 | 32 |
| Total physical steps | 2,048 | 2,048 |

`wider_screen.py` uses these actual region dimensions, not the prior 512-slot region.
The cutoff is 13,107, corresponding to a target minimum distance of at least 13,108.

## What the screen computes

The single-group calculation uses exact expected support-shell counts from
the MDS enumerator and the randomized-symbol law. Their total is exactly 2^256-1.
Conditional on a support, active byte labels are independent uniform nonzero bytes.
The group shuffle makes that support uniform among the 64 regions.

For q at least two, the screen uses the existing uniform-input measure bound

```
beta = 2^512 / (2^32-1)^8.
```

Consequently, each region has Binomial(q,255/256) active packets under the comparison.
Their placement among its 256 slots is without replacement. The ordered matrix
recurrence retains boundary states; regions are not independently reset.
The unchanged local finite-family operator retains the exact shapes of one-step births.
Its positive uniform-density closure remains an upper bound, not an exact transition law.

All evaluations are floating proposals. They have no outward endpoints and
are neither whole-code certificates nor evidence of actual low-weight codewords.
In particular, q=2 uses the uniform envelope, not exact joint outer-shell counts.

## Results and comparison

The combined grids in `wider_v1.json`, `wider_v2.json`, and
`wider_middle_fine_v1.json` cover every integer q from 2 through 256.
Take the largest reported margin separately for each q. The combined weakest
value is -3,696.876641 bits at q=61 and tilt 0.52. The single-group calculation
gives +124.896790 bits, with its support-wise tilts chosen before summation.

The wide q=61 and small q=122 cases have identical active-group fractions.
At the same tilt 0.52 they give -3,696.876641 and -4,022.521743 bits, respectively.
Thus this matched comparison gains 325.645102 bits, far short of the 40-bit target.
The earlier small-outer coarse grid's roughly -4,151-bit value is not its best
available result: `../byte_native_fine_q80_160.json` already reaches about -4,023 bits.

At wide q=61 and tilt 0.52, the all-zero-state path contribution alone has
margin +879.210766 bits. The complete comparison moment is larger by about
4,576 bits. Thus this failure is not explained by the path that never activates
the state. This diagnostic does not identify a particular bad message or prove
that the true outer distribution has the same behavior as its envelope.

## Reproduction and scope

From the repository root:

```text
python -B -m unittest discover -s research/workstreams/packet8_codesign/wider_outer -p test_wider.py -v
python -B research/workstreams/packet8_codesign/wider_outer/wider_screen.py --output <fresh-result.json>
```

The four tests check dimensions, exact support-count mass and envelope domination,
MDS symbol-shell counts, and the ordered placement recurrence on a tiny fixture.
Every result records its full tilt list, map identity, geometry, and source hashes.
All four saved runs verified their source pins on completion. The initial broad
run took about 9 seconds, the extended grid about 28 seconds, and the final
middle refinement about 2 seconds. These are proof-proposal runtimes, not encoder benchmarks.

The single-point control is `small_matched_q122_v1.json`. Its intentionally narrow
tilt list does not provide a useful single-group result; use it only for the matched q=122 comparison.
Likewise, the middle-only refinement's single-group field must not replace the broad run's sparse result.

No authenticated prior sources, implementation files, or production APIs were changed.
The next promising gate is a stronger byte-native inner, not another broad
counting run for this unchanged 16-bit-state design.
