# Construction simplifications for the next Riffle iteration

This note records construction changes that the implementation may use.  None
is part of the frozen `g=4` certificate until a new construction statement and
proof explicitly adopt it.

## 1. Append the graph word

Let `N0=2^21` denote the power-of-two data body.  Encode the complete,
unpunctured data body through the existing permutation and inner map.  Append
the 128 graph coordinates as a separate linear tail.  The resulting length is

```text
N' = N0 + 128.
```

The main permutation still acts on exactly `N0` coordinates, so its fast
power-of-two implementation is unchanged.  The tail needs no large
permutation.  In the transposed encoder, the tail adds 128 field outputs and
the small transposed graph circuit.

This change removes three proof mechanisms:

- puncturing 128 data coordinates;
- placing graph coordinates into the punctured holes; and
- averaging the correlation between punctured data bits and graph bits.

The proof may ignore the nonnegative weight of the appended tail.  It must then
prove that the data body alone exceeds the new distance threshold.  At relative
distance `0.09`, the threshold rises by about 12 bits.  Therefore, removal of
the puncture tax must recover more than the cost of these 12 extra output
weights.  This comparison is a small numerical experiment, not an assumption.

This option best matches the API preference that the user dimension remain a
power of two while the encoded length may be irregular.

## 2. Sample a no-clumping packet layout

Define a clump as two coordinates in one physical `g`-packet that originate
from the same local EBCH block.  A stronger variant also forbids two slots from
the same outer tile or band position.

The construction should sample directly from the allowed layouts.  It should
not sample an unrestricted layout and silently condition afterward.  Direct
sampling gives a clear probability space and avoids an unrecorded factor
`1/Pr[no clump]` in the first-moment bound.

A fast implementation can use `g` lane maps with distinct offsets or distinct
small affine keys.  For packet index `i` and lane `s`, the producer computes the
source from one lane-specific map.  Distinct lane keys enforce the local
no-clumping rule.  The slot map is fused into the existing push address, so the
layout requires no intermediate buffer.

Before fixing the rule, measure two quantities on the current construction:

1. the frequency and type of local-code clumps; and
2. the proof penalty attributable to those clumps.

The weakest rule that removes the measured penalty is preferable.  A stronger
rule can reduce randomness without improving the bound.

## 3. Add a local packet mixer

The mixer is a fixed-width invertible binary circuit on each `g`-packet.  Its
transpose acts on `g` elements of `GF(2^128)` with the same XOR count.  The
mixer is local: it performs no global accumulation and no random memory access.

Three candidate families are worth testing.

### One matching layer

Choose an oriented perfect matching of the packet lanes.  For each matched
pair `(a,b)`, apply `x_b <- x_b + x_a`.  This costs `g/2` field XORs and has
depth one.  A random lane permutation makes the conditional output-weight law
depend only on the input weight.

For `g=8`, the cost is four field XORs per packet.  This is the cheapest probe.

### Permuted local accumulate

Choose a lane permutation `P`.  Apply a prefix-accumulate circuit in the order
defined by `P`, then undo only the lane labeling through fused addresses.  The
circuit costs `g-1` field XORs.  Its transposed circuit is a suffix accumulate
with the same cost.

For `g=8`, the cost is seven field XORs per packet.  This operation is not the
previous global PAP construction.  All dependencies remain inside eight
lanes, and independent packets can execute concurrently.

Uniform `P` acts transitively on inputs of a fixed weight.  Therefore, the
local transition law can be tabulated exactly as a `(g+1)`-class weight kernel.
For `g=8`, direct enumeration of all `8!` lane orders is small.

### Random butterfly

Use `log2(g)` layers of disjoint invertible shears with random switch choices.
This costs `(g/2)log2(g)` field XORs and has logarithmic depth.  For `g=8`, the
cost is 12 field XORs.  This family offers more diffusion but should be tested
only if the cheaper mixers do not give enough proof gain.

## 4. Selection experiment

Evaluate four variants at `g=4` and `g=8`:

```text
baseline
no-clumping layout
no-clumping + one matching layer
no-clumping + permuted local accumulate
```

For each variant, record:

- exact local weight-transition tables;
- the security-margin change at the current hard profiles;
- XORs per field element;
- cycles and memory traffic in the isolated transposed kernel; and
- whether the profile state remains `g+1` Hamming-weight classes.

Adopt a mixer only if its measured proof gain is material and its transposed
throughput remains competitive with Tungsten.  Construction changes do not
replace the current proof task: the existing `g=4` construction should first
receive the corrected Collatz sweep, because that sweep may close without any
change.
