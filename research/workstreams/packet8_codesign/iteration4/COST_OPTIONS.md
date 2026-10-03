# Low-cost map changes for byte packets

Subsequent proof-gate decision: the t32/s16 candidate remains negative at
small support even with the fractional routing proof (-109.11 bits at q=16
with eight-step grouping). Its low operation count alone is not a reason
to implement it. The wider-outer/scalar24 candidate is now the priority
for an implementation gate: its all-occupancy floating bound is positive,
while its runtime remains unknown. See [the iteration result](README.md).

The measured 16-bit-state byte-packet encoder takes 91.120 microseconds at
K=65,536. This note counts additional arithmetic for proposed constructions.
It reports no new encoder timing and does not predict that instruction
counts translate proportionally into time.

The counts use the literal intrinsics in `spin/experiments/packet8_codesign/Fast.cpp`,
the retained shared-parity outer, and the existing tower-field randomizer.
A `Pair` contains two ZMM registers, one for each half of a 128-bit payload.
Applying one byte map to a `Pair` therefore takes two GFNI instructions.
The outer GFNI multiplication instructions and inner affine instructions
need not have identical hardware costs.

## Two-band scaling of the expansion

Work over the AES polynomial basis F=GF(256). The actual24 maps are

```
A(a,b,c)_h = a + h b + h^2 c,       h=0,...,7,
C(x) = (sum x_h, sum h x_h, sum h^2 x_h).
```

The proposed family scales A only: use sigma_h=1 for h<4 and sigma_h=d
for h>=4, where d is a fixed nonzero byte. The map C, byte routing, and
outer code do not change. The original-order step remains

```
y = x + A' s,          s_next = M s + C x.
```

Its input diagonal is the identity, so this finite causal transform is
invertible regardless of C A'. This statement does not apply automatically
to a feedforward variant with a different input diagonal.

For reverse input w and future adjoint state t, the transpose step is

```
xbar = w + C^T t,      t_prev = M^T t + A'^T w.
```

Thus scaling A changes reverse feedback, not reverse emission. In the
following formulas, D_d^T denotes the binary adjoint of multiplication by d.
It is not ordinary multiplication by d in the polynomial basis.

Define the low-band sums L0=w0+w1+w2+w3, Lo=w1+w3, and Lb=w2+w3.
Define the corresponding high-band sums H0=w4+w5+w6+w7, Ho=w5+w7,
and Hb=w6+w7. Compute

```
v0 = D_d^T H0,         vo = D_d^T Ho,         vb = D_d^T Hb,
S0 = L0 + v0,          S1 = Lo + vo,           S2 = Lb + vb.
```

Then the three reverse-feedback coordinates are

```
p0 = S0,
p1 = S1 + D_2^T S2 + D_4^T v0,
p2 = S1 + D_4^T S2 + D_16^T v0.
```

The original unscaled feedback already uses the final four byte maps.
Two-band scaling adds only the first three maps per payload half: six
GFNI instructions per ordinary step. Scaling three separately evaluated
high-band polynomial moments would miss this sharing and cost more.
The current pairwise packet schedule already provides H0 and Hb; it must
retain the low and high odd sums separately.

An arbitrary scale on seven of the eight bytes would instead add seven
maps per payload half. This is fourteen instructions per step. The
two-band family is therefore the preferred small family to screen first.
Choose two or three d values by the separate proof screen, retaining d=1
as the exact unscaled control. No spectrum or distance conclusion follows
from this operation count.

For this particular band partition,

```
C A' = [[0,0,0], [0,0,tau], [0,tau,0]],    tau = 6(1+d) in F.
```

The composition usually is nonzero. It does not invalidate the stated
original-order recurrence or the existing feedback-rank properties.

## Six-product 24-bit updates

The polynomial z^3+z+1 has no root in F, as checked for all 256 bytes.
It therefore defines the field E=F[z]/(z^3+z+1), with 2^24 elements.
Identify the state with a+bz+cz^2 and a sampled scalar with r0+r1z+r2z^2.
Their product needs six F multiplications:

```
p0  = a*r0,                 p1  = b*r1,                 p2  = c*r2,
p01 = (a+b)*(r0+r1),        p02 = (a+c)*(r0+r2),        p12 = (b+c)*(r1+r2),
out = (p0+p12+p1+p2,        p01+p0+p12,                p02+p0+p1).
```

The binary adjoint also needs six byte maps. For incoming coordinates
(t0,t1,t2), compute

```
q0  = D_r0^T(t0+t1+t2),     q1  = D_r1^T(t0+t2),      q2  = D_r2^T(t0),
q01 = D_(r0+r1)^T(t1),      q02 = D_(r0+r2)^T(t2),     q12 = D_(r1+r2)^T(t0+t1),
out = (q0+q01+q02,          q1+q01+q12,               q2+q02+q12).
```

These are binary-transpose formulas, not an assumption that field
multiplication is self-adjoint. The tests compare every binary matrix
column for 40 distinct scalars against an independent polynomial product.

Sampling an independent uniform nonzero scalar at each step gives the
same fixed-nonzero-state action law as a fresh uniform GL(3,F) matrix.
For each fixed nonzero state s, multiplication by s is a bijection on E*.
Thus r*s is uniform on E*. The matrix ensembles differ, but the local
first-moment analysis needs only this conditional action law. Independence
of the new scalar from all earlier setup choices is required.

An implementation-friendly alternative samples the forward update as the
binary adjoint of multiplication by r. This family is also transitive.
For distinct r and r', the difference of their adjoints is the invertible
adjoint of multiplication by r+r'. Hence its images of a fixed nonzero
state are all distinct. The reverse encoder then uses the ordinary
six-product circuit and six GFNI byte multiplications, rather than six
affine maps. Its precomputed coefficients occupy six bytes per step,
instead of six 64-bit affine matrices. This follows the existing outer
randomizer convention, but requires an explicit construction choice.
No sampler or implementation was changed in this task.

## Whole-call arithmetic at K=65,536

There are 2,048 physical steps. The last reverse step omits feedback and
the state update. It still packs the input and applies C^T. The counts below
include that boundary and both payload halves.

| Inner choice | GFNI per ordinary step | Whole inner GFNI count |
| --- | ---: | ---: |
| Measured 16-bit state, GL2 update | 32 | 65,524 |
| 24-bit state, GL3 update | 50 | 102,374 |
| 24-bit state, six-product scalar update | 44 | 90,092 |
| Scalar24 with two-band A scaling | 50 | 102,374 |
| Scalar24 with seven distinct nontrivial byte scales | 58 | 118,750 |

For the measured kernel, the 32 instructions split into 16 for packing,
four for C^T emission, four for A^T feedback, and eight for the GL2 update.
The unscaled scalar24 version uses 16, eight, eight, and twelve, respectively.
No additional routed packet, payload repacking, or outer operation is
required by these inner changes. Extra XORs, state registers, coefficient
loads, and dependency chains remain relevant implementation costs.

## Shorter steps with a quadratic-field update

The literal width32 candidate keeps a 16-bit state and processes h=0,...,3
instead of eight bytes. A and C retain their two-moment formulas. Each
physical step needs one byte map per payload half for C^T emission and
one for A^T feedback. There are twice as many steps, but the total packing
and routed payload volume remains unchanged.

The existing tower E2=F[u]/(u^2+u+0x20) supplies an exact three-product
update. For state a+bu and scalar c+du, compute

```
p0 = a*c,       p1 = b*(0x20*d),       p2 = (a+b)*(c+d),
out = (p0+p1, p0+p2).
```

The coefficients c, 0x20*d, and c+d are prepared during setup. Choosing
the forward update as this multiplication's binary adjoint makes the
reverse update exactly the existing three-product byte kernel. Independent
uniform nonzero scalars give the required fixed-state transitivity by the
same injectivity argument used above. This is not uniform sampling from GL2.

The resulting ordinary step needs 18 GFNI instructions: eight for packing,
two for emission, two for feedback, and six for the update. At K=65,536,
4,096 steps and the omitted final reverse update give 73,720 instructions.
That is only 8,196 more than the measured width64/GL2 inner. Retaining GL2
instead would cost 81,910. The shorter-step scalar candidate has the lowest
added arithmetic among the real construction changes considered here.
It still doubles loop iterations and state-dependency boundaries, so its
operation count alone does not establish an end-to-end timing improvement.

## Reality check for the wider outer

The small outer has 512 groups of 128 input bits. Each of its 16 symbols
uses a GF(2^16) randomizer, requiring three byte products per payload half.
The proposed wider outer has 256 groups of 256 input bits. Its 16 symbols
use GF(2^32) randomizers, requiring nine byte products per half in the
retained quadratic tower circuit.

| Whole-call outer operation | Small outer | Wider outer |
| --- | ---: | ---: |
| GFNI products for symbol randomization | 49,152 | 73,728 |
| GFNI maps in shared RS parity circuit | 30,720 | 30,720 |
| GFNI maps for final output unpacking | 16,384 | 16,384 |
| Total GFNI | 96,256 | 120,832 |
| XORs within symbol randomization | 49,152 | 122,880 |

The wider outer adds 24,576 GFNI multiplications and 73,728 XORs. Its
existing two-half intermediate occupies 8 KiB per group, versus 4 KiB for
the small native-packed outer. There are half as many groups; this does not
make the larger per-group working set irrelevant to register spills or cache
traffic. The total number of routed byte packets remains unchanged.

These counts assume native byte-packed input and the measured 15-GFNI
shared parity circuit. The older K20 wrapper repacks ordinary input and
uses an 18-GFNI parity circuit. Importing that wrapper literally would add
unnecessary work; its parity difference alone adds 6,144 GFNI instructions.
Reuse the field multiplication and parity circuits, not the obsolete
packing interface.

Widening the outer and using scalar24 together raises whole-call GFNI count
from 161,780 to 210,924, a 30.4% increase. It also adds the stated outer XORs
and the inner's extra linear combinations. The 100-microsecond target allows
only 8.880 microseconds beyond the measured baseline. At an assumed 4.5 GHz,
that is 39,960 cycles: 0.813 incremental cycles per added GFNI instruction
before charging other work. This is a budget, not an instruction latency or
a timing estimate. Meeting it would require substantial overlap with work
already limited by memory or other execution resources.

The proof improvement may justify testing the combined candidate, but the
cost audit does not establish that it stays below 100 microseconds. Keep
the wider-outer/16-state and small-outer/scalar24 candidates as single-change
controls. Each adds about 24,576 GFNI instructions, with different XOR,
memory, and proof tradeoffs.

## Comparison with feedforward

For the unscaled maps, the earlier feedforward candidate uses F=I+AC.
Its transpose forms p=A^T w, emits w+C^T(t+p), and updates p+M^T t.
This evaluates A^T and C^T once each, adding only two state-vector XORs
for the 16-bit state or three for the 24-bit state, per payload half.
However, it needs all feedback before its first routed store. Retaining
eight packed packets requires 16 ZMM registers before accounting for state
and temporaries. Repacking instead adds sixteen GFNI packing operations.
The algebraic count is attractive, but register pressure and lost overlap
make it less predictable than the streamed two-band scaling change.

Combining feedforward with A-only scaling requires a fresh invertibility
check and new joint emission/feedback moments. CA=0 from the old maps cannot
be reused. The two-band family was the recommended low-cost scaling screen.
Its first bounded q=16 gates improved the margin by only about 0.2 bits,
so the available evidence does not justify implementing that change. The
shorter width32 step with a quadratic scalar update is a cheaper new
candidate to test in the proof. If a 24-bit candidate earns an implementation
test, use the six-product update as the performance-conscious control.

The subsequent width32 proof screen did not close the smaller-support gap:
q=16 gave -126.06 bits with four-step grouping and -109.11 with eight-step
grouping, at tilt 0.06 and power 0.4. Its q=64 four-step bound was -203.32
bits. The positive q=119 result does not repair those missing classes.
Consequently, width32 is not the primary implementation candidate.
The wider-outer/actual24 combination subsequently reached a complete
floating all-occupancy proposal of approximately 68.89 bits. Its next
decision point is an isolated implementation measurement against the
arithmetic budget above; outward certification remains separate.

Nine portable tests in `test_cost_options.py` verify both fields, ordinary
and adjoint update circuits, the two-band transpose at all 255 nonzero scales,
the CA' formula, and source-level operation totals. All passed. No production source,
frozen numerical input, encoder benchmark, or timing claim was changed.
