# A certified 15-bit restriction

Deleting row10 from the fixed map in `DISJOINT_PAIR_CONSTRUCTION.md` gives
a complete setup-failure margin of **62.04117326862128 bits**. The replay
covers every occupancy q=1,...,512 at 256-bit outward precision.

Retain that document's outer code, routing distribution, coordinate order,
and sixteen functions f0,...,f15. Remove f10=z0*z2+z1*z4 and keep the other
functions in their existing order. Define A15 by evaluating these fifteen
functions on all z in F2^6, and set C15=A15^T.

At step i, the inner receives x_i in F2^64 and state a_i in F2^15. It emits
y_i=x_i+A15*a_i and updates a_(i+1)=M_i*a_i+C15*x_i. Each M_i is sampled
independently from uniform GL(15,2). All setup choices remain mutually
independent. The initial state is zero; state continues across all 2,048
steps and is discarded after the last output.

The retained functions remain independent and satisfy C15*A15=0. Every
consecutive four-coordinate packet still has restriction rank four because
1,z0,z1,z0*z1 remain present. Fresh enumeration gives this expansion spectrum:

| Weight | 0 | 16 | 24 | 28 | 32 | 36 | 40 | 48 | 64 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| States | 1 | 16 | 2,288 | 6,912 | 14,334 | 6,912 | 2,288 | 16 | 1 |

For this ideal setup distribution, the complete replay proves

\[
\Pr[d_{\min}\le 13107]<2^{-62.04117326862128}<2^{-40}.
\]

Thus the binary [131072,65536] code has minimum distance at least 13,108
outside this setup-failure event. The statement does not certify a specific
deterministic setup seed. In particular, GF(2^16) state multipliers do not
transfer to a fifteen-dimensional subspace; this receipt uses GL15 updates.

`reproduce_s15_drop10.py` regenerates the map and every local operator. It
reuses only the rational tilt recipe, not numerical endpoints, from the
completed sixteen-bit replay. It computes the exact q=1 outer support sum,
the q=2 ordered support-pair sum, and a pointwise outer majorant for q>=3.
Every product preserves the state across steps and regions.

The q=1 contribution has margin 62.2857139563 bits; q=2 has 115.0524190862
bits. The weakest q=3,...,128 contribution is q=41 at 66.2361558764 bits.
The final union includes all 512 contributions, summed as exact positive
dyadic numbers and rounded upward once. Its comparison with 2^-40 uses
integer arithmetic.

The receipt is `paired-s15-drop10-whole-p256.json`, SHA-256
`89bed0c4442349a520d54b1aa6e914eac326beb3530123f8273f78fc764b867d`.
All 106 saved source hashes matched at completion and in two independent
receipt audits. The map identity is
`ca6c82af3074fb48e97e4948654e6cc0b4bf0993e7123e1287cc9893786413c5`.

Run a fresh replay from the repository root with an unused output filename:

```text
python -B research/workstreams/k16_codesign_100us/proof/reproduce_s15_drop10.py --output research/workstreams/k16_codesign_100us/proof/new-s15-whole-p256.json
```

The deleted row removes two quadratic channels and one feedback XOR. Those
operation savings do not establish a speed improvement; implementation
timings must decide whether to use this restriction.
