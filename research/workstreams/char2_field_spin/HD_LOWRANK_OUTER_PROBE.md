# Extra moment terms do not give a free stronger mixer

The faster experimental outer uses sixteen RS[16,8] rows per group and
the fixed matrix H2=I+11^T+uu^T between independent nonzero diagonals.
Its sixteen labels are 06+span_GF2(01,0c,50,b0) in the AES representation
of GF256; byte constants here are hexadecimal.
The question here is whether adding one or two moments simplifies the
proof without losing the cheap circuit.

For d=2,3,4, let V_d have row (1,u_i,...,u_i^(d-1)) and set
H_d=I+V_d V_d^T. This changes the construction for d>2. All three
matrices are involutions: V_d^T V_d=0 by the power sums over the affine
coset, hence H_d^2=I. Every column also has full support. Those facts
alone do not retain the existing support enumerator.

The exact local probe enumerated all 65,536 principal subsets, all
14,400 two-by-two minors, and all 255 nonzero ratios for each of the
120 input-coordinate pairs. It completed in 1.828 seconds.

| Moment count d | Singular principal subsets | Zero 2-by-2 minors | Minimum output weight for a two-coordinate input |
|---:|---:|---:|---:|
| 2 | 0 | 0 | 15 |
| 3 | 256 | 56 | 14 |
| 4 | 256 | 116 | 14 |

The principal singularities have nullity one. There is a short example
for d=3. Restrict to the labels A={06,07,ba,bb}, in hexadecimal AES-field
notation. They form an affine GF4 coset. Direct field arithmetic gives

    I_3 + V_A^T V_A = [[1,0,0], [0,1,1], [0,1,1]].

This matrix has nullity one, as does the four-by-four principal matrix
(H_3)_{A,A}. The explicit example is independently checked in the tests.
It does not contradict invertibility of the complete sixteen-by-sixteen
matrix. It prevents transferring the rank2 argument that inverts every
principal restriction and reduces the remaining coordinates to two moments.

Likewise, the weight-14 pair outputs show that more moment terms do not
monotonically improve sparse-input expansion. Neither local defect proves
poor whole-code distance; a different analysis might still succeed.
However, this family adds work and loses properties already available for
rank2. There is no current reason to implement or benchmark it.

The exact results are retained in
`theory/hd_lowrank_outer_probe_exact.json`; the helper is
`theory/hd_lowrank_outer_probe.py`. The probe contains no native kernel,
timing result, distance certificate, or change to the frozen construction.
