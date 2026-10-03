# Wider-outer, 24-bit implementation audit

Follow-up: the [completed checkpoint](README.md) now pairs this implementation
with an outward certificate and a 98.1625-us fresh-seed holdout. The generated
reverse-inner function has no stack spills. The audit below records the
construction review performed before those final checks.

The isolated implementation matches the proposed maps and routing geometry.
This review found no construction-to-kernel mismatch in the files listed below.
The floating proof proposal and the measured implementation remain separate
artifacts; this review does not turn the proposal into an outward certificate.

The review covers `spin/experiments/packet8_wider24`, including setup, both
scalar directions, the SIMD transpose, its retained outer dependencies, and
the correctness and timing harness. No implementation source was changed.

## Dimensions and routing

For input length K, the constructor sets N=2K and G=K/256 outer groups.
Each group contains eight parallel GF(16) RS[16,8] rows and outputs 512 bits.
At each RS position, the eight row symbols form one 32-bit symbol.
After randomization, its four consecutive bytes are four routing packets.

The constructor accepts positive multiples of 2048 within the checked storage
and 32-bit routing limits. Hence G is a multiple of eight, and no physical
inner step crosses a regional boundary. At K=65536, the dimensions are:

| Object | Count |
|---|---:|
| Outer groups | 256 |
| Byte packets per group / regions | 64 |
| Packet slots per region | 256 |
| Eight-packet physical steps per region | 32 |
| Physical steps in the complete inner | 2048 |
| State bits | 24 |

For group g, `columns[g][region]` names the group's packet assigned to that
region. A separate regional permutation assigns one slot to every group.
Thus `route[region*G+slot]` is an inverse route: it gives the outer-side
address read by the forward encoder and written by the transpose.
Both permutation stages use complete Fisher--Yates shuffles; this is not
an old 64-bit route interpreted as pairs of byte packets.

Scratch group stride 516 is measured in 16-byte records. Every group starts
on a 64-byte boundary, and every routed packet occupies eight records.
The four padding records per group are never routed or read by the outer.
The constructor's group bounds also bound N and its allocation sizes.

## Literal inner and its transpose

Let F=GF(256) use polynomial modulus 0x11b. Write the state as (a,b,c), with
a in the least-significant byte. For packet index h in {0,...,7}, define

\[
 (A(a,b,c))_h=a+hb+h^2c,
 \qquad
 Cx=\left(\sum_hx_h,\sum_hhx_h,\sum_hh^2x_h\right).
\]

The forward step emits y=x+As and updates s'=M_e s+Cx. State begins at zero
and persists through all regions; the final state is discarded.
These are exactly the A/C coordinates in `larger_state/maps24.py`.
In the polynomial-basis binary dot product, C is not silently identified
with A's binary transpose.

Let L_r be multiplication by nonzero r in
E=F[z]/(z^3+z+1). The implemented forward update is M_e=L_r^T. Therefore the
reverse step, with future-state adjoint p and output adjoint w, is

\[
 \bar x=w+C^Tp,
 \qquad
 p_{\rm previous}=L_r p+A^Tw.
\]

`Scalar.cpp` applies literal binary columns for both directions. `Fast.cpp`
uses the same equation. Its C-transpose emission directions are

\[
 b+c,\quad D_2^Tb+D_4^Tc,\quad D_4^Tb+D_{16}^Tc,
\]

where D_d multiplies a byte by d. These directions follow from
h=h_0+2h_1+4h_2 and h^2=h_0+4h_1+16h_2.
The feedback uses the raw packed input, not the emitted/routed output.
The raw odd, bit-one, and bit-two packet sums produce A^Tw by the same
three Boolean directions.

The six-product state update implements ordinary multiplication in E.
Its coefficients are r0, r1, r2, r0+r1, r0+r2, r1+r2, in that order.
Both implementations use `updates[epoch]`; reverse evaluation starts with
zero future state at the last epoch. Omitting the final reverse update
at epoch zero is valid because the initial forward state is fixed to zero.

### Update-family transfer

The cubic z^3+z+1 has no root in F and is irreducible. Thus every nonzero
L_r is invertible. For fixed v!=0, the map r -> L_r^T v is injective:
for r!=t, its difference is L_(r+t)^T v, which is nonzero. Its image is
therefore every nonzero 24-bit state exactly once.

Fresh independent uniform nonzero scalars consequently give exactly the
fixed-state transition law used by the GL(3,F) proof model. This is a
transfer of the required action law, not equality of matrix ensembles.
Conditional on any prior state and current input, emission occurs before
the fresh update and is independent of its scalar. The same local comparison
operator therefore applies. Conditioning on routing geometry does not
change this conclusion when inner setup randomness remains independent.

As in the retained benchmark prototypes, `Words` expands deterministic
domain-separated seeds; exact distribution claims refer to ideal independent
setup words, not to information-theoretic independence of that seed expander.

## Outer coordinate convention

For symbol j, let v_c be the encoded RS bit with c=4*lane+bit.
If R is the binary matrix of ordinary GF(2^32) multiplication by its sampled
nonzero scalar, the forward symbol is x_p=sum_c R[c,p]v_c.
Thus the transpose receives ordinary physical coordinates p and computes
vbar_c=sum_p R[c,p]xbar_p.

`outerRows[j][c]` stores precisely row c of R. The scalar forward reads
its columns, and the scalar transpose reads its rows. The native outer
loads four consecutive physical bytes and applies ordinary tower-field
multiplication, so it computes the same R without an input nibble permutation.
Its low output array contains logical lanes 0--3; its high array contains
lanes 4--7. The retained parity/unpack circuit writes those halves at offsets
0 and 128 within the group's 256 output records.

The eight parallel RS rows have the required MDS symbol counts: any eight
symbol positions determine their input, and smaller coordinate restrictions
have the corresponding full rank. This counting fact does not require the
RS nibble basis to equal the tower field's multiplication basis.

Adjoints of nonzero GF(2^32) multipliers are transitive on nonzero binary
symbols by the same difference argument used for the inner. Independent
symbol scalars therefore preserve the proof's fixed-message randomized
symbol distribution. Their joint action on several messages need not match
uniform GL(32,2), and the proof does not require that stronger assertion.

## Layout, tests, and remaining performance question

Packing is the retained byte-plane transform: payload bit b occupies packed
lane 7-b within its byte. The scalar packing and unpacking agree with that
orientation. Ordinary input/output require only 16-byte alignment; SIMD
input loads and final output stores are unaligned. Scratch stores are aligned.
The full transpose consumes every input before writing final output, so
output equal to input is valid with disjoint scratch.

The harness checks the route bijection, one packet per group per region,
scalar agreement, four input/output alignments, padding and output guards,
input preservation, in-place output, route-only output, packing roundtrips,
and the complete forward/transpose dot-product identity. Root reports that
all six K/seed cases pass. The separate nine-test algebra suite also passes,
including the cubic field, six-product multiplication, and binary adjoints.
No additional encoder benchmark was run for this audit.

The normal inner step has 44 GFNI instructions: 16 for input packing, eight
for C-transpose emission, eight for A-transpose feedback, and twelve for
the scalar update. At K16, the skipped last reverse update removes twenty.

The next optimization decision should use the generated assembly. If the
scalar stage spills, stream its six products into three output accumulators
instead of retaining all six products simultaneously. That is a scheduling
change, not a changed field formula. The emission stage also carries three
state pairs and several raw moment pairs; any spill there should be located
before changing the streamed packet order. The outer's fixed 8 KiB temporary
is intentional packed data, not by itself evidence of register spilling.
No speculative kernel change is justified before the fresh-seed holdout and
the assembly inspection requested by root.

## Reviewed source identities

SHA256 values below identify the reviewed isolated sources. Paths in this
table are relative to `spin/experiments/packet8_wider24`.

| File | SHA256 |
|---|---|
| Packet8Wide24.h | `320bab2f396318f0986900cddcfe6cc2f0b6b8e3b0ff5e2e3228a5eca14ec19b` |
| Setup.cpp | `0299f00d361821b9975929a8fa81242ca9e1c3d2b724a64f7d6f3d9f6aefa14e` |
| Scalar.cpp | `3a0bf2d2800005b7ef0ca46a9c92539cf9d1bc6a7089c611281d8484544e2244` |
| Fast.cpp | `ad0712bcbe11e761328ca6bcc0ac75b4b7ada0093d32d03b59550be65437af32` |
| Outer.cpp | `7a133c0cda47dfb02347cbe423fe1f46722b1bf3048706cd676b8f1b032e4c9c` |
| bench.cpp | `68bf0154f3ea9ac4c822b8ee0b88ca1bfa3b8350739161f35584fd9b4033373d` |
| CMakeLists.txt | `dbd8d2fc18eab7009e27f480af58399af43a3b3cabfef8305aef92598f0c9c2f` |
| run_ab.sh | `6d2f69273dd70e2a1176d85edb1ca0988f7bbf6dd4f63b9d713309db71862abd` |
