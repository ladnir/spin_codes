# RS8 setup and certificate interface audit

The public `PacketRs8T64` profile implements the same algebraic maps and
sampling *interfaces* as the K16, K18, and K20 finite certificates. The
information-theoretic certificates sample independent uniform objects;
the executable expands fixed 64-bit seeds through SplitMix. The latter is
a reproducible implementation of the interface, not literal independent
uniform setup. Do not state a statistical distance guarantee for the
seeded implementation without an additional pseudorandomness assumption.

| Object | Certified ideal law | Public implementation |
|---|---|---|
| Local outer | Eight systematic GF(16) RS[8,4] rows per 128-bit group, polynomial basis `0x13` | `spin/src/rs8/SetupScalar.cpp` interpolates at 0..3 and evaluates at 0..7 with that modulus |
| Column maps | For each group and each of eight 32-bit columns, an independent uniform nonzero GF(2^32) multiplier; forward uses its binary adjoint | `prepareOuter` rejects zero scalars from a domain-separated word stream; `Field.h` implements the fixed tower field; transpose uses multiplication and forward uses the adjoint |
| Packet route | Independent uniform permutation of 64 four-bit packets in every group; one packet per region; independent uniform permutation of group slots in each of 64 regions | `paired15/Paired15.cpp` and `certified22/Setup.cpp` build those shuffles with rejection-based divisor sampling |
| Inner | Width 64, fixed expansion A, feedback A-transpose, and independent sampled state updates per physical step; state crosses region boundaries | K<=65536 selects paired s15, with a uniform GL(15,2) update and one fixed coordinate; larger sizes select s22 with uniform GL(22,2) updates. Both samplers reject rank-deficient matrices. The s22 map identity is `e2467b9d5de8403ccd8263f0f55c73df128a8a2dd4c4cc405bc3a8e057fc3581` |

The route, column, and inner seed streams are separately derived from the
public seeds. Under an ideal independent-word source, the rejection
samplers produce the claimed uniform laws. Under the executable's finite
seed space, these choices are correlated and deterministic. This is the
same proof-to-implementation boundary used by the existing finite
certificates; it is not a newly found map mismatch.

The public profile selects paired s15 at K=2^16 and exact s22 at K=2^18
and 2^20, matching the three proof records. The K20 transfer's four
focused arithmetic/enumerator tests passed in this checkout. Its ignored
K16 numerical receipt is absent here, so the K16 test could not start in
this worktree. Running the same source at commit `f0656bf4` in the
canonical checkout, where retained receipts are available, passed all
four K16 certificate tests and all eight K18 certificate tests. No
benchmark was run in this audit.

The [one-group obstruction](ONE_GROUP_OBSTRUCTION.md) uses only the
ideal per-group column and packet-shuffle laws and the inner's causality.
It is independent of the seeded implementation's pseudorandomness.
