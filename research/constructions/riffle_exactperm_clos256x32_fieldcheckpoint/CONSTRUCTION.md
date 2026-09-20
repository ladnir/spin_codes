# ExactPerm-Clos256x32

ExactPerm-Clos256x32 is an evaluator for Riffle ExactPerm FieldCheckpoint v1.
It does not define a new code distribution.

For one transposed region, setup first samples the same uniform permutation
(P\in S_{8192}) used by the existing construction.  Write
(8192=256\cdot32), and partition the input and output positions into 256
groups of size 32.

The permutation (P) induces a bipartite multigraph.  Its left and right
vertices are the input and output groups.  Every permutation edge connects the
group containing its input to the group containing its output.  Every vertex
has degree 32.

Kőnig's line-coloring theorem gives a proper edge coloring with 32 colors.
Each color is a perfect matching between the 256 input and output groups.  The
coloring factors (P) into three stages:

1. 256 input permutations of size 32 route inputs to colors.
2. 32 middle permutations of size 256 route input groups to output groups.
3. 256 output permutations of size 32 route colors to final offsets.

The composition equals (P) exactly.  Therefore, the permutation distribution
and all existing distance-analysis obligations are unchanged.  Any certificate
proved for the uniform region permutation applies without a new structured-
permutation lemma.  Only the evaluator changes.

The route-validation script constructs a seeded uniform permutation of 8192
positions, decomposes it by repeated perfect matchings, and compares the full
output vectors.  The prototype validates schedule construction and exactness;
it does not measure encoder performance.
