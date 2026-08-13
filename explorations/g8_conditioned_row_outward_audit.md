# Audit of the g=8 conditioned-row outward evaluator

The packet-width generalization is mathematically sound for `g=8`.  Let
`t_j` be the frozen fugacity for a packet of weight `j`.  The normalized
64-bit moment at total weight `w` is

\[
R_w(t)=\frac{[x^w]\left(\sum_{j=0}^{g}\binom gj t_jx^j\right)^{64/g}}
{\binom{64}{w}}.
\]

Thus the exact atom polynomial must be raised to `64/g`.  The implementation
does this with `local_block_atoms = BLOCK_BITS // local_group_bits`.  Its
coefficient vector always has degree 64, so the conditioned factors
`R_{w+b}` remain defined for `0 <= w <= 63` and `b in {0,1}`.

No later term acquires a packet-width dependence.  The `42/43/43` bands and
both split spectra describe coordinates of one fixed BCH row.  The puncture
transform and its divisor 42 also concern that row.  The 63 free rows, the
Cauchy split, the 128 normal tiles, the 128 hole tiles, and the graph spectrum
are unchanged.  Packet grouping acts within each 64-bit block and is
orthogonal to these coordinate-level objects.  The only other width-dependent
outputs are the `g+1` affine charges and the reported metadata.

The hardening lifecycle now preserves the selected width.  The parent process
initializes `GROUP_BITS` and `BLOCK_ATOMS` before any cache return.  Worker
processes receive the same width.  The checkpoint key includes `group_bits`,
and the normalization cache is cleared during initialization.  The built-in
full-bijection witness now returns `g+1` zero charges.  These conditions avoid
cross-width cache reuse and stale parent state.

The regression gate should retain four checks.

1. Freeze the complete g=4 outward interval endpoints and require the implicit
   and explicit `group_bits=4` paths to agree exactly.
2. At zero tilt, require the g=8 profile interval to contain `K=2^20`.
   Equivalently, the exact moments are all one, each tile contributes 4096
   bits, and the normalized graph average contributes 524288 bits.
3. Validate exact source masses: each unpunctured split spectrum has mass
   `2^64`, the punctured spectrum has mass `42*2^64`, and the graph spectrum
   has mass `2^24`.  Require one zero word in each source where applicable.
4. Freeze several nonzero g=8 endpoint intervals.  Compare them with an
   independent high-precision evaluation of the same exact dyadic
   fugacities.  Do not require the discovery scalar to lie in the interval:
   the discovery path evaluates a nearby binary64 expression.  A small
   proximity check is diagnostic only.

`scripts/check_packet_group_g8_conditioned_row_outward.py` passes all five
frozen g=8 vectors, the zero-tilt mass identity, and the canonical g=4
asymmetric witness.  Its current endpoint digest is
`3f038cdab9d6ee64867aece01baccab72fc544d69e0a8f3ce1d43d41b15e4de2`.
`scripts/packet_group_conditioned_row_high_precision.py` now supplies the
independent numerical check. It reconstructs the exact-dyadic packet
polynomial, both pair-spectrum sums, and the graph average with 180-digit
Decimal arithmetic. All five high-precision values lie inside their outward
intervals. The reference also verifies the source-spectrum masses before use.
