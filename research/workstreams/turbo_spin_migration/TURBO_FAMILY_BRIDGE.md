# One Turbo architecture at finite and growing lengths

## Replacement criterion

To replace Structured SPIN in the paper, Turbo should be presented as
one outer--route--inner architecture with two outer instantiations, in
the same sense that the current paper uses fixed BCH for its finite
results and growing Golay--BAA for its asymptotic theorem. The finite
and asymptotic Turbo members should share the four-bit packet route,
the width-64 inner recurrence, and the same proof interface. The
outer's size and circuit may differ with the regime.

Let `b` be the outer output length, `k=b/2`, `L` the number of outer
blocks, `N=bL`, and `R=b/4` the number of packet regions. Take `L`
divisible by 16, so a region's `L` four-bit packets form whole 64-bit
inner steps. For each outer block:

1. Encode `k` message bits to `b` bits with the chosen rate-half outer.
2. Split the result into `R` four-bit packets and shuffle those packets
   uniformly into the `R` regions, one packet per region.
3. Independently shuffle the `L` block slots within each region.
4. Serialize the regions and apply the width-64 inner with fixed maps
   `A,C`, a 22-bit state, and independent uniform `GL(22,2)` updates
   at each physical step.

The finite member at `b=256` uses the existing RS8/GF32 outer. At
`K=2^18` and `2^20`, this exactly describes the certified public
`PacketRs8T64` family under its ideal setup law. The tuned `K=2^16`
member uses paired-s15 rather than s22; it is a parameter variation
within the architecture, not the literal s22 instance.

For the first asymptotic member, let `b` grow logarithmically and use
the paper's shared Golay--BAA outer
`Acc_b o tau_2 o Acc_b o tau_1 o H_{24,b}`. This keeps the Turbo route
and s22 inner. The old theorem already certifies this outer's binary
weight spectrum for its own schedule. It does **not** certify this new
composition: four-bit packet occupation, packet labels, and the s22
inner require new bounds. The old constant `39/4` for `b/log2 N` is a
control value, not an approved schedule for this member.

All stages have linear online work if `b=Theta(log N)`: the chosen
Golay--BAA outer has `O(b)` work per block, the packet shuffles touch
each packet a constant number of times, and the fixed-width inner has
`O(1)` work per 64-bit step. This is a complexity statement, not a
performance measurement.

## What a proof must establish

The first proof experiment fixes the known Golay--BAA outer. It then
tests Turbo's route and s22 inner in three regimes:

* **Fixed occupation:** one or finitely many active outer blocks. The
  growing outer must supply enough active packets that late placement
  and pair returns have vanishing total first-moment charge. The
  [fixed-width obstruction](PAIR_RETURN_OBSTRUCTION.md) no longer
  applies verbatim because `R=b/4` grows, but that fact alone is not a
  distance proof.
* **Growing sparse occupation:** control block choices, packet support,
  and repeated state returns uniformly as the active-block count
  grows below a fixed fraction of `L`.
* **Positive occupation:** combine a realized outer-spectrum majorant
  with a packet-route comparison and an s22 output-weight moment.
  The target is a strictly negative exponent at 11% distance for
  every relevant occupancy and outer weight.

The existing BAA certificate counts **binary weight**. Packetization
needs more information: a word of weight `w` may occupy anywhere from
`ceil(w/4)` to `min(w,R)` packets, and its nonzero packet labels need
not be uniform. The initial proof should use the actual packet
distribution or a rigorous worst-case envelope. If that envelope is
too loose, a sampled within-block bit shuffle or packet-label mixer
is an explicit construction change whose cost must be measured in
both finite and asymptotic members before calling them the same
implementation design.

Success means a complete `>11%`, rate-half, `O(N)` theorem for this
growing member **and** retained finite Turbo certificates and
measurements under the same route/inner architecture. A weaker positive
asymptotic distance would show scalability but would not replace the
paper's 11% theorem. Once the Golay control closes, the
[RS8-base growing outer](SCALABLE_RS8_BRIDGE.md) can be evaluated as an
additional performance variant, without making that optimization a
prerequisite for architectural replacement.
