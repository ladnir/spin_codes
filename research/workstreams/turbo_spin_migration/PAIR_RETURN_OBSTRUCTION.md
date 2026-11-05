# A pair-return obstruction for fixed-width Turbo routing

## Claim and assumptions

Fix 64 packet regions and one four-bit packet from each local group in
each region. Let `L` be a multiple of 16, so each region consists of
whole 64-bit physical steps. Within each group, shuffle the 64 packets uniformly. In
each region, independently permute the `L` group slots uniformly. Let
the causal inner process 64-bit steps by

\[
 y_i=x_i+A a_i,\qquad a_{i+1}=M_i a_i+C x_i,\qquad a_0=0,
\]

where the `M_i` are independent uniform elements of `GL(s,2)`, `s` is
fixed, and every nonzero four-bit input supported in one physical packet
slot has nonzero feedback `Cx_i`. Assume the local outer is injective;
its distribution may be arbitrary and may be shared across groups.
For any fixed `delta>0`, the setup probability that the resulting code
has relative distance at most `delta` is bounded below by a positive
constant independent of `L` along these natural sizes. Thus this
fixed-width family cannot have positive relative distance with
probability `1-o(1)` as `L` grows.

This applies to the current RS8/s22 route and inner, and also to a
selected shared RS8 outer or independently conditioned good local
profiles **if their route and GL updates remain as above**. It does
not undermine any finite-size certificate. The lower bound is
astronomically small.
The manuscript already discusses the broader constant-outer-memory
limitation in its Bazzi--Mahdian--Spielman section. The event here is a
specific witness for the Turbo setup law and map, not a claim of a new
general limitation on serial concatenation.

## One-group event

Condition on any local outer realization. Put one fixed nonzero message
in one group and zero elsewhere. Its local output has `v` active
packets for some `1<=v<=64`; injectivity ensures `v>=1`. Write
`m=floor(v/2)`. Choose `m` consecutive pairs of regions
`(0,1),(2,3),...,(2m-2,2m-1)`. If `v` is odd, put the last packet in
region 63. The group's uniform packet shuffle realizes exactly this
support with probability `1/binom(64,v)`. The packet labels can be in
any order.

Let `h=floor(L/32)`. In each pair, require the first packet's group slot
to be among the last `h` slots of its region and the second packet's
slot among the first `h` slots of the next region. If there is an odd
last packet, require its slot to be among the final `h` slots of region
63. Independent regional group permutations give probability
`(h/L)^v` for these placements.

At the first packet of a pair the state was zero, so its outgoing state
is `Cx_i!=0`. All intervening inputs are zero, and invertible updates
preserve a nonzero state. At the second packet, conditional on the
past, `M_i a_i` is uniform over the nonzero `s`-bit states. As
`Cx_i!=0`, the state becomes zero with probability exactly
`1/(2^s-1)`. The `m` return events use distinct fresh updates, hence
their conditional probabilities multiply. After a return, the state
stays zero until the next pair. An odd final packet need not return.

The output is zero outside the intervals from each pair's first to
second *physical step* and, when `v` is odd, from the last packet to
the end. A
physical step emits 16 packets at once, so a pair's interval contains
at most `2h+16` four-bit packet positions, including the remainder of
the second step. An odd terminal interval contains at most `h`.
Thus every output on this event has weight at most

\[
  (8h+64)m+4h(v\bmod2)
  \le256h+2048\le8L+2048\le16L=N/16,
  \qquad N=256L.
\]

The last inequality holds for `L>=256`. For `L>=64`, `h/L>=1/64`.
Conditional on *every* local outer profile, for `L>=256`,
the event therefore has probability at least

\[
 p_v(L)
 =\frac{(h/L)^v}{\binom{64}{v}(2^s-1)^m}
 \ge\frac{1}{\binom{64}{v}64^v(2^s-1)^m}>0.
\]

For the current `s=22`, the coarse bound `p_v(L)>2^-1152` holds
uniformly over `1<=v<=64`. It suffices to refute convergence of the
bad-setup probability to zero at either 10% or 11%. To obtain any
smaller fixed positive relative-weight cutoff, replace `1/32` by a
smaller fixed window fraction; its event probability remains positive
as `L` tends to infinity.

## Check against the implemented inner

The K18/K20 s22 certificate states the physical-step recurrence above,
with independent uniform `GL(22,2)` updates. The public setup creates
one such update per 64-bit epoch. The certified map payload has SHA256
`e2467b9d5de8403ccd8263f0f55c73df128a8a2dd4c4cc405bc3a8e057fc3581`.
Direct binary-rank calculation from its 64 feedback columns gives rank
four for each of the sixteen contiguous four-column packet slots.
Consequently every nonzero packet has `Cx_i!=0`, regardless of its
label or physical slot. The event needs no extra packet-label
randomizer. Reproduce the map check with
`python research/workstreams/turbo_spin_migration/verify_pair_return_map.py`.

## Consequence for the research path

The earlier [late-group event](ONE_GROUP_OBSTRUCTION.md) used rare
low-packet-weight local words and showed bad-setup probability tending
to one for the original independently randomized RS8 outer. Selecting
a better fixed outer removes that event, but the pair-return event
survives. Conditioning independent local profiles on high packet
minimum distance also leaves it intact.

An asymptotic repair must change at least one fixed-width mechanism:
make the number of active local packets grow with length; make the
probability of each pair return vanish with length (for example through
a growing state, with its encoding cost accounted for); or constrain
the route and feedback so close, cancellable pairs cannot occur. The
existing growing-outer theorem takes the first route. We have no
replacement asymptotic proof for a finite RS8 outer.
