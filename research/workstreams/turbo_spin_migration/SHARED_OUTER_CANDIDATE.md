# A selected shared RS8 outer as a scalability candidate

The [late-group obstruction](ONE_GROUP_OBSTRUCTION.md) relies on
independent per-group column maps giving every group a fixed positive
chance of producing a very low packet-weight local word. A possible
next construction samples eight GF(2^32) column multipliers once,
selects a favorable tuple, and reuses it for all groups. Independent
group packet shuffles and regional permutations remain. This is a
candidate, not a distance certificate or an implementation win.
The later [pair-return obstruction](PAIR_RETURN_OBSTRUCTION.md) rules
out an asymptotic high-probability distance theorem for this candidate
with the unchanged fixed-width route and uniform-GL inner. Its scope
is therefore finite-size analysis or a route/inner redesign.

## Exact local selection budget

Let B_v be the number of nonzero messages whose local 256-bit output
has exactly v active four-bit packets before shuffling. Under eight
independent uniform nonzero GF(2^32) multipliers, its expectation C_v is
the coefficient of

\[
 C(z)=\sum_{h=5}^8 A_h
 \left(\frac{(1+15z)^8-1}{2^{32}-1}\right)^h,
\]

where `A_h` counts messages whose eight RS rows have exactly `h`
nonzero 32-bit columns. The rowwise RS code is an `[8,4,5]` MDS code
over the 32-bit alphabet: interpolation acts independently on each of
the eight GF(16) lanes. Thus, with `Q=2^32`,

\[
 A_h=\binom{8}{h}\sum_{j=0}^{h-5}(-1)^j\binom{h}{j}
       (Q^{h-4-j}-1),\qquad 5\le h\le 8.
\]

These counts sum to `Q^4-1=2^128-1`. Exact rational expansion gives
the following cumulative expectations:

| Local packet cutoff d | E[sum_{v<=d} B_v] | -log2 of expectation |
|---:|---:|---:|
| 7 | 2.3226127975e-22 | 71.86667 |
| 8 | 2.9236877748e-20 | 64.89077 |
| 14 | 4.1790252192e-9 | 27.83419 |
| 17 | 4.0901305683e-4 | 11.25557 |
| 18 | 1.6051486513e-2 | 5.96115 |

Markov's inequality therefore proves that a randomly sampled local
profile has packet minimum distance at least 19 with probability at
least `1-0.016051486513 > 0.9839`. In particular, a fixed favorable
eight-scalar tuple exists. This is a statement about the local code
before its packet permutation; it does not prove any global distance.

## What changes and what remains open

Reusing a selected profile removes the independent local-map event
that produced the asymptotic 5/64 upper bound. It leaves the online
RS row circuit, 32-bit multiplications, and route geometry unchanged,
although coefficient reuse might permit a faster kernel. The cost must
be measured rather than inferred.

The finite RS8 certificates average over independent outer maps in
every group. Their product of expected outer measures does not apply
unchanged to one shared selected profile. A proof would need a bound
on the selected profile's realized packet spectrum, followed by a
fresh first-moment calculation with the group shuffles. Packet minimum
distance 19 alone does not control multiplicities or inner cancellation.
Also, sampling one profile and merely hoping it is good leaves a
constant setup-failure probability; an asymptotic 1-o(1) theorem needs
a fixed certified profile or a selection law that excludes bad tuples.

Efficiently certifying one profile is a separate problem. A low-packet-
weight search can be formulated as binary rank tests on the subcode
shortened outside each candidate packet support, exploiting the RS
column structure. No complete search or explicit tuple is recorded yet.
Because the pair-return obstruction survives every fixed profile, such
a search is useful only if a fixed outer improves finite-size setup or
is paired with a redesigned route/inner. It is no longer a standalone
asymptotic-proof gate.
