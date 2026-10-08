# Smaller packets: a cheaper GF256 inner candidate

The subsequent [whole-code checkpoint](GF8_PACKET_CERTIFICATE.md) closes
the repaired cubic width8/state4 profile with GRS[128,83] at K20. The
screens below preserve the design search and its higher-rate gaps.

The width-eight, four-state inner is the most promising arithmetic
candidate in this bounded screen. A repaired cubic map now improves its
return-weight bounds and admits a nine-GFNI paired-map circuit. The
conservative GRS[128,83] profile has complete directed middle and dense
component certificates, recorded [separately](GF8_PACKET_MIDDLE_CERTIFICATE.md).
The earlier normalized-Cauchy screens remain below as controls.

This file records middle-weight screens, not a complete code certificate.
The separate sparse-row screen still has an unresolved class for
dimension 84. Positive samples do not supply interval coverage, and an
instruction count does not establish encoding time.

## Map interface and transferred argument

Fix emission `A : F^s -> F^t` and injection `C : F^t -> F^s`, where
`F=GF256`. Both matrices must be superregular: every square minor is
nonzero. Require `CA=0`, and choose a fixed nonzero scalar eta. The two
maps need not be transposes.

At each packet, independent nonzero input-coordinate multipliers produce
X. Independently sampled nonzero state-coordinate multipliers produce R
from the incoming state. The packet computes

\[
Y=X+\eta AR,\qquad Q'=R+CX.
\]

The source, the one regional routing stage, and all setup independence
requirements remain those of the GF256 construction. The new map family
is recorded separately by the map investigation. The following transfer
does not assume that a diagonal samples a uniform invertible matrix.

The generic return-output minimum is

\[
\max(1,t-s+1-w)
\]

for an input packet of weight w returning a nonzero incoming state to
zero. Superregular emission gives at least `t-s+1` nonzero emitted
coordinates; the input can cancel at most w. The minimum is at least one
because `CA=0` makes the full packet block invertible. Independently
proved map-specific minima may strengthen this bound.

For the joint return/output-zero cover, fix an incoming actual state
support of size h. Let i selected output zeros lie on active input
coordinates, and let z lie outside the input support. Then `z<h` by
superregularity. A possible return requires `h+w>s`; the return equations
have rank s. Selected identity pivots and inactive emission rows have
joint rank at least `i+z`. Finally, invertibility gives joint rank at least
`h+w-v`, where `v=t-i-z` counts unselected output coordinates. Therefore
one may charge

\[
\max(s,i+z,w-t+i+2z+1).
\]

These are three lower bounds on the same equation matrix, not independent
events. The implementation takes the minimum of the complete pure-return
and joint-return covers; it does not multiply them. State zero-flat bases
are counted before the current independent diagonals are exposed.

The source comparison uses the existing modal-conditioning factor. The
feasible-support polynomial is used only when its source tilt satisfies
the previously checked convexity condition. Otherwise the established
relaxed polynomial is retained.

## Width-eight, state-four middle frontier

The normalized half-block maps have return-weight bounds
`[5,5,4,2,1,1,1,1]` for input weights one through eight. With those bounds,
the K20 screens give:

| Outer dimension | Requested rate | Input density | Sampled margin |
|---:|---:|---:|---:|
| 88 | 128/187 = 68.449% | .10 | -23243.597 bits |
| 88 | 128/187 | .14 | -25323.593 bits |
| 86 | 128/191 = 67.016% | .12 | -10122.441 bits |
| 85 | 128/193 = 66.321% | .03 | +64.767 bits |
| 85 | 128/193 | .10 | +499.333 bits |
| 85 | 128/193 | .12 | -237.021 bits |
| 84 | 32/49 = 65.306% | .025 | +1241.160 bits |
| 84 | 32/49 | .03 | +2390.980 bits |
| 84 | 32/49 | .10 | +8584.114 bits |
| 84 | 32/49 | .11 | +8447.542 bits |
| 84 | 32/49 | .12 | +9796.942 bits |
| 84 | 32/49 | .24 | +486678.523 bits |

All these optimizations converged. The twenty dimension-84 densities were
`.025,.03,.035,.04,.05,.06,.07,.08,.09,.10,.11,.12,.13,.14,.15,.16,.18,.20,.22,.24`.
Every sample was positive. Dimension 80, a 5/8 nominal-rate control, also
passed the six-point coarse grid.

The dimension-85 failure is narrow enough that a sharper transfer could
matter. However, the sparse-row obstruction at dimension 84 must be
addressed before pursuing that higher-rate middle refinement. These
floating upper-bound failures are not evidence of actual low-distance
words.

## Wider-packet controls

Widening the packet without increasing the state does not help the
current proof. Width 32 with state eight and dimension 87 has a
converged margin of -49566 bits at density `.10`; dimension 88 gives
-56848 bits. Both return to positive margins at higher densities.

Width 32 with state sixteen and the analogous half-block return bound
passes six middle samples at dimension 88. Its margins range from about
3148 bits at `.03` to 310455 bits at `.24`. The `.03` optimization reached
its fifty-iteration limit; its positive candidate witness is usable for
further screening, but optimality is not asserted. The other five runs
converged. This profile has much larger dense map cost, so the extra
proof slack does not yet justify an implementation.

The generic width-32/state-16 cover at dimension 96 fails by approximately
14930, 41277, and 71568 bits at densities `.03`, `.10`, and `.14`. The
packet-width sweep did not identify a useful 3/4-rate candidate.

## Cost gate and the older asymmetric map

A straightforward four-record SIMD packing costs 19 GFNI512 instructions
per eight-symbol packet: sixteen for the two maps and three for input
and state diagonals. Eta can be absorbed into the fixed emission map.
The analogous width-16/state-8 dense path costs 70 per sixteen symbols.
At dimension 84, the extra output records reduce the net inner-instruction
saving to about 43%, still a substantial potential gain. Loads, routing,
packing, and outer work remain unmeasured for this changed construction.

Both maps may admit a cheaper shared circulant circuit. That is a separate
operation-count and scheduling question; it is not assumed in these
proof screens.

The older embedded-GF4 asymmetric map has a singleton return weight of
eight, but its injection matrix fails several two-by-two minors. Its
singleton bound alone cannot justify using the present transfer.

The repaired map keeps `A[u,j]=(u+v_j)^3+1` and replaces the injection by
`C[j,u]=(u+v_j)^3+3`. The map audit checks all 7904 square coordinate minors
of A and C across eight phases, verifies `CA=0`, and obtains exact return
minima `[8,6,5,4,3,2,2,1]`. Both matrices retain five-product circulant
circuits, giving a paired-map arithmetic count of nine GFNI512
instructions per packet including the two diagonal stages and eta.
The count includes an explicit record-transpose schedule, but no native
timing is claimed.

These stronger return bounds materially improve the middle frontier:

| Outer dimension | Input density | Repaired-map sampled margin |
|---:|---:|---:|
| 85 | .03 | +1545.345 bits |
| 85 | .10 | +23938.581 bits |
| 85 | .12 | +28897.450 bits |
| 86 | .035 | +234.719 bits |
| 86 | .04 | +1375.620 bits |
| 87 | .045 | -857.008 bits |
| 87 | .05 | +59.169 bits |
| 87 | .055 | +1006.968 bits |
| 87 | .10 | +7776.817 bits |
| 88 | .10 | -223.080 bits |
| 88 | .14 | -922.339 bits |

All reported optimizations converge. Dimension 85 passes the nine-point
grid from `.03` through `.24`. Dimensions 86 and 87 need progressively
larger sparse-class coverage before these middle bounds can help. Because
dimension 84 still fails a sparse sampled class at the 40-bit target, the
first directed campaign conservatively targets dimension 83.

## Reproduction and checks

```text
python -B research/workstreams/char2_field_spin/theory/gf8_packet_width_screen.py --width 8 --state 4 --dimension 84 85 86 --density .03 .06 .10 .11 .12 .13 .14 .20 .24 --return-minima 5 5 4 2 1 1 1 1
python -B research/workstreams/char2_field_spin/theory/gf8_packet_width_screen.py --width 8 --state 4 --dimension 83 85 86 87 88 --density .03 .06 .10 .11 .12 .13 .14 .20 .24 --return-minima 8 6 5 4 3 2 2 1
python -B research/workstreams/char2_field_spin/theory/gf8_packet_width_screen.py --width 32 --state 8 --dimension 87 88
python -B research/workstreams/char2_field_spin/theory/gf8_packet_width_screen.py --width 32 --state 16 --dimension 88 --return-minima 17 17 17 17 17 17 17 17 16 14 12 10 8 6 4 2
python -B research/workstreams/char2_field_spin/theory/gf8_packet_width_tests.py
```

All five tests pass. They include exact agreement with the previously
certified width-16/state-8 local cover, exact local-sum checks of the
batched floating aggregate for five geometries, and 512 numerical joint
rank checks. The rank checks supplement the preceding argument; they do
not replace its quantifiers. The tests also check `CA=0` for each explicit
polynomial map geometry and finite source/route accounting.

The next decision should combine a complete sparse/middle/dense proof
screen for a conservative outer dimension with the shared-circuit cost
audit. Native implementation work is premature until those two gates
agree on one profile.
