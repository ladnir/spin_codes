# Physical t32/s8 screen

The t32/s8 candidates pass the occupancy-one screen, but the present
intermediate-occupancy bounds do not close the 10% distance target.
No encoder was implemented or benchmarked for these candidates.

The outer and routing remain the retained K16 RS16 construction: 512 groups,
128 message bits per group, 64 four-bit packets per group, and 64 regions.
The new physical inner emits 32 coordinates per step and retains eight state
coordinates. Each update is independently uniform GL8. Initial state is zero,
state persists across steps and regions, and there is no final flush.

The expansion rows begin with the constant and five coordinate functions on
the 32 five-bit points. The first screen appends `x0*x1` and one other
quadratic monomial. Feedback is the transpose of the expansion. All maps
have full state rank and rank-four packet restrictions. Products of two
rows have degree at most four, so their sums over five variables vanish;
therefore `C A = 0`.

Fresh 256-bit outward occupancy-one bounds give 42.17301995671 bits when
the second quadratic shares a variable with `(0,1)`. Disjoint second pairs
give 42.19173423094 bits. These bounds include all choices of the active
group and all its nonzero messages, at output-weight cutoff 13,107.

The tail proposal composes four physical steps into one 128-bit macro.
Each composition uses exact hypergeometric placement and retains state.
This is a composition of four independently refreshed 32-bit steps, not
a single 128-bit inner step. The regional calculation uses 16 macros and
retains state through all 64 regions. Its arithmetic is logarithmic floating
point; proposal values are not certificate endpoints.

At occupancy 84, refined tilts .24 through .36 give -1,156.10 bits for
the disjoint monomial pair and -873.04 bits for the shared-variable pair.
These are failures of the present upper bound, not demonstrated bad codes.

A stronger eight-row expansion uses the two quadratic forms
`x0*x1 + x2*x3` and `x0*x2 + x1*x4`. Its complete fresh weight spectrum is
`{0:1, 12:48, 16:158, 20:48, 32:1}`. Its occupancy-one bound reaches
42.19241505595 bits. Nevertheless, occupancy 84 remains at -709.74 bits
on the tested grid, with best tilt .33.

Reproduction:

```text
python -B research/workstreams/k16_codesign_100us/t32/screen_t32.py
python -B research/workstreams/k16_codesign_100us/t32/screen_t32.py --tail --pair 2 3 --q-min 84 --q-max 84 --tilts .24 .27 .30 .33 .36
python -B research/workstreams/k16_codesign_100us/t32/screen_t32.py --rank4
python -B research/workstreams/k16_codesign_100us/t32/screen_t32.py --tail --rank4 --q-min 84 --q-max 84 --tilts .24 .27 .30 .33 .36
```

The search stopped after this bounded negative screen. The next useful step
is the separate paired-quadratic t64/s16 candidate, whose stronger state
refresh retains the original 16-bit return denominator.
