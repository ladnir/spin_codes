# Eight-bit packet screen

The eight-bit packet candidate passes the occupancy-one bound but leaves a
large gap in the present intermediate-occupancy bound. It is not a certified
replacement for the retained four-bit packet construction.

Each of 512 groups encodes 128 bits with four GF16 RS[16,8] rows. Independent
uniform GL16 maps randomize the 16 aligned symbols. Each symbol supplies two
eight-bit packets. Independent group shuffles assign packets to 32 regions;
independent regional shuffles order the 512 packets within each region.

The inner has 64 output coordinates and a 16-bit persistent state. At each
step it emits `x + A a` and sets the next state to `M a + C x`. Each `M` is
independent uniform GL16. State starts at zero and persists across all steps
and regions. There is no final flush.

Adjacent eight-coordinate blocks of the retained expansion have feedback
rank seven. The candidate swaps coordinates 7/15, 23/31, 39/47, and 55/63.
Apply those swaps to the expansion coordinates and take `C = A^T` afterward.
Every resulting eight-coordinate packet has feedback rank eight. The swaps
preserve the global expansion spectrum and `C A = 0`.

`screen_packet8.py` freshly enumerates all 65,536 expansion images and all
2,040 possible single-packet births. Its exact birth census retains both
input weight and the expansion weight of the next state. For a nonzero
entering state, the refreshed state is uniform nonzero. Thus a single active
packet returns to zero with probability `1/65535`, independently of emission
weight. The other destinations have density at most the same tilted emission
mass relative to uniform nonzero state.

The regional calculation averages all without-replacement packet positions.
The support calculation averages all ordered supports among 32 regions.
Both retain the state between factors. The calculation uses exact expected
RS support counts and sums over every nonzero message in one active group.
Fresh 256-bit outward arithmetic gives **52.58019371235 bits** at the weight
cutoff 13,107, using per-support minima from tilts .00256 and .00512.
Only occupancy one is covered by this result.

`screen_tail.py` extends the local operators to every occupancy from zero to
eight. A Walsh census computes the weighted birth classes. Optional capped
density representations use outward dyadic Walsh inversion. For nonzero
states, dropping the requirement `C x != 0` gives a valid return upper bound.
The subsequent regional and occupancy calculations use logarithmic floating
arithmetic. Their outputs are proposals, not certificate endpoints.

The pointwise outer majorant is unchanged: an active group's expected output
measure is at most `2^256/(2^16-1)^8` times uniform bits. Consequently the
regional mixture uses `Binomial(q,255/256)` active packets. A screen through
occupancy 128, with tilts .00512 through 1.6, leaves a **-4,323.26-bit** bound
at occupancy 128. The best tested tilt there is .4096. This is a gap in the
bound, not evidence of a low-distance code.
An independent 256-bit outward regional replay confirms that component at
**-4,323.26118859182 bits**.

The local rank audit helps explain why the single-packet improvement does not
settle larger occupancies. Among the 28 packet pairs, four restrictions have
rank 11, eight have rank 12, and sixteen have rank 13. More extensive coordinate
interleaving could improve these ranks, but its encoding cost is unmeasured.
A diagnostic uniformly shuffled the 64 coordinates using Python's seed 42.
Its pair ranks were 14 once, 15 eighteen times, and 16 nine times. Nevertheless,
the occupancy-128 floating bound improved only to -4,303.62 bits. Improving
pair ranks alone did not close this bound.

Reproduce from the repository root:

```text
python -B research/workstreams/k16_codesign_100us/packet8/screen_packet8.py
python -B research/workstreams/k16_codesign_100us/packet8/screen_tail.py --q-max 128 --tilts .00512 .01024 .0256 .0512 .1024 .2048 .4096 .8192 1.6
python -B research/workstreams/k16_codesign_100us/packet8/screen_tail.py --outward-q 128 --tilts .4096
python -B -m unittest discover -s research/workstreams/k16_codesign_100us/packet8 -p "test_*.py" -v
```

No benchmark or remote workload was run. The next useful proof experiment
would strengthen the return bound at the failed occupancy before implementing
an encoder. A separate exact-map optimization can reduce outer
packing cost while retaining four-bit packets, so that implementation work
does not depend on closing this candidate's proof.
