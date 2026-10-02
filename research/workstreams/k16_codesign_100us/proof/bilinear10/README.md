# Bilinear t64/s10 bounded screen

The explicit bilinear map passes all algebra checks and gives a fresh q=1
bound of 53.65471457271895 bits. The tested middle-occupancy bounds do not
certify the target. No full certificate or low-weight counterexample is
claimed. This screen changes neither pinned proof modules nor kernel code.

## Construction and checked algebra

The 64 evaluation points are the bit vectors `z=(z0,...,z5)`, in integer
order. Write `X=z0+z2*a+z4*a^2` and `Y=z1+z3*a+z5*a^2` in GF8, where
`a^3+a+1=0`. The ten expansion rows are the constant, six linear coordinate
functions, and the three coordinates of `X*Y`. Explicitly, the quadratics
are:

```
q0 = z0*z1 + z2*z5 + z3*z4
q1 = z0*z3 + z1*z2 + z2*z5 + z3*z4 + z4*z5
q2 = z0*z5 + z2*z3 + z1*z4 + z4*z5
```

All sums are binary. The feedback is `C=A^T`. Every physical step emits
`x+A*s`, then updates the state to `M*s+C*x` for a fresh independent uniform
`M` in GL10. Initial state is zero; it persists across steps and regions.
There is no final flush. The small RS16 outer, 512 groups, 64 four-bit packet
regions, and ideal independent routing law are unchanged.

Fresh enumeration verifies:

- expansion and feedback rank 10, and `CA=0`;
- rank 4 on each of the sixteen consecutive four-coordinate packets;
- polar rank 6 for each of the seven nonzero quadratic combinations;
- an injective census of all 1024 states with spectrum
  `{0:1, 28:448, 32:126, 36:448, 64:1}`.

The proof wrapper uses two chronological physical t64 steps, with separate
fresh updates, per t128 macro. Direct regional placement of 32 physical
steps agrees with 16 such macros for occupancies 0, 1, and 2.

## Partial numerical results

`refined-screen-p256.json` contains the fresh 256-bit q=1 calculation. It
uses exact RS expected shell counts and separate minima over valid tilts
for each support size. Its resulting margin is:

```
53.65471457271895218078049693461055786442703488283882765119966355432590544007 bits
```

The same record contains floating critical-tail proposals. Selected
middle-occupancy objectives were then recomputed with 256-bit outward
arithmetic in `middle-outward-p256.json`:

| Active groups q | Best tested tilt | Outward margin, bits |
| ---: | ---: | ---: |
| 80 | 0.30 | -860.28388725303885 |
| 96 | 0.35 | -942.75866916505610 |
| 128 | 0.50 | -977.05905884125104 |

The exact trial set was `0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.80, 1.00`.
A negative margin means this upper bound on the expected bad-word count
exceeds one. It is vacuous for certification; it does not establish a bad
codeword or prove that no other proof method can succeed.

The initial `critical-screen-p256.json` is an interrupted coarse-screen
checkpoint. Its final floating trial at 0.8192 triggered the existing
underflow guard. The earlier entries are retained, but the two completed
records above are the primary results. High-tilt q=1 entries in the outward
tail record are not the optimized q=1 result.

## Reproduction

From the repository root:

```
C:/Python314/python.exe -B research/workstreams/k16_codesign_100us/proof/bilinear10/test_bilinear.py -v
C:/Python314/python.exe -B research/workstreams/k16_codesign_100us/proof/bilinear10/screen.py --tilts .003 .004 .005 .006 .008 .3 .35 .4 .45 .5 .6 --tail-q 64 80 96 112 128 --output NEW-REFINED-RECEIPT.json
C:/Python314/python.exe -B research/workstreams/k16_codesign_100us/proof/bilinear10/screen.py --tilts .3 .35 .4 .45 .5 .6 .8 1.0 --tail-q --outward-q 80 96 128 --output NEW-OUTWARD-RECEIPT.json
```

All four tests pass. Every screen requires a fresh output path and records
source hashes. No benchmark or remote workload was run.

Next, retain the certified s15/s16 alternatives. Do not prioritize a new
s10 kernel unless a stronger middle-occupancy proof route offers a concrete
reason to overcome the observed gap.
