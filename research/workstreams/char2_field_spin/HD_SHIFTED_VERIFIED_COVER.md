# Pinned composition of observed directed replays

`theory/hd_shifted_verified_cover.py` checks a manifest of previously
observed 256-bit and 384-bit replay artifacts. It does not infer arithmetic
validity from hashes or from agreement between two saved files.

The caller supplies the expected manifest SHA-256 separately. This trust
anchor authorizes the manifest's artifact pins and its record of observed
fresh replays. The checker verifies every pinned byte sequence before use.
It then recomputes the coverage and probability union.

The cached path checks:

- the GF256 Toeplitz4 construction, independent pair swaps, scalar-code
  scope, source dimension, output length, and forbidden weight;
- exact rational witnesses, explicit outer-density mode, and matching
  rational ceilings at both precisions;
- every integer weight of a potential or sparse band, without interpolation;
- each accepted leaf of a box plan and its deterministic witness identity;
- disjoint assigned bands and complete coverage of weights 1 through 8192;
- the exact count of feasible support profiles and a fresh directed sum.

Potential-band unions are recomputed from the per-weight rational ceilings.
This can be slightly weaker than the saved band union, which sums unrounded
probabilities. The saved band union is not substituted for missing rows.
Restricting a box receipt retains its whole probability bound; mass is never
prorated by volume or by a profile count.

The unsigned tail uses only the established unsigned-envelope transfer or
the route-independent inverse-image bound. Fixed-oriented receipts do not
transfer. The terminal calculation is cheap and is recomputed during cache
checking. The stronger Toeplitz density remains an explicit potential-band
choice; the old sparse and box paths retain their inherited density.

## Manifest interface

The module docstring specifies the schema. Its top-level fields are
`schema_version: 1`, the exact `SWAP_CONSTRUCTION` string,
`trust_basis: externally-pinned-observed-fresh-256-and-384-bit-replays`, and
`entries`. Each entry names an assigned weight range and both receipt pins.
A potential entry also pins its exact proposal and selects `outer_measure`.

Checkpoint exports retain accepted sparse receipts or accepted box leaves.
Dense and tail entries pin their source plans and the replay implementation.
The checker may reconstruct deterministic witness metadata from those
inputs. Its output explicitly records that reconstruction. This step does
not recompute coefficient arithmetic or turn metadata into arithmetic
evidence. Summary-only box caches remain insufficient.
The reconstruction is bound to the reviewed replay implementation's exact
hash; a different same-basename source file cannot authorize it.

For the two historical potential receipts that lack a proposal hash, a
manifest may explicitly authorize `legacy_missing_proposal_pin: true`.
The receipt bytes are still pinned externally, and every exact witness
field must match the pinned proposal. This exception never accepts a
present but incorrect proposal hash.

Run the inexpensive cached composition with:

```text
python -B hd_shifted_verified_cover.py --manifest MANIFEST.json --expected-manifest-sha256 REVIEWED_SHA256 --output VERIFIED.json
```

For a full arithmetic reproduction, add `--fresh`. This invokes the
established `replay_swap_cover` separately at 256 and 384 bits, with an
acceptance threshold of -40 rather than the older -80 search threshold.
This mode can take substantially longer. It uses no saved numerical bound.
Neither mode certifies an incomplete manifest.
The output option creates a new file and refuses to overwrite an existing
verification receipt before beginning work.

## Independent composition check

The reviewed manifest `hd_shifted_closed_cover_manifest.json` has SHA-256
`accb040d153fca5833cf60efdfb93e1f61eeebca80ccd1cb4f1858524eb8d743`.
Its 37 assigned intervals partition every weight from 1 through 8192.
An independent count gives 11465830401 feasible support profiles.

A separate aggregation read the pinned receipt rows and leaves directly,
without using saved band unions or the ledger's summation routine. All
10084 rational ceilings agree between the two precisions. At 512-bit
working precision, max-shifted positive summation encloses their union by

    log2 U = -55.09273850822352338434219682710286647634309918042...

with enclosure radius below `2.18e-150`. Thus this independent sum is
strictly below -40. Receipt, proposal, and selected implementation hashes
also match the authorized manifest. The full geometric cache check is a
separate step; this calculation does not infer arithmetic validity from
the matching hashes or dual ceilings.

The subsequently checked 256-bit and 384-bit ledger outputs agree on all
37 entry summaries and the global rational ceiling
`-55092738508223/1000000000000`. They contain 6059 accepted boxes and zero
unresolved feasible profiles. The 1671 remaining geometric boxes contain
no feasible profile; they are not coverage gaps. Both outputs retain the
explicit trusted-cache qualification and the ideal-setup scope above.

## Final-gap reuse: weights 4068 through 4095

The earlier plan for 4096..4351 retained a useful FLOAT witness bank.
Each accepted-box screen value determines its witness's affine intercept
after subtracting the known histogram term and profile-count term.
The recovered values agree across original leaves within 1.9e-9 bits.
This recovery is a proposal mechanism, not a numerical certificate.

Extending the r=4096 boundary rectangles downward and subdividing failing
boxes required no new operator evaluations or witness fits. Over the full
4051..4095 window, only the singleton profiles (4052,2026,0) and
(4054,2027,0) remained unproved by that frozen bank. The retained proposal
therefore starts at 4068; a separate whole-weight bank handles the prefix.

The resulting 4068..4095 plan has 981 disjoint leaves and 158 witnesses.
Fresh replays at both precisions cover all 116696664 feasible profiles.
Every leaf's rational ceiling agrees, and both directed unions satisfy

    log2 U <= -15988140130039 / 200000000000
            = -79.940700650195.

The complete artifacts are:

- `hd_shifted_swap_4068_4095_frozen_plan3.json`:
  `5d1b62ac6c49cf5c81e04db8b518d3e3fb5e6589cbcc33290661eaf56718ffbe`.
- `hd_shifted_swap_4068_4095_directed256.json`:
  `cfde5173152678d5d84c51cd361ab9f416d8d9ca0608d866fbba8eb1028416d1`.
- `hd_shifted_swap_4068_4095_directed384.json`:
  `e09b697b4a8af6a8d5053643e6d972a5d8479018424483137097fbfc5059d749`.

`hd_shifted_4068_replay.py` reproduces both arithmetic passes and checks
every leaf ceiling. It refuses to overwrite frozen outputs. The plan's
directed replay uses threshold -60; an old default of -80 is an unnecessary
per-leaf acceptance restriction, not the final 40-bit target.

This band alone is not a whole-code certificate. The
[completed manifest](HD_SHIFTED_CLOSED_CERTIFICATE.md) now combines it
with the remaining observed exports. Both aggregation precisions give
the exact whole-code ceiling `-55092738508223/1000000000000` and zero
unresolved feasible profiles. The authorized manifest SHA-256 is
`accb040d153fca5833cf60efdfb93e1f61eeebca80ccd1cb4f1858524eb8d743`.

## Independent outer-density check

For a nonzero target group supported on h columns, invertibility forces
the four RS[64,32] source rows to vanish outside those columns. If h<33,
no source exists. Otherwise each shortened row has dimension h-32.
Conditioning on mixers at any h-32 target columns determines at most one
source group. Exactly 32 independent column mixers remain.

Each remaining nonzero input hits its prescribed output with probability
at most `1/(255*256^3)`. A zero input makes the event impossible. Thus the
expected counting measure at each nonzero output is at most that atom
bound to the 32nd power. Relative to uniform GF256^256, the density is
`256^160/255^32`.

The atom bound uses surjectivity of the unconditioned Toeplitz action on
each fixed nonzero input and the nonsingular fraction `1-1/q`, proved by
[Kaltofen and Lobo](https://kaltofen.math.ncsu.edu/bibliography/96/KaLo96_issac.pdf).

One may fix the independently sampled inner route, pair swaps, and refresh
maps before this argument. They do not change the outer counting measure.
Independent group setups give the product density bound. The scalar
GF256 theorem also applies to its scalar extension: a nonzero component
codeword's support is contained in the complete output support. No lane
union or additional projective factor is required.

These claims concern the ideal setup law. The native Mersenne Twister
setup remains a reproducible research implementation, not independent
uniform randomness for a deterministic-seed distance theorem.
