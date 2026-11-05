# Sparse supports with independently swapped pairs

Exact coefficient extraction removes the positive-alias penalty for small
active-group supports. The implementation covers every support at each
requested weight, not selected matching profiles. It changes the proof
calculation, not the encoder.

The construction is `toeplitz-gf256-h4-alternating-swapped-pairs-k20-rate-half`.
It has 8192 outer groups, input dimension 1048576, and output length 2097152
over GF(256). The forbidden output weight is 734003. The setup uses the
independent outer and refresh maps of the existing Toeplitz certificate,
independent uniform regional pair permutations, and fresh independent fair
swaps within every routed pair. This is an ideal-randomness statement, not
a guarantee supplied by the native reproducible PRNG.

## Regional coefficients

Fix a nonempty indexed support of weight `r`. Let `a,b` count its full
pairs under the even and odd matchings. In a region with `f` full pairs,
there are `s=r-2f` singleton pairs and `4096-s-f` empty pairs.

Fix rational witnesses `0<z<1` and `0<alpha<=1`. For a local arrangement
`gamma` of 32 pairs, let `A_gamma` be its established nonnegative 2-by-2
packet cap before the entrywise power. It includes the diagonal density
factor and the improved empty-packet quiet cap. Define

    B[s,f] = 2^(-s) sum_gamma A_gamma^[alpha],

where the sum contains every labeled arrangement with `s` singleton and
`f` full pairs, including all singleton orientations. Each orientation has
probability `2^(-s)` after the independent swaps. The power precedes the
average; the implementation never raises the averaged raw cap to `alpha`.

Each half-step contains 16 pairs. Its oriented counts `(u,v,f)` have exact
multiplicity `binom(16,u) binom(16-u,v) binom(16-u-v,f)` and occupancies
`(u+f,v+f)`. Combining two halves gives the four block occupancies.
Exchanging the halves supplies multiplicity two only when they differ.
As an independent count check,

    sum_gamma 1 = 2^s binom(32,s) binom(32-s,f).

Retain the existing good-good improvement: good steps have at most one
singleton. Write `G,D` for the good and remaining parts of `B`. Compose
the raw matrices of two good steps before applying `alpha`, obtaining `C`.
Good caps are invariant under their singleton's orientation, so `C` is
unchanged by the swap average. The ordered two-step polynomial is

    T = C + GD + DG + DD.

All coefficients are nonnegative. The regional majorant is

    R(r,f) = [X^(r-2f) Y^f] T^64
             / (binom(4096,r-2f) binom(4096-r+2f,f)).

The denominator counts unsigned arrangements under a uniform pair
permutation. Coefficient extraction therefore retains sampling without
replacement. It does not make packet occupancies independent. At `alpha=1`,
`T=B^2`, and the implementation directly computes `B^128`.

## Truncation and directed arithmetic

To obtain every coefficient of weight at most `R`, keep only bidegrees
`s+2f<=R` after each multiplication. Discarded terms cannot later enter this
set because every degree is nonnegative. This loses no target coefficient.

The implementation encodes `(s,f)` by `s+(2R+1)f` in an Arb polynomial.
Both operands are already truncated. Hence their summed singleton degree
is at most `2R`, and no carry can turn excess singleton degree into full-pair
degree. After each ordered matrix product, the implementation retains only
the valid bidegrees and takes Arb upper endpoints. Positivity preserves
the coefficientwise majorant. It introduces no positive coefficient aliases.

The rational marker `R/4096` weights both variables to control coefficient
magnitudes. Its exact target monomial cancels in the regional normalizer.
Packet caps are exact rational formulas, rounded upward to 60 mantissa
bits before Arb exponentiation. FLOAT arithmetic proposes witnesses only.

## Support union and fixed witnesses

Let `N(r,a,b)` be the indexed-cycle support-count upper bound from
`hd_shifted_pairs_count.profile_count_upper`. It already includes every
even/odd split and does not identify rotations. Put

    beta = 256^256 / 255^128,  e0=(1,0),  one=(1,1)^T.

The established outer comparison and projective normalization give

    U_r <= (beta^r / (255 z^734003))^alpha
           sum_(a,b) N(r,a,b) e0 (R(r,a) R(r,b))^128 one.

The region order is unchanged. Counts stay outside the fractional power;
the single projective factor stays inside it. The sum includes every
nonzero-count matching profile for this `r`, and `r=0` is excluded.

The module freezes this rational schedule:

| Active groups | Weight marker `z` | `alpha` |
| --- | --- | --- |
| 1..8 | `1-r/1200` | `1` |
| 9..17 | `1-r/1200` | `2/5` |
| 18..32 | `1-r/1400` | `1/5` |
| 33..64 | `1-r/1400` | `1/10` |

Both 256-bit and 384-bit directed replays cover all 64 weights. Every
per-weight rational ceiling agrees across the two precisions. The weights
are disjoint, and their directed aggregate is

    log2 U_[1,64] <= -107393558356229 / 1000000000000
                   = -107.393558356229.

This domain contains 23936 feasible unsigned triples `(r,a,b)`, including
every even/odd split. Selected per-weight ceilings are:

| `r` | Directed log2 upper ceiling |
| --- | ---: |
| 1 | -107.393558356230 |
| 2 | -222.005473608528 |
| 8 | -986.579921280187 |
| 16 | -687.951888330612 |
| 32 | -690.057289231772 |
| 48 | -362.432638452113 |
| 64 | -398.943133093601 |

The complete dual-precision pass took 193.9 seconds locally. This is proof
replay time, not an encoder benchmark. Eight focused tests pass.

The frozen `weight_receipt`, `family_witness`, and `family_receipt` APIs
remain scoped to `1<=r<=64`. The batched API below has a separate resource
guard. Every receipt identifies the swap construction and sets
`whole_code_certificate=False`; neither API claims whole-code coverage.

## Reusing regional coefficients across weights

Fix one rational witness `(z,alpha)`, a coefficient cutoff `R<=256`, and
a nonempty contiguous requested interval `[first,last]` inside `[1,R]`.
`batch_receipt` computes `regional_directed(z,alpha,R)` once. For every
requested `r`, it then evaluates the complete support union `U_r` above,
including its own `beta^r` factor and every nonzero-count `(a,b)` profile.
Finally it sums the directed per-weight rational ceilings. The intervals
of different receipts must be disjoint before adding them to a coverage
ledger without adjusting overlaps.

The shared computation is valid because the retained bidegrees form a
down-set. A term discarded at cutoff `R` cannot contribute to any requested
coefficient. Increasing the cutoff from `r` to `R` changes neither that
exact coefficient nor its normalization. The marker changes to `R/4096`,
but each coefficient and its normalizer contain the same factor
`(R/4096)^(r-f)`, which cancels. Directed rounding may change numerical
slack; it does not change the claimed majorant. No interpolation between
weights or additional independence assumption is used.

Only `[first,last]` is claimed in the receipt, even if `R` is larger.
The batch API rejects zero weight, reversed intervals, a target beyond
the cutoff, and `R>256`. The upper cutoff is a bounded resource policy,
not a claim that the coefficient argument fails above 256. The frozen
1..64 witness schedule is unchanged.

The separate literal `BATCH_SCHEDULE` closes 65..128 with four witnesses.
All batches use `R=128`, coefficient marker `1/32`, and `alpha=1/10`:

| Requested weights | `z` | Feasible triples | Directed union log2 ceiling |
| --- | --- | ---: | ---: |
| 65..72 | `333/350` | 9808 | -374.252761257220 |
| 73..90 | `659/700` | 31116 | -224.654677819326 |
| 91..110 | `13/14` | 52180 | -270.555735995112 |
| 111..128 | `32/35` | 66000 | -283.900775505865 |

The per-weight ceilings and all four batch ceilings agree at 256 and 384
bits. The directed sum of these disjoint batches is

    log2 U_[65,128] <= -8986187112773 / 40000000000
                     = -224.654677819325.

This covers all 159104 feasible triples in the weight band, not samples of
the even/odd matching profiles. Both precision passes together took 127.51
seconds locally; the longest single batch took 18.10 seconds. These are
proof-replay times, not encoder benchmarks. Independent read-only review
confirmed down-set truncation, carry-free encoding, marker cancellation,
per-weight normalization, and the claimed interval scope.

### Cutoff-256 extension and boundary repair

The extension retains `alpha=1/10` and uses `R=256`, coefficient marker
`1/16`, and four predetermined weight-marker centers 144, 176, 208, and
240 in the formula `z=1-center/1400`. The initial equal-width partition
failed the 40-bit target only in its last band, 225..256: the directed
256-bit union ceiling was -24.099841133022. Weight 225 dominated, with
ceiling -24.099841133583; weight 226 already had ceiling -55.360239494309.
This was an inadequate assignment of weights to witnesses, not evidence
that the construction fails. The prescribed check stopped before replaying
that unsuccessful partition at 384 bits.

The bounded repair assigns weight 225 to the preceding witness `149/175`,
where its per-weight ceiling is -492.099018280752, and leaves 226..256 with
`29/35`. No rational witness changes. The frozen `BATCH_SCHEDULE_256` is:

| Requested weights | `z` | Feasible triples | Directed union log2 ceiling |
| --- | --- | ---: | ---: |
| 129..160 | `157/175` | 171200 | -113.557692540288 |
| 161..192 | `153/175` | 254144 | -111.971464058234 |
| 193..225 | `149/175` | 366241 | -81.942941581917 |
| 226..256 | `29/35` | 456415 | -55.360239487290 |

All 128 per-weight ceilings and all four repaired batch ceilings agree
exactly between the 256-bit calculations and complete 384-bit replays.
This includes every unchanged per-weight ceiling in 226..256. The new
domain contains all 1248000 feasible triples, and its directed union is

    log2 U_[129,256] <= -11072047894587 / 200000000000
                      = -55.360239472935.

Adding the disjoint preserved 1..64 and 65..128 receipts gives

    log2 U_[1,256] <= -27680119736467 / 500000000000
                    = -55.360239472934.

The combined sparse domain contains 1431040 feasible triples. This meets
the 40-bit target for the stated family; weights above 256 remain outside
this result. It is not a whole-code certificate.

The initial four 256-bit runs took 368.74 seconds, the one-weight repair
check took 71.87 seconds, and the four full 384-bit replays took 287.42
seconds. These are proof-replay times. Independent read-only review also
confirmed that the cutoff-256 extension preserves the down-set, radix,
and marker arguments.

## Reproduction and checks

From `research/workstreams/char2_field_spin/theory`:

```text
python -B -m unittest hd_shifted_swap_sparse_tests -v
python -B hd_shifted_swap_sparse.py --family 64 --precision 256
python -B hd_shifted_swap_sparse.py --family 64 --precision 384
python -B hd_shifted_swap_sparse.py --batch 65 72 --maximum 128 --z 333/350 --alpha 1/10 --precision 256
python -B hd_shifted_swap_sparse.py --batch 73 90 --maximum 128 --z 659/700 --alpha 1/10 --precision 256
python -B hd_shifted_swap_sparse.py --batch 91 110 --maximum 128 --z 13/14 --alpha 1/10 --precision 256
python -B hd_shifted_swap_sparse.py --batch 111 128 --maximum 128 --z 32/35 --alpha 1/10 --precision 256
python -B hd_shifted_swap_sparse.py --batch 129 160 --maximum 256 --z 157/175 --alpha 1/10 --precision 256
python -B hd_shifted_swap_sparse.py --batch 161 192 --maximum 256 --z 153/175 --alpha 1/10 --precision 256
python -B hd_shifted_swap_sparse.py --batch 193 225 --maximum 256 --z 149/175 --alpha 1/10 --precision 256
python -B hd_shifted_swap_sparse.py --batch 226 256 --maximum 256 --z 29/35 --alpha 1/10 --precision 256
```

Repeat each batch command with `--precision 384` for the second directed
pass. Programmatic callers can iterate `BATCH_SCHEDULE` and invoke
`batch_receipt(first,last,z,alpha,maximum=128,precision=...)`, or iterate
`BATCH_SCHEDULE_256` with `maximum=256`. Every batch returns its full list
of per-weight receipts as well as the range union.

Focused tests compare small coefficients with an independently implemented
oriented sum, check exact geometry counts, and test noncommuting products.
They also test an explicit radix-carry hazard, regional multinomial and
marker cancellation, indexed support counts, scope rejection, and directed
precision agreement. An independent read-only audit found no gap in the
truncation, chronology, normalization, or support union.

The batch tests also compare directed per-weight ceilings with the original
single-weight API at different cutoffs and markers, verify exactly one
regional call, and check that receipt scope excludes unrequested weights.
Both fixed batch schedules are tested for disjoint and complete coverage.
A cutoff-256 boundary test checks the degree index set and
requested-range metadata without repeating the large polynomial transfer.
Thirteen focused tests pass.

The maximum univariate degrees are 4128, 16448, and 65664 for `R=64,128,256`.
At `R>=128`, the good-good geometry set has reached its fixed maximum;
only coefficient transfer and the final profile sums keep growing.
These sizes suggest batching narrow intervals, not replaying every weight
with a distinct witness. The construction and its ideal setup law do not
change when regional calculations are reused.

The next proof step should integrate these disjoint sparse receipts into
the construction-specific ledger, then use the scalable run/potential
method for larger weights. The degree growth and measured replay cost do
not justify an automatic increase of the sparse cutoff beyond 256.
