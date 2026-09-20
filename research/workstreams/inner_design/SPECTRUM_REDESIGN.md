# Removing the all-one expansion state without increasing state size

The cheap mixer leaves a fixed nonzero state unchanged with probability 1/2.
That creates a difficult weighted tail when the expansion image contains the
all-one word. Merely tightening arithmetic cannot remove that state.

`remove_constant.py` constructs a different RM2 subspace with the same
dimensions t=128 and s=19. Its exact audited properties are:

| Property | Original map | New map |
|---|---:|---:|
| Rank / state bits | 19 | 19 |
| Minimum nonzero image weight | 48 | 48 |
| Maximum nonzero image weight | 128 | 80 |
| Minimum kernel weight | 6 | 5 |
| Distinct, nonzero B columns | yes | yes |
| B=A^T and BA=0 | yes | yes |

These map properties alone are not a SPIN distance certificate. The subsequent
all-occupancy computation now passes in outward arithmetic, and this exact map
has been implemented and timed. See [BALANCED_RESULT.md](BALANCED_RESULT.md).
Original-map transvection timings remain a separate experiment.

## Construction and exact audit

The retained RM2 subspace contains the constant function, seven linear
functions, and eleven selected quadratic forms. Extend its quadratic space
by the form with mask `0x93`. Exhaustive enumeration of the resulting
20-dimensional parent confirms minimum nonzero weight 48. A preliminary
rank-two screening found 49 admissible extension cosets; the full enumeration,
not that screen, is the audit of the selected parent.

Project each parent evaluation column by dropping its constant coordinate.
Shift all resulting 19-bit columns by `0x17`. This selects a hyperplane of the
parent that does not contain the constant function. The shift avoids every
singleton and triple syndrome of the projected evaluation columns, excluding
kernel words of weights one and three. Distinct columns exclude weight two.
Even-weight kernel words remain in the parent's kernel, excluding weight four.

The final generator has rank 19 and the following exact nonzero spectrum:

| Weight | Multiplicity |
|---:|---:|
| 48 | 5,166 |
| 56 | 110,288 |
| 64 | 293,455 |
| 72 | 110,128 |
| 80 | 5,250 |

All 2^19 states are enumerated. MacWilliams inversion gives an integral kernel
spectrum with minimum weight five and the correct total mass 2^109. Direct
row inner products check BA=0. The map and source hashes are recorded in
[NO_CONSTANT_MAP.json](NO_CONSTANT_MAP.json).

The lack of an all-one image also strengthens the dense syndrome analysis:
every nonzero difference of image words has weight between 48 and 80.
This restricts the possible intersections in its Fourier bound. It does not
justify reusing the original numeric certificate without recomputation.

## Reproduce

```text
python workstreams/inner_design/remove_constant.py
python workstreams/inner_design/general_occupancies.py --maximum 128 --dense --map workstreams/inner_design/NO_CONSTANT_MAP.json
```

The second command is a full-range binary64 diagnostic using the old cover,
not outward verification. The successful refined covers are retained separately.
Replay their full certificate with:

```text
python workstreams/inner_design/certify_no_constant.py --verify
```
