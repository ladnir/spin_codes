# Feedback bounds from compressed character sums

The existing exact feedback census covers six occupied windows. Above
that cutoff, several transitions use weaker bounds. This candidate
computes exact zero-feedback and expansion-class counts through ten
windows without reconstructing every feedback distribution. Its atom
bound is an upper bound, not an exact maximum.

Fix a shape h of nonzero four-bit packets. Its packets occupy distinct
windows, with the same uniform window and lane-mask choices as the
construction. Let n_h(t) count choices whose feedback is state t, and
let D_h be the total number of choices. There are S = 2^19 states.
Write the unnormalized Walsh transform as

    n_hat_h(a) = sum_t n_h(t) (-1)^(a dot t).

The existing character-polynomial routine computes these coefficients.
It groups characters by their histograms of four-bit restrictions across
the 32 windows. All characters in a group H have the same coefficient
f_h(H). Let m_H be that group's size.

## Exact counts and a maximum-atom bound

Walsh inversion at zero gives the exact integer

    n_h(0) = (1/S) sum_H m_H f_h(H).

For an expansion weight v, define A_v = {t : weight(E(t)) = v}.
Only nonzero expansion weights are used, so source zero is absent.
Let g_v(H) be the sum, over characters in H, of the Walsh transform of
the indicator of A_v. Character orthogonality gives the exact count

    sum_(t in A_v) n_h(t) = (1/S) sum_H f_h(H) g_v(H).

Finally, dropping signs in inversion gives, for every target t,

    n_h(t) <= (1/S) sum_H m_H |f_h(H)|.

Round the right side upward to an integer and cap it by D_h. This
bounds both zero and nonzero atoms. Reusing it as a bound for nonzero
atoms is conservative; it must not be described as their exact peak.

The spectral atom bound is weak for small shapes. The adapter retains
the exact peaks through six windows. For longer shapes, it intersects
the spectral bound with the conditioned-subset bound from
[MASS_DENSITY.md](MASS_DENSITY.md). The integer ceiling is applied after
multiplying the conditioned probability bound by D_h. The exact zero
and expansion-class counts are retained unchanged.

## Arithmetic and local integration

`candidates/spectral_feedback.py` uses the guarded signed-int64 character
polynomials through degree ten. The largest counting denominator is
1438322789449728000, below 2^63. Each partial polynomial update is a
signed count of a subset of the full shape assignments, so the same
denominator bounds its magnitude.

An indicator transform has magnitude at most S. Compressing at most S
such values therefore fits within S^2 = 2^38. The compression uses
integer indexed addition, not floating-point weighted bin counting.
Every subsequent dot product uses unbounded Python integers. The routine
checks exact divisibility, nonnegativity, the trivial character, and that
the zero and nonzero-class counts sum to D_h.

The small-state tests compare these formulas with complete distributions,
including a case whose dot products exceed signed-int64 range. The
production-map run through ten windows passes independent comparisons
with all fourteen one- and two-window shapes. The integration run also
passes comparisons with all 209 shapes in the existing full census
through six windows.

The existing full-feedback refinement uses zero counts exactly, including
in complements such as 1-Pr[BX=0]. It uses peaks only as upper bounds,
and expansion-class counts in positive moment sums. Thus records with
exact zero/class counts and certified peak bounds satisfy that interface.
The adapter keeps all input penalties paired with their reciprocal outer
counts and takes maxima over complete shape sets.

`candidates/mass_verify.py --spectral-feedback 10` applies these bounds
before the mass-based density and zero-return alternatives. Its default
mode still covers only selected homogeneous support vectors. No complete
occupancy certificate currently uses this extension.

```sh
python -B research/workstreams/permutation_locality/independent_rows/candidates/spectral_feedback.py --maximum 10
python -B research/workstreams/permutation_locality/independent_rows/candidates/mass_verify.py --groups 80 --supports 192 200 208 --tilts .064 .068 --spectral-feedback 10
```

## Selected occupancy-eighty results

The command above passes its local integer checks and selected-point
outward replay at 192-bit precision. Each row concerns only the event
where all 80 active groups have the stated union support. It includes
all compatible outer words and group locations. Logarithmic uppers are
rounded upward.

| Union support per group | Selected tilt | Outward log2 upper |
|---|---:|---:|
| 192 | .068 | +223.496836 |
| 200 | .068 | +256.803777 |
| 208 | .068 | +158.901375 |

All three remain vacuous as probability bounds. Relative to the best
previous mass-based results, they reduce the logarithmic uppers by about
282, 333, and 374 bits. This is evidence that the six-window cutoff
caused substantial loss, not a certificate for these events or the
complete occupancy.

Combining this extension with the fixed-mixture optimizer at support 200
gives outward log2 upper +256.190414, rounded upward. It gains only
about 0.61 bits. The selected fractions are one for the zero column at
local occupancies 6, 7, 8, 11, 12, 13, and 14; they are one for the
density column at 6 through 12. All other fractions are zero.
Improving the local feedback information has therefore been more useful
than further mixture tuning in these tests.

```sh
python -B research/workstreams/permutation_locality/independent_rows/candidates/mass_optimize.py --groups 80 --support 200 --tilt .068 --maximum 16 --iterations 30 --passes 2 --spectral-feedback 10
```

## Exact peaks by limb inversion

The opt-in `candidates/exact_feedback.py` computes the complete feedback
distribution through ten windows, including its exact nonzero peak.
It uses the same character polynomials and construction. Unlike the
compressed method, it reconstructs all 2^19 state counts for each shape.

A direct signed-int64 Walsh transform can overflow before division by
S, even when every final count fits. The new routine first checks a
bound on every butterfly. If that bound is too large, it writes each
signed coefficient as

    f = 2^32 h + l,    0 <= l < 2^32.

Linearity gives W(f) = 2^32 W(h) + W(l). Both limb transforms fit
signed-int64 under separately checked bounds. Since S divides 2^32,
integrality requires W(l) to be divisible by S. The reconstruction is

    n = (2^32/S) W(h) + W(l)/S.

The routine checks the multiplication and addition bounds before
reconstruction. It then checks nonnegative counts, the total using
unbounded integers, parity, the zero count, and selected characters.
Known short-census records are compared exactly. Small-state tests also
cover signed coefficients whose unsplit inverse exceeds 64-bit range.

`--exact-feedback 10` selects this route in `candidates/mass_verify.py`
or `candidates/mass_optimize.py`. It is mutually exclusive with the
compressed spectral option. The longer exact census has now passed all
1000 shapes, including exact agreement with the previous 209 short-shape
records. No complete occupancy result uses it yet. Neither candidate
changes the production verifier or writes its operators to the baseline
cache.

```sh
python -B research/workstreams/permutation_locality/independent_rows/candidates/mass_verify.py --groups 80 --supports 192 200 208 --tilts .064 .068 .072 .076 --exact-feedback 10
```

At occupancies seven through ten, the worst all-target atom probabilities
are respectively 5/560976, 11/1753050, 49/9349600, and 5/1075204.
Each maximum occurs at the shape consisting only of weight-four packets.
These are exact local probabilities, not full-code failure bounds.

The command above gives the following selected-point results at 192-bit
precision. Log2 uppers are rounded upward.

| Union support per group | Selected tilt | Outward log2 upper |
|---|---:|---:|
| 192 | .068 | +181.416331 |
| 200 | .072 | +210.043130 |
| 208 | .072 | +60.370886 |

The exact peaks improve the best compressed results but still do not
close these events. The four-window joint-output refinement in
[JOINT_FOUR.md](JOINT_FOUR.md) and the density extension in
[DENSITY_EXTENSION.md](DENSITY_EXTENSION.md) target different remaining
losses. The four-window replay removes about another 40--45 bits at
these points; the density extension remains under evaluation.
