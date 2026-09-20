# Separate constant rows in the dense range

The full-occupancy calculation is sensitive to the all-one BCH codeword.
An auxiliary Bernoulli row assigns positive probability to defects in that
codeword. The resulting moment bound can be poor even when the actual
all-one input has deterministic, large output weight. We therefore keep
all-one rows exact and apply the Bernoulli comparison only to other rows.

This note uses the construction and setup probability space in `README.md`:
BCH [256,128], RM2Sub t128_s15, 8192 rows, 256 regions, and cutoff 209716.
Multipliers are fresh independent nonzero field elements at every epoch.
They are shared across messages after setup. The results below bound the
expected number of bad nonzero messages over this setup randomness.

## Counts and a bound for one composition

Call a BCH row *ordinary* if it is neither zero nor all-one. Let d count
ordinary rows and h count all-one rows. There are 8192-d-h zero rows, and
the nonzero occupancy is Q=d+h. A composition means one pair (d,h).

Fix z in (0,1). Let R_j(z) be the nonnegative 3-by-3 envelope for one
uniformly permuted region of weight j, as defined in
`SHARED_RANGE_1024.md`. Its state classes retain activation separately
from the near-uniform live state. Write e_Z=(1,0,0) and let 1 be the column
vector of three ones. Matrix inequalities below are entrywise.

Partition the possible ordinary weights into the twelve bands in
`occupation_three.BANDS[:-1]`. For each band g, choose p_g in (0,1), and
let U_w be a certified upper bound on the number of BCH words of weight w.
Define

\[
\Gamma_g=\max_{w\in g}
 \frac{U_w}{\binom{256}{w}p_g^w(1-p_g)^{256-w}},
\qquad \rho_g\geq\Gamma_g^{1/256}.
\]

The counting measure of a permuted ordinary row in band g is dominated by
Gamma_g times the product Bernoulli(p_g) measure. Every all-one row adds
exactly one bit to every region. For this fixed composition, initialize
K_j^(0)=R_(h+j) for 0<=j<=d, and recursively define

\[
K_j^{(r+1)}=
 \max_g\rho_g\bigl((1-p_g)K_j^{(r)}+p_gK_{j+1}^{(r)}\bigr),
\qquad 0\le j<d-r.
\]

For every fixed sequence of ordinary bands, induction bounds its region
matrix by K_0^(d), including its share of the density cost. All matrices
are nonnegative, so multiplication preserves this domination. There are
12^d band sequences. The expected number of bad messages in the composition
is therefore at most

\[
U_{d,h}:=\binom{8192}{d}\binom{8192-d}{h}\,12^d z^{-209716}
 e_Z(K_0^{(d)})^{256}\mathbf1. \tag{1}
\]

For d=0, use K_0^(0)=R_h without a density cost. Witnesses p_g and z may
differ between compositions: the first moment is summed after each bound.
Checking pure bands alone does not prove (1) for mixed bands.

`certify_constant_split_points.py` evaluates (1) with outward bounds.
Its producer uses 256-bit Arb inputs and final operations, plus an upward
binary64 recurrence with a separate exponent for each degree. The replay
uses 512-bit inputs. `fast_replay_constant_split.py` also removes redundant
coefficient pairs using an exact rational upper hull. This pruning retains
the factor 12^d; the number of retained hull lines is not a band count.

## A whole class from a region-weight floor

Many compositions can be bounded together without any BCH shell caps.
Fix integers D and b, and consider messages with d<=D and h>=b. Every
region has weight at least b, regardless of its ordinary rows. Define

\[
E:=\max_{b\le j\le8192}R_j,
\qquad
C_D:=\sum_{d=0}^{\min(D,8192-b)}
 \binom{8192}{d}\,2^{8192-d}(2^{128}-2)^d.
\]

For each eligible message, condition on the row permutations. Its region
matrices are bounded by E, so its moment is at most e_Z E^256 1.
The bound is uniform in those permutations. For each ordinary row there
are exactly 2^128-2 choices. Enlarge the count of constant assignments to
2^(8192-d) only after applying the uniform bound to eligible messages.
The combined first moment of this class is at most

\[
C_D z^{-209716}e_ZE^{256}\mathbf1. \tag{2}
\]

Some assignments counted by C_D do not satisfy h>=b. They enlarge the
count; the argument does not claim that their moments satisfy this bound.

## Verified coverage added

Both of these outward certificates passed their 512-bit replays:

| Scope | Combined contribution | Receipt in `generated/` |
|---|---:|---|
| d<=1200 and h>=513 | <2^-5390 | `region_floor_d1200_h513_outward.json` |
| d=8192 and h=0 | <2^-9382 | `constant_split_ordinary8192_outward.json` |

The first uses (2), with z=exp(-exp(-1.4)). The second uses (1), with
z=exp(-exp(0.8)), and certified BCH caps through the snapshot in its screen.
The diagnostic margins are 5390.919941821194 and 9382.329058863646 bits.
These are margins for the indicated entire class or composition, not
per-message margins.

An earlier, weaker certificate covers Q>=1025 and d<=512 with contribution
<2^-67182. Its 512-bit replay also passed. That class is contained in the
new region-floor class and need not be added again.

An improved-cap shared calculation additionally certifies **every occupancy
Q1025--Q1251**, including all mixtures of zero, ordinary, and all-one rows.
Their combined contribution is <2^-64 (64.88651381072486 diagnostic bits).
It reuses the Q1024 witness at tilt -18 with the stronger exact caps and
the original thirteen-band assignment count. Each recurrence depth is
evaluated, rather than inferred from an endpoint.

The receipt `exactcap_shared_q1025_q1280_outward.json` retains all 256
attempted occupancies. Q1252--Q1280 do not pass and are excluded from the
ledger. The frozen producer's embedded replay compares different key
domains: the stored caps omit weight zero, but the reconstructed dictionary
includes it. `replay_exactcap_shared_range.py` compares the nonzero domains
explicitly and recomputes all 256 rows. This 512-bit replay passed; neither
the original producer nor its receipt was changed.

A further shared calculation certifies every Q1252--Q1655, with combined
contribution <2^-79 (79.94601511536166 diagnostic bits). The witnesses were
optimized at Q1280,1536,1792,2048; the adaptive recurrence still includes
all thirteen bands. The producer and 512-bit replay both passed for every
stored row of `full_exactcaps_q1252_q1792_outward.json`. Its nonpassing
Q1656--Q1792 rows remain stored but are excluded from the ledger.

Together with Q1--Q1655, the certified union still has contribution <2^-49.
The region-floor class overlaps some low occupancies; adding both bounds
is safe. This does **not** cover all compositions at Q8192, nor any complete
occupancy above Q1655. `verify_constant_split_coverage.py` checks the
partial ledger and reports the uncovered compositions.

Replay commands:

```powershell
python -B workstreams/bch_rm2sub_bridge/certify_region_floor_class.py --ordinary-limit 1200 --floor 513 --tilt -14 --tag d1200_h513 --verify
python -B workstreams/bch_rm2sub_bridge/fast_replay_constant_split.py --screen constant_split_ordinary8192_exactcaps_screen.json --certificate constant_split_ordinary8192_outward.json
python -B workstreams/bch_rm2sub_bridge/replay_exactcap_shared_range.py
python -B workstreams/bch_rm2sub_bridge/certify_full_exactcaps_range.py --screen full_exactcaps_middle01_screen.json --lower 1252 --upper 1792 --tag q1252_q1792 --verify
python -B workstreams/bch_rm2sub_bridge/verify_constant_split_coverage.py
```

## Remaining work and failed diagnostics

The remaining compositions satisfy Q>1655 and either d>1200 or h<513,
apart from the verified composition (8192,0). The middle mixtures at
Q8192 remain unresolved. Tested witnesses for (2048,6144), (4096,4096),
and (6144,2048) gave vacuous bounds. This is a limitation of the tested
envelopes, not a counterexample to the code.

The direct iid-epoch identity agrees with the region-polynomial calculation
at the full-occupancy pure-band probes. Exponential-mode fits, including
an interior-slope refinement, were too loose. A separate single-exponential
envelope preserved total BCH weight across regions but also gave vacuous
bounds in the first probes (`screen_weight_conservation.py`). None of
those floating screens adds certified coverage.

The Q2048 single-shell probes also remain vacuous for some weights, even
without any band-assignment cost. For example, the best tested witness
at tilt -8 loses about 15,352 bits at weight 96. Thus broad bands and
adaptive maximization alone do not explain the observed gap. Several
middle shell caps are still older bounds; a serial exact sweep is filling
those gaps. The setup and map parameters are unchanged.

`screen_shared_mode_slopes.py` tries a different fit for the fixed-weight
mode formula in `DENSE_RANGE_ALTERNATIVES.md`. All matrix entries on a
segment use the same rate, chosen from the zero-to-zero entry. Every
entry receives its own covering coefficient. This avoids choosing growing
rates solely because an activation entry initially rises. The choice of
rate is a discovery heuristic; a future certificate must recompute every
entrywise covering inequality. No mode-fit screen counts as a certificate.

The next objective is shared coverage of ranges of (d,h), with tighter
control of ordinary-band mixtures where current witnesses fail. Exact
shell caps and composition-dependent tilts remain available inputs; they
do not by themselves establish range coverage.
