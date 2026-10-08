# Regional designs: exact rejection through period sixty-four

The first screening targets are rate 3/4 at 10% symbol distance, or 30% distance at rate 1/2, with competitive encoding cost.
An intermediate rate may still be useful if its complete proof and encoding cost improve the existing half-rate baseline; 3/4 is not a user-imposed hard cutoff.
Forty bits of whole-code setup-failure margin is the initial proof target.
This checkpoint does not establish either complete guarantee.

The analysis rules out the inexpensive regional designs below by actual failure-probability lower bounds.
The period-64 modification passes the one-row and paired-window gates but fails a thirty-seven-row window event at K20.
Its exact failure lower bound exceeds 2^-38, so it cannot meet the 2^-40 target.

## A cheap deterministic repair for one active row

Work first over F=GF(256), with its GF16 subfield U.
Let A0 have rows (1,u,u^2), indexed by u in U, and put C0=A0^T.
Then C0 A0=0, and every nonzero state emits at least fourteen nonzero symbols through A0.

Assign each of sixteen outer-coordinate regions a distinct additive coset beta+U.
Define

    P_beta = [[1,beta,beta^2], [0,1,0], [0,0,1]],
    C_beta = P_beta^T C0,       A_beta = A0 P_beta^T.

The injection columns are (1,u+beta,(u+beta)^2)^T.
The emission map deliberately differs from C_beta^T.
Both maps retain the needed MDS properties, and C_beta A_beta=0.

The global recurrence is Y=X+A_beta Q and Q'=Q+C_beta X.
In local coordinates S=P_beta^T Q, it becomes

    Y=X+A0 S,       S'=S+C0 X.

The boundary change from beta to gamma is S <- P_(beta+gamma)^T S.
It requires only two GF256 state products per region boundary.
There is no refresh inside a region.
The transpose uses the transposed boundary map; the scalar reference checks this adjoint identity.

For one nonzero outer row, every injection lies in a distinct region and therefore at a distinct Vandermonde point.
Any three such columns are independent over F and over every compatible scalar extension.
Starting from zero, the state cannot first return to zero before the fourth injection.
Such a return has at least two complete intermediate regions between its first and last injections.
Throughout those regions, the state before and after the single active packet is nonzero.

Write L=16m for the regional length.
Each complete intermediate region contributes at least 14m-1 output symbols: its active input can cancel at most one emitted symbol.
If the state never returns, outer distance at least three leaves at least two complete regions after the first injection.
Thus every one-row word has output weight at least

    2(14m-1) = (7/64)N-2,       N=256m.

This includes partial first and last regions and discarded final state.
It exceeds 10% at the tested sizes.
The scalar tests include a four-injection return followed by a final-region injection at the boundary packets.

### Two rows invalidate the cheap version

Take the fixed GRS[16,12] outer, independent nonzero GF256 input diagonals, and the sixteen-label inner above.
Two source rows can have the same lane in each coordinate of a minimum support S of size five.
If their diagonal ratios are constant across S, scalar multiples of the unique minimum-support outer word cancel pairwise.
Each active region then starts and ends with zero state.
Its output lies between the two routed inputs.

This event gives exact setup-failure lower bounds:

| Nominal size | Rounded K | Failure lower bound |
|---:|---:|---:|
| K16 | 65664 | 6.6327710e-6 |
| K18 | 262272 | 1.0788694e-4 |
| K20 | 1048704 | 1.7252443e-3 |

The displayed decimals round downward.
They are lower bounds on actual ensemble failure, rather than failures of a proof envelope.
The full finite calculation sums row-pair and support events, then subtracts upper bounds for every pairwise intersection.
An independent audit checked the support-overlap counts, diagonal-ratio ranks, permutation probabilities, and gap count.

## A stronger outer still needs a longer injection period

The next outer is fixed systematic GRS[32,24] over GF256, with one independent permutation in each of its thirty-two regions.
It has rate 3/4 and minimum distance nine.
Choose independent uniform GL(4,256) state orientations once per region, retaining the fixed width-16 GF16 kernel inside each region.
Input diagonals are independent uniform nonzero GF256 elements.
All orientations, diagonals, and permutations are independent setup choices.

The once-per-region refresh suppresses a one-row return with probability 1/(256^4-1).
It does not distinguish repeated injection labels within that region.
This leaves a stronger two-row obstruction that requires no favorable diagonals.

Suppose two rows have matching labels in exactly seventeen outer-coordinate regions S.
The shortened GRS[32,24] space on S has dimension nine.
Two such words give eighteen source dimensions, while equality of their diagonally scaled coordinates imposes only seventeen equations.
A nonzero cancelling pair therefore exists for every choice of diagonals and regional orientations.
Both rows are nonzero, since equality would force the other row to vanish if either vanished.

Let the injection period be B, a multiple of sixteen, and write L=Bm.
Conditioned on a matching label, the two copy indices are a uniformly sampled ordered pair of distinct elements of {0,...,m-1}.
Their positive gap d has probability 2(m-d)/(m(m-1)).
If the seventeen gaps sum to at most floor((e-17)/B), the output weight is at most e, where e=ceil(N/10)-1.
Indeed, each paired region contributes at most Bd+1 symbols.

The cumulative gap count is evaluated exactly from

    W(x) = sum_(d=1)^(m-1) (m-d)x^d,
    (1-x)^2 W(x) = (m-1)x - mx^2 + x^(m+1).

Requiring exactly seventeen matches gives a unique set S for each row pair.
For intersections of two row-pair events, a bivariate binomial calculation bounds the probability that both match counts reach seventeen.
It uses the exact per-region shared-row or disjoint-pair permutation law and drops the gap restrictions.
The first two Bonferroni terms then give:

| Period B | Nominal size | Rounded K | Failure lower bound |
|---:|---:|---:|---:|
| 16 | K16 | 65664 | 6.8736015e-9 |
| 16 | K18 | 262272 | 1.3202456e-7 |
| 16 | K20 | 1048704 | 1.9696510e-6 |
| 32 | K18 | 262656 | 1.5562446e-12 |
| 32 | K20 | 1049088 | 2.7468742e-11 |

Each displayed decimal rounds downward and exceeds 2^-40.
Thus periods sixteen and thirty-two are rejected at the target sizes.
The independent audit also verified this rank witness and the exact bivariate calculation.

## Period sixty-four: narrow gates passed before rejection

For the period-64 design, let A0 evaluate monomials of degrees zero through three on U, and put C0=A0^T.
Choose a,b in GF256 such that U, a+U, b+U, and a+b+U are distinct cosets.
Cycle the packet offsets through 0,a,a+b,b.
This produces sixty-four distinct injection points before repeating.

Let T_beta be the translation matrix satisfying

    T_beta (1,u,u^2,u^3)^T = (1,u+beta,(u+beta)^2,(u+beta)^3)^T.

For packet offset beta, use C_beta=T_beta C0 and A_beta=A0 T_beta^-1.
All these maps satisfy C_beta A_beta=0.
An independent regional orientation V multiplies C_beta on the left and A_beta on the right by V^-1.
The global state starts at zero, carries across every boundary, and is discarded at the end.
There is exactly one global routing stage: the thirty-two independent regional permutations.

The local basis changes between packets use only T_a or T_b, because the offsets follow the stated cycle.
For a four-symbol state, the two-shear reference implements each translation with four fixed GF256 products and four XORs.
The regional GL(4,256) change occurs only thirty-two times, not at every packet.
Within each packet, A0 and C0 use fixed GF16 coefficients.
No operation requires multiplication by a general element of GF(2^m), for m in {8,32,64,128}.

The [GRS32_24_COST.md](GRS32_24_COST.md) circuit computes the outer parity transpose with 34 fixed GF256 products and 114 Boolean instructions per four source rows.
It also requires systematic additions, transposes, loads, and stores as recorded there.
The period-64 inner still adds the input diagonal, fixed A0/C0 work, packet translations, routing, and output traffic.
These are operation estimates; competitive total runtime has not been established.

### Exact one-row bound

Fix an entire one-row word, its route, and all input diagonals, then expose the regional orientations in order.
Before an active region, a nonzero incoming global state is determined by earlier orientations.
The fresh regional orientation makes its local coordinates uniform over the nonzero four-dimensional state vectors.
For the single active input in that region, the conditional return probability is at most p=1/(256^4-1).
No uniform-state assertion is made at intermediate packets.

For a word of weight w with r returns, at least w-2r-1 active regions are neither a birth nor a return region.
Each such region contributes at least 13(L/16-1) output symbols, even after dropping its active packet entirely.
At all three rounded sizes, a bad word permits at most three such regions.
It therefore requires r_w=max(0,floor((w-3)/2)) returns.
Return positions cannot be adjacent in active-region order or occupy the first position.
There are binom(w-r_w,r_w) possible selected return sets.

Let A_w be the number of GRS[32,24] words on a prescribed support of size w.
The exact one-row union bound is

    L sum_(w=9)^32 binom(32,w) A_w/(255)
        * binom(w-r_w,r_w) / (256^4-1)^r_w.

The Fraction replay with L rounded to a multiple of sixty-four gives:

| Nominal size | K | N | One-row margin |
|---:|---:|---:|---:|
| K16 | 66048 | 88064 | 45.5514 bits |
| K18 | 262656 | 350208 | 43.5598 bits |
| K20 | 1049088 | 1398784 | 41.5619 bits |

Weight ten dominates this bound.
These are exact rational bounds covering all one-row words, not the complete code.
The independent audit recomputed the return counting and projective MDS multiplicities.

For period sixty-four, the exact seventeen-match event family has a K20 union upper bound of about 2^-51.7809.
Summing the tested paired-window family over every support weight gives an upper bound below 7.202230e-16, about 2^-50.3024.
This family requires matching labels at every active coordinate and bounding windows that themselves fit inside the bad-weight budget.
It does not include all bad two-row words or all multirow words.

## Why ordinary word counting also fails

At period-64 K20, fix twelve source rows and the final 4368 positions in each region.
The probability that all their coordinates land in these windows is

    P_window = (binom(4368,12)/binom(43712,12))^32.

The source space has dimension 24*12.
Imposing zero state after the first thirty-one regions costs at most 4*31 dimensions.
The remaining kernel therefore has dimension at least 164.
Every kernel word has output confined to the thirty-two windows, whose total length is 139776, below the bad-weight threshold 139878.
The final state is discarded, so the last region adds no return constraint.

For this fixed row set, the expected number of bad projective words is at least

    P_window * (256^164-1)/255 > 2^27.

The exact value has binary logarithm about 27.3511, although log2(P_window) is about -1276.6546.
This is a lower bound on an expected word count, not a lower bound on setup failure.
It prevents an ordinary full-word first moment from proving the desired failure margin.

A successful argument for a revised construction must bound the existence of a bad subspace before averaging rare route events.
The following regional interface remains available for that analysis.
For one fixed region, write Cbar X for the sum of its packet injections, Abar R for its emissions from a fixed incoming state, and B for its deterministic zero-initial-state scan.
Then

    Y = B X + Abar R,       Q_next = V(R+Cbar X),
    B^2=I,                 Cbar B=Cbar,       B Abar=Abar.

These identities follow from C_beta A_beta=0 packet by packet, including changing packet offsets.
The independent orientation supplies a uniform nonzero R only at region entry.
The uniform-state transfer from a per-packet refresh cannot be reused inside the region.
The next calculation counts route events directly and rejects the period-64 design despite its narrow gates.

## A thirty-seven-row event rejects state four

Fix the K20 parameters L=43712, outer GRS[32,24], and four field elements of persistent state.
The encoder is linear and causal, with invertible identity-feedthrough packet maps, homogeneous state updates, and zero initial state.
The claimed family, including all choices of its diagonals, packet phases, and regional orientations, has these properties.
The only randomness used below is the thirty-two independent uniform regional permutations.

Use three packet-aligned windows of length W=15520 in each region, starting at 0, 14096, and 28192.
Adjacent windows intersect in I=1424 positions; the first and third windows are disjoint.
A box selects nine regions and one window in each selected region.
There are

    M = binom(32,9) * 3^9 = 552084530400

boxes.
Its occupancy X counts source-row labels routed into all nine selected windows.

If X is at least thirty-seven, choose any thirty-seven such rows.
The outer space shortened to the selected nine coordinates has dimension one in each row.
These rows therefore supply thirty-seven source dimensions.
Impose zero state immediately after each of the nine selected windows.
These homogeneous linear constraints cost at most 4*9=36 dimensions, leaving a nonzero source word.
Induction over the windows shows that its output vanishes outside the selected windows: no input or incoming state occurs in any intervening interval.
Invertibility ensures that the output word is nonzero.
Its weight is at most 9W=139680, below the bad-weight threshold e=139878.
This implication holds for every realization of the diagonals, packet phases, and regional orientations; no favorable field values are needed.

The occupancy distribution permits exact alternating bounds without counting one event repeatedly for each thirty-seven-row subset.
Write (n)_j=n(n-1)...(n-j+1).
For a fixed box,

    mu_j = E[binom(X,j)] = (W)_j^9 / (j! (L)_j^8).

The alternating tail identity is

    Pr[X >= 37] = sum_(j=37)^W (-1)^(j-37) binom(j-1,36) mu_j.

Truncating after j=68 gives a lower bound; truncating after j=69 gives an upper bound.
Both are evaluated as exact rational numbers.

To bound intersections of distinct box events, use degree m=29:

    Pr[X >= 37 and Y >= 37]
        <= E[binom(X,m) binom(Y,m)] / binom(37,m)^2.

Consider two m-subsets of row labels sharing c labels.
Their number, as an ordered pair, is binom(L,m) binom(m,c) binom(L-m,m-c).
For one shared region, put D_c=(L)_(2m-c).
The probabilities for equal, adjacent, and disjoint windows satisfy

    s_c = (W)_(2m-c)/D_c,
    a_c <= I^c W^(2(m-c))/D_c,
    d_c = 1_(c=0) (W)_m^2/D_c.

The adjacent-window bound permits position collisions, so it is an upper bound on injective routing.
Among the nine ordered choices of two windows, three are equal, four are adjacent, and two are disjoint.
Thus each shared region contributes at most 3s_c+4a_c+2d_c when window choices are summed.
Each unshared region contributes 3(W)_m/(L)_m.

Group ordered coordinate-set pairs by their intersection size r.
There are binom(32,9) binom(9,r) binom(23,9-r) such pairs.
For r=9, subtract (3s_c)^9 before summing: those terms represent identical boxes.
Finally divide by two for unordered box pairs and by binom(37,29)^2 for the moment normalization.
Every relaxation increases this second Bonferroni term.

The exact replay obtains the following displayed approximations:

| Quantity | Value |
|---|---:|
| Sum of box-event probability lower bounds | 6.387983463554045e-12 |
| Sum of distinct box-pair intersection upper bounds | 2.886945594218527e-16 |
| Difference: actual setup-failure lower bound | 6.387694768994623e-12 |

The rational comparison verifies that the difference exceeds 2^-38, and hence 2^-40.
An independent implementation reproduces all three rational quantities.
Small exhaustive tests check the alternating tails, injective regional assignments, and the full box-pair sum.
This is a setup-failure lower bound, unlike the earlier lower bound on an expected word count.
Changing the injection period alone cannot remove this event at the stated geometry and state dimension.

The next construction change should be screened against this occupancy event before further kernel work.
Adding one outer parity coordinate, such as GRS[33,24], is a possible knob: its rate 8/11 falls below the first screening target but remains substantially above the half-rate baseline.
It is only a proposal, with neither a complete distance certificate nor a runtime claim.
Changing the state dimension also changes the one-row return bound and must be checked separately.
Neither proposal has been certified here.

## Reproduction

The references are `rate_design_coset.py`, `rate_design_region_refresh.py`, `rate_design_collision.py`, and `rate_design_window_boxes.py` under `theory`.
Eighteen exact or finite reference tests pass:

```sh
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p rate_design_tests.py -v
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p rate_design_window_tests.py -v
python -B research/workstreams/char2_field_spin/theory/rate_design_window_boxes.py
```

The references check GF256 coset algebra, global/local equivalence, the adjoint, endpoint packets, gap counts, permutation laws, diagonal-ratio ranks, Bonferroni tails, one-row return counting, and both window obstructions.
The floating regional-moment search remains a separate one-row screen and is not used for the displayed exact margin.
No native code, remote commands, benchmarks, or commits were added in this checkpoint.
