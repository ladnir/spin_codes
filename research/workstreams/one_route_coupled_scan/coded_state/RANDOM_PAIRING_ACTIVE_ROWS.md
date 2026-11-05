# Retain the exact active outer-row set

The capped count improves when it distinguishes inactive outer rows from
active rows. An outside crossing requires a nonzero input, so it cannot
land on an inactive row. This constraint was lost in the deficit moment.
The construction, setup distribution, and accepted pair arithmetic remain
unchanged.

This continuation extends restricted coverage from8192 to16384 swaps.
It does not close the whole code. In particular, a sum over all active
row sets still fails at moderate turnover in the longer retained suffix.

## Conditional bound

Use the support histories and enlarged union U from
[the capped proof](RANDOM_PAIRING_CAPPED.md). There are M=8192 outer rows,
R suffix regions, L=16R updates, and T=L-1 transitions after activation.
Each state after activation has weight256. Let B count births, hence also
deaths. The enlarged union has size A=L*256+B+256. The field has q=2^32.

Fix a history, route, and matchings. Write J for the outside crossing
coordinates. Fix the exact set S of nonzero outer rows. Let b_i and j_i
count row i's positions in U and J. A compatible message requires:

- every position of J belongs to a row in S;
- every row in S has b_i+j_i>=129 covered coordinates.

On this event, shortening the RS rows gives dimension
d_S=sum_{i in S}(b_i+j_i-128). Therefore

    d_S-|J| = sum_{i in S} b_i -128|S|.

For z selected zero-output coordinates and C possible selections, the
conditional bad-message probability is at most

    1[compatibility] min(1, 2 C q^[sum_{i in S} b_i-128|S|-B-z]).

The count includes all messages on those rows. The cap follows that union,
before averaging the shared route and matchings. Discarding the second
compatibility condition gives a valid upper bound, including when its
formal dimension exponent is negative: the actual event is then empty.

## Conserved route moment and sum over all row sets

As in the capped proof, compare external crossing endpoints with independent
fair coins at density cost rho^T, where log2 rho is approximately4.82645.
Represent regional U memberships by independent Bernoulli indicators with
parameters p_j, conditioned on their actual column sizes. The conditioning
cost is at most (M+1)^R.

Choose 1<r<=q and t>0. Apply min(1,x)<=x^theta with theta=log_q r.
Insert t^(A-sum_i b_i), which equals one on the conditioning event.
For an active row, after dropping its coverage constraint, the row moment
is product_j(1-p_j+p_j*r/t). For an inactive row, no outside-U coin may
be one; its moment is product_j((1-p_j)/2+p_j/t).

The arithmetic-geometric mean inequality bounds each product using
p=A/(RM). Define

    a=(1-p+p*r/t)^R,
    b=((1-p)/2+p/t)^R.

Summing over all subsets S, including the harmless empty-set overcount,
gives the bound

    H rho^T (M+1)^R
      r^[(1+log2 C)/32-B-z] t^A (r^-128*a+b)^M.             (1)

Here H counts the support histories in the chosen family. Histories and
row sets are fixed before averaging setup, then unioned. This is not a
claim about a message chosen independently of setup only.

The FLOAT helper also supports a coverage marker v>=1. It replaces a by
v^-129*((p*r/t+(1-p)/2)*v+(1-p)/2)^R. All selected all-row-set optima below
use v=1. The directed replay therefore does not need that refinement.

## Directed restricted replay

At distance46%, W=floor(.46*2^21)=964689. For a turnover interval [l,h],
use z=max(0,L*256+l-W) selected zeros and C<=binom(T*256-l,z).
Bound the history count by a positive marker for all B<=h. Bound each row
factor at the worse endpoint p, the positive t^A at its worse endpoint,
and the death charge using B>=l. These choices give an upper bound on
every history in the interval, without relying on FLOAT optimization.

`certify_random_pairing_active_rows.py` uses exact rational values and
integer-power logarithm enclosures. Nine intervals per suffix cover every
integer B from0 through16384. Fixed weights and worst cell results are:

| R | r | t | Worst directed log2 bound |
| --- | --- | --- | --- |
| 192 | 8 | 4 | -404016 |
| 256 | 16 | 12 | -35799 |

Each cell is below2^-120. The union of18 cells is below2^-115; the replay
conservatively quotes2^-113. This covers all support masks, swap placements,
active row sets, and messages in these two half-full-state suffix families.
It does not cover changing state weights, repeated resets, other suffixes,
or extension to GF(2^128). Earlier restricted certificates remain intact.

## The remaining gap

Selected exact-B FLOAT proposals from (1), in log2 units:

| B | R192 | R256 |
| --- | --- | --- |
| 16384 | -422762 | -57792 |
| 32768 | -279306 | +54115 |
| 65536 | -63543 | +200645 |
| 131072 | +208414 | +297971 |
| 262144 | -286606 | -215200 |

These are upper-bound screens, not certificates or exhibited bad words.
Seven sampled row counts all passed at R256,B32768, but their sum over
every row set fails. The sampled counts must not be reported as coverage.
Large turnover eventually forces enough nonzero output; the remaining
problem lies between nearly constant supports and that regime.

The next structural test should retain constraints at births and deaths.
For example, a birth whose incoming pair is empty requires a nonzero
outer input. Equation(1) currently permits that birth on an inactive row.
An active-inactive pair also forces nonzero inputs at certain deaths and
zero outputs, by superregularity. Test these constraints jointly with S
before tuning more numerical weights or changing the construction.

## Validation and execution

`test_random_pairing_active_rows.py` passes175 dimension identities,
12 exact row-moment checks, and12 exact conserved-routing comparisons.
The tests enumerate small regional routes and outside-U coins, including
the sum over all active row sets and both t<1 and t>1.

Runs were serial under the20GiB disk guard, without raw output files.
An initial finer interval replay exceeded its25s budget and was stopped
by its own check. A coarser nine-cell replay completed in7.80s; its bounds
remain directed. Exact tests took0.008s. No native benchmark, library or
paper changes, commits, or promotion. The goal remains active.
