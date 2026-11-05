# Limits of two MDS refinements in the current actual-state bound

2026-10-31. Two precise limitations explain why the next proof step should
not be another output-MDS marker inside the current fractional moment.
The combined deterministic output budget is redundant on the selected
failing profile. At the current output marker and fractional exponent,
selected-zero MDS filters can improve the untruncated counting expression
by less than three bits globally. An independent read-only audit passed
both limitations and the optional-coefficient derivation below. These
are limitations of the bounds, not a failure of the fixed construction.

Use the actual-state definitions from
[the actual-state profile note](QUARTET_ACTUAL_STATE_PROFILES.md):
Ls=sum|S_t|, factual=f, F=sum|Gamma_t(S_(t-1)) union S_t|,
G=F-f, D the coordinate deaths, and a=|S_T|. In particular
D<=F-Ls and a<=512. These are actual-state identities and inequalities;
no availability-only feasibility condition is used below.

## Even the combined deterministic output budget is redundant

In a quartet with nonzero incoming state, let r be its incoming support
size, d its coordinate deaths, b its coordinate births, and m=d+b the
number of changed coordinates. Every changed coordinate has nonzero
output. Let j be the number of positions outside the eligible-row set H,
and let q_out count changed positions outside H.

On H's complement the outer input is zero, so the output equals
(A+delta I)Q_prev. Superregularity therefore forces at least
max(0,j-r+1) nonzero outputs there. Combining this fact with the changed
coordinates gives the lower output budget

    beta = m + max(0,j-r+1-q_out).                         (1)

This is a lower bound on the actual output weight, not an expression for
that weight. Since j-q_out counts outside-H positions that did not change,

    j-q_out <= 4-m,
    beta <= max(m,5-r).

Moreover r+b is the size of the union of the old and new supports, so
r+b<=4. Hence

    m-(5-r) = d+(r+b)-5 <= d-1,
    beta <= 5-r+max(0,d-1).                               (2)

Unhit quartets contribute exactly f forced-birth outputs to this budget.
Let Q_death count live quartets containing at least one coordinate death.
There are G/4 live quartets and their incoming support sizes sum to Ls-a.
Summing (2) therefore gives

    B_joint <= 5F/4-Ls-f/4+a+D-Q_death
            <= 9F/4-2Ls-f/4+512.                         (3)

At the profile

    (h,R,Ls,f,F)=(4096,786432,786432,4,1048576),

the last quantity is 786943, strictly below W=943718. Thus every history
in this exact profile already satisfies B_joint<=W. Even the combination
of forced coordinate changes and the outside-H output-MDS lower bound
excludes no history in the profile.

This statement does not say that every such history has output weight at
most W. The lower budget (1) is weaker than the exact cumulative output
cutoff. A marker rho^(W-B_joint), rho>=1, cannot strengthen the exact
profile sum when (3) is already below W. Any numerical improvement in a
larger relaxed transfer would concern that relaxation, not a new exclusion
within this profile.

## Selected-zero filters gain less than three bits at the current markers

Fix an actual support history and the H labels. Write i for its number
of coordinate survivals. In the existing weighted projective count, each
survival contributes either a selected zero output, of weight 1, or a
nonzero output, of weight qx. Their unrestricted factor is

    (1+qx)^i.

Consider an A+delta I output-zero filter, or the
[mixed-minor selected-zero filter](QUARTET_MIXED_MINORS.md), added to this
same untruncated weighted expression. The already imposed A-minor rule
requires

    |H^c intersect (complement of new support)| <= r-1

in every live quartet. When no surviving output is selected to be zero,
the A+delta I rule only tests outside positions where both states are
zero; these are a subset of the displayed set. The mixed-minor rule at
an empty selected-zero set is exactly the already imposed A rule.
Consequently the all-nonzero survivor-output term remains feasible in
both combinatorial filters.

The filtered optional-output polynomial therefore retains a term of
weight (qx)^i. Before taking a fractional power, the ratio between the
unfiltered and filtered expressions is at most

    (1+1/(qx))^i.                                        (4)

At q=2^32, x=2^-16 and theta=1/16, positivity and i<=N=2^21 give the
following uniform ratio after the theta power and the history sum:

    (1+2^-16)^(N/16)
      = (1+1/65536)^131072 < (11/4)^2 = 121/16 < 8.     (5)

The finite rational series estimate used in the actual-state note gives
the first strict inequality. Corresponding regional bounds produced by
the same positive operations, quota conditioning and common-basis
entrywise p minima preserve this uniform comparison. Arbitrarily chosen
independent upper bounds need not preserve it. The improvement in the
current counting expression is therefore less than three bits. A
separate mixed-minor setup-exception allowance does not improve this
comparison.

The scope of (5) matters. It concerns selected-output-zero exclusions
that retain the all-nonzero term inside the current weighted sum. It
does not bound the gain from extracting the exact output coefficient,
truncating at total weight W before taking a fractional power, or proving
additional amplitude/rank constraints. In particular, the all-nonzero
term has output weight Ls+D; an exact cutoff removes it when Ls+D>W.

## Exact optional-output coefficients give a concrete next interface

For a fixed support history, total coordinate births are b=D+a and
survivals are i=Ls-D-a. The b+D changed coordinates force nonzero outputs.
Set

    m = min(i,W-2D-a).

If m<0 the output event is impossible. Otherwise the same projective
rank argument, summing only admissible optional outputs, gives

    P_S <= min(1,(gamma/q) q^(b-R)
                    sum_(k=0)^m C(i,k) q^k).             (6)

Each term preceding the last is at most i/q times its successor. Thus
the sum is at most C(i,m)q^m/(1-N/q). Since
gamma/(1-N/q) <= (1024/989)(2048/2047)<2, define

    Delta = max(0,Ls+D-W).

For a feasible output support, (6) implies the terminal-uniform bound

    P_S <= min(1, 2^-31 q^(Ls-R-Delta)
                          C(Ls-D,Delta)).                (7)

No fractional exponent has been applied in (6)-(7). The binomial
coefficient uses the worst terminal value a=0; this is an upper bound,
not a terminal reset. The inequality is strict whenever the displayed
uncapped bound is at most one.

At the selected profile R=Ls=786432, the threshold is
D_0=W-Ls=157286. For D<=D_0, (7) is already below 2^-31. For
D=D_0+Delta>D_0 it becomes

    P_S < 2^[-31-32Delta+log2 C(629146-Delta,Delta)].      (8)

For example, put n=629146-Delta and k=Delta. At Delta>=32768, n<=20k.
The binomial theorem and the integer inequality 17^10<2*16^10 give

    C(n,k) <= (17/16)^n 16^k < 2^(6k).

Thus the exponent is below -31-26Delta<=-851999 throughout this death
range. This is a hand-off for
an exact death-conditioned count, not a claim that its remaining history
census has been bounded.

For a fixed D and terminal a with Delta>0 and m>0, the optimized weighted
Chernoff marker instead satisfies qx=m/Delta. At Delta=32768 and a=0,
this is qx=563610/32768, whereas the current marker has qx=65536. Exact
binomial extraction saves only logarithmic factors over an appropriately
optimized Chernoff bound for those fixed statistics. Indeed m is a mode
of Bin(i,m/i), so its probability is at least 1/(i+1): the optimized
coefficient envelope differs from C(i,m) by at most a factor i+1. The
large potential gain comes from retaining D before selecting x and taking
a fractional power; coefficient extraction alone does not supply an
exponential gain.

The full distance proof remains open. These derivations justify stopping
the two weak MDS refinements at the current markers and examining exact
output coefficients or genuinely stronger joint rank information next.
