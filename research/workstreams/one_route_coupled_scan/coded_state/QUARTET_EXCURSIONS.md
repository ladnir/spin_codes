# Leaving the band and returning through zero

2026-10-26. Follow-up to [the composable band](QUARTET_HISTORY_BAND.md).
The construction and all previous checkpoints are unchanged.

## Two larger restricted families close

Retain the same K=2^20, rate-one-half construction over J=GF(2^32), with
RS[256,128] rows, regional routing, and the accepted field-quartet inner.
Setup is sampled once and fixed for all messages. Let W=943718, the bad
output threshold at 45% distance. Both statements below include every exact
inactive-row count r in [1850,2100], every such row set, all possible first
activation times, and every nonzero J-message in the specified family.

Call a transition regular if both adjacent states have weights in [216,232]
and exactly 24--40 coordinates leave their support. Every other transition
after the first activation is exceptional. This is a classification of
candidate histories, not an independence assumption about realized states.

1. **At most 25 exceptional transitions.** The initial nonzero state may
   have any weight. Exceptional transitions may have any incoming and next
   weights, including full resets, reactivation, and termination at zero.
   The failure probability for this family is at most

       2^(-3901196831485731/43980465111040) < 2^-88.70.

   Positions and types of all exceptional transitions are included in the
   union. A departure followed by a return generally consumes two slots.
   Each additional zero-to-zero transition also consumes a slot.

2. **Arbitrarily many short reset episodes.** Every nonzero state has weight
   in [216,232]. Adjacent nonzero states have 24--40 deaths per transition.
   After first activation, each run of zero states has length at most three.
   The number of these runs has no separate cap beyond the 4096-update stream.
   Termination at any allowed zero age is included. This family has failure
   probability at most

       2^(-13146227725709835/2199023255552) < 2^-5978.

These are two different families. Neither statement permits arbitrary
excursions together with arbitrarily many short resets. Their union can be
bounded by adding the two displayed bounds. Other row counts and the full
history space remain open. As before, these restricted support-history
statements are over J; only a whole-code J-distance bound transfers directly
to the extension-field encoder.

## A cheap envelope for all state weights

Use the positive transfer T and potential f(w)=X^w ell^(theta w) from the
preceding note. The local polynomial and all proof markers remain frozen.
Let G retain only regular transitions. The existing integer coefficient
calculation, extended by zero rows outside the band, gives

    Gf <= beta f,
    log2(beta) <= -11954019317728937/43980465111040
               = -271.802930859... .

To bound the full transfer F, including every state weight and turnover,
there is no need to build a 513-cubed coefficient table. For a canonical
incoming quartet support U of size a, let V be the next support, T_0 its
inactive-row input labels, b=|V minus U|, and d=|U minus V|. Let c be the
covered inactive-label count and P(ell) the same permitted-zero polynomial.
The local row sum after multiplying by f(next)/f(incoming) is

    L_a ell^(-theta a),

    L_a = sum_(V,T_0) p^|T_0| (1-p)^(4-|T_0|)
          (P(ell)/q^c)^theta X^b [X^-1 (q ell)^(-theta)]^d,

where q=2^32. All fractional powers are rounded upward using exact integer
comparisons. For any actual partition, the product over its 128 quartets
is at most the 128th power of the largest local row sum. Averaging the
partition cannot increase that upper bound. Therefore

    Ff <= eta f,
    eta = (max_(a=0,...,4) L_a ell^(-theta a))^128.

Fortieth powers compare the five rational/algebraic expressions exactly.
The maximum occurs at a=0, with

    L_0=(1+(1-p)X)^4,
    log2(eta) <= 3645960473881/68719476736 = 53.055707741... .

The directed local log bounds for a=0,...,4 are approximately
0.41450, -4.21878, -3.13669, -2.20279, and -1.27163. The empty incoming
quartet, not a particularly expensive nonzero quartet, controls this envelope.
These are weighted counting bounds, not physical reset probabilities.

## Paying for arbitrary exceptional transitions

For a suffix of n transitions with at most h exceptional positions, expand
the transfer into regular and exceptional edges. Replace an exceptional
operator by F; this only increases the count. Positivity gives, for every
ordering of the operators,

    G^(...) F^(...) ... f <= beta^(n-j) eta^j f

when j factors are exceptional. No commutation or state independence is used.
The complete count is consequently bounded by

    sum_(j=0)^h binom(n,j) beta^(n-j) eta^j.

Here 2048<=n<=4095 and h<=32. Since eta>=beta, the terms increase through
h. The replay uses the conservative upper bound

    (h+1) binom(4095,h) beta^(n-h) eta^h.

It also allows every initial nonzero weight and every terminal weight. Their
combined boundary factor is at most

    (1+(1-p)X)^512 / f(512).

The preceding outer prefactor, row-set union, and exact zero-prefix quota
cancellation remain present. All 2048 activation times are checked, and
binomial log-concavity bounds all 251 row counts by their endpoints. The
25-exception result includes these unions. At 26, this particular bound
becomes positive; that is a bound failure, not a code counterexample.

## Separating reset, reactivation, and waiting

The full-transfer envelope treats a reset as costly as the worst arbitrary
step. Its actual weight is much smaller. For a live quartet with a incoming
nonzeros and empty next support, local compatibility is exactly |T_0|<a.
Put r_0=1 and, for a>0,

    r_a=sum_(t=0)^(a-1) binom(4,t) p^t(1-p)^(4-t) q^(-theta t).

The reset row sum for incoming weight w is

    k_reset(w) = [z^w](sum_(a=0)^4 binom(4,a) r_a z^a)^128
                 / binom(512,w).

This coefficient is evaluated as a product in the Bernstein basis. Each
stage is a hypergeometric convex combination. A 512-bit fixed-point integer
vector, exact binomial denominators, and upward rounding preserve tiny
coefficients without underflow. The vector has only 513 entries.

After including the transfer and potential factors, define

    b = max_(w=216)^232 k_reset(w) q^(-theta w) X^-w ell^(-2 theta w),
    a = sum_(w=216)^232 binom(512,w) ((1-p)X)^w.

Here b bounds a band-to-zero edge and a bounds a zero-to-band edge. The
directed logarithms are at most -834.09978 and -311.10508, respectively.
The zero-to-zero edge has weight one. Use the simpler certified constants

    g=2^-271 >= beta,  a_bar=2^-311 >= a,  b_bar=2^-834 >= b,
    lambda=2^-270.

The following finite automaton counts every short-reset pattern. Its states
are the live band and zero ages 0,1,2. A reset enters age 0; a further zero
state increments the age. Every zero age may return to the band. There is
no age-2-to-age-3 transition. Thus it permits up to three consecutive zero
states, rather than three zero-to-zero transitions.

Give the band potential one and age i the potential

    psi_i=a_bar sum_(k=1)^(3-i) lambda^-k.

Every zero-age row satisfies the lambda-potential inequality exactly. The
band row satisfies it because

    g/lambda + a_bar b_bar (lambda^-2+lambda^-3+lambda^-4)
      = 1/2 + 2^-605 + 2^-335 + 2^-65 < 1.

The combined potential is f(w) on band states and psi_i on zero ages. Its
minimum is still f(232): even psi_2=2^-41 is larger. Therefore the same
free-terminal bound covers termination in any zero age. The original outer
prefactor, activation union, and row-count union then give the displayed
2^-5978 result.

All episodes refer to one shared outer message space and one fixed setup.
The proof does not resample outer rows, restart the row-set union, or assume
independent messages at reactivation. Outer shortening retains input zeros
throughout every quiet interval. Only the finite transfer representation
changes.

## Validation, interpretation, and next test

The source-only replay is

    python -B research/workstreams/one_route_coupled_scan/coded_state/quartet_excursions.py

It uses no optimizer and asserts both frozen final exponents. It also replays
the previous band bound unchanged. Checks include 40 exact noncommuting-path
bounds, 125 local occupancy products, 72 Bernstein coefficient cells, all 64
local reset masks, four zero-age potential rows, and 32 terminal-age bounds.
The prior frozen pair-defect helper's SHA256 remains unchanged.

Serial runs took about 11 seconds, with a 60-second watchdog, bounded arrays,
and at least 20 GiB free. No generated files, native benchmark, construction
change, library/paper edit, commit, or push was made.

The useful distinction is now explicit: reset and reactivation pay enough
in this weighted count to absorb several quiet updates. Short resets are
not an unresolved obstruction within the stated band. An additional zero
age fails this particular small-potential comparison, not the code.

Next retain regional inactive-row quotas and the RS coverage requirements
jointly with longer quiet intervals. Also cover other nonzero weights and
row counts. A global quiet loop of weight one cannot simply inherit the
live-band decay. Prefer a regional continuation bound that accounts for
the surviving outer message space over further local mixer changes or an
unbounded sweep of exceptional-step budgets.
