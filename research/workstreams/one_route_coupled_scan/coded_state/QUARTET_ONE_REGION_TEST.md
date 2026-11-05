# Bounded full-region test: composition remains unverified

2026-10-31. The construction is unchanged. The target remains 45% distance
and 40 bits of margin at K=2^20. This test does not close that target.

The question is whether the two-update analysis can be composed through
one real region without dropping population transitions or replacing its
inactivity quota by a global average. A region has 16 updates of width512
and exactly2100 original-inactive inputs. The seven pending types are
Z,V0,V1,B,L0,L1,L2, as in the [local transfer](QUARTET_CLOSED_TRANSFER.md).
Internal region boundaries do not apply a terminal test or reset these types.

## Operator tested

For incoming population n and outgoing population m, let R(n,m) be the
positive closed-set transfer through all16 updates, with the exact regional
quota and its binomial normalization. This is a moment operator, not a
transition probability or a code-failure probability. Write

    Phi_h(n) = product_i h_i^(n_i),
    B_h(n) = sum_m R(n,m) Phi_h(m) / Phi_h(n).

The retained weights have base-two logarithms

    (0, -.022766994810343107, -.03148991732888965,
     -.038897622883545944, -1.9738559013803745,
     -1.4394120075915136, -1.4317249625558666).

For the exact lower witness below, h is the positive binary64 vector
constructed by the helper, interpreted as exact rational numbers. The
moment weights are a=2^(707/1024), d=2^(-1905/1024), and y=1/8.
There is no new failure allowance and no new assertion about minimum distance.

## Full-width diagnostic and why it was rejected

`quartet_one_region_probe.py` samples the exact positive transfer by
sequential importance sampling, including the complete sampled population.
It never projects the population into a small neighborhood.

Within an update, suppose n labels remain. The probability of an incoming
quartet histogram m is product_i binom(n_i,m_i)/binom(sum_i n_i,4).
If R inactive inputs remain among V unprocessed regional positions, the
probability of each particular inactive subset of cardinality j in the
next quartet is

    binom(V-4,R-j) / binom(V,R).

It is NOT the probability of cardinality j, which has an extra binom(4,j).
The kernel already counts choices of coordinate labels. Multiplying this
probability into each local atom and decrementing R retains the exact
regional quota. The outgoing population becomes the next update's input.
At a region boundary only R resets. No terminal test is applied.

Particles sample atoms in proportion to their local positive weights and
are resampled after each quartet. In exact arithmetic, the product of
mean incremental weights is an unbiased normalizing-constant estimator.
Its logarithm is not unbiased, and no finite estimate is an upper bound.

The estimates failed two useful controls:

| Incoming population | Particles / seed | Region-one log2 estimate |
|---|---|---:|
| (104,176,72,32,16,16,96) |128 /1| -1220.60 |
| same |1024 /2| -856.91 |
| all Z |1024 /3| -753.07 |
| all Z |1024 /4| -862.77 |
| all Z, flat h |512 /5| -2185.76 |

These negative values are NOT evidence of contraction. Increasing the
particle count does not stabilize them. More decisively, exact lower
bounds below contradict their apparent scale. Local importance sampling
misses rare, heavily weighted complete histories. It cannot be used to
decide whether the full-width operator composes usefully.

## A rigorous full-width transition with growth

Start with all512 coordinates of type Z. Fix a suffix of L coordinates,
aligned to updates, with at least three updates. Given the original
inactive set I, pay exactly P=I intersect suffix. All suffix inputs are
then allowed. The canonical support selects the entire suffix and ends
with all512 coordinates of type L2.

This paid set is closed. Every payment is selected and thus touched;
an inactive coordinate outside the suffix is disallowed and cannot be
touched by a feasible support. The canonical score is at least lambda L,
where lambda=707/2560: selecting only the full suffix is already feasible,
and it incurs no terminal death. Additional canonical runs cannot lower
the maximizing score. This is a subset of the actual transfer, not an
extra terminal test on the operator.

Consequently the single matrix entry satisfies

    R(all Z,all L2) Phi_h(all L2)/Phi_h(all Z)
      >= h_L2^512 a^L
         * sum_c binom(L,c) binom(8192-L,2100-c) y^c
         / binom(8192,2100).

The right side keeps the exact quota and is independent of the quartet
partitions. Different suffix lengths are not summed: their paid sets can
overlap. No global paid-count cap is imposed on this regional operator;
the eventual capped whole-stream proof must retain that additional constraint.

At L=4096 (eight suffix updates), rational arithmetic proves that the
right side is at least **2^75**. The helper rounds a downward to a dyadic
rational and checks a_lower^1024 <= 2^707 exactly. It evaluates the entire
hypergeometric sum as an integer over a common denominator. The displayed
integer exponent is verified by exact rational power-of-two comparisons.

Thus a one-region norm bound below1 with the retained h is impossible.
This does NOT rule out a useful whole-code bound: some regional growth
can fit the global threshold, and this particular transition cannot
repeat from its own source. Its outgoing population is all L2, not all Z.
It does not establish a positive cycle or a low-distance codeword.

A simpler independent control uses only the empty paid set. For every
original inactive set and every partition sequence, its terminal-accepting
canonical path is included in the regional operator. All completed runs
have nonnegative score because 3*707>1905. With J inactive positions in
the final update, its weight is at least (2/5)^(512-J), using a>=8/5 and
every h_i>=1/4. Averaging the exact hypergeometric J gives B_h(all Z)>=2^-442.
With flat h, the same argument gives B_1(all Z)>=1. Even these weaker
controls reject the particle estimates.

## Complete small-instance composition

An exact integer computation uses width4, all16 updates, quota16 out of64
positions, and rational weights a=13/8,d=9/32,y=1/8. It retains every
reachable endpoint population and carries them through a second region.
The weights are a small-instance test fixture, not a full-width certificate.

- One region: log2 row sum lies in [12,13), display12.6465572831.
- Two regions with carried populations: log2 row sum lies in [18,19),
  display18.1429876176, not twice the first result.
- There are75 endpoint populations and a peak of1115 population/quota states.
- Arithmetic uses a common integer scale; no floating rounding enters
  the exact operator calculation.

The particle probe also underestimates this small control: at1024 particles
it reports8.52 and10.07 instead of12.65 and18.14. The exact calculation,
not agreement between particle runs, is the reference.

Another718 normalization checks and the small negative controls for a
global quota and an internal terminal test pass. There are1232 exact
small suffix checks of closure, score, and the all-L2 endpoint.

## Decision and next test

The algebraic interface composes in the complete small calculation.
The full-width diagnostic does not establish a composable upper bound.
Reject the particle estimates rather than interpreting them as successful
closure. There is no construction counterexample or basis for redesign
from this test alone.

Next use deterministic upper bounds on transitions between boundary
population cells. Begin with the all-Z to high-selected transition and
its continuation, keeping the quota counter and accounting for every exit.
Use the directed suffix entry and the complete small operator as mandatory
controls. A reported one-region upper bound below either control is wrong.
Do not repeat local particle sampling or discard exits from a finite set.

All runs were file-free on Peach with180-second deadlines and at most1024
particles. No native benchmark, construction, production, paper, or frozen
helper changes were made. Two early diagnostic runs briefly overlapped
after a command returned a live session; both completed. Subsequent runs
were serial. Their runtimes are not performance measurements.
