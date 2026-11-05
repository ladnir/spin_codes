# A whole-column bound for zero-coverage components

2026-10-30. This strengthens the read17 calculation in
[the component audit](QUARTET_COMPONENTS.md). It keeps the construction,
regional quotas, and free terminal boundary unchanged. It is a restricted
geometric certificate, not a whole-code distance certificate.

Fix arbitrary quartet partitions and an inactive-row set T of size2100.
Let S be a state support, A its total size, C its inactive-input coverage,
and D its number of deaths. Write

    F(S)=lambda*A-C-tau*D, lambda=1/4, tau=37/50.

The support may depend on the entire sampled route. We bound its maximum
over all S with C=0, then union-bound over all choices of T.

## Group complete coordinate sequences before applying Holder

A coordinate (t,j) is free if its own input is active and all four inputs
in its quartet at update t+1 are active. At the last update, only its own
input must be active. Exactly these coordinates may be selected when C=0.
Optimization separates into512 vertical coordinate sequences.

Set alpha=3lambda-tau=1/100. For column j, let W3_j and W4_j count windows
of three and four consecutive free coordinates, and set

    g_j=alpha*W3_j+(lambda-alpha)*W4_j.

A nonterminal free run of length l contributes max(0,lambda*l-tau), exactly
its contribution to g_j. A terminal run contributes lambda*l; the difference
is at most tau. Consequently

    max_(C(S)=0) F(S) <=512*tau+sum_j g_j.                 (1)

Temporarily replace the actual quotas by independent active-input bits with
probability p=1-2100/8192. An input bit at (t,k) affects its own column k
and the preceding free-status bits of the four columns in its quartet.
Column k is already among those four. Thus each input bit is read by at
most FOUR whole-column functions g_j, regardless of how the partitions
change between updates. The read4 Holder inequality yields, for b>=1,

    E b^(sum_j g_j) <= product_j (E b^(4g_j))^(1/4)
                    = M(b)^128.                         (2)

This groups overlapping windows before the inequality. Applying Holder
separately to every window loses the temporal dependence and gives read17.
No independence between columns is assumed in (2).

## Compute the one-column moment

For a fixed column, let A_t indicate that its own input is active, and let
V_t indicate that all four inputs of its current quartet are active. Under
the temporary product measure, the pairs (A_t,V_t) are independent across
updates and have probabilities

    (0,0):1-p, (1,0):p-p^4, (1,1):p^4.

The free bit is z_t=A_t*V_(t+1), except z_4096=A_4096. Keep the current A_t
and the current free-run length capped at three. This gives eight states.
A true free bit contributes weight b^(4alpha) when extending length two,
weight b^(4lambda) when extending a length at least three, and weight one
otherwise. A false bit resets the run length. The terminal step uses A_4096
without an additional quartet condition. The final total is M(b).

The directed recurrence rounds every nonnegative transition upward. State
values have80 fractional bits and coefficients64 fractional bits; Python
integers avoid overflow. Rational upper roots and directed logarithms
complete the tail computation. Exhaustive enumeration of1092 category paths
through six updates checks the recurrence. At three updates it also checks
the identity M=1+(b^(4alpha)-1)*p^9, including the terminal boundary.

## Return to the actual quotas

Both alpha and lambda-alpha are nonnegative. Hence b^(sum_j g_j), expanded
as a polynomial in active-input indicators, has nonnegative coefficients.
Reduce repeated indicator powers using x^k=x. For any d distinct positions
within one region, the actual probability that all are active is

    binom(8192-d,2100)/binom(8192,2100) <= p^d.

Regions are independent. Comparing each monomial therefore upper-bounds the
actual expectation by the product-measure expectation used in (2), without
a conditioning penalty. Quartet partitions were arbitrary throughout, so
the bound also holds after averaging their actual distribution.

Combining (1), (2), and the union over binom(8192,2100) row sets gives

    log2 Pr[exists T,S: |T|=2100, C_T(S)=0, F(S)>B]
      <= -(B-512*tau)*log2(b)
         +128*log2 M(b)+log2 binom(8192,2100).             (3)

The frozen choice b=2,B=15877 yields the directed upper bound

    -1670147885140353/27487790694400 < -60.7596.

`quartet_run_components.py --columns` reproduces the tests and certificate.
It writes no files, requires20GiB free, and has a30-second watchdog.

## What remains

The selected interaction components have additive A,C,D and score. The union
of all components with C=0 is itself a zero-coverage support. Therefore,
outside the event in (3), those components contribute at most15877 for
every maximizing support and every row set of size2100. To obtain the
sufficient whole-support offset18171 by this decomposition, it remains to
bound the sum of positive scores of components with C>0 by2294.

The existing fixed-geometry sample gives447.30 for that latter sum. This is
not an ensemble bound and is not used in (3). Components need not be
independent. Other inactive-row counts also remain required by the full
45%/40-bit goal. No additional allowance is inserted into the whole-code
failure budget until a complete argument uses this restricted certificate.

## Retain the field cutoff through a penalized objective

The [filtered-count obstruction](QUARTET_FILTERED_COUNT_OBSTRUCTION.md)
uses histories with C+D=290432, above the surviving field cutoff255463.
This suggests retaining the cutoff before selecting a geometric maximizer.
The following deterministic implication is valid without constraining the
maximizer itself. For eta>=0, define

    F_eta(S)=F(S)-eta*(C+D)
            =A/4-(1+eta)C-(37/50+eta)D,
    B_eta=18171-eta*255463.

If max_S F_eta(S)<=B_eta, then every actual history with C+D<=255463
satisfies F(S)<=18171. Indeed, add eta*(C+D) to the asserted bound.
Thus existing field exclusions can be used without assuming that geometric
maximization preserves a codeword's turnover. Canonical filters must be
derived for this new objective if it is used in a filtered count.

For C=0, the column argument changes only the death weight tau. For general
tau>0, put m=ceil(tau/lambda), alpha=m*lambda-tau. Replace W3,W4 by Wm,W(m+1)
and cap the hidden run state at m. The same identity and read4 argument hold.
At b=2, four bounded checks give:

| eta | Required B_eta | 60-bit zero-coverage bound | Remainder for positive coverage |
|---:|---:|---:|---:|
| 0 | 18171 | 15877 | 2294 |
| 1/250 | 17149.148 | 15678 | 1471.148 |
| 1/100 | 15616.37 | 15384 | 232.37 |
| 1/50 | 13061.74 | 15095 | -2033.26 |

These are directed sufficient offsets at one fixed moment parameter, not
optimal offsets. The old subcount's excess over B_eta is

    108.6-34969*eta.

At eta=1/250 it is negative31.276, so that subcount no longer consists of
violations. This does not prove that its weighted first-moment contribution
is small; its expected multiplicity remains large. The stronger penalties
also rapidly consume the available zero-coverage allowance. Therefore this
test establishes a valid proof interface and its tradeoff, not a successful
replacement certificate. It does not justify another large transfer run.

`quartet_run_components.py --penalty` checks3276 additional exact category
paths and12288 generalized run identities, then replays this table. Runs are
serial, file-free, and bounded by the existing watchdog. The next useful
test is a bound for positive-coverage components, either under the original
score or the small penalty; tightening the large penalty's scalar marker
alone does not address the missing component count.
