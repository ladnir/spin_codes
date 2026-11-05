# Search the actual route and pairing geometry

The [rank criterion](RANDOM_PAIRING_BIRTH_DIAGNOSTIC.md) gives a sufficient
test for a genuine distance obstruction. This pass applies that test to
sampled full-size routing and matchings. It also tests a joint geometric
moment. Neither test closes the proof or finds a full-size obstruction.
The construction and accepted pair arithmetic remain unchanged.

## Search objective

Use K=2^20, N=2^21, R=256, M=8192, and512 state coordinates. Sample each
regional route uniformly and each state matching independently and
uniformly, with seed20261027. The samples are fixed during the search.
No field matrices are needed for the coefficient-independent rank test.

For proposed half-full state masks, let V contain the coordinates touched
by a live incoming pair or the next mask. Let c_i count row i's positions
in V, and define d=sum_i(c_i-128)_+. Let J count crossing endpoints outside
the consecutive-state union, and D count deaths. Choose

    z=max(0, |union of consecutive supports|-964689).

The tested gap is d-J-D-z. A positive gap, with enough intersection
positions to select z zeros, would establish a nonzero message of weight
at most964689. A nonpositive gap is inconclusive, not a distance bound.

To search efficiently, temporarily select1974 target rows T. For a fixed
T, minimize its coverage in V plus D. The parity deficit satisfies

    sum_i(128-c_i)_+ >= 128|T|-sum_(i in T)c_i.

The search alternates this objective with selecting the1974 least-covered
rows. It keeps each state mask half-full and constant within blocks of
updates. Within a block, it tries swaps between active and inactive
coordinates. The exact cost includes the pairing interaction between the
two swapped coordinates; summing single-flip costs without this correction
would give incorrect descent steps. Independent small tests check it.

## Selected runs

The constant-start run used block lengths64,16,4, six sweeps, and at most
four accepted swaps per block per sweep. It completed in18.34s. The best
rank gap remained the initial -83887, at zero turnover. That result alone
would be a weak search of the intermediate-turnover region.

A second run therefore used independently sampled half-full masks in
four-update blocks. Its25s budget ended during sweep3. The best/final
history had:

| Quantity | Value |
| --- | --- |
| Deaths/swaps D | 89408 |
| Outside crossings J | 412590 |
| Shortened outer dimension d | 501998 |
| Required selected output zeros z | 173295 |
| Rank gap | -173295 |
| Minimum row coverage | 135 |
| Rows with coverage below128 | 0 |
| Mean coverage of the1974 least-covered rows | 155.096 |

No positive gap was found. Every row in that selected history remains
above the parity threshold, so its parity deficit is zero. This contrasts
with the low-coverage rows dominating the failed fractional moment.
The search is bounded, uses one setup seed, and need not reach a global
optimum. It neither proves that the bad geometry is impossible nor
certifies the actual distance.

## Joint matching and route moment

Fix a target row set T of size s and a support history with B deaths.
Write C_T for the coverage of T in V. At constant state weight256,

    |V|=A+B+J, where A=4096*256.

An exact local polynomial retains J and B jointly. For a current support
of size w=256 and k crossing pairs, introduce markers x for deaths, u for
next-state weight, and h for outside crossings. The polynomial is

    (x+u)^w (h+u)^k (1+u)^(256-k).

An incoming active coordinate either dies or stays. A crossed inactive
coordinate either becomes a birth or remains an outside crossing. An
uncrossed inactive coordinate either becomes a birth or stays zero.
The coefficient of u^256 therefore sums the exact local support choices.
Average over the previously derived P_w(k), then bound this coefficient
using any positive u. Raising the resulting bound to4095 counts the
remaining steps; initial activation has binom(512,256) choices. A positive
x then bounds the coefficient for total B.

For routing, set p=s/M and introduce a positive conservation marker t.
For 0<v<=1, define

    a0=1-p+p*t, a1=1-p+p*t*v, h=a1/a0.

Conditioning Bernoulli target labels on exactly s labels per region gives

    E_route[v^C_T | V] <= (M+1)^R t^(-Rs) a0^(RM) h^|V|.

Thus the exact matching polynomial can be composed with the routing bound
without replacing crossing endpoints by independent coins. The helper
uses positive markers for both coefficient extractions; its numerical
optimizer proposes weights, not a certified optimum.

For a positive rank gap, C_T must be at most
A-128M-z+128s-1. The structural screen includes the union over all
binom(M,s) target sets and all histories with B deaths. A separate field
screen uses the capped fractional-moment count and its output-zero choice
factor; it restricts v>=2^-32. For the field screen, each fixed s is only
one term of a sum over deficit-row subsets, not an exact active-message
class or a whole-code bound.

Selected FLOAT log2 bounds:

| R | B | Target rows s | Positive-rank-gap event | Field-count term |
| --- | --- | --- | --- | --- |
| 256 | 131072 | 1974 | +178752 | +309953 |
| 256 | 32768 | 1200 | -4285 | +68931 |
| 192 | 131072 | 4096 | +233217 | +233222 |

The hardest selected geometry remains unresolved. The negative structural
screen in the middle row is not a directed replay and does not establish
distance: absence of a positive dimension gap is insufficient for that.

## Next discriminating test: a cut-based lower bound

For fixed target rows, coverage plus deaths is a submodular cost in the
binary state-support variables. A target coordinate contributes the OR of
its two incoming pair coordinates and its outgoing coordinate. A death
contributes x_prev*(1-x_next).

This admits a directed cut representation. Represent each OR by an
auxiliary node of unit cost and implications from its three variables.
Represent each death by a unit-capacity edge from the previous state
variable to the next. A negative unary price -mu*x favors activity.

For a window with n_state variables, every half-full history has exactly
n_state/2 active variables. Consequently

    min_all_assignments(cost-mu*activity) + mu*n_state/2

is a lower bound on the cost of every admissible history in that window.
Boundary terms can be dropped for a valid lower bound. Summing disjoint
windows gives a lower bound for the chosen target set without relying on
local-search convergence. Small windows keep memory bounded.

This cut computation has not yet been implemented or run. It would test
whether the observed gap is an optimization artifact. Bounds for one fixed
target set would still not cover every target set or establish the full
distance proof; that later step must retain the setup probability and
field constraints.

## Validation and preservation

Tests pass2240 exact swap-delta checks and140 descent checks, eight exact
joint-polynomial checks, and30 exact conserved-target moment checks.
The search recomputes the shortened dimension, crossings, deaths, and
rank identity from complete masks rather than trusting incremental costs.
The full-size runs were serial, guarded by20GiB free space, and produced
no raw files. The second search checked its25s budget and returned after
25.14s; the moment screen took1.25s. No native encoding timing, construction
change, paper edit, commit, or promotion. Prior certificates remain intact.
