# Realized-route obstruction from shortened dimensions

Two deterministically replayed routes admit a nonzero output below 10% weight
for **every fixed assignment of the four-state inner coefficients**. The
witnesses use different outer-coordinate supports in different source rows.
Their shortened dimensions exceed the cost of forcing zero state after each
selected window.

These are deterministic failures of the two realized routes. They do not by
themselves establish a lower bound on failure probability for uniformly random
regional permutations. The bounded search also does not compute a global
maximum of the dimension objective.

## Construction and obstruction criterion

Consider \(L=47680\) independent GRS[32,22] source rows. A fixed regional route
maps the rows to positions independently in each of the 32 outer-coordinate
regions. The stream length is \(N=32L=1525760\). Each region boundary is a
packet boundary for packet width sixteen.

In region \(j\), choose a contiguous window \(U_j\) whose start and length are
multiples of sixteen. Fix these windows and the realized route. For row \(i\),
let \(S_i\subseteq\{1,\ldots,32\}\) contain the outer coordinates routed into
their selected windows, and put \(c_i=|S_i|\). Different rows may have different
sets \(S_i\).

The subspace of one GRS[32,22] row supported within \(S_i\) has dimension
\(\max(0,c_i-10)\). The row spaces are independent. Thus the source subspace
whose routed input is confined to the windows has dimension

\[
 D:=\sum_{i=1}^{L}\max(0,c_i-10).
\tag{1}
\]

Independent nonzero input scalars do not change this dimension or the support
restriction. Fix any assignment of all inner coefficients. The inner must be
linear and zero-preserving, have four field symbols of state, and start in
state zero. The Cauchy diagonal-refresh construction has these properties.
Write \(\operatorname{Enc}_{\theta}\) for the resulting encoder with this
coefficient assignment \(\theta\) and the fixed route.

Impose the four zero-state equations immediately after each window. These are
at most \(4\cdot32=128\) homogeneous equations on the source subspace in (1).
Their common kernel has dimension at least \(D-128\).

If \(D>128\), choose a nonzero source in this kernel. Before the first window,
both the input and state vanish. After each window, the imposed state is zero;
the input remains zero until the next window. Induction along the stream
therefore confines the output to \(\bigcup_jU_j\). The causal identity
feedthrough makes the inner injective, so this output is nonzero. Consequently,

\[
 D>128\quad\Longrightarrow\quad
 \exists u\ne0:\quad
 0<\operatorname{wt}(\operatorname{Enc}_{\theta}(u))
   \le\sum_{j=1}^{32}|U_j|.
\tag{2}
\]

Equivalently, for each fixed coefficient assignment there exists a bad source;
the source may depend on that assignment. The proof does not assert one source
that works for every assignment simultaneously. It also does not require the
128 exit equations to be independent.

The strict bad-output threshold is

\[
 e=\lceil N/10\rceil-1=152575.
\]

The initial equal windows have length
\(16\lfloor e/(32\cdot16)\rfloor=4752\), for total length 152064. Using
\(\lfloor N/10\rfloor\) instead would give length 4768 and total 152576,
which does not meet this strict threshold. Only admissible windows are used
in the recorded witnesses.

## Exact objective and bounded search

The implementation is
[`window_dimension_search.py`](theory/window_dimension_search.py). A route
array maps positions to row labels. Its windows determine every \(c_i\)
exactly; the search uses no floating-point objective.

For a coordinate-ascent step, remove one region's current window from the hit
counts. Adding row \(i\) then increases (1) exactly when its remaining count
is at least ten. Prefix sums of these indicators along the region's
permutation score every packet-aligned start. The search takes a strictly
improving start and retains the current start on a tie. It reconstructs the
exact dimension after every sweep.

Each route first receives sixteen deterministic-seed restarts. Each restart
has an initial ascent and three perturbations of four region starts. Every
ascent has a limit of 24 sweeps. These settings give 64 bounded local searches
per route.

Starting from the best equal-window result, 512 proposals transfer length
between two regions. Transfer sizes are 16, 64, 256, or 1024 symbols. Each
proposal keeps all windows nonempty, preserves the total length, and runs at
most twelve coordinate-polishing sweeps. The search accepts only
nondecreasing exact dimension; ties may change the layout. This remains a
bounded heuristic, not an optimization certificate.

Routes use NumPy PCG64 with the route seeds below. The start-search seed is
193711, and the length-exchange seed is 9173. The recorded generator version
is NumPy 2.4.4. Full route hashes detect incompatible replay behavior.

## Realized witnesses and padding

All listed windows have total length 152064. The equal-window searches found
no dimension above 128; the length exchanges did.

| Route seed | Equal-window dimension | After length exchanges | After excluding final 18 rows | Kernel dimension at least |
|---:|---:|---:|---:|---:|
| 0 | 109 | 123 | not needed | no witness |
| 1 | 111 | 136 | 136 | 8 |
| 2 | 110 | 133 | 132 | 4 |

The kernel lower bounds use the conservative 128 exit equations and the
tail-excluded dimensions. The seed-one window lengths range from 2240 to
6416. The seed-two lengths range from 2464 to 7936. Both sets remain
packet-aligned and inside their respective regions.

The padded geometry has \(22L=1048960\) source symbols, whereas the nominal
requested size is \(2^{20}=1048576\). To avoid relying on row-major padding,
set the last eighteen complete source rows to zero. The remaining complete
rows use

\[
 22(47680-18)=1048564\le2^{20}
\]

source symbols. The seed-one dimension remains 136. For seed two, excluded
row 47669 had eleven hits, so its one-dimensional contribution is removed;
the remaining dimension is 132. This restriction still leaves a nonzero
kernel in both cases. A different padding layout must be checked against
its actual linear source constraints; the full padded-code obstruction does
not need a padding convention.

The literal compact fixtures are
[`window_dimension_witnesses.py`](theory/window_dimension_witnesses.py).
They contain the 32 starts and lengths, positive-dimension row counts, and
expected hashes and dimensions. They do not contain the full permutations.
Their verifier regenerates the routes and independently counts hits using
concatenation and `bincount`, without running the optimizer. Both fixtures
passed this replay.

An independent reviewer also regenerated both routes and counted membership
through inverse permutations, without using the fixture verifier's counting
method. That replay confirmed the tail-excluded dimensions 136 and 132 and
the surviving dimensions eight and four.

## Verification and scope

The nine tests in
[`window_dimension_search_tests.py`](theory/window_dimension_search_tests.py)
cover exact hit counts, prefix-score identities, monotonicity, tie handling,
bounded restart replay, and length-budget preservation. A tiny exhaustive
search checks coordinatewise optima without assuming they are global optima.
A separate GRS[4,3] example over GF5 constructs shortened row bases and checks
the exit-kernel confinement argument for twelve arbitrary linear inner
assignments. The frozen large fixtures are also verified without search.

```sh
python -B research/workstreams/char2_field_spin/theory/window_dimension_search_tests.py
python -B research/workstreams/char2_field_spin/theory/window_dimension_witnesses.py
```

The original bounded search can be replayed with:

```sh
python -B research/workstreams/char2_field_spin/theory/window_dimension_search.py --route-seeds 0 1 2 --restarts 16 --perturbations 3 --max-sweeps 24 --length-exchanges 512 --exclude-tail-rows 18
```

The obstruction depends on the routed shortened space and the four-symbol
state dimension, not on the particular state refresh. Replacing diagonal
refresh with arbitrary linear four-by-four state matrices cannot repair these
two route realizations. The same dimension argument applies over scalar
extensions of GF256.

No native encoder, benchmark, remote run, or distribution-level failure bound
was produced. The next proof gate is to quantify the frequency of this
dimension excess under ideal regional permutations, or change the relevant
outer/state/routing geometry. Additional coefficient-only rank estimates
cannot eliminate these realized witnesses.
