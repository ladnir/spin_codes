# Extending the distance certificate to K = 2^30

We have certified occupancies Q=1,...,32 at K=2^30. Their combined
first-moment bound is less than 2^-40.51009417. This is not yet a full
distance certificate: occupancies 33 through 8,388,608 remain uncertified.

The objective is a full failure bound below 2^-40 for the same fixed
BCH [256,128] outer and selected RM2Sub t64_s20 map. The setup distribution
remains the fresh, independent nonzero multiplier model in
`T64_S20_FULL_CLOSURE.md`. Row permutations and region permutations are
independent; the state persists across regions. We change only the size:

- message length K=2^30;
- number of outer rows L=K/128=8,388,608;
- output length N=256L=2^31;
- bad-output cutoff H=floor(N/10)=214,748,364.

For each Q, let Z_Q count messages with exactly Q nonzero outer rows whose
output weight is at most H. All expectations and probabilities below are
over the shared random setup. A full certificate needs bounds U_Q with

\[
\Pr[\exists x\ne0:\operatorname{wt}(E_\theta(x))\le H]
\le \sum_{Q=1}^{L}\mathbb E[Z_Q]
\le \sum_{Q=1}^{L}U_Q < 2^{-40}.
\]

This target concerns distance under the stated setup, not the complete SPIN
protocol or a pseudorandom replacement for the multiplier stream.

## Exact partial certificate

The new receipts use the existing certified BCH spectrum constraints, not
the modeled spectrum. No complete BCH-256 spectrum is assumed.

| Covered class | Margin of its certified upper bound |
|---|---:|
| Q1 | 40.51009419 bits |
| Q2 | 67.55814558 bits |
| Q3 | 105.91788896 bits |
| Q7 | 261.43580945 bits |
| Q8 | 300.09152107 bits |
| Q32 | 937.72715636 bits |
| Sum over every Q1..Q32 | 40.51009418 bits |

Every integer Q1..Q32 is covered. The table displays selected rows only.
The exact rational sum and remaining budget are recorded in
`generated/k30_partial_coverage_v1.json`. The uncertified occupancies must
fit within the remaining budget

\[
2^{-40}-\sum_{Q=1}^{32}U_Q
\approx 0.297823399823\,2^{-40}
\approx 2^{-41.74747098}.
\]

These decimals summarize rational inequalities; they are not the arithmetic
used by the ledger. Q1 consumes almost all of the covered contribution.

`certify_k30_q1.py` evaluates four exact rational tilt witnesses with
256-bit Arb. Its 512-bit replay uses polynomial matrix powering instead
of paired degree-zero/degree-one matrices. It checks all 92 shell
coefficients and the reused BCH weighted inequality. The outer inequality
is reused by coefficient domination, not by substituting an old probability.

`k30_sparse.py` certifies Q2..Q7. The selected B-kernel has minimum distance
8, so every nonempty epoch input in these classes has nonzero syndrome.
`k30_kernel.py` includes zero-syndrome inputs and certifies Q8..Q32.
Both use degree-truncated region polynomials, positive matrix powers,
and the existing adaptive bound over 13 outer-weight bands. Their replays
use 512-bit Arb and reverse the direction of binary epoch powering.
The replays share the mathematical transfer formulas; they are not
independent implementations of the entire proof.

All producers verify the 804-file migration manifest before using the
frozen map, outer caps, and weighted certificate. The partial audit checks
source hashes, replay provenance, parameter agreement, exact range sums,
and disjoint coverage. It does not replace numerical replay.

## Zero-syndrome extension

The transfer distinguishes zero state Z, arbitrary nonzero state D,
uniform nonzero state U, and nonzero state with density at most kappa
relative to uniform, denoted L. Here M=2^20-1 and kappa=M/(M-1).

Fix an epoch input weight j and let X be uniform among weight-j vectors.
Write beta_j=Pr[BX=0], obtained exactly from the B-kernel spectrum.
For z in (0,1), let m_j be the expectation of z^wt(X+Aq) when q is
uniform nonzero and independent of X. Hypergeometric averaging of the
exact A spectrum gives m_j. Let A_w^inner count its nonzero words of
weight w. Define the pointwise bound

\[
d_j(z):=z^{\min_{w:A_w^{\rm inner}>0}|w-j|}.
\]

For the entering classes D,U,L, upper bounds on the unrestricted emission
moment are respectively d_j,m_j,kappa*m_j. For any such bound v, use

\[
b_0:=\min(\beta_j d_j,v),\qquad
b_1:=\min((1-\beta_j)d_j,v).
\]

The weighted zero-syndrome mass is at most b_0; the weighted nonzero-syndrome
mass is at most b_1. This does not assume independence between syndrome
and emitted weight. The two upper bounds may overcount their sum.

For an entering nonzero state, the zero-syndrome branch refreshes the next
state to U. The nonzero-syndrome branch reaches Z with probability 1/M.
Its surviving nonzero state has density at most kappa relative to uniform.
Thus the outgoing masses to U,Z,L are bounded by b_0,b_1/M,b_1(M-1)/M.
The Z row has exact masses beta_j*z^j and (1-beta_j)*z^j to Z and D.
At j=0, the implementation uses the sharper existing refresh transfer.

Exact toy tests enumerate all input weights of a self-orthogonal [8,4]
map. They check arbitrary states, the uniform law, and every uniform law
with one nonzero state excluded. Additional tests check polynomial powering
and agreement with the sparse transfer below the kernel distance.

## A scalable route for dense occupancies

The small-Q polynomial calculation alone is not a practical way to cover
millions of occupancies. A comparison with independent input bits gives
a bound whose matrix size does not depend on Q.

First, let J be the sum of L independent Bernoulli variables, with possibly
different probabilities and mean Lr. Then, for every j,

\[
\Pr[J=j]\le (L+1){L\choose j}r^j(1-r)^{L-j}.
\]

For 0<r<1, the arithmetic-geometric mean inequality bounds the moment
generating function of J by that of Binomial(L,r). Chernoff's argument
therefore gives exp(-L D(j/L || r)) as a point-mass upper bound.
The binomial point mass is at least that quantity divided by L+1:
under Binomial(L,j/L), j is a mode, whose mass is at least 1/(L+1).
The deterministic endpoint cases r=0 and r=1 satisfy the inequality too.

After uniformly shuffling the L bits, all vectors of a given weight have
equal probability. The same factor L+1 consequently bounds their density
relative to independent Bernoulli(r) bits. Applying the comparison to all
256 independent regions costs (L+1)^256.

To connect this fact to the outer code, fix a Bernoulli probability p_g
for each outer-weight band g. Let U_w be the existing certified cap on
the number of outer words of weight w, and define

\[
\Gamma_g:=\max_{w\in g}
\frac{U_w}{{256\choose w}p_g^w(1-p_g)^{256-w}}.
\]

The all-one band uses p_g=1 and Gamma_g=1. The probability comparison
for this deterministic band is interpreted directly, without division
by zero. Every ordinary band uses 0<p_g<1.

For a fixed assignment of bands to Q active rows, the product Bernoulli
reference dominates the summed row distributions after multiplication
by the product of the Gamma_g values. Its count in each region has the
Poisson-binomial distribution above; inactive rows contribute zero bits.
Let x=Q/L and let r be the average reference bit density over all L rows.

Let h be the upper concave hull of the points (p_g,log Gamma_g).
Concavity bounds the log of the product cost by Q*h(r/x). This covers
every band assignment, including mixed assignments. There are at most
13^Q assignments and binomial(L,Q) choices of active rows.

Let T(r,z) be the nonnegative epoch transfer averaged over independent
Bernoulli(r) input bits. The existing three-state transfer in
`tightened_occupancy.py` is sufficient here. With e_Z selecting its initial
zero state and 1 denoting the all-ones terminal column, the resulting bound is

\[
U_Q\le {L\choose Q}13^Q(L+1)^{256}
\sup_{r\in[x\min_g p_g,x\max_g p_g]}
\left\{e^{Qh(r/x)}
\inf_{\lambda>0}
e^{\lambda H}e_Z T(r,e^{-\lambda})^{N/64}\mathbf1\right\}.
\]

Each actual band assignment has one fixed r. It may use its own tilt
witness, which explains the order of the supremum and infimum. The
comparison changes only the bounding distribution, not the encoder.

### What the dense screen does and does not establish

`k30_dense_unrestricted_screen.py` evaluates this expression on sampled
occupancies, sampled mean densities, and a finite tilt grid. It uses
binary64 arithmetic, so it supplies discovery evidence only.
At each tested occupancy it chooses one set of row probabilities before
taking the worst sampled density. The retained receipt is
`generated/k30_dense_unrestricted_screen_v1.json`.

Every tested power of two from Q=32,768 through Q=8,388,608 has a positive
diagnostic margin. Selected worst-grid margins are approximately:

| Occupancy | Diagnostic margin, not certified bits |
|---|---:|
| 32,768 | 17,764 |
| 131,072 | 384,335 |
| 1,048,576 | 7,819,258 |
| 8,388,608 | 3,950,131 |

The grid does not bound the supremum between sampled densities or cover
the intervening occupancies. These numbers cannot enter the partial ledger.
The same screen fails to close its sampled points Q=1,024 through 16,384.

Two preceding screens are retained as source-level experiments.
Midpoint row probabilities give poor bounds. Tilting those probabilities
helps, but splitting the nonzero emission moment into separate kernel and
nonkernel upper bounds loses too much at dense occupancies. The existing
three-state transfer retains the unrestricted emission moment and avoids
that particular loss. The original midpoint screen also has an unresolved
nonfinite diagnostic at the deterministic all-one endpoint; it is not an
audit input. The retained final screen evaluates that endpoint exactly.

## Reproduction and next work

Restore the inputs using `MIGRATION.md`. Run the following sequentially,
using fresh output names because receipts are write-once:

```text
python -B workstreams/bch_rm2sub_bridge/test_k30_certificates.py
python -B workstreams/bch_rm2sub_bridge/certify_k30_q1.py --output workstreams/bch_rm2sub_bridge/generated/k30_q1_new.json
python -B workstreams/bch_rm2sub_bridge/certify_k30_q1.py --output workstreams/bch_rm2sub_bridge/generated/k30_q1_new.json --verify
python -B workstreams/bch_rm2sub_bridge/k30_sparse.py --output workstreams/bch_rm2sub_bridge/generated/k30_sparse_new.json
python -B workstreams/bch_rm2sub_bridge/k30_sparse.py --output workstreams/bch_rm2sub_bridge/generated/k30_sparse_new.json --verify
python -B workstreams/bch_rm2sub_bridge/k30_kernel.py --output workstreams/bch_rm2sub_bridge/generated/k30_kernel_new.json
python -B workstreams/bch_rm2sub_bridge/k30_kernel.py --output workstreams/bch_rm2sub_bridge/generated/k30_kernel_new.json --verify
```

The replay commands create adjacent `_replay.json` receipts. Change the
`_new` suffix when repeating these commands on the same checkout.
`audit_k30_partial.py --output <fresh-path>` validates the retained v1
receipts and reconstructs their exact partial sum.

Next on the certificate track: convert the dense expression into outward
interval bounds over both Q and r, then bridge the small/intermediate-Q
gap. The current small-Q method can extend beyond Q32, but its cost grows
with the polynomial degree. A joint input/output tilt is another candidate
for reducing the loss from replacing fixed row weights by independent bits.

Keep the estimator track separate. None of these new calculations assumes
a heuristic BCH spectrum, measures the true failure probability, or removes
the need to assess the low-shell model in `DUAL_TRACK_PROGRESS.md`.
