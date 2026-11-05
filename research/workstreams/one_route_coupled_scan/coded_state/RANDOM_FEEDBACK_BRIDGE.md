# Random feedback isolates a compact-state obstruction

Follow-up: the [256-state investigation](IDEAL_RESET_LADDER.md) leaves the
running-sum candidate open and adds certified ideal MDS-reset controls.
It distinguishes insufficient memory, weak reset algebra, and outer/route
concentration. The exclusions below remain unchanged.

The user requested an incremental path from a provable random inner to an
efficient structured inner. The first bridge replaces only the feedback
randomness. It removes the hole-enumeration problem, but an independent
dimension argument rules out the desired >49.1%/40-bit guarantee for the
64-symbol state under the current route distribution. The same obstruction
also applies to the structured inner. This is a genuine setup-failure lower
bound, not a failed proof envelope.

The structured-inner notes and their restricted certificates are preserved.
They do not claim a complete distance bound and are not contradicted by
the new bad trajectories. Do not continue trying to assemble them into
a whole-code certificate at the excluded target.

## 1. The first bridge and its exact inverse distribution

Keep F=GF(2^128), K=2^20, N=2K, the fixed RS[256,128] outer, and the 256
independent uniform regional permutations. Each region has M=8192 symbols.
Keep batch size b=128, state size s=64, compression B, and the independent
nonzero GF256 diagonal scalings R_t. Sample independent uniform matrices
T_t in F^(128 by 64), once during setup, independently of the route and R.

    Y_t = X_t + T_t Q_(t-1),
    Q_t = Q_(t-1) + B R_t Y_t,       Q_0=0.

Emit all N symbols and leave the terminal state free. The inner is linear
and invertible for every setup: reconstruct the state from Y, then X from
Y and T. For every fixed nonzero Q, T_t Q is uniform on F^128.

Fix an output y and all R. Its state sequence is independent of T. Let A
be the union of batches with nonzero incoming state. The inverse input is
uniform on y+F^A. For the routed outer C, MDS shortening and complementation
give

    Pr_T[inverse(y) in C] <= q^(d_C(A)-|A|)
                         = q^(-K+d_C(A^c)).             (1)

There are no holes inside a live batch. The distinction from the previous
inner is full feedback randomness in (1), not merely a uniform next state
for a fixed forward input. The previous inner already had the latter.

For example, if all zero-incoming-state batches lie in at most 128 regions,
then d_C(A^c)=0 for every route. Every outer block has at most 128 coordinates
there and cannot support a nonzero shortened RS word. Union over all
outputs of weight <=W=floor(.491N)=1029701 gives

    failure <= 2^N q^(W-K) = 2^-318848.                 (2)

No union over masks is needed: (1) has the SAME route-uniform bound q^-K
for every candidate in this class. Equation (2) includes arbitrary outputs
in a full-length live episode, not just 32 exceptional intervening batches.
This confirms that stronger feedback randomness removes the earlier hole
obstruction. It does not cover all zero-state patterns.

## 2. A deterministic reset-supported subcode

The following argument works for any causal linear inner with an s-symbol
state for which zero input and zero state produce zero output and state.
It does not depend on how much randomness its feedback uses.

Fix the route and a union U of c disjoint, batch-aligned scan intervals.
Restrict outer codewords to U. Their dimension over F is

    D_U = sum_(i=1)^8192 (t_i-128)_+,

where t_i is the number of block-i coordinates in U. Impose that the inner
state is zero at the end of each interval. These are at most s*c homogeneous
linear equations on that SAME shortened message space. If D_U>s*c, a
nonzero input remains. Its output is zero before, between, and after the
chosen intervals. Injectivity makes its output nonzero, so

    D_U>s*c  implies  d_min <= |U|.                    (3)

The final reset is an imposed constraint used to zero the suffix, not an
assumption that encoding requires a zero terminal state. Rank dependencies
only strengthen the argument. The witness may depend on the entire setup,
as it may in a minimum-distance claim.

## 3. An exact K20 obstruction under the actual route

There are 64 batches per region. Pair regions (0,1),(2,3),...,(254,255).
In each pair choose a suffix of the first region and a prefix of the
second; they meet at the regional boundary and form one interval.
Use 32 selected batches in 108 regions and 31 in the other 148. These
choices define c=128 disjoint intervals, with

    |U| = 128(108*32+148*31) = 1029632 < W.

For one outer block, independent regional permutations give independent
Bernoulli inclusions with probabilities 32/64 or 31/64. Write

    p_j = [z^j] (32+32z)^108 (33+31z)^148 / 64^256.

Then, exactly,

    E_route[D_U] = 8192 sum_(j=129)^256 (j-128) p_j
                 = 17728.325879... .                  (4)

Linearity of expectation suffices. Different outer blocks need NOT have
independent occupancies. Also D_U<=|U|/2=514816 pointwise. For R=s*c below
the expectation, these facts give

    Pr_route[D_U>R] >= (E[D_U]-R)/(514816-R).            (5)

Equation (3) holds for every inner realization on those routes. Thus (5)
lower-bounds the complete setup failure probability, regardless of the
feedback randomness or the compression coefficients.

| State dimension | Reset-rank cap | Rigorous failure lower bound, rounded for display |
|---:|---:|---:|
| 64 | 8192 | >0.0188 |
| 128 | 16384 | >0.00269 |

The exact checks also exclude every state dimension <=138 at the target
failure probability 2^-40. State139 is merely not excluded by this particular
expectation bound; it is NOT certified. Likewise this argument does not say
every route is bad, or exclude a separately specified rejection-sampled
setup distribution. Such a setup change would need its own analysis.

The resulting witnesses have relative weight at most 1029632/2097152,
approximately 49.09668%. This is below the target, despite being close to it.
It explains why a sequence of restricted proofs could not close the whole
64-state ensemble. The missing trajectories can repeatedly turn the state
off across many separated regional intervals.

## 4. The next controlled bridge

Do not thin the feedback randomness yet. The first change must retain more
memory. Merely allocating a larger vector with Q'=Q+B Y and fixed B cannot
help beyond rank(B): starting from zero, all reachable states lie in im(B).
With 128 input symbols per batch that effective dimension is at most128.

A simple next proof model uses b=s=256 and

    Y_t = X_t + T_t Q_(t-1),
    Q_t = Q_(t-1) + Y_t,

with independent uniform 256-by-256 T_t. The outer and regional route are
unchanged; each region now contains32 batches. The state is fully reachable,
nonzero quiet outputs are impossible, and the inverse law (1) still holds.
This is a proposed ideal-randomness checkpoint, NOT a selected fast kernel
or a complete certificate. No performance is inferred from its operation
count. Later, a multi-bank state could preserve a smaller batch, but that
requires a new transition proof rather than pretending all state coordinates
refresh in one step.

The first bounded screen combines forward reset probabilities with inverse
output enumeration. For c episodes occupying v batches, at least c-1
episodes reset before the free terminal. For a fixed input, uniform feedback
and an onto compression give a factor q^(-s(c-1)). In the inverse direction,
an episode of length L with output weight w has reset rank at least
min(s,ceil(w/L)): some output batch has that many independent compression
columns. With p=1-beta, 0<=beta<=1, its weighted exponent is at most

    max(L*p-1, b*L*p-s) <= b*L*p-min(s,1+2(b-1)*p).

Here a closed episode has L>=2. The free terminal episode receives no reset
charge. Nonzero quiet batches, when permitted in the earlier bridge, each
contribute at most b*p-s. These inequalities give an output-count interface,
not a claim of independent reset constraints on overlapping episodes.

The screen also retains that c intervals cut at most 2c partial regions.
It maximizes the one-block route moment over the feasible number of full
regions, equalizing only the partial inclusion probabilities. Selected
256-state points pass and others remain positive. For example, v4096,c128
still has a positive sampled bound. This is a limitation of that diagnostic,
not a new exclusion of the 256-state construction.

Next complete or decisively reject this wider, uniformly randomized bridge
before imposing structured feedback or adding checkpoint spacing. Keep the
already certified window construction as the complete-proof control.

## Checks and reproduction

```text
python -B research/workstreams/one_route_coupled_scan/coded_state/test_random_feedback_bridge.py
python -B research/workstreams/one_route_coupled_scan/coded_state/compact_reset_obstruction.py
python -B research/workstreams/one_route_coupled_scan/coded_state/random_feedback_bridge.py
python -B research/workstreams/one_route_coupled_scan/coded_state/uniform_feedback_screen.py
```

Exact tests cover48 fixed-input feedback laws, all4096 three-batch GF4
outputs with672 weighted-count comparisons, all1296 toy route restrictions,
and128 small outer/inner fixtures with explicit reset-supported witnesses.
The K20 obstruction uses exact polynomial coefficients and Fractions, not
simulation. Floating values are display only in that checker. The separate
screens are FLOAT diagnostics and cannot accept a certificate. Runs are
serial and bounded; no generated data or encoder benchmarks are produced.
