# Cell-aware column generation for factorized witnesses

Independent inner/outer mixtures reduce the discrete bank from all cross-pairs
to two marginal banks. Continuous column generation can extend those banks at
the cells where the current finite minimax problem is weakest.

## Restricted cell problem

Fix one leaf with exact vertices `v_1,...,v_m`. For each inner component `i`
and outer component `j`, write

```text
I_i(a)=c_i-u_i*a,                 O_j(a)=d_j-w_j*a.
```

Let `H(a)` be the packet-profile orbit normalization. For restricted eligible
banks `A` and `B`, solve

```text
minimize    U
subject to  sum_i alpha_i I_i(v_k) + sum_j beta_j O_j(v_k) - H(v_k) <= U
            for k=1,...,m,
            alpha_i>=0,  sum_i alpha_i=1,
            beta_j>=0,   sum_j beta_j=1.
```

Let `y_k>=0` be the dual multiplier for vertex `k`, normalized so that
`sum_k y_k=1`. Define

```text
mu_I=min_(i in A) sum_k y_k I_i(v_k),
mu_O=min_(j in B) sum_k y_k O_j(v_k).
```

Strong LP duality gives

```text
U=mu_I+mu_O-sum_k y_k H(v_k).
```

The two minima are the dual values associated with the two marginal-sum
constraints. They need not equal averages under the primal weights.

## One pricing point

Define the dual barycenter

```text
abar=sum_k y_k v_k.
```

Every component is affine. Therefore

```text
sum_k y_k I_i(v_k)=I_i(abar),
sum_k y_k O_j(v_k)=O_j(abar).
```

The inner pricing oracle minimizes `I_theta(abar)` over valid inner parameters
`theta`. The outer oracle independently minimizes `O_phi(abar)` over valid
outer parameters `phi`. Their reduced costs are

```text
rho_I=I_theta(abar)-mu_I,          rho_O=O_phi(abar)-mu_O.
```

A candidate violates the restricted dual when its reduced cost is negative.
Discovery should require `rho<-tau_price`, where `tau_price` exceeds solver
and oracle noise. Insert inner and outer candidates independently. A useful
candidate on one side does not require a simultaneous candidate on the other.

This pricing rule is cell-aware. The dual assigns mass only to vertices that
support the current cell optimum. Uniform vertex averages and the cell's
geometric centroid generally solve a different problem.

## The common normalization

The normalization contributes the fixed dual term

```text
sum_k y_k H(v_k).
```

In general, this term is not `H(abar)`. The Gamma extension of `H` is
nonlinear. Neither pricing oracle may replace the dual term by `H(abar)` when
computing a reduced cost.

An existing inner tuner may report

```text
P_theta(abar)=I_theta(abar)-H_Gamma(abar).
```

Minimizing `P_theta(abar)` still selects the same `theta`, because
`H_Gamma(abar)` is independent of `theta`. Pricing must then recover
`I_theta(abar)` by adding that same Gamma value. The outer oracle already
returns its unnormalized affine part. The combined vertex evaluator subtracts
`H(v_k)` once, after adding the two marginal values.

## Fractional pricing profiles

The barycenter is nonnegative, has total mass `M`, and lies in the leaf. It
need not be an integer profile. This causes no proof problem. Tuning only
chooses deterministic parameters; every accepted parameter row is later
evaluated at the leaf's proof vertices.

The analytic objectives extend naturally to real profiles:

- the conditioned-row outer objective depends on the profile only through an
  affine dot product;
- the inner constant depends on its pole, fugacities, transfer kernel, and
  Collatz vector, while its profile charge is affine;
- `log Gamma(a_k+1)` supplies a finite discovery normalization for `a_k>=0`.

The current code has two interface hazards.

1. `refine_inner_sparse_fugacities` and `refine_inner_full_fugacities` cast the
   supplied profile to `int64`. They must not receive `abar`; truncation changes
   its mass, support, and anchor. A pricing-specific tuner should retain a
   binary64 profile and optimize the unnormalized affine objective.
2. `evaluate_conditioned_outer` does not require integral entries, but it tests
   the binary64 sum for exact equality with `M`. A pricing wrapper should build
   one coordinate as `M` minus the others, or bypass this check through a
   dedicated affine objective.

`tune_inner_density` and the shaped-density initializer retain floating profile
entries. They are suitable bounded starts. They do not turn their local search
into a global pricing oracle.

Within an exact-support shard, every vertex has positive coordinates on the
shard support. Thus `abar` has the same positive support. The tuner may set
variables to zero only outside that support. Final eligibility remains a
component-level verifier check.

## Discovery iteration

For each selected leaf, use the following bounded loop.

1. Solve the restricted factorized LP and validate its primal and dual
   residuals.
2. Form `abar` from the returned vertex dual.
3. Run the inner oracle from active inner rows and fixed generic starts.
4. Run the outer oracle from active outer rows and fixed generic starts.
5. Freeze each returned binary64 parameter row and compute its complete column
   on all leaf vertices.
6. Recompute its price from that complete column. Reject a candidate unless
   the measured reduced cost is below `-tau_price`.
7. Deduplicate parameter rows by canonical parameter digest. Add accepted rows
   to the global component bank.
8. Re-solve the leaf LP. Stop at the round cap or after two consecutive rounds
   with no accepted column.

The producer should also reprice every newly added global component on other
selected leaves before calling another continuous oracle. One cell's pricing
column may close several neighboring cells.

A local nonlinear optimizer cannot prove that no negative reduced-cost column
exists. The no-column stop is only a discovery stop. It must not be recorded
as pricing optimality or used as a completeness premise.

## Certificate boundary

The dual vector, barycenter, reduced costs, optimizer traces, and column
generation stop are diagnostics. A certificate contains only:

- content-addressed inner and outer parameter rows;
- two fixed sparse rational marginals per certified leaf;
- the exact leaf geometry and count; and
- the declared single normalization rule.

The independent verifier reconstructs every component interval from its
parameters. It evaluates both rational marginals at every exact vertex and
subtracts one outward normalization interval. It need not reproduce the
continuous tuning or establish global pricing optimality. If rationalization
or outward evaluation loses the diagnostic improvement, the producer retains
the prior selector or leaves the cell unresolved.

## Bounded first experiment

Use the eight active leaves with the largest current union contributions.
Deduplicate leaves by their original `h2` parent before filling the eighth
slot. Run at most three pricing rounds per leaf. Each round emits at most one
inner and one outer candidate, for at most 48 candidate rows before global
deduplication. Reprice every accepted row across all eight leaves.

Use `tau_price=64` bits for this first diagnostic. Keep the old selector as an
explicit fallback, so no replayed leaf may regress. Do not interpret failure
to cross this tolerance as evidence that the continuous family is exhausted.

Proceed to a larger wave only if all correctness conditions and the payoff
condition hold.

Correctness conditions:

1. Every affine barycenter identity agrees with the full vertex-column price
   within `1e-6` bit before the `64`-bit acceptance margin.
2. Every accepted candidate is support-eligible and has a stable canonical
   parameter digest.
3. Rational marginal replay and outward replay never exceed the retained old
   bound on any of the 181 leaves.
4. At least one new inner and one new outer row receive positive rational
   weight. Otherwise classify the experiment as one-sided, even if useful.

The payoff condition is:

```text
at least four of the eight leaves improve outward by 16384 bits,
and either the 181-leaf expanded diagnostic or its worst contribution
decreases by at least 32768 bits.
```

If correctness passes but payoff fails, retain any useful columns and stop
continuous pricing. The next run should compare geometry refinement against a
broader component family. If both pass, price one additional finite wave
before changing the proof geometry.
