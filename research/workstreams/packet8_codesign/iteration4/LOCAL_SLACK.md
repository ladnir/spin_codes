# Local filling cannot explain the remaining actual24 q16 gap

The actual24 fixed-eight-step proposal gives about 26.7524 bits at q=16,
weight tilt .06, and fractional power .4. Reaching 40 bits requires another
13.2476 bits. The local relaxation that fills forbidden output destinations
cannot provide that improvement while retaining the same fractional coordinates.

The analytic ceiling is less than 1.5 bits. A cached fixed-four-step diagnostic
improves by only 0.41521 bits even with an optimistic lower comparison.
This conclusion concerns local filling, not the outer envelope or the loss
from expanding hidden paths at fractional block boundaries.

## Physical law and the existing positive comparison

Fix an entering state a!=0, a physical potential support S, and z=exp(-theta).
There are j=|S| potential bytes. Under the uniform outer comparison their
labels are independent uniform bytes, including zero. The emitted block is
Y=X+Aa. The next state is Ma+CX, where M is fresh and sends every fixed
nonzero state uniformly over the L=2^24-1 nonzero states.

Define the weighted emission and syndrome masses

```
T = E[z^wt(X+Aa)],
N_u = E[z^wt(X+Aa) * 1{CX=u}],
b0 = E[z^wt(X+Aa) * 1{X=0}].
```

The exact weighted mass at a target state u is (T-N_u)/L. At a nonzero
target, the current comparison uses T/L. Its zero-target coefficient is
(T-b0)/L: potential-label thinning preserves the old rule that zero input
cannot return a nonzero source to state zero.

These statements concern nonzero source states. For a zero source, the
weighted zero law and the normalized potential-birth family are retained
exactly. No correction below changes those rows. It also leaves the entire
j=0 operator unchanged.

## Atom bound after potential-label thinning

Weighting by z^wt(X+Aa) preserves independence between potential bytes.
For an offset v in one byte, the normalizer is

```
sum_{x in GF256} z^wt(x+v) = (1+z)^8.
```

Consequently each tilted byte atom is at most delta=(1+z)^(-8). This is a
bound for full potential labels. It does not incorrectly condition every
potential byte to be nonzero.

For the actual24 map, every restriction of C to one, two, or three byte
positions has rank 8, 16, or 24. Fixing all other bytes therefore determines
at most one assignment to any r=min(j,3) selected bytes for a specified
syndrome. Hence

```
N_u / T <= delta^min(j,3),                    j>=1.
```

The inequality holds for each fixed potential support and entering state.
Averaging supports and any normalized nonzero source family preserves it.
Thus the exact mass at every nonzero target is at least

```
(1-delta^min(j,3)) * T/L.
```

## False zero returns are much smaller

At j<=3, injectivity of the complete support restriction gives CX=0 only
when X=0. Therefore N_0=b0 and the existing return coefficient is exact.
This does not say that returns are impossible: nonzero input can cancel
the independently refreshed nonzero state. Removing these returns would
delete genuine transitions.

For j>=4, put n=N_0/T and p0=b0/T. The fraction incorrectly added to the
zero-return coefficient is

```
(N_0-b0)/(T-b0) = (n-p0)/(1-p0) <= n <= delta^3.
```

The first inequality follows from 0<=p0<=n<=1. This argument accounts for
the all-zero input already excluded by the comparison. The exact return
mass is therefore at least (1-delta^3) times the filled return coefficient.

An alternative active-byte proof gives eta^3, where
eta=1/((1+z)^8-1). The potential-byte argument above is both simpler and
slightly sharper; the diagnostic uses delta^3.

## A lower comparison and its scope

Use the same potential-birth coordinates as the existing proposal. Form a
nonnegative shadow family S_j from the current family T_j as follows:

- Preserve the zero-source row and the entire j=0 operator.
- For j>=1, multiply each nonzero-source U coefficient by
  1-delta^min(j,3).
- For j>=4, multiply each nonzero-source zero coefficient by 1-delta^3.
  Preserve those coefficients at j<=3.

Let E map coordinate coefficients to their represented state measures, and
let K_j denote the exact weighted physical transition averaged at potential
occupancy j. The local inequalities give, in row-operator convention,

```
S_j E <= E K_j <= T_j E.
```

All measures and transitions are nonnegative, so these inequalities compose
chronologically. They also survive averaging exact placements. The shadow
is a lower comparison, not an upper bound that could certify distance.

At occupancy q, all regions together contain exactly Rq potential slots,
where R=32. There can be at most Rq nonempty physical steps. Since S_j is
entrywise at least (1-delta)T_j for j>0, local filling changes the unpowered
conditional comparison moment by at most

```
-Rq log2(1-delta) bits.
```

For any fixed deterministic partition into fractional blocks, multiply the
local lower factors inside each block and then raise that product to alpha.
Entrywise monotonicity gives the corresponding same-coordinate ceiling

```
-alpha Rq log2(1-delta) bits.
```

At theta=.06, q=16 and alpha=.4, these estimates are respectively
3.663942466 and 1.465576986 bits. The latter applies to the saved g8 operator
without repeating its tuple census.

There is also a conservative exact check. Since exp(-3/50)>=47/50,
delta<=(50/97)^8. Integer arithmetic verifies

```
(1-(50/97)^8)^2048 > 2^(-15).
```

This proves the fractional ceiling is strictly below 1.5 bits, without
depending on rounded evaluations of exp or log. At most Rq/4=128 steps can
have j>=4. False-return subtraction alone can improve by at most

```
-alpha*128*log2(1-delta^3) = 0.000008948021 bits.
```

An exact rational estimate using -ln(1-x)<=x/(1-x) and 1/ln2<3/2 proves
this quantity is below 0.00001 bits.

These ceilings isolate the physical filling relaxation. A different basis
or longer fractional blocks may reduce path-expansion loss while preserving
the unpowered moment. Such improvements are not excluded by the sandwich.
Neither the sandwich nor its shadow gives a lower bound on actual code failure:
the beta-uniform outer comparison remains a majorant.

## Cached bounded diagnostic

`local_slack.py` authenticates `iteration3/cache24_v2/theta_3_50.json`, thins
potential labels, and applies the frozen potential-birth rebase. It computes
only the fixed-four-step q16 proposals. It runs no large-state census and
does not rebuild a g8 operator.

| Operator used for the diagnostic | Margin | Gain over unchanged |
|---|---:|---:|
| Unchanged comparison | 17.3541827921 bits | 0 |
| False-return shadow only | 17.3541828094 bits | 0.0000000173 bits |
| Full local shadow | 17.7693926083 bits | 0.4152098162 bits |
| Delete every nonzero-to-zero return | 253.5189699128 bits | 236.1647871207 bits |

The last three rows are not certificate operators. The last row is especially
misleading as a construction prediction: it removes real cancellation events
at j=1,2,3, where the current return coefficient is already exact.

`local_slack_q16_v1.json` records the authenticated inputs, source hashes,
analytic checks, and all-log regional evaluations. Four tiny tests independently
check the sandwich in a GF4 three-coordinate fixture, preservation of zero-input
rules, the rational ceilings, and scalar/adjoint transitivity. All four pass.

## Transferring GL3 updates to field scalars

Identify the same 24-bit state space with GF(2^24) through a fixed binary-linear
basis. Sample c uniformly from its nonzero elements, independently at each
update and independently of the routing and outer setup. Multiplication by c
fixes zero and sends each fixed a!=0 uniformly over all nonzero states:
c maps to ca bijectively.

Uniform GL(3,GF256) has exactly the same action law on one fixed state.
Condition on any fixed message, outer setup, route, and preceding execution.
The fresh update sees a determined entering state. Both update families have
the same conditional next-state law and hence the same full emitted-word law,
by induction over the physical steps.

This is equality of the law for each fixed message, not equality of random
linear-map ensembles or joint laws for several messages. It is sufficient
here. Conditional on routing geometry J, the first moment sums fixed-message
probabilities. Those summands transfer individually, so the conditional
expectation bound and min(1, conditional expectation) clipping still transfer.
No independence between different messages is used.

The binary-adjoint scalar family is also nonzero-transitive. Let L_c denote
the binary matrix of multiplication by c. For fixed a!=0, the map
c maps to L_c^T a is linear and injective: a nonzero difference c-d makes
L_(c-d)^T invertible. Including c=0 gives a bijection of all 2^24 states;
excluding it gives the required nonzero law. This handles the transpose kernel
without silently assuming the multiplication matrix is symmetric.

The polynomial X^3+X+1 is irreducible over GF256. Its roots lie in GF8 and
not in GF2, while GF8 intersects GF256 in GF2. A cubic without a root is
irreducible. Thus the proposed AES-byte cubic representation is suitable.
Literal implementation verification must still check its coordinate basis,
adjoint orientation, nonzero-uniform sampling, and freshness at each update.

## Recommendation

Do not pursue an expensive exact excluded-atom census to recover the 13.25-bit
q16 deficit. The available local slack is too small. Prioritize outer/routing
geometry, map changes, or an explicitly tighter fractional representation.
The field-scalar update substitution is proof-compatible for this first-moment
method, but its implementation and timing remain separate checks.
