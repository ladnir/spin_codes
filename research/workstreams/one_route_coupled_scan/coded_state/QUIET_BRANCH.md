# Cheap column scaling closes the entirely quiet branch

The post-diagonal inner still permits nonzero inputs in ker(B) while its
state is zero. This note changes only the inner compression coefficients.
It closes every low-weight trajectory with no live incoming state, including
possible activation on the last batch. Mixed live/quiet trajectories remain
open; this is not a whole-code distance certificate.

## Construction and scope

Keep F=GF(2^128), E=GF256, K=2^20, N=2K, and the fixed RS[256,128] outer.
Keep the 256 independent uniform regional permutations. There are M=8192
outer blocks and n=16384 batches, each containing 128 consecutive symbols.
Each batch lies inside one region. Use A and B from [DESIGN.md](DESIGN.md),
with a 64-symbol state.

Independently sample every diagonal entry of D_t from F and every diagonal
entry of R_t from E minus zero. Sample all coefficients once during setup,
independently of the route, and fix them for every message. Define

    Y_t = X_t + D_t A Q_(t-1),
    Q_t = Q_(t-1) + B R_t Y_t,       Q_0=0.

The final state is unrestricted. This is a separate variant of
[POST_DIAGONAL.md](POST_DIAGONAL.md), not a reinterpretation of its fixed B.
It adds one GF256 scaling per output symbol, with 2 MiB of raw coefficients.
It adds no full-F products. Prepared storage and native runtime are unmeasured.

For every fixed nonzero R_t, any 64 columns of B R_t are independent and
its kernel has minimum distance 65. A live incoming state therefore still
gives an exactly uniform next state over F^64. The earlier sparse and
regional proofs use precisely these properties, not B A=I. Consequently
their combined <2^-46 bound at >49.1% transfers for all messages using at
most 64 outer blocks, even after fixing every R_t.

## Descend the quiet subcode to E

Here an entirely quiet trajectory means Q_0 through Q_(n-1) are zero;
Q_n may be nonzero. Its output is Y=X, and its constraints are

    B R_t X_t = 0,        1 <= t < n.

Fix the route and all R_t. These constraints and the outer are defined
over E. Their solution code over F is the scalar extension of their
solution code over E. Decompose an F-valued codeword in an E basis of F.
Every component is an E-valued solution, and each component's support is
contained in the original support. At least one component is nonzero.
Conversely, every E solution is an F solution. Thus these two solution
codes have the same minimum symbol distance, pointwise for each setup.

This descent applies only to the quiet subcode. Live trajectories involve
the full-F entries of D_t and do not inherit this argument.

## Charge the nonzero quiet batches

Fix a nonzero E-valued batch x before sampling its R. If wt(x)<=64, then
B R x cannot vanish. Otherwise select 64 nonzero coordinates of x and
expose every R entry except those coordinates. Their B columns are
independent, so at most one choice of those 64 nonzero scalars solves
B R x=0. Therefore

    Pr_R[B R x=0] <= 255^-64 <= 256^-kappa,
    kappa = 1279/20 = 63.95.

The last inequality follows from the exact integer comparison
255^160 >= 2^1279. Independent batches give a factor 256^(-kappa*j)
for j prescribed nonzero quiet batches of a fixed E-valued outer word.

Let J contain its nonzero batches among 1,...,n-1, and put j=|J|.
Always include the last batch in the allowed coordinate set U, whether
or not its input is nonzero. Then |U|=128(j+1). A nonzero outer word uses
at least 129 regions, so j>=128. Every nonzero quiet batch has weight at
least 65. Weight at most W=floor(.491N)=1029701 therefore implies

    128 <= j <= floor(W/65) = 15841.

For a fixed route, write t_i for the coordinates of outer block i in U.
The dimension of the outer shortened to U is

    d(U) = sum_i max(0,t_i-128).

There are at most 256^d(U) E words supported in U. Restrict to words
nonzero on every batch in J before applying the j local probability
factors; this restriction only reduces that count. The probability that
one such nonzero quiet word exists is at most

    1[d(U)>0] min(1, 256^(d(U)-kappa*j)).

This union is over all messages, not a fixed-message distance claim.

## Average the actual regional route

For 1<tau<=256, the preceding expression is at most

    tau^(-kappa*j) 1[d(U)>0] tau^d(U).

Set k=j+1, v=min(k,n-k), and let Z have binomial distribution
Bin(256,v/n). Define H=E[tau^max(0,Z-128)]. The retained regional-moment
argument bounds the moment by H^M on the smaller of U and its complement.
It uses negative association between labels within each region and convex
averaging of the 256 independent regional inclusion probabilities. Regions
with probability zero or one are allowed; this is not an independent-label
approximation. The exact MDS complement identity is

    d(U) = |U|-N/2+d(U^c).

For k<=n/2, keep the positive-dimension correction

    E[1[d(U)>0] tau^d(U)]
      <= min(H^M, (H^M-1)/(1-1/tau)).

For k>n/2, use tau^(128*k-K) H^M instead. There are binom(n-1,j)
choices of J. These expressions give a complete one-dimensional union
bound, with no assumptions about the arrangement of its chosen batches.

## Finite verification

[quiet_certificate.py](quiet_certificate.py) covers all 15714 values of j
using four boxes. Floating point proposes boxes and rational markers only.
The final bound uses exact rational binomial moments and outward-rounded
arithmetic at 80 and 120 decimal digits. Both replays pass; the worst
box upper bound on log2 of a summand is -11327.1404373998... .
We claim only the preselected threshold of 2^-80 per summand. Since the
number of summands is less than 2^14, the entire quiet branch fails at
>49.1% distance with probability at most 2^-66 over route and R setup.
The bound is independent of D setup.

[test_quiet_branch.py](test_quiet_branch.py) independently checks 15 local
quiet-probability fixtures and 432 uniform-refresh fixtures over GF4.
It explicitly embeds GF4 into GF16 and exhaustively compares a nontrivial
quiet subcode over both fields: minimum weight 6 in both, with 15 and 255
nonzero words respectively. The last batch is unconstrained. It also
checks the integer coefficient comparison and coverage count.

Run serially from the repository root:

```text
python -B research/workstreams/one_route_coupled_scan/coded_state/quiet_certificate.py
python -B research/workstreams/one_route_coupled_scan/coded_state/test_quiet_branch.py
```

## Next discriminating step

The two closed classes are small source support and no live incoming state.
Neither covers large-source messages with both live episodes and quiet
nonzero batches. The latter batches satisfy small-field kernel constraints;
the live episodes introduce full-field coefficients. A valid combined proof
must retain the common message space and their shared constraints. It cannot
multiply separately averaged bounds or descend the entire encoder to E.

First bound the dimension left by a prescribed collection of quiet-batch
constraints, then combine that dimension with the live-episode calculation.
The new [mixed rank interface](MIXED_RANK.md) does this with a common
independent-coordinate set and a random-diagonal rank bound. Its dense-class
diagnosis shows that this forward count alone cannot complete the proof,
even with ideal quiet rank; output-side counting is the next step.
If this loses too much, test a structural change to compression rather than
repeatedly tuning the same loose counting bound. The fixed-B candidate,
two-tap candidate, and certified window control remain preserved.
