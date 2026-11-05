# Keep the actual feedback support in the routing bound

The containing-interval bound failed on a full-length live episode, even
when every other output batch was zero. Keeping the feedback holes repairs
that loss, provided their output constraints and mask counts are retained.
The construction in [ALTERNATING_MAPS.md](ALTERNATING_MAPS.md) is unchanged.

The result below closes one full-length trajectory class at **>49.1%
distance with failure probability <2^-61**. It is not a whole-code
certificate. Arbitrary intervening outputs and multiple episodes remain open.

## 1. Invert a fixed output before averaging the route

Use F=GF(2^128), q=|F|, K=2^20, N=2K, M=8192 independent RS[256,128]
outer blocks, and n=16384 inner batches. Fix all nonzero compression
scalings R. The 256 regional permutations and all full-field diagonal
entries D remain independent random setup choices.

For a fixed output y, reconstruct Q from

    Q_0=0,   Q_t=Q_(t-1)+B R_t y_t.

Define S(y,R) as the coordinate support of A_(epsilon_t) Q_(t-1)
across the stream. In particular, S depends on y and R, but not on the
route or D. Over D, the inverse input is uniform on y+F^S. For a fixed
route pi, let C_pi be its routed outer and d_pi(S) its shortened dimension.
The affine intersection is empty or has size q^d_pi(S). Therefore

    Pr_D[inverse(y) in C_pi] <= q^(d_pi(S)-|S|).          (1)

Each outer block has one coordinate in each region. If t_i counts its
coordinates in S, MDS shortening gives d_pi(S)=sum_i(t_i-128)_+.
Because the total occupancy is exactly |S|,

    d_pi(S)-|S| = -K+d_pi(S^c).                          (2)

Apply this identity before relaxing the regional sampling. Keeping the
complement removes the spurious contribution from all feedback holes
inside an otherwise live interval.

## 2. Regional moment for the actual mask

Fix S and let a_j=|S^c intersect region j| for j=1,...,256.
For tau>=1, define

    G(a,tau) = E[tau^((Z-128)_+)],
    Z=sum_j Bernoulli(a_j/M),

where these 256 Bernoulli variables are independent. Then

    E_pi[tau^d_pi(S^c)] <= G(a,tau)^M.                  (3)

This does not assume independent occupancies of different outer blocks.
In one region, their inclusion indicators form a uniform subset of size a_j.
Conditioning on other regions, expand the product of increasing block
factors as a polynomial with nonnegative coefficients in those indicators.
Every degree-r inclusion probability is (a_j)_r/(M)_r <= (a_j/M)^r.
Replacing one region by independent indicators thus increases the moment.
Repeat across regions; the resulting outer-block factors are independent.
This proves (3), using the retained
[regional comparison](../STRATIFIED_FRACTIONAL_BOUND.md).

Equation (3) preserves every regional count. A convenient uniform bound
uses only the total h_c=|S^c|. Convexity of tau^((z-128)_+) permits averaging
the 256 Bernoulli parameters at fixed sum:

    G(a,tau) <= G_bin(h_c,tau),
    G_bin(h_c,tau)=E[tau^((Bin(256,h_c/N)-128)_+)].        (4)

The equality in (2) is exact even though (3) and (4) are upper bounds.
No assumption about the number of partially occupied regions is needed.

For each exact mask S, suppose at most C_S output candidates belong to a
specified class. Union over its outputs first, conditional on the route.
For 1<tau<=q and alpha=log_q(tau), (1)--(3) imply

    Pr[there is a bad output with this mask S | R]
        <= C_S^alpha tau^(-K) G(a,tau)^M.               (5)

Then sum (5) over masks. Mask multiplicities are outside alpha. Applying
a single fractional power to varying masks while using one mask's route
moment would be invalid. The following class counts all its masks explicitly.

## 3. The full-length class

Consider outputs satisfying both conditions:

- Q_1,...,Q_(n-1) are nonzero; Q_n is unrestricted.
- y_2=y_4=...=y_(n-2)=0.

Thus the episode activates in batch 1 and stays live through the last
incoming state. Pair feedback observations (2,3),(4,5),...,(n-2,n-1).
Within each pair, the zero intervening output preserves the state. The
two maps evaluate one degree-<64 polynomial at 256 distinct points.
Each pair has at most 63 feedback holes. The unpaired final batch has
at most 63 holes among its 128 coordinates.

There are P=8191 pairs and one unpaired observation. Write H for their
total holes. Then

    0<=H<=516096,   |S^c|=128+H.                        (6)

The extra 128 coordinates are the activation batch, whose incoming state
is zero. No union over interval locations is needed: this class occupies
the entire stream and fixes the zero-output parity.

Define P_l(x)=sum_(h=0)^63 binom(l,h)x^h. The number of hole masks with
total H is at most

    L_H=[x^H] P_256(x)^8191 P_128(x).                   (7)

This counts exact feedback masks; some counted masks are unrealizable.
Two useful upper bounds are binom(N-128,H), and, for every rational x>0,

    L_H <= x^(-H) P_256(x)^8191 P_128(x).               (8)

At H=516096 the coefficient equals binom(256,63)^8191 binom(128,63).

## 4. Charge the outputs that create the holes

There are F_b=8193 potentially nonzero output batches: the 8192 odd
batches and the final batch. Let U=128 F_b=1048704. The final output has
no state constraint, since the terminal state is free.

Each paired state is created by the immediately preceding odd output;
the final unpaired state is created by output n-1. These are 8192 distinct
fresh output blocks. Fix their supports and all proposed hole positions.
For an observation with h holes and a preceding output of weight w,
the matrix on that output's fresh variables is A_H B_W R_W.
Its rank is at least

    r(h,w)=max(0,h+min(64,w)-64).                       (9)

Indeed A_H has row rank h, B_W has column rank min(64,w), and R_W is
invertible. For a pair, A_H uses distinct rows from the combined 256-point
evaluation map. Earlier constraints do not contain the fresh variables.
Selecting triangular minors therefore adds (9) across observations.
The incoming state is never charged as independent randomness.

For every 0<=h<=63 and 0<=w<=128, the remaining local dimension satisfies

    w-r(h,w) <=128-h-(64/127)(128-w).                   (10)

For w>=64, the left side is w-h. For w<64, it is min(w,64-h);
the tight slope among these cases occurs at h=63,w=1. The unconstrained
last output also satisfies (10) with h=0 as an upper bound on its dimension.

At a fixed total output weight w_total<=W=floor(.491N)=1029701,
the total kernel dimension is at most

    min(w_total, U-H-ceil(64(U-W)/127)).                 (11)

Let v_H=min(W,U-H-ceil(64(U-W)/127)). Summing every output support size
gives the following bound, uniform over the exact hole mask and fixed R:

    C_H <= (U+1) binom(U,max(v_H,U/2)) q^v_H.            (12)

To see this, bound the summand of weight w by binom(U,w)q^min(w,v_H).
For w<=v_H the q-weighted binomial terms increase because q>U.
For w>=v_H use the maximum binomial coefficient from v_H through W,
bounded by binom(U,max(v_H,U/2)). U is even. There are at most U+1 terms.
This count may include outputs without the required liveness or exact holes;
relaxing those conditions only enlarges C_H.

## 5. Complete restricted-class bound

Combine (4)--(8) and (12), fixing one tau independently for each H:

    failure <= sum_(H=0)^516096
        L_H C_H^(log_q tau) tau^(-K)
        G_bin(128+H,tau)^8192.                         (13)

Every output amplitude, hole position, output support, and regional route
in the class is included. The argument holds for each fixed nonzero R,
so averaging R changes nothing. It covers witnesses chosen after setup.

The directed verifier checks all 516097 integer H values through 37
nonoverlapping intervals, at both 80 and 120 decimal digits. Each H term
is at most 2^-80. As 516097<2^19, the full class has failure below 2^-61.
The weakest interval upper bound is approximately -232.10 in log2 units;
the certificate deliberately retains the conservative -80 threshold.

For each interval [a,b], the verifier uses the dimension at a, the largest
binomial-support bound across that interval, and the route moment at b.
Positive polynomial bounds in (8) use a or b according to the tilt's sign.
The binomial moment is formed with exact integer numerators. All subsequent
arithmetic uses outward Decimal endpoints; logarithms are enclosed by the
neighbors of their correctly rounded values. Floating arithmetic chooses
candidate tau values only. It cannot accept a bound.

## 6. Checks and next step

Run serially from the repository root:

```text
python -B research/workstreams/one_route_coupled_scan/coded_state/test_feedback_support.py
python -B research/workstreams/one_route_coupled_scan/coded_state/feedback_support_certificate.py --precision 80
python -B research/workstreams/one_route_coupled_scan/coded_state/feedback_support_certificate.py --precision 120
```

Tests enumerate 4096 incidence matrices and 768 exact moment comparisons.
They also check 512 multitime fresh-column rank fixtures, every scalar
tradeoff through state64, the small integer output-volume domain, directed
rational operations, and rejection on immediate timeout. No generated
files, benchmarks, or construction changes are involved.

The next extension should allow nonzero intervening outputs. Then the two
observations can have different state polynomials; their holes need not
form one root set of size at most 63. Retain the new actual-support route
moment and fresh-output accounting, but replace that pair enumerator with
the joint live-to-live constraints. Sparse one-half updates have a retained
[local interface](HALF_UPDATES.md); dense updates and multiple episodes
still need a complete count. Do not extrapolate (13) to those trajectories.
