# A uniform bound for nearly full outer support

The ordinary rank transfer already controls every message whose routed outer
word has at most one percent zero coordinates. No lookahead augmentation,
outer enumeration, or numerical matrix search is needed for this band.
This is a partial distance result, not a certificate for all messages.

## Statement and probability space

Use the coupled scan from [CONSTRUCTION.md](CONSTRUCTION.md), with width four
or eight, K=2^20, N=2K, initial state zero, and unrestricted terminal state.
Fix any injective F-linear outer map of dimension K and any coordinate
permutation before sampling the inner coefficients. In particular, this
includes every realization of the proposed outer and structured route.

At each step and coordinate, independently sample a coefficient from
theta+GF(256). The resulting setup is fixed for every message. Let V(m) be
the routed outer word and Y(m) the complete encoder output. Then

    Pr_inner[there exists m!=0 such that
                 wt(V(m)) >= N-floor(N/100)
             and wt(Y(m)) <= floor(7N/20)]
        < 2^(-167780).

The bound is uniform over the fixed outer and route. It therefore also holds
after averaging their actual setup distributions. The message may depend on
the complete setup; this is not a fixed-message statement.

The proof uses the common-quotient rank accounting in
[RANK_CONTINUATION.md](RANK_CONTINUATION.md), but not its lookahead improvement.

## A constant positive vector for full input support

Fix 1<=tau<=256 and 0<z<=1/tau. Write

    a = 1/tau,
    f = z+(1+z)/tau = z+a+az.

For one step, let I, H, J, and H' denote the input, current-state, output,
and next-state support masks. The ordinary rank charge is

    p = |H minus (J intersection H')|.

Each admitted branch receives weight z^|J| tau^(-p). Preliminary
deterministic rank gains are omitted, as in the existing relaxed transfer.
Sum over all output and next-state masks to obtain a row sum.

For coordinate j, write h=1[j in H] and k=1[j in R(H)]. The following table
gives its row-sum factor. The next-state bit is summed out. The column
"input zero" applies when j is outside I.

| (h,k) | Input nonzero | Input zero |
|---|---|---|
| (0,0) | z | 1 |
| (0,1) | 2z | 1 |
| (1,0) | z+a | z |
| (1,1) | f | z(1+a) |

These factors follow directly from the two necessary support inclusions

    I symmetric_difference H subset J subset I union H,
    J symmetric_difference R(H) subset H' subset J union R(H).

The coordinate product includes the forbidden empty-input extinction branch
if that branch satisfies the mask conditions. Keeping it only enlarges the
row sum, so the product remains a valid upper bound for the actual transfer.

First suppose I is the full mask. Since R permutes coordinates, the numbers
of coordinate types (0,1) and (1,0) are equal. Pair their factors. Because
z<=a,

    (2z)(z+a) <= (z+a)^2 <= f^2,

and the (0,0) factor satisfies z<=f. Every full-input row sum is therefore
at most f^w. Equivalently, the all-ones positive vector satisfies

    M_full(z) 1 <= f^w 1.

This argument does not depend on the ordering of the state cycle.

## Charge each missing input coordinate by at most z^-1

For each row of the table, its input-zero factor is at most z^-1 times its
input-nonzero factor. For the first two rows this is immediate. For the last
two rows the input-zero factor is already no larger than the input-nonzero
factor, and z<=1.

Consequently, for an arbitrary input mask I,

    M_I(z) 1 <= f^w z^(-(w-|I|)) 1.

Now fix an exact nonempty routed outer support S, with h=N-|S| missing
coordinates. Shorten the message space to S and denote its dimension by D_S.
If D_S=0, it contributes no nonzero message. Otherwise D_S<=K. Require every
coordinate in S to be nonzero, retaining these requirements in the same
common quotient throughout the scan.

Multiplying the row bounds over N/w steps gives

    e_empty^T (product_t M_(I_t)(z)) 1 <= f^N z^(-h).

There is no terminal-state factor: the final all-ones vector already sums
every terminal support. The rank-event bound and the output-weight marker
therefore give

    Pr_inner[there exists a nonzero message with exact support S
             and output weight at most e]
        <= tau^(D_S-1) z^(-e-h) f^N
        <= tau^(K-1) z^(-e-h) f^N.                         (1)

The input forms need not be independent or uniformly distributed. The
ordinary rank lemma is uniform over those forms and the incoming quotient.

## Sum every support in the band

Put H_max=floor(N/100). There are at most binomial(N,h) supports missing h
coordinates, including supports that the fixed outer cannot realize. Their
inclusion is a conservative union bound. For any 0<x<=1,

    sum_(h=0)^H_max binomial(N,h) z^(-h)
        <= x^(-H_max) (1+x/z)^N.                           (2)

Indeed, x^h>=x^H_max for the included terms. Multiplying each term by
x^h/x^H_max and then adding the remaining binomial terms proves (2).

Combining (1) and (2), the probability in the statement is at most

    U = tau^(K-1) z^(-e) f^N x^(-H_max) (1+x/z)^N,

where e=floor(7N/20). Choose the exact rational values

    tau = 256,
    z = 7/(13*257) = 7/3341,
    f = 5/832,
    x = z/99.

They satisfy 0<x<=1 and 0<z<=1/tau. Using K=N/2,
e<=7N/20, and H_max<=N/100 gives

    U <= tau^(-1) B^(N/100),

    B = tau^50 z^(-35) f^100 (100/99)^100 (99/z).

The required finite comparison is the exact rational inequality

    B < 1/256.

It can be checked without logarithms, floating-point arithmetic, or interval
matrix calculations. After clearing denominators, it is precisely

    256^51 * 3341^36 * 500^100
        < 7^36 * 832^100 * 99^99.

It follows that

    U < 2^(-8-8N/100)
      = 2^(-167780.16)
      < 2^(-167780)

at N=2^21. Nonintegral exponents cause no rounding assumption: equivalently,
raise the positive bound for U to the 100th power before comparing.

## Scope and next obstruction

This proof covers all exact outer supports with at most floor(N/100) zeros,
not just the single full-support mask. It is conditional on any fixed outer
and route, and it makes no timing claim.

It does not control messages with more missing outer coordinates. In
particular, it does not contradict the identity-outer example supported in
only 33 regions. That example lies far outside this band.

The next discriminating task is to control the remaining sparse and
intermediate supports using their actual shortened dimensions and outer
distribution. A large augmented-state calculation is unnecessary for this
already-closed dense band.
