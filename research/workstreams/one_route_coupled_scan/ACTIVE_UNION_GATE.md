# A coarse complete gate for active unions

Status: an analytic complete-count envelope with one exactly bounded
active-union band, not a whole-code distance certificate. This bound is
for the actual one-diagonal inner and packet route. It does not reuse the
two-diagonal reset probability. Root owns numerical evaluation.

The guarded regional moment retains detailed support counts. A simpler
first gate can instead replace every random inverse support by its common
whole-batch activity union. This loses information, but permits all output
weights and support masks within one excursion layout to be summed before
fractionalization. The resulting complete sum has only two indices.

## Fix one excursion layout

Use the construction in [One-diagonal feedback](SINGLE_DIAGONAL.md), with
field size q, batch/state width b, N=bn coordinates, and a free terminal.
The outer has M=N/B MDS [B,k] rows, k=B/2, grouped into G=M/p packets.
Require p to divide b and b to divide M; put L=M/b. Thus each region has
L complete batches and each activity union contains whole packets.

An excursion activates at a batch with zero entering state and nonzero
output. It continues until the first return to zero, or to the stream's
end. Closed excursions have at least two batches and at least b+1 nonzero
output symbols, by the global MDS reset equations. Consecutive excursions
may abut. All but possibly the last excursion are closed.

Fix an ordered layout with c excursions occupying v batches. Let U be its
whole-coordinate activity union. Delete each activation batch to obtain
the set of live entering-state batches. Write

    u=bv,             a=b(v-c),
    w0=(b+1)(c-1)+1,  e'=min(e,u),
    Q=sqrt(q),
    V(u,e,c)=sum_(w=w0)^e' binom(u,w) Q^w.                    (1)

An empty sum is zero. Here e is the bad-output-weight threshold. The layout
and its sets are fixed before the route and inner randomness are averaged;
every output is eventually covered by the union over layouts.

For an exact state-support mask A and output-support mask W of total weight
w, let E be the dimension bound for its linearized output assignments.
The amortized rank lemma in [Delayed syndrome rank](SINGLE_DIAGONAL_RANK.md)
gives

    E-|A| <= (w-a+b)/2.                                      (2)

There are at most 2^a choices of A: its activation masks are empty and it
is supported on the a live coordinates. Counting all such choices includes
impossible masks. There are at most binom(u,w) choices of W. For every
route, D(A)<=D(U), where D denotes outer shortened dimension.

The inverse membership law and projective counting therefore imply,
conditionally on the route,

    Pr_inner[an eligible nonzero output in this layout]
      <= 1[D(U)>0] min(1, C(u,a,e,c) q^D(U)),

    C(u,a,e,c)=2^a q^((b-a)/2) V(u,e,c)/(q-1).                (3)

This sums all masks and all eligible weights before the next step.
Different masks now share the same D(U) factor. That replacement, not an
independence assumption, permits the entire factor C to enter a fractional
power. The guard follows because the nonzero inverse input belongs to the
outer shortened to U. No condition D(A)>0 is imposed.

## Average the guarded whole-packet moment

Choose 1<tau<=q and put nu=log_q(tau), t=tau^p, and rho=v/n. The actual
whole-packet route has

    D(U)=p sum_(g=1)^G (T_g(U)-k)_+.

In particular, its positive values are at least p. Pointwise,

    tau^D(U) 1[D(U)>0]
      <= t/(t-1) (tau^D(U)-1).                              (4)

The regional without-replacement comparison and convex Bernoulli
equalization give

    E_route[tau^D(U)] <= H(rho,t)^G,
    H(rho,t)=E[t^((Z-k)_+)],       Z distributed as Bin(B,rho).

The function z -> t^((z-k)_+) is convex for t>=1. This is the unguarded
moment, so no convexity assertion about the marked J from
[Guarded regional moments](GUARDED_REGIONAL_MOMENT.md) is needed.
Subtracting the exact constant one in (4) yields

    F(v,c;tau)=C(u,a,e,c)^nu
                   t/(t-1) [H(v/n,t)^G-1].                  (5)

This bounds the layout's failure probability; it may be capped at one.
The layout count remains outside the fractional power. There is no
renormalization by the probability of D(U)>0.

For stable positive evaluation, write

    H=1+x,
    x=sum_(i=k+1)^B binom(B,i)rho^i(1-rho)^(B-i)(t^(i-k)-1),
    H^G-1=x sum_(r=0)^(G-1)(1+x)^r.                          (6)

Both expressions have nonnegative terms. A floating subtraction of two
near-one approximations is not a directed certificate.

## Count every layout

The feasible integer indices satisfy

    1<=c,       2c-1<=v<=n,       w0<=e.                     (7)

For fixed v,c the number of layouts is at most

    J(v,c)=binom(v-c,c-1) binom(n-v+c,c).                     (8)

Assign minimum lengths two to the first c-1 intervals and one to the
last. Distributing the remaining v-(2c-1) batches among c intervals gives
the first binomial coefficient. Distributing n-v quiet batches among the
c+1 prefix, interior, and suffix gaps gives the second. Zero interior
gaps allow abutting intervals. Allowing a free last interval to end before
the stream ends only adds impossible layouts. Every actual output remains
covered, with no artificial terminal reset charge.

As a check, summing (8) over v gives binom(n+1,2c), by composing all
remaining interval and gap lengths together.

Some indices exclude every route before applying the averaged bound.
An interval of length l batches touches at most

    1+ceil((l-1)/L)=2+floor((l-2)/L)

regions. The formula also holds for l=1. Summing over intervals, even if
their touched regions overlap, proves

    touched(U) <= min(B, 2c+floor((v-2c)/L)).                (9)

The easier min(B,2c+floor(v/L)) is also safe. If either upper bound is at
most k, then D(U)=0 for every route, so the whole term is zero. Apply this
test before the homogeneous binomial moment; that moment can otherwise
assign positive mass to a deterministically excluded class.

Thus a complete failure bound is

    sum_(v,c satisfying (7), not excluded by (9))
      J(v,c) min(1, inf_(1<tau<=q) F(v,c;tau)).              (10)

A finite list of verified markers may replace the infimum. Restricting
the v range gives a bound for that entire active-union size range, not a
sample of profiles. It is valid to seek a sparse-U closure first while
leaving the remaining indices explicitly open.

## Volume bounds and the next evaluation

For upper bounds one may omit the lower cutoff w0 in (1). The resulting
truncated weighted-binomial sum is increasing in u and e. For any 0<z<=1,

    V(u,e,c) <= z^(-e') (1+Qz)^u.                          (11)

If 0<e'<u and e'<=Qu/(1+Q), choose z=e'/[Q(u-e')]. This gives

    log V <= e' log Q+u h(e'/u),
    h(x)=-x log x-(1-x)log(1-x).

Otherwise use log V<=u log(1+Q). In particular, Q=2^16 at q=2^32 is
smaller than many K20 values of u. Do not import a volume formula whose
proof assumed its alphabet weight was at least u. The displayed
condition is the one needed here. The simpler full-space bound remains
valid in every case.

At a fixed marker, H(v/n,t) and its guarded excess increase with v.
For q>=4, the factor 2^a q^(-a/2) decreases with a=b(v-c).
These facts permit conservative rectangular bounds, but do not make the
entire summand monotone in v or c. Bound the two binomial factors in (8)
with their actual ranges rather than assuming the complete term is
monotone. All sums must retain (7), the guard, and the free terminal.

Equation (10) is deliberately coarse. Equalizing all B regional densities
can lose many orders of magnitude, and replacing D(A) by D(U) can lose
the dense-profile gains. Failure of this upper bound is not a low-weight
codeword. Preserve the more informative guarded regional operator and
the deterministic regional cap as complementary tools. No whole-code
closure follows from the partial checkpoint below.

## Exact partial checkpoint

At K=2^20, N=2^21, q=2^32, B=64, p=4, and b=128, root's serial exact
checker proves

    Pr[an eligible nonzero output with weight <=734003,
       8096<=v<=8192, and 1<=c<=32] < 2^(-255).              (12)

This covers all 3104 integer (v,c) pairs in that rectangle, every marked
excursion placement, all exact masks, and every eligible amplitude. The
probability is over the actual packet route and the independently sampled
inner coefficients, fixed once for all messages. All other (v,c) classes
remain outside this checkpoint.

The calculation uses only one rational marker, tau=33/32, with t=tau^4.
Here n=16384, G=8192, and rho=v/n<=1/2. A crude volume bound is enough:

    V(u,e,c) <= 2^u Q^e.

Dropping the projective saving 1/(q-1) from C and using a=u-bc gives

    log_2 C <= 16e-14u+15bc+16b
             <= 16*734003-14*(128*8096)+15*128*32+2048
              = -2700496
             <= -32*84390.

Therefore C^nu<=tau^(-84390), without approximating nu=log_q(tau).
For every pair in the rectangle, the layout count obeys

    J(v,c) <= binom(8191,31) binom(8320,32).

The first top argument is at most 8192-1, and the second is at most
16384-8096+32. Both lower indices are below their respective binomial
modes, so replacing them by 31 and 32 increases the bound. Put

    H0=2^(-64) sum_(j=0)^64 binom(64,j)t^((j-32)_+).

Equation (10), monotonically bounding H at rho=1/2, is at most

    3104 binom(8191,31) binom(8320,32)
      tau^(-84390) t/(t-1) (H0^8192-1).                    (13)

Every factor in (13) is rational. The helper
`theory/check_active_union_band.py` forms its integer numerator and
denominator, retains the exact subtraction in H0^8192-1, and verifies the
strict comparison with 2^-255. Root's serial run passed. An independent
read-only source audit checked the amplitude rounding, layout monotonicity,
guard factor, and final integer comparison; it ran no additional jobs.

This is a complete band bound, not a sample of support patterns. It is
still a partial result: a separate bounded floating diagnostic of 101
selected pairs found mostly uninformative positive bounds. Those probes
neither certify their surrounding classes nor exhibit bad codewords.
The homogeneous moment has not closed the full two-dimensional sum.

The next sparse gate should retain regional activity densities when
summing excursion layouts, using the guarded regional weights before
global Bernoulli equalization. Preserve (12) as a checked band; do not
extrapolate its margin to smaller unions or more excursions.
