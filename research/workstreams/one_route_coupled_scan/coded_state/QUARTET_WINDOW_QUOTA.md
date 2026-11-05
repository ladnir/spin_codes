# Keep the unpaid-input count in the window moment

2026-10-30. The construction and full45%/40-bit objective remain unchanged.
The [window-cover formulation](QUARTET_WINDOW_COVER.md) eliminated support
histories. This continuation identifies a second obstruction: replacing
the remaining input quota by independent bits admits impossible dense
streams. Sharpening temporal dependence alone cannot repair that step.

The new results are an obstruction to a class of moment calculations, an
exact comparison retaining the unpaid count, and a checked finite transfer.
They are not a whole-code distance certificate.

Continuation: [the regional-quota audit](QUARTET_REGIONAL_QUOTAS.md) now
rules out even exact global fixed-count coefficient extraction for this
comparison. A long all-active suffix has excessive moment mass under the
uniform-subset law but is impossible under the actual regional quotas.
The exact five-type kernel remains valid; its continuation must retain
the regional counts, not just their sum.

## Why an exact Bernoulli moment would still fail

Retain N=2097152, R=537600 original inactive input positions, and the middle
covered count c=111572. Put u=R-c and p=1-R/N=1523/2048. The window value
H_(lambda,tau) is the maximum of lambda*A-tau*D over allowed supports.
For an all-active input stream it equals lambda*N, including the free terminal.
It is nonnegative for every stream.

The earlier field multipliers mu,nu,eta>=0 give

    lambda=39/512+mu+nu, tau=381/512+eta,
    B=(1+mu+eta)c+713983mu+743225nu-255463eta-111127.

Every remaining low-output message requires H>B. In this entire multiplier
family,

    B/lambda <= Lmax=713983+c=825555.                    (1)

Indeed, c-111127<= (39/512)Lmax,743225<=Lmax, and
eta(c-255463)<=0. These are exact rational comparisons, not an optimization
over a finite parameter grid.

For a paid-set marker t>0, the independent-bit comparison used active
probability p_t=1-(1-p)/(1+t). Let M(p_t,b) be the EXACT whole-stream moment
E[b^H] under those independent bits and the actual quartet partitions.
The proposed first-moment expression, even before the row-set union, is

    U=(1+t)^R t^-c b^-B M(p_t,b), b>=1.                  (2)

Since H>=0 and the all-active event has probability p_t^N,

    M(p_t,b)>=max(1,b^(lambda*N) p_t^N).

When B>0, interpolate these two lower bounds with weight B/(lambda*N).
Equation (1) ensures this weight lies in[0,1]. We obtain

    b^-B M(p_t,b)>=p_t^(B/lambda)>=p^Lmax.

If B<=0, the same lower bound follows from b^-B M>=1. Also
(1+t)^R t^-c>=binom(R,c), by positivity of the binomial expansion. Thus
EVERY choice in this proof-parameter family satisfies

    U>=binom(R,c) p^825555,
    log2 U>=11904170035564781/274877906944
             >43307.1183.                              (3)

This remains true with exact spatial and temporal dependence inside M.
It rules out closing (2) by a sharper read-degree bound or different scalar
markers. It does not rule out other field constraints, exact quota counts,
capacity certificates, or the code itself.

There is a simpler problem at the previously selected point t=1/6,
lambda=707/2560,B=165556,log2 b=3. The actual unconditioned paid-set
polynomial includes P=I, whose size is R rather than c. This one term yields

    t^(R-c) b^(lambda*N-B) >2^(665912/5)=2^133182.4.

The inequality uses log2(6)<13/5, checked by6^5<2^13. This is not a term of
the fixed-c coefficient. It shows why coefficient extraction cannot be
silently replaced by a single positive evaluation at that marker.

The negative independent-column diagnostics in the preceding note therefore
do not describe a nearly completed proof. Their apparent slack disappears
once actual shared-input dependence permits the all-active event in the
relaxed probability space. That event is impossible when u>0 is retained.

## Compare the actual unpaid set with a uniform set of the same size

For a fixed original row set, each of Breg=256 regions independently has an
r=2100-subset I_j of its m=8192 positions. Conditional on |P|=c, selecting
P uniformly from I is equivalent to selecting U=I outside P uniformly with
size u=R-c. For a fixed u-subset U, write u_j=|U intersect region j|. Then

    Pr_nested[U]
      = product_j binom(m-u_j,r-u_j)
         / [binom(m,r)^Breg binom(R,u)],                 (4)

with value zero when some u_j>r. The number of supersets I of U proves (4).
The distribution is NOT uniform among global u-subsets.

Set f(v)=binom(m-v,r-v). Its successive ratio is

    f(v+1)/f(v)=(r-v)/(m-v),

which decreases in v. Hence the product in (4), for fixed sum u_j=u, is
maximized when the u_j differ by at most one. Write u=Breg*h+e,0<=e<Breg.
For the uniform global u-subset law,

    Pr_nested[U] <= D(u) Pr_uniform[U],
    D(u)=binom(N,u) f(h)^(Breg-e) f(h+1)^e
           / [binom(m,r)^Breg binom(R,u)].               (5)

This pointwise comparison applies to every nonnegative function of U and
the independent quartet partitions. It retains the exact unpaid count.
Directed log2 bounds for D are:

| c | u | Upper bound on log2 D(u) |
|---:|---:|---:|
| 27476 | 510124 | 272185290389085/549755813888 <495.103 |
| 111572 | 426028 | 135987214195495/549755813888 <247.360 |
| 255463 | 282137 | 30306365478071/274877906944 <110.254 |

Thus the exact fixed-c first moment can be upper-bounded by

    binom(R,c) D(u) E_uniform_u[b^H].                    (6)

The row-set union and threshold factor still have to be included. Equation
(6) is a proof interface, not a claimed small probability. The comparison
is verified on2673 exact small unpaid subsets, including normalization and
attainment of the balanced maximum.

## One binomial conditioning bound also loses too much

It is tempting to evaluate E_uniform_u[b^H] by choosing independent active
bits of probability s and dividing the unconditioned moment M(s,b) by
Pr[Binomial(N,s)=N-u]. That is an upper bound, but a lower bound on this
EXPRESSION proves that it cannot close the current middle class either.

Drop D(u)>=1 and the row-set union to make the expression smaller. The
same interpolation gives b^-B M(s,b)>=s^Lmax. Therefore the expression is
at least

    binom(R,c) / binom(N,u)
      * s^(-(N-Lmax-u)) (1-s)^(-u)
      >= binom(R,c) binom(N-Lmax,u) / binom(N,u).         (7)

The second inequality is the binomial point-mass bound. It holds for every
s in(0,1), every b>=1, and every current nonnegative field multiplier. Its
directed log2 lower bound is

    42687892632741515/1099511627776 >38824.4121.

Thus keeping the quota only as a denominator after an unrestricted moment
is insufficient. The count must constrain the moment calculation itself,
or the proof must avoid this first-moment expression. Exact fixed-count
coefficients are not ruled out by (7).

## An exact transfer retaining temporal overlap and the count

For the middle profile,2lambda<=tau<=3lambda. Put alpha=3lambda-tau.
Let A_t,j denote the active-input bit and V_t,j the indicator that its entire
quartet at update t is active. The allowed-state bit is

    z_t,j=A_t,j V_(t+1),j,
    z_final,j=A_final,j.

After processing free bits through t-1, retain A_t,j and the preceding
free-run length capped at three. Only five types are reachable:

    0: A_t,j=0, run length0;
    1: A_t,j=1, run length0;
    2: A_t,j=1, run length1;
    3: A_t,j=1, run length2;
    4: A_t,j=1, run length at least3.

If V_t,j=0 then the preceding free bit is zero, so an inactive current
input cannot have a positive preceding free-run length. This removes the
three unreachable states from the earlier eight-state column description.
Record the five type counts n=(n_0,...,n_4), totaling the width w.

At the next update let l be the active-input count and k the number of
positions in wholly active quartets. The number of masks for a fixed
partition with these statistics is

    J_(l,k)=[z^l y^k] ((1+z)^4-z^4+y^4 z^4)^(w/4).     (8)

Averaging the independent uniform quartet partition makes those k positions
a uniform k-subset, conditional on l,k. If v_i of them have old type i,
their intersection probability is product_i binom(n_i,v_i)/binom(w,k).
The new type counts are exactly

    n'=(w-l, l-k+v_0, v_1, v_2, v_3+v_4),              (9)

and this transition has moment weight b^(alpha*v_3+lambda*v_4). Its raw
coefficient is J_(l,k) times that hypergeometric probability. Keep l as
the increment of an explicit active-input counter. Summing all transitions
without moment weights gives2^w, not a fresh fixed-quota probability law.

Initialize with counts(w-l,l,0,0,0), coefficient binom(w,l), and counter l.
After all remaining updates, the exact free-terminal factor is

    b^[lambda*n_1+2lambda*n_2+3lambda*n_3+(lambda+tau)*n_4]. (10)

Extract the coefficient with total active count N-u and divide by
binom(N,u). This yields E_uniform_u[b^H]. The extraction is INSIDE the
calculation, unlike the failed scalar conditioning bound in (7). No
independence is asserted between successive window starts or columns.

The helper checks the kernel against26880 explicit partition/input cases.
It checks all65536 four-coordinate/four-update input streams and17 fixed-
count coefficients. Unweighted transfers return binom(N,a) for every a.
An eight-coordinate/five-update run retains multiple overlapping windows
and reaches5101 states. These exact toy checks use lambda=1,tau=5/2,b=4;
they validate the transfer, not the large-code parameters or distance.

## Scope and next step

`quartet_quota_audit.py` reproduces (3),(5),(7) with directed logs and exact
integers; its run takes about20.33 seconds. `quartet_window_transfer.py`
checks (8)-(10) in under one second. All runs are serial and remote, have
90-second watchdogs, require20GiB free, and write no files. No construction,
library, paper, frozen helper, or whole-code failure-budget changes were made.

The five-type transfer at width512 has more than2.9 billion possible type
histograms before adding the active-count coordinate. It has NOT been run,
and brute-force scaling is not proposed. The next task is a bounded,
quota-aware envelope or block reduction of this exact transfer, first at
c=111572. It must retain the count strongly enough to avoid (3) and (7).
The other covered counts, row counts, and all original goal requirements
remain open. This is a proof-method correction, not evidence for redesign.
