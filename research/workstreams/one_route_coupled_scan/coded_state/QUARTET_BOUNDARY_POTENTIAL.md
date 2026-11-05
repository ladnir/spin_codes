# Keep the boundary state in the regional moment

2026-10-30. The construction and full45%/40-bit goal are unchanged.
The [regional-quota transfer](QUARTET_REGIONAL_QUOTAS.md) is valid, but its
free-terminal block relaxation cannot close the middle class. This note
proves that obstruction and derives a smaller upper-bound interface that
retains boundary state. The latter is not yet a numerical certificate.

Continuation: [the history-conditioned audit](QUARTET_HISTORY_POTENTIAL.md)
rules out the stationary scalar formula below for every positive potential.
It retains input-count history instead and replaces the arithmetic-mean
step with a type-separated coefficient bound. Exact local tests and a
bounded full-width diagnostic are recorded there; the full goal is open.

## Free terminal states are not an affordable regional relaxation

Keep c=111572, r=2100 inactive positions per8192-position region, and256
regions. Assign435 payments to44 regions and436 payments to212 regions.
For each region require every unpaid position to occur before its last six
512-coordinate updates. Conditional on the regional payment count c_j,
the unpaid set is uniform of size r-c_j. Its weighted contribution is

    g(c_j)=binom(r,c_j) binom(5120,r-c_j)/binom(8192,r-c_j).

These requirements respect the actual regional quotas. If every region is
given a free terminal state, its six-update suffix permits full state
support with no deaths. The total relaxed score is at least lambda*786432.

This exceeds B throughout the current onset-three multiplier family. To
see this, put a=39/512, t0=381/512, L0=713983+c, M=255463-c, and
K0=a*L0-(c-111127). The current field weights satisfy

    B=lambda*L0-K0-nu*(L0-743225)-eta*M,
    tau=t0+eta,  lambda=a+mu+nu,
    mu,nu,eta>=0.

For tau>=2lambda, eta>=max(0,2lambda-t0). Since K0>0 and t0*M-K0>=0,
considering lambda on either side of t0/2 gives

    B/lambda <= L0-2K0/t0 =250598845/381 <786432.

Thus, for every b>=1, the EXACT free-terminal regional product expression
after multiplying by b^-B is at least

    binom(256,212) g(435)^44 g(436)^212.

Its directed log2 lower bound is

    66206280811877073/1099511627776 >60214.26162.

The binomial factor counts disjoint assignments of payment counts, not
overlapping suffix choices. This obstruction uses neither a read-degree
bound nor a scalar paid marker. Improving those steps cannot repair the
one-region free-terminal product within this multiplier family.

The real stream has only ONE free terminal state. At an internal boundary,
a selected last state coordinate also requires its next quartet's inputs
to be active. To obtain a support guaranteed feasible under the suffix
event, delete each internal region's final selected update and charge the
preceding death. At lambda=707/2560,tau=381/512 this changes the displayed
suffix support's score from217190.4 to83978.4; B=165556 lies between them.
These are scores of specified supports, not upper bounds on the true
geometric maximum and not codeword witnesses.

Two-region free-terminal blocks also fail at these selected weights.
Fix44 regions with435 payments and six-update suffixes,104 with436
payments and six-update suffixes, and108 with436 payments and seven-update
suffixes. Fix their locations; do not count nested suffix events as disjoint.
Within each two-region block delete the last selected update of its first
region and charge its death. The guaranteed block score is165594.4>B.
The exact weighted event count has directed log2 lower bound

    914033532750807/34359738368 >26601.87696.

This does not exclude larger blocks or a boundary-preserving calculation.
It does exclude proceeding with the proposed one-region product, or merely
doubling its block length at the selected field weights.

## A coefficientwise bound for a product potential

The exact transfer has five type counts n=(n_0,...,n_4), summing to w=512.
Enumerating every histogram is impractical. Instead fix positive numbers
h_0,...,h_4 and use the potential phi_h(n)=product_i h_i^n_i. This section
bounds one transfer step while retaining its exact active-input count l.

Assume2lambda<=tau<=3lambda and put alpha=3lambda-tau. Set

    d=(h_1,h_2,h_3,h_4*b^alpha,h_4*b^lambda),
    s(n)=sum_i n_i*d_i/w.

Let T_l include the raw allowed-mask counts, the average over the uniform
quartet partition, and the moment reward, restricted to l active inputs.
If k positions belong to wholly active quartets, the exact kernel gives

    (T_l phi_h)(n)
      =sum_k J_(l,k) h_0^(w-l) h_1^(l-k)
         * e_k(d repeated according to n)/binom(w,k).   (1)

Here J_(l,k) is the previously defined quartet occupation count, and e_k
is the elementary symmetric polynomial. Averaging two arguments of e_k
with their sum fixed increases their product and cannot decrease e_k.
Repeated averaging therefore proves

    e_k(d repeated according to n)/binom(w,k) <= s(n)^k.

Substitution in (1) gives the COEFFICIENTWISE envelope

    (T_l phi_h)(n) <= Q_l(s(n)),
    Q_l(s)=[x^l] ((h_0+x*h_1)^4-(x*h_1)^4+(x*s)^4)^(w/4). (2)

The quartet polynomial has nonnegative coefficients, including s^4 as
its degree-four coefficient. No count is replaced by independent inputs.
Equation (2) applies at width512 without enumerating its histograms.

## Reducing verification to one dimension

For a second positive vector g, let psi_g(s) be the lower convex envelope
of the five points (d_i,ln(g_i)), on[min_i d_i,max_i d_i]. Repeated d_i
values use the lowest ordinate. For every histogram,

    sum_i (n_i/w)*ln(g_i) >= psi_g(s(n)).

Consequently any upper bound gamma_l satisfying

    gamma_l >= sup_s Q_l(s)*exp(-w*psi_g(s))             (3)

proves T_l phi_h <= gamma_l phi_g for ALL histograms. The relaxation to
real convex combinations in psi only enlarges the supremum. Verifying
(3) is a one-dimensional problem on at most four linear hull segments,
rather than a table of billions of histograms. Directed evaluation is
still required; sampling s or using a local optimizer is not a certificate.

A stationary choice g=h gives a particularly simple full-stream interface.
Let Gamma(x)=sum_(l=0)^w gamma_l*x^l. The exact final factor has type weights

    e=(1,b^lambda,b^(2lambda),b^(3lambda),b^(lambda+tau)).

It is bounded by phi_h times(max_i e_i/h_i)^w. The initial type count is
(w,0,0,0,0). Thus define C=(h_0*max_i e_i/h_i)^w. It is charged once, not
once per region. Applying the coefficientwise inequality at every update
and the exact quota normalization at every regional boundary gives

    Z_c(b) <= C * [y^c] P(y)^256,
    P(y)=sum_(j=0)^r y^j binom(m-r+j,j)/binom(m,r)
          * [x^(m-r+j)] Gamma(x)^16,  m=8192.           (4)

The state potential telescopes across regional boundaries. Only the
regional input counter is reset; no new free terminal state is introduced.
The paid coefficient remains explicit. Formula (4) is an upper-bound
interface, not a claim that a useful h or small gamma has been found.
Time-dependent vectors can give sharper bounds, but have not been tested.

As a normalization check, b=1 and h=(1,1,1,1,1) give
gamma_l=binom(w,l), C=1, P(y)=(1+y)^r, and Z_c<=binom(256r,c), with equality.

## Checks and next step

The initial floating-point regional read-four screen had log2 expressions
328301.86,351276.71,482491.59,678317.54,1167113.27,1722853.89 for
log2(b)=1,2,3,4,6,8. Those were non-directed diagnostics, not certificates
or universal lower bounds. The exact suffix argument above independently
proves failure of the underlying free-terminal regional product.

`quartet_quota_audit.py` reproduces both new directed obstructions in its
24.22-second full audit. `quartet_regional_transfer.py` checks108 exact
potential coefficients against the direct kernel and an independent quartet
polynomial expansion. All satisfy (1)-(2). Its prior12544 quota cases,
65536 boundary cases, width-eight partition checks, and count normalization
still pass. Total transfer audit0.34 seconds. Runs were serial, bounded,
remote, and file-free; no full-width histogram table was launched.

Next test a bounded set of stationary potentials in (3)-(4) at c=111572.
Keep both the regional quota and the boundary potential. A promising
floating-point candidate needs directed supremum and coefficient bounds;
a failed stationary envelope is not evidence that the actual code fails.
All other covered counts and row counts remain required for the full goal.
