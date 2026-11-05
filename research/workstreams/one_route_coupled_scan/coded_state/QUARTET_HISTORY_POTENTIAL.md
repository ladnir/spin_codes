# Couple the potential to input history

2026-10-30. The full45%/40-bit goal and construction remain unchanged.
The [stationary potential](QUARTET_BOUNDARY_POTENTIAL.md) has a structural
loss even before numerical optimization. A count-conditioned replacement
passes exact local tests. A sharper coefficient bound also gives a useful
full-width diagnostic, but no whole-code or covered-count certificate.

## A stationary scalar factor forgets which state can precede a count

Let phi be ANY strictly positive function on the finite five-type histogram
space. Suppose T_l phi<=gamma_l phi for every histogram, with the same phi
for all input counts. For b>=1, each row sum of T_l is at least binom(w,l):
there are that many input masks and every incremental moment weight is at
least one. Evaluate the inequality at a histogram minimizing phi to obtain

    gamma_l >= binom(w,l).

Also, the all-type-four histogram is a self-transition under a fully active
input, with weight b^(lambda*w). Hence gamma_w>=b^(lambda*w).

In the scalar polynomial Gamma(x)^16, fix the last six updates of every
region to be fully active. Use the second bound on these six factors and
the first bound on the other ten. After exact quota normalization, the
same balanced payment-count family as the boundary audit contributes at
least2^60214.26162 after b^-B. The initial/terminal factor C is at least one,
because the terminal value at the all-type-zero state is one.

Thus the stationary scalar formula cannot close the middle class for ANY
positive potential, throughout the previous tau>=2lambda multiplier family.
This is stronger than failure of a particular product potential or of the
Maclaurin envelope. The maximizing incoming state can change independently
between successive scalar factors, recreating the favorable run histories
that the real transfer must charge.

## Retain the preceding count and consecutive full updates

After processing an input of weight a, n_0=w-a exactly. Keep a as a state
label. Also keep f, the number of consecutive fully active input updates,
capped at four. Use f=0 when a<w. There are w+4 labels, not one label per
histogram. A fully active update shifts every type one position toward
type four. Consequently:

- For a<w, permitted histograms have n_0=w-a.
- For a=w and1<=f<=4, n_0 through n_(f-1) are zero.
- At f=4, every coordinate is type four.

The next label is (l,0) if l<w and (w,min(4,f+1)) if l=w. These constraints
are preserved by the exact kernel. They are necessary history constraints;
they do not assert that every permitted histogram is jointly reachable
along every count sequence.

For a positive product potential phi_h, choose gamma_(a,f;l) to bound
T_l phi_h / phi_h only over the permitted incoming histograms. The scalar
calculation must now propagate the label, instead of multiplying unrelated
worst-case factors. Regional quota counters and the global paid count are
still required. The label is NOT reset at regional boundaries.

Remembering only a is insufficient to charge a long full-input burst's
exit. For example, use toy parameters w=4,lambda=1,tau=5/2,b=4 and
h=(1,1,4,16,32). On input weights(0,4,4,4,4,4,4,0):

| Calculation | Moment value |
|---|---:|
| Exact, and count/full-run history bound | 2^20 |
| Remember only the previous active count | 2^40 |
| Stationary scalar factors, before terminal factor | 2^48 |

The state after several full updates is old, not newly activated. The
full-run flag preserves that distinction. All1024 full/empty sequences of
length ten are exact under this toy potential, including free termination.
This does not establish exactness for partially active updates.

## Keep the types separate in the symmetric coefficient

The full-width screen shows another loss: replacing the type mixture by
its arithmetic mean can be expensive. Use a coefficient bound before
maximizing over the permitted types instead. Retain the previous notation
d=(h_1,h_2,h_3,h_4*b^alpha,h_4*b^lambda), alpha=3lambda-tau.
For every z>0,

    e_k(d repeated according to n)
       <= z^-k product_i (1+d_i*z)^n_i.

For a<w, fix n_0=w-a and let S={1,2,3,4}. For a=w, let S={f,...,4}.
The exact kernel therefore gives the uniform coefficientwise bound

    (T_l phi_h)(n)/phi_h(n) <= U_(a,f;l),
    U_(a,f;l) = sum_k J_(l,k) h_0^(w-l) h_1^(l-k)/binom(w,k)
       * inf_(z>0) z^-k ((1+d_0*z)/h_0)^(w-a)
           * (max_(i in S) (1+d_i*z)/h_i)^a.            (1)

Set gamma_(a,f;l)=U_(a,f;l); no optimal gamma needs to be computed.
Each positive z gives an upper bound. Taking
the infimum preserves it. Choosing a different z for each k is permitted
because the nonnegative k-terms are bounded separately.

The maximization over types is a maximum of at most four affine functions
of z. Its logarithm, and the complete objective in ln(z), are convex.
A positive evaluation point gives an upper formula even if optimization
has not converged. Floating-point rounding still needs certification.
This bound retains k and the input count; it does not identify independent
columns or discard the regional quotas.

For completeness, define a (w+4)-state matrix M(x) whose entry from label
(a,f) to its successor for input weight l is gamma_(a,f;l)*x^l. Apply the
regional normalization coefficientwise:

    R(y)=sum_(j=0)^r y^j binom(m-r+j,j)/binom(m,r)
                      * [x^(m-r+j)] M(x)^16.

The initial row vector has value h_0^w at label(0,0) and zero elsewhere.
The upper bound for Z_c is that vector, times[y^c]R(y)^256, times
a terminal vector majorizing the actual terminal factor divided by phi_h.
For h=(1,1,b^lambda,b^(2lambda),b^tau), that terminal ratio is exactly
b^(lambda*a). For a general h with h_0=1, the maximum over permitted
histograms gives a valid terminal vector. This formula defines the required
count/history calculation; no large matrix-polynomial allocation was made.

## Full-width diagnostic and what it does not prove

The screen fixes an input-count path with44 regions containing fifteen
408-count updates and one407-count update, and212 regions containing sixteen
408-count updates. It has exactly c=111572 payments under the regional
normalization. All incoming labels after the first update are partial,
so the full-run flag alone does not help this path.

Two potential families were tested: the full-run potential above, and the
positive right Perron vector of the single-coordinate transition at
active probability408/512 and full-quartet probability(408)_4/(512)_4.
This vector is only a candidate potential; independence between coordinates
is NOT used to certify the full-width kernel.

At log2(b)=4 with the Perron candidate, the floating-point fixed-path
log2 expression, including its quota factors, threshold, and row-set union,
has the following diagnostics:

| Local treatment | Fixed-path diagnostic |
|---|---:|
| Sampled lower bound on the previous Maclaurin envelope | +375024.36 |
| Sampled lower bound using exact kernel ratios | -33128.59 |
| New type-separated Cauchy upper formula | -9333.94 |

The first number already defeats that Maclaurin envelope at this candidate.
The second is a LOWER bound on a worst-case factor, not a certificate.
The third bounds only one fixed count path, in non-directed arithmetic.
Summing all count paths, all covered counts, and all row counts remains
open. In particular, -9333.94 is NOT a security margin. It shows that the
new coefficient bound survives a test the previous envelope fails.

## Reproduction and next test

`quartet_history_potential.py` checks15 stationary coefficient floors,
58 exact type-separated Cauchy coefficients (widths four and eight),65536
partial input streams grouped into625 count sequences, and1024 full/empty
sequences. All count-history bounds dominate direct enumeration. Runtime
0.18 seconds. `quartet_history_screen.py` computes integer quartet counts,
then evaluates the full-width diagnostic in floating point in2.78 seconds.
Both run serially and remotely with90-second watchdogs and no files.

Next sum a bounded neighborhood of the407/408 input-count path using the
type-separated envelope, keeping regional totals. This tests whether the
negative fixed-path value survives count-path entropy. If it does, develop
directed bounds and cover the count tails rather than treating that
neighborhood as the whole class. Preserve all previous exclusions and the
unchanged full45%/40-bit objective.
