# A full-width two-update bound and a diversity potential

The [tail diagnostic](RANDOM_PAIRING_FLOW_TAIL.md) ruled out one separable
one-step proof template. This continuation tests a two-update envelope at
the actual width512. The resulting bound is valid but its numerical screen
does not close the selected1974-row class. No encoder component changes.

## Two updates as a pairing sum

Use the IID target-label model for this envelope. Exact regional counts
will be restored by conditioning, with that cost included in the screen.
Let n=2m, let z>=1 and y>0 mark loss and target-label count, and choose
positive endpoint weights v(a), a in {-2,-1,0,1,2}.

For a pair with incoming auxiliary states a,b and j target labels, write
u(a,b,j),w(a,b,j) for its outgoing states and ell(a,b,j) for its loss.
These are the causal-flow quantities, not encoder field values. Define

    W_ab(j) = binom(2,j) p^j (1-p)^(2-j) z^ell y^j
              v(u) v(w) / (v(a) v(b)),
    S_ab = sum_j W_ab(j).

For a fixed intermediate state vector, the second update's moment ratio
is the uniform perfect-matching average of product S_ab. This is a
pairing polynomial; independence of partner states is not assumed.

Choose a symmetric matrix C that is positive definite and entrywise at
least S. In the diagnostic, C=S+dI for a small rational d. Exact Schur
complements verify positive definiteness. A centered real Gaussian vector
G with covariance C satisfies the pairing identity for its even moments.
Thus replacing S by C bounds the second matching average by

    E_G[product_i G_(intermediate state i)] / (2m-1)!!.

This replacement is justified term by term in the pairing polynomial.
The product inside the Gaussian expectation need not be nonnegative.

Now fix the first matching. Each of its m pairs contributes the quadratic
polynomial

    Q_ab(G) = sum_j W_ab(j) G_u G_w.

Define the probability vector

    q_ab(c) = sum_j W_ab(j) (1[u=c]+1[w=c]) / (2 S_ab).

The elementary inequality |xy|<=(x^2+y^2)/2 gives

    |Q_ab(G)| <= S_ab sum_c q_ab(c) G_c^2.

Apply this inequality to every first-layer pair and then apply Holder's
inequality with m equal exponents. The two-update ratio, for every input
histogram and every first matching, is at most

    B_2 = max_(a,b) S_ab^m
          E_G[(sum_c q_ab(c) G_c^2)^m] / (2m-1)!!.

Averaging the first matching preserves the bound. This is an all-histogram
envelope at arbitrary even n, not a toy-to-full-size extrapolation.
The numerical screen below is floating point, so it is not a directed
finite certificate even though the bounding formula is justified.

## Evaluate the Gaussian moment without enumerating histograms

Put A=diag(q) C and define H(t)=det(I-tA)^(-1/2). If h_k=[t^k]H(t), then

    h_0=1,
    k h_k = (1/2) sum_(i=1)^k tr(A^i) h_(k-i).

The normalized Gaussian moment in B_2 is

    h_m / ([t^m](1-t)^(-1/2)).

This requires only5-by-5 matrix powers and a scalar recurrence through
m=256. No large histogram transition matrix is built.

For a homogeneous incoming state a, first-layer pairs have equal
outgoing states. The same coefficient formula with C=S gives the exact
two-update moment from that homogeneous start. Here the identity is
formal algebra on pairing polynomials; it does not require S to be a
real Gaussian covariance. Do not use this formal identity to justify
the absolute-value or Holder bounds when S is not positive semidefinite.

## What the screen says

The helper `random_pairing_flow_block.py` uses z=27/25, y=7/6 and positive
endpoint weights rounded to rational values from the one-step spectral
control. The diagonal addition d=0.005 makes C positive definite. At
n512, the complete selected-class union expression is approximately
2^7689, versus2^6026 for repeating the one-step bound. The additional
Gaussian majorant loses more than the block gains at these markers.

This alone would not rule out a better majorant or endpoint weights.
Therefore the helper also checks a necessary condition that does not use
the Gaussian majorant: the exact two-update moments from homogeneous
starts. For fixed z,y, define the positive homogeneous operator

    F_a(v) = E[ z^(two-update loss) y^(two-update tags)
                product_i v(final state i) | all initial states a ]^(1/n).

The minimum and maximum of F_a(v)/v(a) bracket the best homogeneous-start
envelope. For the lower bound, compare any candidate vector with a scalar
multiple of v at a coordinate attaining the maximal ratio; monotonicity
and homogeneity give the inequality. Iteration finds a narrow floating
bracket, not a directed certificate.

A bounded nine-point marker screen, tuning endpoint weights at each point,
still has an optimistic exponent about+3664 bits at its best sampled point
(ln z,ln y)=(0.10,0.20). This omits nonnegative boundary costs and includes
the row-set union and regional conditioning. It is not an exhaustive
optimization over markers or a theorem excluding all two-update methods.
It does indicate that polishing this particular uniform block bound is
not the best next step.

The weakness is concrete: a scalar envelope must tolerate homogeneous
input histograms at every block boundary. It does not retain diversity
created in previous blocks. That is a relaxation in the proof, not an
actual reset or an operation in the code.

## Carry diversity across updates

A candidate nonseparable potential is

    V_kappa(h) = (product_a v(a)^h_a) exp(kappa sum_a h_a^2/n),
    kappa > 0.

The quadratic term assigns a larger value to concentrated histograms.
It can decrease when a homogeneous state produces several auxiliary
states, retaining that gain across update boundaries. No supermartingale
claim is made yet.

There is an exact initial check under the regional model. Start with all
n coordinates in state a, and condition on m_reg remaining positions
containing r target labels. Equal-state pairs produce equal outgoing
states b_a(j)=clamp(a+2j-1,-2,2). For two distinct first-layer pairs, put

    P_equal = sum_(j,k: b_a(j)=b_a(k)) binom(2,j) binom(2,k)
       (r)_(j+k) (m_reg-r)_(4-j-k) / (m_reg)_4.

Then

    E[sum_c h_next(c)^2] = 2n + n(n-2) P_equal.

The first term counts the two equal endpoints within each pair. The
second counts coincidences across distinct pairs. This is sampling
without replacement, not an IID substitute.

At n512, m_reg8192, r1974 and kappa=0.01, the expected change in the
quadratic log-potential is about-2.708 for initial states-2,-1,0 and
-2.491 for initial states1,2. These negative drifts only motivate testing
the potential. They do not bound its exponential moment; mixed histograms
may increase the same potential.

The identity

    exp(kappa ||h||^2/n)
      = E_Z exp(sqrt(2kappa/n) Z dot h),  Z~N(0,I_5),

expresses the candidate as a mixture of product potentials. A future bound
must retain the histogram-dependent normalization of that mixture. Taking
a separate all-histogram maximum for every component would discard the
intended benefit and need not escape the separable obstruction.

## Validation and next step

Exact rational checks cover all70 histograms at width4 for two updates,
15 independent fourth-moment identities, and5 regional collision-drift
identities. The selected full-width block screen took1.636s. The completed
nine-point necessary-condition screen took15.263s, with a22s iteration
deadline and a25s total bound. Its first invocation printed all nine
points but failed in a misplaced reporting call; that scope error was
corrected and the whole screen reran successfully. No failed-run output
is used as a substitute for the successful rerun.

The final exact-only suite took0.094s. Runs were serial, Python-B and
file-free, with the20GiB disk guard checked in each invocation. No native
encoder timing, construction change, paper edit, commit, or promotion.

Next test the conditional exponential moment of V_kappa on mixed
histograms, retaining actual uniform matching and regional counts. Reject
it if a concrete mixed profile defeats the needed rate. The current
results do not close the structural tail, the field-cancellation count,
or whole-code distance.
