# Dense occupation intervals without a full region coefficient table

This diagnostic method bounds a whole interval of occupations using a
positive coefficient inequality. It retains the activation-aware epoch
transfer and the original uniform permutations. Its purpose is to make
dense evaluation possible when the number L of outer rows is too large
for a table of all region coefficients. A coarse counting majorant or a
wide interval can leave substantial slack; this method is not an outward
certificate.

Fix a deterministic outer spectrum, or simultaneous shell caps for one
reused constituent. Let B be the constituent length, A_w the nonzero shell
multiplicities or caps, and Q the number of nonzero outer rows. After its
coordinate permutation, one row has counting measure

    nu(x) = A_w / choose(B,w),  w = weight(x) > 0,
    nu(0) = 0.

For 0<p<1 let mu_p be the product Bernoulli(p) law and define

    Gamma(p) = max_w A_w / (choose(B,w) p^w (1-p)^(B-w)).

Then nu <= Gamma(p) mu_p pointwise. Applying this inequality to Q row
positions contributes Gamma(p)^Q. It does not resample the outer code.

Let T_j(z) be the existing three-state epoch envelope for a uniform
weight-j input of length t. There are E=L/t epochs in each region. Under
the auxiliary row law, a region contains Q eligible positions chosen
uniformly from its L positions, each carrying an independent Bernoulli(p)
bit. Its transfer is

    K_Q = [u^Q] S(u,p)^E / choose(L,Q),

where

    S(u,p) = sum_a choose(t,a) u^a
               sum_(j<=a) choose(a,j) p^j (1-p)^(a-j) T_j(z).

This follows by counting eligible positions per epoch. Conditional on the
number j of actual ones, their support is uniform, including the kernel
cases already represented by T_j. It is also exactly
sum_j choose(Q,j) p^j (1-p)^(Q-j) R_j(z), with the previously defined
fixed-weight region matrices R_j.

For any positive xi, all coefficients are nonnegative, so entrywise

    K_Q <= S(xi,p)^E / (choose(L,Q) xi^Q).

The binomial identity gives a small matrix representation:

    S(xi,p) = (1+xi)^t M(theta),
    theta = p xi/(1+xi),
    M(theta) = sum_j choose(t,j) theta^j (1-theta)^(t-j) T_j(z).

Since state continues across B regions, choose(L,Q) original row choices
and the counting-measure domination yield

    U_Q = z^(-H) e_Z M(theta)^(BL/t) 1
          (1+xi)^(BL) Gamma(p)^Q xi^(-BQ) choose(L,Q)^(1-B).

Write eta=log(xi), lambda=-log(z), and
J=log(e_Z M(theta)^(BL/t) 1). Then

    log U_Q = lambda H + J + BL log(1+exp(eta))
              + Q (log Gamma(p) - B eta)
              - (B-1) log choose(L,Q).

For a fixed witness (lambda,p,eta), this expression is convex in integer Q:
the binomial sequence is log-concave, B>=1, and the remaining Q-dependent
term is affine. Therefore every integer lo<=Q<=hi satisfies

    log U_Q <= max(log U_lo, log U_hi),

and the interval union is bounded by

    (hi-lo+1) max(U_lo,U_hi).

Different disjoint intervals may choose different witnesses. Their bounds
must be summed, and their union must cover every intended occupation. An
endpoint check is sufficient here because of the explicit convexity
argument; checking endpoints of an unrelated bound would not suffice.

`dense_occupation_ranges.py` implements this calculation. Its tests compare
the eligible-position bound with every relevant exact region-mixture term
on a toy map, check the serialized matrix power, and check every integer
in toy intervals against the endpoint envelope. The numerical implementation
uses nearest binary64, including log-binomial evaluation. Outward replay
would need directed arithmetic for all these steps.

For random outers, Gamma uses the same simultaneous integer caps and setup
event as the sparse calculation. The setup-failure term is paid once in
the final union. This method does not justify multiplying expected spectra.
