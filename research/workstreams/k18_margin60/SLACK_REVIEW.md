# Proof tightening for the current K18 construction

This review keeps the certified t64/s22 maps, four parallel GF16 RS[16,8] rows,
and the ideal independent uniform GL setup distribution. It introduces no
construction change and supplies no new numerical certificate.

The retained whole bound has 42.090840 bits of margin. Its q=1 component has
60.751047 bits; q=2 has 126.546358 bits. Keeping those endpoints requires the
new q>=3 tail to exceed approximately 61.30106 bits for a 60-bit whole union.
The retained intermediate tail has 42.090844 bits. Its weakest individual
occupancy is q=89, at 44.841583 bits.

## Cheapest valid choices

`k16_design/packet_rs_state_ladder.py::run_tail` fixes the activity argument of
`kernel_t64.local_operators` to 1/2. This argument only selects a complete
newborn-state representation in `s16_closure/kernel_birth_density.py::refine_local`.
It does not change the packet input law. Each activity in [0,1] therefore
gives another valid witness. At q=89, the expected physical-slot activity under
the uniform outer comparison is 15q/(16L)=1335/32768, where L=2048.
An initial screen can compare activities .02, .04, .08, .16, and .5.

Each newborn-row candidate represents the same tilted state measure using
uniform-density and expansion-weight coordinates. Choose one complete row
per physical occupancy. Entrywise minima across different candidates can
discard both parts of that measure and are invalid. Taking the minimum of
complete scalar bounds for a fixed q is valid.

`k16_design/packet_rs_k20_dense.py::physical_candidates` caches all existing
outward row candidates at a fixed tilt. The continuation potential in
`select_physical` is computed from the baseline matrices before row selection.
Recomputing that potential after installing the selected rows offers another
small heuristic refinement. Every iterate remains a valid collection of rows;
its final regional scalar bound must still be evaluated.

`kernel_birth_density.thresholds` considers only six quantiles of the positive
state masses. Additional caps are valid, but their numerical benefit is unknown.
This is a lower-priority knob because it needs the state-law arrays.

## Fractional conditioning on routing

The following extension can use the same physical operators if ordinary
witness tuning is insufficient. It changes the probability bound.

Fix a set S of q active message groups. Let E_S be the event that some message
with exactly those active groups has output weight at most D=52428. All
probabilities below are over the ideal setup. Let J record the number of
slots assigned to S in each chronological physical step, including slots
whose packet value will be zero. Each region has exactly q such slots.

Condition on a feasible value of J. The outer maps and inner GL updates retain
their independent distributions. Within each step, the assigned slots form a
uniform subset of the specified size. These subsets are independent across
disjoint steps, conditioned on their counts. Regional permutations are
independent between regions.

For lambda>0 and z=exp(-lambda), let L_k(z) be a valid physical envelope for k nonzero
packets with independent uniform nonzero four-bit labels. Define

\[
 P_j(z):=\sum_{k=0}^j {j\choose k}(15/16)^k(1/16)^{j-k}L_k(z).
\]

Thus P_j includes the zero values of j uniform four-bit packets. The existing
pointwise outer comparison has beta=2^256/(2^16-1)^8. It remains valid after
fixing routing, because the outer setup is independent of routing. If e_0
selects the initial zero state and 1 sums all terminal envelope coordinates,
the conditional expected count is at most

\[
 A(J):=\beta^q z^{-D}\,e_0\prod_t P_{J_t}(z)\,\mathbf1.
\]

Markov's inequality gives Pr[E_S | J] <= min(1,A(J)). For 0<alpha<=1,
min(1,x)<=x^alpha. Partition the chronological matrices into blocks and
multiply within each block first. For any resulting nonnegative matrices Q_b,
subadditivity of x^alpha gives

\[
 \left(e_0\prod_b Q_b\mathbf1\right)^\alpha
 \le e_0\prod_b Q_b^{\circ\alpha}\mathbf1.
\]

Here the superscript circle denotes entrywise exponentiation. This inequality
expands scalar matrix paths. It does not require envelope
coordinates to be disjoint probability states.

For two physical steps per block, define the block family

\[
 T_j:=\sum_{a=\max(0,j-16)}^{\min(16,j)}
 \frac{{16\choose a}{16\choose j-a}}{{32\choose j}}
 (P_aP_{j-a})^{\circ\alpha}.
\]

The hypergeometric weights remain outside the fractional power. Exact regional
placement of q potential slots among 2048 slots gives

\[
 R_q:=\frac{[x^q](\sum_{j=0}^{32}{32\choose j}T_jx^j)^{64}}
                 {{2048\choose q}}.
\]

The contribution of all q-group subsets is consequently bounded by

\[
 {2048\choose q}\,\beta^{q\alpha}z^{-D\alpha}
 e_0R_q^{64}\mathbf1.
\]

There is no additional outer binomial mixture after this placement: P_j already
includes packet values of zero. The group-subset union remains outside alpha.
The full matrix power retains the state between regions.

At alpha=1, thinning and hypergeometric placement recover the existing
uniform-input bound for the same selected local rows. This supplies a useful
regression check. Numerical evaluation needs nonnegative outward endpoints,
including the fractional powers, before it can improve a certified margin.

The same conditioning argument appears in
`packet8_codesign/iteration3/fractional_gate.py` and
`packet8_codesign/iteration3/fine_grouped_gate.py`. Their byte-native numerical
results are not endpoints for the current four-bit construction.

Recommended order: refine the ordinary tilt and complete-row choices; screen
the two-step fractional bound near q=89 if needed; certify and union every
affected occupancy only after a useful numerical witness is found.
