# Count-path entropy and fresh coordinates

This continues [the history-conditioned potential](QUARTET_HISTORY_POTENTIAL.md).
The construction, once-sampled setup, and full 45%/40-bit goal are unchanged.
The new statements concern a positive transfer used to upper-bound the
geometric moment. Numerical values below are floating-point diagnostics,
not directed certificates or whole-code security margins.

## Why retaining the input count was insufficient

The five-type histogram records the current input activity and the length
of the preceding free run. Type0 has inactive input. Type1 has active input
but no preceding free run. Types2,3,4 have preceding run lengths1,2,>=3.
The previous envelope fixes type0 and maximizes over all four active types.
It can therefore make too many coordinates old after every update.

At the diagnostic class r=2100,c=111572, use
lambda=707/2560, tau=381/512, log2(b)=4, and threshold B=165556.
Keep the regional quotas:44 regions have6527 allowed inputs and212 have6528.
These are only selected payment counts; other allocations remain open.

The count-band calculation sums every input-count path in its stated band,
and every placement of the44 lower-total regions. Its best tested envelope
with the refined point-mass correction gives the following log2 expressions:

| Allowed input counts per update | FLOAT upper expression |
| --- | ---: |
| 400..416 | -1705.04 |
| 396..420 | +15501.82 |
| 392..424 | +33746.05 |

Thus a negative near-flat path does not settle count-path entropy.
Independently adapting the original five-type potential to each input count
performs worse in the tested narrow bands. Its wider underflowing diagnostic
was discarded, not interpreted as a zero contribution.

## A coefficient bound with a point-mass factor

For positive d_i and z, let independent auxiliary Bernoulli variables have
probabilities p_i=d_i z/(1+d_i z). Then

    e_k(d) = z^(-k) product_i(1+d_i z) Pr[sum_i Bernoulli(p_i)=k].

These Bernoulli variables represent a polynomial coefficient. They are not
assumed to be the actual routed input bits. Write V=sum_i p_i(1-p_i).
Fourier inversion and the absolute characteristic function imply

    Pr[sum=k] <= (1/(2pi)) integral[-pi,pi] exp(-2V sin^2(t/2)) dt.

Using sin(|t|/2)>=|t|/pi gives min(1,sqrt(pi/(8V))). For
0<theta<pi and c_theta=1-theta^2/12>0, splitting the integral at theta gives

    Pr[sum=k] <= 1/sqrt(2pi V c_theta)
                + exp(-V theta^2 c_theta/2).

Indeed, on the central interval sin^2(t/2)>=t^2 c_theta/4.
On its complement sin^2(t/2)>=theta^2 c_theta/4. Both displayed
upper bounds decrease with V, so a uniform variance lower bound suffices.
The implementation takes the minimum with1 and the coarse bound, at
theta in {1/2,3/4,1,5/4}. Directed evaluation remains a separate obligation.

## Retain the actual fresh count

Let w=512. Before an update, let a be the active-input count and r_1 the
type1 count. Thus n0=w-a and n2+n3+n4=a-r_1. Use the positive potential
phi(n)=product_i h_i^n_i, normalized by h0=1, and define

    d=(h1,h2,h3,h4 b^(3lambda-tau),h4 b^lambda).

At the next update let l be the active-input count and k the number of
positions in fully active quartets. The multiplicity J_(l,k) is the
coefficient [x^l y^k] of ((1+x)^4-x^4+y^4 x^4)^(w/4).
Conditional on l,k, the k positions form a uniform subset of the w positions.
Let v0 be the number selected from the known type0 positions. The outgoing
fresh count is EXACTLY

    r_1' = l-k+v0.

Keeping only the preceding k gives r_1>=a-k, which is weaker. The following
transfer keeps r_1' instead and sums over k and v0 jointly.

For each admissible k,v0, the potential-weighted term is

    J_(l,k) / binom(w,k) * binom(w-a,v0) * d0^v0 * h1^(l-k)
    * [x^(k-v0)] product_(i=1..4)(1+d_i x)^n_i
      / product_(i=1..4)h_i^n_i.

Group incoming fresh counts into an interval [rlo,rhi]. For any z>0, put
L_i=(1+d_i z)/h_i and M=max(L2,L3,L4). The remaining coefficient is at most

    z^(-(k-v0)) max_{r in {rlo,rhi}} L1^r M^(a-r).

The maximum over r occurs at an endpoint because its logarithm is affine.
Importantly, this maximum is taken BEFORE choosing z. Optimizing separately
at the endpoints would not justify an upper bound for the interior.

For the point-mass factor let v_i=d_i z/(1+d_i z)^2. A uniform variance
lower bound is

    min_{r in {rlo,rhi}} [r v1 + (a-r) min(v2,v3,v4)].

Sum the terms into the bucket containing l-k+v0. This gives a nonnegative
matrix between incoming and outgoing fresh-count buckets for every a,l.
Products of these matrices retain the hidden count across regional
boundaries. The actual initial state is all type1 after the first input;
only the final update receives the terminal factor. No regional reset or
independent new state bits are granted.

The original potential telescopes as before. The bucket matrix is a uniform
bound over every histogram in the bucket; it is not an assertion that the
bucket labels alone determine the true state distribution.

## Hard-path test and a larger family

The previous scalar envelope selects the regional words

    low  = (424 repeated7,407,392 repeated7,408),
    high = (424 repeated7,408,392 repeated7,408).

For44 low words followed by212 high words, the previous refined scalar
bound gives +21839.24. Retaining only the preceding fully active quartet
count improves it to +14421.31. Five constant-potential centers404..420
do not close it; the best tested center is408. This is a failed bound,
not a code counterexample or a proof that every potential fails.

Retaining the actual fresh count gives:

| Fresh-count bucket width | FLOAT log2 expression |
| --- | ---: |
| 32 | +872.50 |
| 16 | -2261.97 |
| 8 | -3609.80 |

Next cyclically rotate the first15 entries of each word, leaving the final
408 unchanged. There are15 distinct words of each kind. Sum all rotations
independently in all256 regions, and all choices of which44 regions use
the lower total. A regional matrix product, followed by a coefficient DP
in the number of low regions, sums this family while carrying the fresh
count. There are15^256 binom(256,44) input-count paths, with log2 count1165.71.
At bucket width8, the FLOAT expression is **-3133.02380**.

This demonstrates that the new bound survives a nontrivial count-path sum.
It does not cover other regional words, other payment allocations, other
covered counts, or other outer-row counts. A negative floating expression
for this restricted family must not be reported as a whole-code margin.

## Checks, resources, and next test

- `quartet_history_potential.py`:8800 exact outgoing fresh-count coefficient
  comparisons over every width4/8 histogram pass. The previous2129 freshness
  tests,58 coefficient tests,15 stationary floors,65536 input masks, and1024
  full/empty sequences still pass.
- `quartet_fresh_count_screen.py`:64 exact small matrix paths verify the
  count/hidden-state DP;30 regional words have the required distinctness
  and totals. The full-width rotation result replays identically.
- `quartet_count_band.py` retains regional totals and sums bounded bands.
  `quartet_freshness_screen.py` preserves the weaker preceding-k diagnostic.
- Numerical replays use Peach, one thread, a90s watchdog, no files, and
  bounded arrays. The final exact audit takes about1.08s; the rotation
  screen about7.57s. These are proof-computation times, not encoding timings.

Next extend the fresh-count transfer beyond selected regional words.
First test all count paths in a bounded band while retaining the regional
total, then bound the omitted count tails. Only after this diagnostic
survives should a directed replay and full class coverage be attempted.
Preserve all older proofs; no construction or failure-budget change follows.
