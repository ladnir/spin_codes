# Temporal selector capacity without a state-history union

This note strengthens the quartet selector argument for the unchanged
construction. It adds the number of coordinate edges between consecutive
selected quartet layers. The resulting setup event controls every matching
actual state history at once. It does not establish the full 45%/40-bit
distance certificate.

Use M=8192 outer rows, T=4096 updates of width n=512, N=nT, and
J=N/4=524288 quartets. Fix a candidate eligible-row set H of size h.
At each quartet node let j be its number of H-labelled input slots.
For r in {1,2,3}, define

    P_r = the union of quartet nodes with j>=r,
    Z_r = sum over all quartet nodes of (j-r)_+.

Join a quartet at update t to a quartet at update t+1 once for each
coordinate position they share. Each node has at most eight incident
temporal edges. Let E_r count edges whose two endpoints belong to P_r.
There is no edge before the first or after the last actual update.

## A deterministic capacity inequality

For an actual support history let U_t be its live quartet union, let
B_t be its forced births outside U_t, and use the established statistics
Ls, f, F, G=F-f and C_out=R-128h. Every forced birth belongs to H.
The selected live nodes therefore have H-label sum G-C_out and count
G/4. Put

    z_profile = (1-r/4)G-C_out.

The selector slack has the exact expression

    Delta = Z_r-z_profile
          = sum_unselected (j-r)_+ + sum_selected (r-j)_+.

At most Delta selected nodes lie outside P_r. Consequently the number
E_selected of temporal edges between selected nodes satisfies

    E_selected <= E_r+8 Delta.                         (1)

Write a=|S_T| and f_T=|B_T|. Every state position before the last update
that is not a forced birth gives a selected-selected edge. Hence

    E_selected >= Ls-a-f+f_T >= Ls-512-f.

Combining these inequalities gives the universal setup constraint

    E_r+8Z_r >= Ls-512-f+8 z_profile.                   (2)

This argument does not enumerate support histories or assume that temporal
edges are independent. The terminal state is free.

## A setup-only moment

First label input slots independently eligible with probability p, where
0<p<1. Conditioning each region to have exactly h eligible slots restores
the actual uniform regional routing for fixed H. Put

    pi_h(p) = binom(M,h) p^h (1-p)^(M-h).

The conditioning costs pi_h(p)^(-256). For beta>1 define

    q0 = sum_(j<r) binom(4,j) p^j (1-p)^(4-j),
    q1 = sum_(j>=r) binom(4,j) p^j (1-p)^(4-j)
                                  beta^(8(j-r)).

After weighting by beta^(8Z_r), the mass for |P_r|/4=b in one update is

    w_b = binom(128,b) q0^(128-b) q1^b.

Conditional on b, the set P_r is a uniform 4b-subset. These weighted sets
are independent at distinct updates before conditioning the regional
quotas. If the preceding population is 4a, the overlap moment is

    H_ab = sum_j binom(4a,j) binom(4b,j)
                    (beta-1)^j / binom(512,j).

Thus the setup moment equals e0 K^T 1 for K=H diag(w). The artificial
initial population is empty and adds no temporal edge. The final vector
1 retains every terminal population. The floating diagnostic
[quartet_selector_temporal.py](quartet_selector_temporal.py) evaluates
this positive 129-state expression; it is not a rounding certificate.

## A rational scalar supersolution

An exact bound needs no matrix powers. Choose z>=1 and lambda satisfying

    lambda >= q0+q1 z,
    lambda >= q0/z+q1 beta^4.                            (3)

Log-convexity of q0+q1 z beta^c for 0<=c<=4 gives

    (q0+q1 z beta^c)^4 <= lambda^4 z^c.                 (4)

All five integer values of c can also be verified directly with rational
arithmetic. Define

    F(X)=q0(1+X)^4+q1 z(1+beta X)^4.

For v_a=z^a, the exact row identity is

    (Kv)_a = [X^(4a)] F(X)^128 / binom(512,4a).

Applying (4) coefficient by coefficient proves Kv<=lambda^128 v.
Since v_0=1 and 1<=v, positive composition yields

    E beta^(E_r+8Z_r) <= lambda^J.                      (5)

Equivalently, for any fixed current partition, its quartet occupancies c
obey the same local inequality. Averaging partitions supplies the row
identity above but is not needed to justify the scalar majorant.

After quota conditioning and the one H-set union, every score threshold B
has the bound

    Pr[there exists H of size h with E_r+8Z_r>=B]
      <= binom(M,h) pi_h(p)^(-256) lambda^J beta^(-B). (6)

There is no state, profile, death-count or marker union in this fixed-tail
bound. Choosing the best valid markers for one threshold only improves
an upper bound on the same event.

## Exact checker and initial classes

[quartet_selector_temporal_certificate.py](quartet_selector_temporal_certificate.py)
checks (4) exactly and evaluates (6) with the directed rational logarithm
helpers from [quartet_selector_capacity.py](quartet_selector_capacity.py).
It forms no production matrix and needs no floating witness. Its default
rational parameters are

    p=5/16, r=2, beta=12/11, z=77/64,
    lambda=10613108988597669204677/8760324352986371325952.

The serial exact replay passed all 90 local rational comparisons and
90 independent small-row comparisons. Directed logarithms verify the
production setup moment bound

    E beta^(E_2+8Z_2) <= 2^145119.

The resulting setup classes and example scores are:

| h | Certified score threshold for a 2^-92 tail | Example score | Integer upper exponent at the example score |
| --- | ---: | ---: | ---: |
| 3584 | 2062541 | 2155602 | -11774 |
| 3072 | 1444966 | 1492050 | -6002 |

Every actual profile at the specified h whose score in (2) reaches the
certified threshold belongs to that entire setup class. The examples
used to identify the comparisons are respectively

    (h,R,Ls,f,F)=(3584,819200,845414,4,1048576),
    (h,R,Ls,f,F)=(3072,819200,837222,4,1015808).

All death counts are covered simultaneously. A profile union allowance
must not be charged again. The exact checker and this derivation passed
independent static review; no production matrix or floating witness is
part of the certificate.

The callable `temporal_scalar_threshold(h)` evaluates the same fixed
rational markers for any 1<=h<8192. One certified threshold for each h
at r=2 costs less than 2^-79 in total. Even if all three r values are
used, at most 8191*3 such 2^-92 classes cost less than 2^-77. Marker
choices for the same h,r can be minimized before choosing its threshold,
because they bound the same nested setup tails.

The fixed-marker thresholds have been applied to the existing
diagnostic grid after the selector-capacity cuts. They exclude another
727 points, including 377 previously FLOAT-open points. There remain
23979 FLOAT-open points. The worst sampled tuple is now
(h,R,Ls,f,F)=(2560,786432,779879,4,950272).
A finite grid cannot prove coverage of its complement, and a surviving
tuple is not evidence of a bad codeword. The next global inequality
retains the selected-quartet count as well as temporal overlap.
