# Regional coverage bins before the quota minimum

This note records a bounded comparison for the unchanged quartet
construction. It uses the actual-state kernel from
[the joint MDS count](QUARTET_ACTUAL_MDS_JOINT.md). The full 45%/40-bit
certificate remains open. The reported evaluations are floating upper-bound
diagnostics; their positive values are not codeword witnesses.

The implemented control is
[quartet_actual_region_bins.py](quartet_actual_region_bins.py).
It retains a regional coverage interval before taking the eligible-quota
minimum. It does not extract an exact global coverage coefficient.

## Comparison and interface

Fix a complete candidate support history before averaging partitions and
eligible labels. In one region let F_reg=G_reg+f_reg be the total covered
positions, with 0<=F_reg<=8192. The existing markers z and nu count G_reg
and f_reg. Replacing both by z*beta and nu*beta therefore multiplies every
history by beta^F_reg.

For the Bernoulli eligible-label tilt p, define

    Q_p(beta) = K_p(z*beta,nu*beta)^16 - Z,
    Z = e_0 e_0^T,
    pi_h(p) = binom(8192,h) p^h (1-p)^(8192-h).

Here Q_p(beta) is formed by positive composition, with Z kept separately.
Its (0,0) entry includes nonempty regional excursions and is retained.
The exact conditional operator for a nonempty coverage bin I=[a,b] obeys
the entrywise bounds

    Q_I <= beta^(-a) Q_p(beta)/pi_h(p),   beta>=1,
    Q_I <= beta^(-b) Q_p(beta)/pi_h(p),   beta<=1.       (1)

Indeed, the appropriate multiplier times beta^F_reg is at least one on
the bin, and the quota indicator can then be dropped from the positive
sum. Every candidate in (1) acts in the same physical-population basis.
For bins partitioning all integers 1 through 8192, set

    B_I = entrywise min_(p,beta) of the right sides of (1),
    B = entrywise min(sum_I B_I,
                      entrywise min_p Q_p(1)/pi_h(p)). (2)

Then Z+B dominates the exact conditional regional operator. There is no
additional bin union factor. Global composition uses (Z+B)^256 with
initial row e_0 and free endpoint vector 1. The usual global profile
marker and field corrections remain external. There is no separate
routing factor or external quota divisor.

The public `ActualRegionBinComparison` accumulator accepts scaled,
already quota-adjusted Q matrices carrying beta^F_reg. Its `add` method
updates the per-bin entrywise minima; `finish` sums and caps them as in
(2). `actual_region_binned_kernel` constructs the candidates using the
existing MDS evaluator. The default uses 32 bins of width 256, three p
values, and three beta values, so the bins share nine kernel builds.

The unsplit cap is also retained in the reported moment bound. This
matters for finite Collatz iteration, which can choose different witness
vectors for the two matrices. Optional upward positive floors retain
structural zeros; they do not provide directed-rounding certification.

## Checks and bounded controls

The serial replay passed 625 exact bin-comparison inequalities, two
empty/excursion checks, and one free-endpoint check. The reference count
enumerates width-four old/new support and eligible-label masks, records
both regional coverage and eligible count, and conditions on five
eligible labels among eight positions. The replay also passed 25 floating
comparisons against the rational matrix and 50 floor-comparison checks.
The implementation and probability bookkeeping passed independent static
review.

Both controls below use the fixed profile

    (h,R,Ls,f,F) = (4096,786432,786432,4,1048576),

with x=2^-16, theta=1/16, s=2, z=1, nu=2^-18, eta=1/4,
t=(qx^2/(1+qx))^theta, and p in {1/5,1/3,1/2}. Each retains all 32
coverage bins and uses the unchanged fixed-profile correction.

| Coverage tilt log2(beta) | Kernel builds | Result |
| --- | ---: | --- |
| {-1,0,1} | 9 | Zero entrywise gain; the profile bound remains approximately +675044.0045 bits. |
| {-3/4,-3/8,-1/16,0} | 12 | Zero entrywise gain; the profile bound remains approximately +675044.0045 bits. These transition tilts were derived from the recorded per-p growth intercepts. |

Neither completed control improves the bound. This does not prove that
every coverage partition or tilt family is ineffective. Separately, the
three-p exact-f=4 control returned
approximately +674889.164 bits, only about 155 bits below the earlier
bound. This agrees with the
[persistent birth-free family](QUARTET_LOW_BIRTH_MOMENT_OBSTRUCTION.md):
low birth count alone does not remove the dominant relaxed histories.

Both regional-bin controls are now parked. A separately audited
selector-capacity argument bounds the probability of any matching
hard-profile footprint by 2^-4094, simultaneously over all selected
quartet sets and all death counts. Its hypergeometric comparison removes
the regional quota penalty. It closes this profile without a
state-history union. The positive moment diagnostics above remain method
limitations, not evidence of a bad codeword. The full-code certificate
remains open; the next structural step is the selector argument for the
other eligible-row counts and profiles.
