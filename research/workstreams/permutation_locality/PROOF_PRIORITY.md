# Close a roughly 5 ms design

## Current objective

The user accepts the measured 1.88--1.90x speedup. A strict 2x speedup is no
longer the gate. Prioritize a complete distance certificate for any nearby
design that takes roughly 5 ms at K=2^20 on 128-bit elements. Preserve the
10% relative-distance target and at least 40 bits of setup-failure margin.
Setup is precomputed and excluded from encoding time.

## What is established

The one-column, four-row grouped permutation with IMT(128,19), one update,
takes 5.078410 and 5.086134 ms on the two recorded route seeds. Matched
production controls take 9.652769 and 9.576125 ms. These measurements use
the full transposed encoder, not just the permutation.

For this same distribution, outward certificates cover all messages with
one through thirteen active groups, including all ranks, supports, and
group locations. Their summed failure contribution is below 3.188222e-15
(48.1562 bits). There are 2048 groups in total. No bound on the complete
remaining occupancy range has been established. In particular, thirteen
covered occupancies must not be described as nearly complete coverage.

The latest full fourteen-group binary64 search, with pair-conditioned
tails, stopped at log2 upper -19.8612 after 150 splits and 1128 retained
leaves. It did not pass to outward replay. Conditioning tail probabilities
on the entering expansion-weight class gained about 1.45 bits at one
diagnostic point. A separate restricted-density refinement gained about
1.07 bits. Neither result is a certificate or a low-distance witness.

## A nearby alternative

Two independent transvections per step preserve the permutation and fixed
expansion/feedback maps but change the inner. They reduce the probability
of retaining a fixed nonzero state from the lazy part 1/2 to 1/4. The
remaining part refreshes uniformly over nonzero states.

The two-update implementation now measures 5.40--5.43 ms. A complete
fourteen-group cover passes outward replay at 152.28 bits. See TWO_UPDATES.md
for the rescaling argument, checks, timings, and exact certificate scope.
The one-update certificates are not transferred. Full-code closure remains
open, but the measured overhead makes this a competitive proof-first path.

## Closure-first work order

1. Completed the two-update timing gate, with dense-reference, adjoint,
   and sanitizer checks. The added cost is 6.5--7.0%; retain this variant
   for the proof-first investigation, without another speedup campaign.
2. Establish a route to the entire occupancy range before investing in
   successive small-occupancy refinements. The current support-box cover
   is not yet a scalable proof for all 2048 occupancies. Separate sparse,
   intermediate, and dense regimes, with explicit overlap and failure
   budgets. Do not presume that the dense regime closes automatically.
3. For a competitive two-update implementation, audit the transformed
   state envelope and test full occupancy covers, including mixed supports.
   Otherwise retain the one-update implementation and target its shared
   cancellation/density bounds, rather than isolated support points.
4. Close every regime, replay the final bounds with outward arithmetic,
   sum the failure contributions, and match the certified setup and maps
   to the measured implementation. Only then claim a full certificate.

Production defaults and the paper remain unchanged. See ONE_COLUMN.md for
the distribution and matched timings, and CDF_COVER.md for certificate scope.
