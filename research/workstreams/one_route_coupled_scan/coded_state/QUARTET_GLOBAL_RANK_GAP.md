# Global rank gap and forced activation

The construction is unchanged. The full K=2^20, rate-1/2, 45%-distance,
40-bit target remains open. This checkpoint adds a proved global class
and records complete-history upper bounds that still fail.

The probability space is the once-sampled accepted setup in
[the inverse-rank reset](QUARTET_GLOBAL_RANK_RESET.md). The field is
J=GF(2^32), q=2^32, N=2^21, W=943718. The accepted-density factor satisfies
gamma<=1024/989. Every argument retains the free terminal state.

## Rank need only exceed the available state positions

Fix the route and quartet partitions. For a candidate covered footprint
V=(V_1,...,V_T), let R(V) be its RS parity rank. Let A_t contain positions
of V_t whose next quartet lies wholly in V_(t+1); set A_T=V_T. Define
L=sum_t |A_t|. These are the positions where a compatible state can be
nonzero, including the terminal positions.

The grouped-footprint argument has summed all compatible state and output
supports. At output marker x=1, it gives

    Pr_field[nonzero word with V | route,partitions]
      <= (gamma/q) q^(-R(V)) (q+2)^L.                 (1)

The event R(V)-L>=65537 therefore has total probability at most

    2^N (gamma/q) q^(-65537) (1+2/q)^N
      = 2^-64 gamma (1+2/q)^N
      <= 2^-64 * 1048576/1011747
      < 2^-63.                                      (2)

The rational majorant uses (1+2/q)^N<=1/(1-2N/q)=1024/1023. There are at
most 2^N footprints. No additional union over state supports, output
supports, eligible-row sets, or quartet partitions belongs in (2).
Condition on the geometry, apply (1), then average.

This additional class covers low absolute rank when L is smaller. It
does not replace the previous absolute-rank theorem. The remaining event
must satisfy R-L<=65536. Two independent audits agree; the helper's
`rankgap_certificate` checks the rational arithmetic exactly.

## Impose the remaining event before enlarging its count

Fix a candidate eligible-row set H of size h. On the actual classification,
R=128h+C_out, where C_out counts covered positions outside H. For eta>=1,
multiply the remaining-event count by eta^(65536+L-R), which is at least
one on this event. In the transfer this gives

    global factor eta^(65536-128h),
    m -> eta m,  b -> eta b,
    r = p + (1-p)/(a eta),  a=q^theta.                (3)

Here p is the independent eligible-label probability used before exact
regional-quota conditioning. Scale both m and b after computing them;
do not recompute b from the scaled m. This marks available positions,
not total coverage. Keep eta fixed across the history. Do not inherit
the ordinary outside-degree marker's different normalization.

The earlier minority exclusions were also inserted directly into the
130-state transfer. For C_H=|V intersect H|, their remaining event implies

    m_V-32C_H <= 7217-4096h,
    m_V+11L-32C_H <= 7217-4096h+11W.

Two nonnegative log markers implement these inequalities. Minority uses
the actual min(v,512-v) for each outgoing footprint, including the final
one. The eligible-count normalization is compensated explicitly. These
substitutions passed independent audit and 25 exact small comparisons,
but did not improve the hard tested bounds.

## A propagation error in the earlier relaxation

This is a proof overcount, not an implementation defect. The old
130-state kernel permitted coverage 1->512->1. Only one quartet in the
first expansion can meet the old covered position. The other 127 full
quartets contain 508 forced births. Their next quartets must be fully
covered, making the immediate collapse to 1 impossible.

Three valid alternatives were tested without changing the construction:

- A 513-state count charges every fresh birth's eligible-row cost, but
  drops forced propagation. It loses more than it gains on the tested
  middle classes. It remains a diagnostic, not the preferred bound.
- A 307-state factorization retains a guaranteed number of fresh full
  quartets. Incoming coverage bins have upper endpoints 16,64,127,512.
  If u<=U, at least (k-U)_+ of k full new quartets are fresh. Select that
  many actual fresh quartets with auxiliary randomness of total weight
  one, independently of new routing labels. No binomial penalty is
  needed. Their positions remain forced in the next update. This
  preserves every previous forced-position restriction.
- A separate 129-state relaxation retains every canonical fresh birth
  and its propagation, while permitting extra fully forced quartets.
  This is the most compact current global control.

For the last alternative, z_f,z_o mark incoming forced/optional coverage,
and y_f,y_o mark outgoing coverage. Its local polynomial is

    y_o^4 [(1+m z_f+b z_o)^4-1]
       +(1+z_o)^4 (1+y_f)^4.                         (4)

An optional full quartet must meet incoming coverage. All other outgoing
covered positions are forced. Relative to the canonical kernel, (4)
adds the nonnegative term y_f^4[(1+z_o)^4-1]. This permits a fully forced
quartet with incoming optional coverage, whose state must then be zero.
No canonical trajectory is discarded. The relaxation is not claimed to
be pointwise smaller than every previous relaxation.

For (3), apply eta^L to the canonical kernel first, then add the positive
extra term. The extra annotations need not track geometric L exactly;
they still dominate the canonically weighted count. Terminal factors
must use the same scaled m,b.

## The new 129-state factorization

For i,k in 0,...,128, let chi mark coverage. Before additional markers,
r=p+(1-p)/a. The reduced matrix is

    C_(i,k) = binom(128,i)(r chi)^(4i)/binom(512,4i)
       * [z^(4i)] [ (1+m p chi+bz)^4-1 ]^k
                         (1+z)^(512-4k).            (5)

All coefficients are nonnegative. Binomial forced-position counts cancel
the physical multinomial denominator, leaving degree-512 univariate
polynomials. The terminal vector is

    g_i = binom(128,i)(r chi b)^(4i)
                       (1+m p chi)^(512-4i).        (6)

Exactly 4096 updates contribute e_0 C^4095 g: 255 regular regions C^16
followed by D=C^15 g. At zero coverage, Z=e_0 e_0^T and D_0=e_0.
Positive recurrences compute C^16-Z and D-e_0 without cancellation.

Each region has 8192 routing labels. For pi=Pr[Bin(8192,p)=h], conditioning
charges pi^-1 only on the positive parts. The physical factor A and zero
operator are common across p. Entrywise minima of the quota-adjusted
regional matrices and final vectors are therefore valid.

For those positive bounds Q,d, define A_j=e_0 Q^j e_0 and B_j=e_0 Q^j d.
The weight of exactly k occupied regions is
[x^k] A(x)^(256-k)(1+x B(x)). Sum k>=129. This treats the last region
separately and retains every endpoint.

## Results and limits

All full-width numbers below are FLOAT diagnostics, not certificates.
The new 129-state control uses 9 quota tilts.

| Eligible rows h | theta | log2 z | Unrestricted | Best tested rank-gap marker |
|---:|---:|---:|---:|---:|
| 512 | 1/8 | -1/4 | +55024 | +55024 |
| 2048 | 1/16 | -1/2 | +151203 | +151203 |
| 4096 | 1/16 | -1 | +454173 | +411663 at log2 eta=1/2 |

Positive upper bounds do not exhibit low-distance words or a structural
failure of the code. The propagation fix and rank-gap theorem do not
close the global remainder by themselves.

All helper checks pass: 261 prior factorization entries, 8 prior composed
paths, 25 eligible substitutions, 6 renewal comparisons, 25 minority
comparisons, 106 fresh-birth entries, 261 coarse-forcing entries, 13 new
factorization entries, 8 new composed paths, terminal/region identities,
and the rational rank-gap bound. The exact width-8 negative control
rejects the impossible 1->8->1 trajectory.

Jobs were serial, single-threaded, file-free on Peach, with 175/180-second
and 2 GiB guards. No benchmark, construction, production, paper, commit,
or push change was made. The frozen pair helper is unchanged.

Next isolate the joint rank/available-position classes driving the global
bound. Distinguish a limitation of fractional-moment grouping from a
remaining geometric overcount. Preserve the invariant-plateau control;
do not resume local seven-type refinements or promote a selected negative
class into a full certificate.
