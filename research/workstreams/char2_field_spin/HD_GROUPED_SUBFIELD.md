# Grouped outer with byte-field randomization

## What changed

The [grouped-outer study](HD_GROUPED_OUTER.md) retained GF65536 randomness.
Its small RS constituents save arithmetic, but the extra random diagonal
can erase that saving. This candidate instead samples all outer and inner
diagonal entries in GF256. Each random scalar multiplication then needs one
GFNI byte multiply, without the second product, byte swap, and XOR of the
GF65536 kernel. This changes the setup distribution and needs its own proof.

At K=2^20, N=2^21, the new complete certificate gives at least35% symbol
distance with72.3975395430 bits of whole-code setup-failure margin. Directed
256-bit search and an independent384-bit replay agree. The first matched
measurements put the four-row candidate around 4.8 ms, essentially tied with
the old control. These shared-host measurements do not establish a speedup.

## Construction

Let F=GF256. A concrete group contains four RS[64,32] constituent rows over
F, for128 input symbols in total. Collect the four outputs at each column
and multiply by D2 H4 D1. Here H4 interpolates nodes0,...,3 at4,...,7, and
the diagonals have independent uniform nonzero F entries. Each group and
column has independent diagonals. Splitting the four outputs gives256
coordinates per group. Coordinate4*j+a enters region4*j+a.

The8192 groups give256 regions, each independently and uniformly permuted.
The inner processes64-symbol packets with16-symbol state Q, initially zero:

    Y = X + A Q,
    Q_next = D2 H16 D1 (Q + C X).

A is the direct sum of four RS[16,4] emission maps and C=A^T. H16
interpolates nodes0,...,15 at16,...,31. All inner diagonals are independent
uniform nonzero F entries, independent of the outer and route. Setup is
sampled once and reused; the final state is unrestricted. As before, CA=0
and the systematic recurrence is invertible.

The same outer proof allows g=2,8,16 constituents of length256/g and
dimension128/g, with the matching g-coordinate column mixer. These are
alternative implementations of the same256-region proof interface.

The native encoder applies the scalar code to128-bit records. It can be
viewed as sixteen F lanes or eight GF65536 lanes. To justify the latter,
expand any nonzero GF65536 message in a basis over F. At least one component
message is nonzero, and its encoded support is contained in the full output
support. Thus extending scalars preserves the guaranteed symbol distance.
There is no additional union bound over lanes.

## Outer measure and the sparse-case problem

For every fixed256-symbol output, the grouped-outer membership probability
is at most255^(-128), by the conditional-information-set argument in the
grouped-outer note. Hence the expected nonzero-message measure is dominated
by beta times uniform F^256, with beta=256^256/255^128. Independent groups
retain the regional proof's product comparison.

For the inner, write d=q^16-1 and b=d/(q-1)^16, with q=256. The established
local operator can charge b only on transitions to nonzero state:

    [[K,       b(T_j-K)],
     [Jbar/d,  b V]].

Even this bound wastes substantial mass during long runs with no potential
input positions. The initial one-active-group screen gave log2 bound+86.33.
The improvement below changes only that proof bound, not the construction.

## An inverse-refresh bound for empty packets

Fix a packet whose potential input set is empty, so X=0. Suppose the incoming
weighted measure is at most B/d at each nonzero state. For a fixed nonzero
next state y, the outgoing weighted density is at most

    (B/d) E_T[z^wt(A T^-1 y)],       0<z<1.

This follows by changing variables from the previous state to T^-1 y in
the weighted transition sum. The current refresh T is independent of the
incoming measure. Because H16 is totally nonsingular, H16^-1 is also totally
nonsingular, by the complementary-minor identity. Moreover,
T^-1=D1^-1 H16^-1 D2^-1 has independent nonzero diagonals.

Here is an exact way to evaluate a sufficient upper bound. For state size s,
let P[a,b] be the probability that T^-1 maps a fixed weight-a vector to a
weight-b vector. Every b-element output support has the same probability.
Indeed, [I;H^-1] restricted to the a active input coordinates is MDS; counting
vectors with prescribed nonzero input and output supports uses only a,b,s.
Put t=a+b-s. The count for each prescribed output support is zero for t<=0;
otherwise it is

    M(a,b) = (q-1) sum_{i=0}^{t-1} (-1)^i binom(a+b-1,i) q^(t-1-i).

Consequently P[a,b]=binom(s,b) M(a,b)/(q-1)^a. These rows normalize exactly.
The final inverse diagonal also makes values uniform conditional on support,
although the support-only bound below does not need that stronger property.

Partition the s state coordinates into blocks of h coordinates. Write
L=64/(s/h), and let c count nonzero state blocks. Each active MDS emission
block contributes at least L-h+1 output symbols. Conditional on state weight
b, support is uniform, so define

    C_b(z) = sum_c n[b,c] z^(c*(L-h+1)) / binom(s,b),
    n[b,c] = [x^b y^c] (1+y*((1+x)^h-1))^(s/h).

Then L0(z)=max_a sum_b P[a,b] C_b(z) bounds the inverse expectation for
every nonzero y. For the empty packet, replace the old live-to-live entry
by min(b V,L0(z)). The zero-to-zero entry remains1; both cross entries
remain0. Other packet profiles retain the existing local operator.

This removes the unnecessary repeated density inflation in empty stretches.
It does not assume that the refresh is uniform, that successive codewords
are independent, or that a nonempty potential set always contains a nonzero
symbol. The numerical implementation evaluates the support counts exactly
as integers and the bound as a rational before upward rounding.

## Certificate and validation

The inverse bound changes the one-active-group FLOAT estimate from+86.33 to
-114.07. The directed value for the pinned witness is at most-113.8732.
The complete certificate uses37 rational witnesses and covers every active
group count1,...,8192. It also includes a directed dense tail starting7665;
overlap is harmless because the driver sums the bounds conservatively.

The final log2 union bound is

    -72.39753954305973881403725601752170114922055371676681560409900212211082.

The384-bit replay encloses this value with radius below1.42e-107. It excludes
weights0,...,734003, hence proves minimum distance at least734004>0.35*N.
No density loss is added afterward: each packet's actual refresh bound is
already inside the fractional moment. The claim concerns ideal independent
setup randomness, not the reproducible MT seed used by the native tests.

Exact small-field tests enumerate every mixer input/output atom and reproduce
the inverse weighted bound. They also check the weight kernel's stationary
distribution, normalization, and agreement with the old proof interface.
The existing35% proof regression tests continue to pass.

The native implementation uses a distinct byte-refresh setup type and
compile-time dispatch. It preserves the fixed-width SIMD schedule, padded
regional layout, prefetching, and streaming stores. Release checks cover
all four outer group widths, both local ranks, alignment, guards, and aliasing.
Independent K20 dense checks for the certified h4 path give:

| Outer grouping | K20 checksum |
|---|---:|
|Four RS[64,32] rows|10655490065220255253|
|Eight RS[32,16] rows|6472899548487601650|

The same small-size matrix of group widths, local ranks, alignment cases,
and aliasing checks also passes ASAN/UBSAN.

The cheaper h2 inner is tested arithmetically but is not covered by this
certificate. The current GF256 proof is for35%, not40%. At40%, even the
all-groups dense-tail envelope is positive (about+359679 log2 units), so the
existing40% certificate cannot simply be transferred to this smaller field.
That is a limitation of this proof, not a counterexample to the code.

## Reproduction and next decision

Run `python -B hd_grouped_subfield_witnesses.py --precision 384` in the
workstream's theory directory. The pinned replay never performs a numerical
search. The exact inverse-refresh kernel is in `hd_grouped_quiet.py`.

For the native candidate, build `spin_char2_hd_grouped_forward` with
`HDGROUP_SUBFIELD=1` and `HDGROUP_FUSE_LAST=1`. Run `--check` and `--large-check`.
Measure `--bench 501 full 20 nt 4 4` and the g8 alternative, bracketed by
the frozen GF65536/h2 RS256 control. Setup and allocation are outside the
timed interval. No production profile or paper timing has changed.

The validated candidate binary has SHA256
`995b76399f71bc2c15d8dbeed43fe64427adb6e5a86983b483da261988d16071`;
the preserved35% control has SHA256
`81bc8f93fe7f1dc1380a4c91ab9d5679975e32b9b6a8ba668d3a60a85b722785`.
All24 targeted Python regression tests and the rebuilt older native encoder
checks pass. The local kernel probe exits77 (unsupported CPU), and the
configured AWS connection times out. The initial timing pause ended when
the user authorized isolated-core measurements alongside the unrelated job.

## CPU 2 comparison, 2026-10-08

All our measurements ran serially on Peach CPU 2, using both benchmark
locks and its existing 4.5 GHz setting. The unrelated single-threaded
`bre-bench` remained pinned to CPU 0. These cores share L3 cache, so this
is a shared-host comparison, not an uncontended benchmark.

Each full-encoder invocation measured 501 calls after three warmups.
Both encoders use K=2^20 input records, N=2^21 output records, 128-bit
records, precomputed setup, and streaming output stores. Input generation
is outside the timed interval. We bracketed candidates with frozen-control
runs and reversed the g4/g8 order. Full-output checksums matched the oracles.

| GF256 candidate | Full-encoder median, ms | Evidence |
|---|---:|---|
| g2, fused | 4.890 | One invocation |
| g4, fused | 4.675--4.988 | Six invocations; median of medians 4.789 |
| g4, unfused | 4.843 | One invocation |
| g8, fused | 8.672--8.764 | Two invocations, reversed order |
| g8, unfused | 8.748 | One invocation |
| g16, unfused | 9.255 | One invocation |

The frozen RS256/GF65536/h2 control medians ranged from 4.734 to 4.808 ms
throughout the session. In the final alternating sequence, g4 medians were
4.942, 4.728, and 4.851 ms, bracketed by controls at 4.759--4.779 ms.
The candidate's variability prevents a reliable small-speedup claim.
The much larger g8 regression appears in both fused and unfused schedules;
its cause has not yet been profiled.

Standalone phase measurements used 301 calls:

| Encoder | Outer, ms | Routed inner, ms |
|---|---:|---:|
| Frozen RS256/GF65536/h2 | 2.484 | 2.307 |
| Grouped g4/GF256/h4, fused | 2.021 | 2.613 |
| Grouped g8/GF256/h4, fused | 6.082 | 2.677 |

These are separate diagnostic runs, not an exact decomposition of the full
median. The g4 outer saves about 19%, but the certified h4 inner costs about
13% more than the old h2 inner. Arithmetic savings alone therefore do not
establish a full-encoder win.

The unfused GF256 build also passes the small-size check matrix and both
independent K20 checks. Its SHA256 is
`fb2bd4190a2853608096377148fc556f5e1b42b721f3d7d4cb361e691ca60de1`.
Disabling fusion changes only the evaluation schedule, not the sampled map.

At this checkpoint, retain the old implementation as the performance control.
The next useful
optimization target is the certified h4 routed inner; g8 additionally needs
an outer memory/code-generation diagnosis. The 35% / 72.3975-bit certificate
remains a separate successful result. No new proof search is required for
exact-map implementation changes.

## Subsequent compact-key and Toeplitz comparison

The [outer continuation](HD_GROUPED_OUTER.md#invertible-toeplitz-column-mixers)
improves this checkpoint. Compact byte coefficients preserve the DHD map
and reduce its prepared outer setup from 14 to 3.5 MiB. Five g4 full-encoder
run medians span 4.443--4.584 ms, with median 4.492 ms. An independently
sampled invertible Toeplitz column mixer retains the same conservative
72.3975-bit certificate and uses a nine-product g4 circuit. With compact
keys, its five run medians span 4.369--4.503 ms, with median 4.387 ms.
The eight bracketed frozen-control medians span 4.745--4.792 ms.

The same shared-host caveat applies. That checkpoint takes about
4.4--4.5 ms at K20, rate 1/2, 35% distance; its extra benefit
over compact DHD is smaller than the coefficient-layout improvement. The
proof, negative screens, hashes, correctness checks, and timing protocol
are recorded in the continuation. No library or paper profile is changed.

The subsequent [exact-map inner optimization](HD_GROUPED_OUTER.md#exact-map-inner-optimization)
reaches about 4.0 ms with the same code, output checksum, and certificate.
It uses short cosets and a persistent striped state, and computes feedback
from the output using CA=0. The final matched inner timings are 2.040 ms
versus 2.540 ms for the preceding schedule. The full validation and timing
record is in that continuation; production remains unchanged.
