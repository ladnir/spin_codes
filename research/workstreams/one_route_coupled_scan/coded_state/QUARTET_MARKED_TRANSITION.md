# Joint test of the hole-birth constraint

2026-10-30. This tests the marked local kernel from
[the component audit](QUARTET_COMPONENTS.md). The construction remains fixed.
The full45%/40-bit certificate is still open. The new joint bounds are valid
but positive; the marginal birth probabilities cannot close the calculation
by themselves.

## Regional test and its normalization

The test keeps state weight256 and44 deaths at each of16 updates. At every
step, the incoming marked set consists of the preceding44 hole births.
The region contains exactly2100 inactive labels. This is a conditional
regional continuation family, not an entire trajectory starting at zero.

Let P(x,z,y,w,h) be the128-quartet local product before coefficient
extraction, with markers for incoming selected positions, inherited hole
births, outgoing selected positions, outgoing deaths, and inactive labels.
Each quartet term also carries b^(l/2-c-(37/50)d), where l,c,d are its
inactive-label, covered-label, and death counts. The exponent accumulates
the averaged-reference excess G for this fixed weight profile.

The per-step incoming normalization is

    binom(512,256)binom(256,44).

The regional quota contributes the separate divisor binom(8192,2100).
The helper uses positive coefficient bounds at rational markers. It
multiplies the sixteen per-step bounds before extracting the regional
quota. The arithmetic is directed; a floating-point Newton solve only
proposes rational markers. An untight coefficient bound remains an upper
bound and is not presented as the exact moment.

With the new birth gate alone, directed log2 regional upper bounds are
1176.039 at log2(b)=4 and2225.891 at log2(b)=8. Adding the earlier
canonical quartet-erasure constraints reduces them to1086.490 and2217.868.
These are not directly comparable to earlier bounds that summed every
death count. None is a probability of a low-distance codeword.

## A stronger deletion rule

Fix a maximizing hole support H and a quartet at update t+1. Let E be
its e>0 newly created holes from update t. The birth gate already implies
that all four parents are holes. Let f count the members of E that are
selected, rather than holes, at update t+1. Let k count inactive inputs
in this quartet whose outgoing coordinates remain holes.

Delete E from H. No own input at update t was uncovered by a member of E,
because that coordinate was not a hole at update t-1. Deletion therefore
loses exactly k uncovered inputs, all at update t+1. It removes e holes
and reduces the number of hole births by f. Consequently the change in
the hole score J is

    lambda*e+tau*f-k.

Maximality requires

    k >= lambda*e+tau*f.                                  (1)

No strict inequality or tie-breaking assumption is needed. This condition
can be combined with the existing canonical rules for the selected
maximizer. At lambda=1/4,tau=37/50 it is implemented as100k>=25e+74f.
It is a condition on maximizing geometric witnesses, not arbitrary codewords.

The original kernel has265 aggregated local terms. The birth budget alone
leaves242; the old canonical rules alone leave160; both leave153.
The fully filtered directed regional upper bounds are

    log2(b)=4: 585065420629073/549755813888 = about1064.228,
    log2(b)=8: 1212336584967415/549755813888 = about2205.227.

The stronger local condition helps modestly but does not resolve the gap.

## Count prefixes, without a score tilt

The helper also introduces a covered-label marker v_C<=1 and bounds the
sum of histories with total regional coverage C at most a specified cap.
The other four per-update statistics and the exact regional quota remain
the same. The directed log2 count bounds are:

| Coverage cap | Upper bound |
|---:|---:|
| 530 | 1066.425 |
| 472 | 819.771 |
| 432 | 640.374 |

These are upper bounds on expected COUNTS of continuations. The cap432 is
motivated by the true score of a repeated constant profile: if all256
regions had A=4096,D=704,C=432, their total score would be18186.24, slightly
above18171. Such a repeated bulk profile does not by itself satisfy the
empty initial state and terminal birth gate. No full violating history is
asserted. Cap472 instead probes the looser averaged-reference threshold;
it must not be confused with the true constant-profile score threshold.

The remaining positive bound might include coefficient-extraction slack
and extensive multiplicity of overlapping histories. The calculations do
not prove a large first moment for actual violations. An exact subcount
or sharper extraction is needed to distinguish these possibilities.

## A bounded check of the row-set assumption

The sufficient geometric inequality must cover every inactive-row set,
not just the fixed set used in the earlier min-cut diagnostics. The C++
helper now permits up to three alternating exact support optimizations
and row-set best responses. For a fixed support S, its best row set of
size2100 consists of the2100 smallest per-row coverage counts.

On the existing slope0.25 geometry, the first response keeps precisely the
original set: C remains1238, the score remains3988.44, and no rows change.
The run stops immediately rather than repeating an unchanged state.
This is only a local best-response check. It does not bound the worst
row set or the ensemble. Its support also lies below the remaining
field/activity cutoff, so it is not a witness in the unresolved code class.

## Verification and next step

`quartet_marked_transition.py` verifies all four local-filter combinations
against independent colored-mask enumeration:20736 cases each. A separate
two-quartet calculation enumerates all35 partitions of eight coordinates
and the outgoing/label masks. The four weighted continuation counts are
162/245,279/490,237/490,447/980; the normalized polynomial matches exactly.
The numbers are expected counts, not automatically probabilities.

The exact optimal-triple toy test also verifies(1), including nonvacuous
hole births. The full-size row-response run retains flow feasibility and
an independent support recount. All work was serial, remote, and bounded.
No raw data were saved; temporary C++ builds were removed. The frozen
prior helper, construction, paper, and library remain unchanged.

Do not yet allocate or run the full131841-state transfer. First audit a
fixed-type subcount at the coverage432 profile, or extract more of that
profile exactly. The purpose is to determine whether the remaining
positive bound reflects large history multiplicity or primarily the
coefficient upper envelope. A positive subcount would diagnose a counting
method, not establish a bad codeword or justify changing the construction.
Full r=2100 closure and all other unresolved row counts remain required.
