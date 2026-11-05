# XOR2: one single-row group and sixteen two-row groups

The complete directed calculation for this 17-group cohort gives a
70.3541-bit bound at 35% distance. It covers all 81 ordered source matching
classes and all constituent-pair assignments. The 256/384-bit replays agree
on all 324 cell ceilings; independent 512- and 1024-bit rounded-cell sums
agree. No encoder, setup distribution, or production profile changed.

The candidate is the eight-RS[32,16] outer with independent diagonal--XOR--
diagonal column randomization and the paired-and-swapped route described in
[the main XOR2 record](HD_XOR_PAIRED_GATE.md). Its matched full encoding time
remains about 3.46 ms. The proof here concerns ideal independent setup
randomness, not a benchmark seed.

## Class and count

Fix the 17 active outer groups. One group has one nonzero RS constituent;
each of the other sixteen has two. Every constituent row set remains fixed
across all 32 RS columns. Write `(a,b)` for the numbers of full source pairs
under the two alternating route matchings. There are 81 possible ordered
pairs; their counts come from the exact cyclic source-count helper.

Let `mu` be the largest multiplicity of a constituent pair among the sixteen
two-row groups. The four disjoint bins are 1--8, 9--14, 15, and 16. Their
exact counts sum to `28^16`. Each cell is multiplied once by its source
count, 17 choices of the single-row group, 8 choices of that group's row,
and its constituent-pair assignment count. No factor is repeated per RS
column.

In one column, let `one` indicate whether the single-row constituent is
nonzero, and let `s,d` count the two-row groups with one and two nonzero
inputs. Their source-support multiplicity is

    C(1,one) C(16,s) C(16-s,d) 2^s.

The existing resource bound applies with `single=one+s` and `double=d`.
The extra single-row group only contributes a missing-region count; it
cannot increase the multiplicity of a cancelled two-row edge. Consequently
the `mu` constraints remain valid. An independent audit and an exhaustive
one-singleton/one-pair test check this interface, including all labelled
row positions and cancellation branches.

There are 33 active RS words, each with minimum support 17. A positive
support polynomial retains total source support at least `33*17=561`
across its 32-column power and divides by `255^(33*16)=255^528`.
The region operators remain in chronological order. Cancellation branches
are averaged before the common potential quotient. The final bound uses
one global projective factor `1/255` and one endpoint-potential factor.

## Witness refinement and validation

The inherited witnesses were tuned for seventeen two-row groups. Applied
to this different cohort, their worst tested corner gave only 7.20 bits.
A bounded 32-evaluation search changed only the Chernoff mark and common
potential for that corner. The new literal is

    z = 4487948476156805 / 4503599627370496,
    log2(potential) = -8814329 / 131072.

That corner improves to 262.55 bits. Testing the new literal alongside
each inherited literal, without further fitting, gives 70.3765 bits on
the four source corners and 70.3541 bits on the complete 324-cell union.
The worst remaining cell is `(a,b)=(0,8), mu=16`; it retains its inherited
witness. Different cells are disjoint message classes, so choosing a fixed
witness per class introduces no witness-bank multiplicity.

The new directed helper recomputes all counts and bounds. It reads only
the proposed rational witnesses, not saved floating-point numerical values.
Its verifier requires the complete 81-by-4 partition, compares every
256/384-bit rational ceiling, and sums rounded cell bounds independently
at 512 bits. Both replays have finished. Their raw aggregate ceiling is
`-2814165921867/40000000000`; the rounded-cell union ceiling is
`-35177074023337/500000000000`, or -70.354148046674. A separate 1024-bit
sum reproduces the latter ceiling. This covers the stated cohort, not
all nonzero messages of the candidate.

Twenty-six related unit tests pass, including exact toy support-polynomial
enumeration, reduction to the homogeneous two-row polynomial, fixed-source
chronology, assignment counts, and the extra-singleton resource constraints.
The independent directed-wrapper audit found no mathematical discrepancy.

## Reproduction and frozen inputs

Run from `theory`:

```text
python hd_xor_sparse_mixed_refine.py --output hd_xor_sparse_mixed_refine32_corners_float.json
python hd_xor_sparse_mixed_refine.py --bank-from hd_xor_sparse_mixed_refine32_corners_float.json --all-sources --output hd_xor_sparse_mixed_refine32_all_sources_float.json
python hd_xor_sparse_mixed_certificate.py --witnesses hd_xor_sparse_mixed_refine32_all_sources_float.json --precision 256 --output hd_xor_sparse_mixed_full_receipt256.json
python hd_xor_sparse_mixed_certificate.py --witnesses hd_xor_sparse_mixed_refine32_all_sources_float.json --precision 384 --output hd_xor_sparse_mixed_full_receipt384.json
python hd_xor_sparse_mixed_certificate.py --witnesses hd_xor_sparse_mixed_refine32_all_sources_float.json --verify hd_xor_sparse_mixed_full_receipt256.json hd_xor_sparse_mixed_full_receipt384.json --output hd_xor_sparse_mixed_full_dual_verified.json
```

Scripts refuse to overwrite an existing artifact. Choose new output names
for independent replays; preserve the frozen evidence.

| Input | SHA256 |
|---|---|
| `hd_xor_sparse_mixed_refine.py` | `11c90f5e0cce144f8ec7ea6743f16657b44cd995a3ecfcf5120567cd2d07f97e` |
| `hd_xor_sparse_mixed_refine32_all_sources_float.json` | `10a7116e0376f203fee0256f524ce1c3614cbc358f872dba0afb863d9d831b75` |
| `hd_xor_sparse_mixed_certificate.py` | `d3ff58009fcbcb4ab1618b96f816f176bfeb2642d4126e9fb2ab7484a48164d0` |
| `hd_xor_sparse_mixed_full_receipt256.json` | `fe89c05e39bcb6f48cda7a687a20c2006736f09880adff5b11159da4ad0ec887` |
| `hd_xor_sparse_mixed_full_receipt384.json` | `7f0c8c43bcc9ac63d0f9e312301fa354c14e3250400db4f4c7c33b7ec56075e4` |
| `hd_xor_sparse_mixed_full_dual_verified.json` | `e78f2b3320fb073187915608ce64725917463b5820f1db12f6c7c95bb1631366` |

Next for this retained candidate: address other sparse/dense cohort mixtures
using a shared count bound rather than another isolated corner fit.
