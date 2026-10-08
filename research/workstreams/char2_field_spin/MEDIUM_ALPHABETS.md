# Medium alphabets at 10% distance

The current evidence favors GF(2^8) as the performance baseline. Full
GF(2^32) arithmetic works in four lanes per 128-bit record, but its measured
coefficient pass is about four times slower at K16/K18. Its additional
randomness barely changes the sampled classes blocking higher rate or
smaller state in the existing construction. This is a first comparison,
not a claim that every GF32 construction is slower.

The objective remains higher rate at 10% symbol distance, one regional
routing stage, and competitive full encoding cost. Another half-rate/10%
scalar extension of binary SPIN is not the objective. Production and the
completed GF256 certificates are unchanged.

## Alphabet versus application record

Every timed record contains 128 payload bits. The GF256 kernel treats it
as sixteen field coordinates; the GF32 kernels treat it as four. At each
record position one sampled coefficient is shared by all its coordinates.
Thus all lanes use the same scalar code, rather than independently sampled
codes. A nonzero record word has a nonzero component word, whose support
is contained in the record support. A scalar distance guarantee therefore
gives the same record distance without a union over lanes.

At relative distance .10, the asymptotic GV rate benchmark is approximately
84.14%, 88.53%, and 89.63% for alphabets of 8, 32, and 128 bits. These are
evaluations of 1-H_q(.10), an existence bound, not rate upper limits or
finite certificates for this construction. See
[Guruswami's GV notes](https://www.cs.cmu.edu/~venkatg/teaching/codingtheory/notes/notes2.pdf).
Moving from 32 to 128 bits buys little in this benchmark. The difference
between 8 and 32 bits is larger, but must pay for the arithmetic.

## Implemented arithmetic controls

The new [header](../../../spin/experiments/char2_field_spin/MediumField32.h)
and [frontend](../../../spin/experiments/char2_field_spin/medium_field.cpp)
implement four isolated diagonal passes. They do not implement a new full
encoder. Setup chooses a nonzero coefficient per record; the timed call
reads records and coefficients, multiplies, and writes records.

| Control | Coefficient set | Implementation |
|---|---|---|
| `gf8` | All 255 nonzero GF256 scalars | Bytewise GFNI |
| `poly32` | All nonzero GF32 scalars | Four VPCLMUL products per 512-bit vector, then reduction |
| `tower32` | All nonzero GF32 scalars | Three variable GFNI products plus two fixed reduction products and shuffles |
| `poly32-small` | 255 nonzero degree-at-most-seven polynomials in GF32 | Two affine GFNI maps, shifts, and reduction |

The native polynomial is X^32+X^7+X^3+X^2+1. Its irreducibility is checked
by Rabin's criterion. The tower uses AES GF256 and basis [1,u,v,uv], with
u^2+u+delta=0, v^2+v+delta*u=0, and delta=0x20. The trace of delta to GF2 is
one, proving both quadratic extensions irreducible by trace transitivity.

In the tower, multiplication by an embedded GF256 coefficient is exactly
bytewise GFNI. Consequently the existing fixed outer and Cauchy maps can
remain bytewise without conversions. That statement does not hold for
the native polynomial representation: it needs an explicit embedding or
different fixed maps. The two GF32 kernels use different representations,
so their raw output bytes need not agree.

The small-polynomial control is not uniform GF32 randomness. A fresh
pivot drawn from its nonzero coefficient set has the elementary atom
bound 1/255, not 1/(2^32-1). Its speed cannot be combined with the full-GF32
proof screen.

## Serial native measurements

Peach Ryzen 7950X, CPU2, GCC15.2, Release -O3, AVX512/GFNI/VPCLMUL,
`-mtune=znver4 -fno-lto`. All benchmarks ran serially, after native
correctness checks. Each launch uses three warmups and 101 samples. At
4096 records a sample averages 31 calls; at K16 it averages three; larger
sizes use one. Allocation, coefficient generation, and checksum are
outside the timed call. Ordinary vector storage is used; no explicit
huge-page request is made. Coefficients are precomputed, not generated
on the timed path.

The ranges below contain the medians from seeds 113 and 7919, not confidence
intervals. Each pass transfers 32 payload bytes per record plus one or four
coefficient bytes. The K20 working set is large enough for cache and memory
effects to obscure the arithmetic ratio.

| 128-bit records | GF256 | GF32 VPCLMUL | GF32 tower | Small polynomial coefficients |
|---:|---:|---:|---:|---:|
| 4096 | 0.969--0.975 us | 5.257--5.258 us | 5.069--5.350 us | 2.141--2.164 us |
| 2^16 | 21.126--21.340 us | 84.358--84.461 us | 84.672--85.243 us | 35.450--35.630 us |
| 2^18 | 84.517--85.239 us | 337.660--337.671 us | 340.115--341.327 us | 145.241--145.562 us |
| 2^20 | 546.450--602.255 us | 1428.417--1431.522 us | 1441.070--1442.944 us | 723.030--772.161 us |

This is not an encoder timing table. In particular, it does not include
the outer, persistent-state recurrence, or route. A redesign that removes
enough of those costs could still favor GF32. These kernels have not
exhausted representation-specific optimization.

Release executable SHA256:
`c785f425c57904f729bcc5f8f1784db981afc830edde1f86281dc96f6da59ca9`.
The timed loop retains out-of-line kernel calls; output checksums consume
the results. Release and ASAN/UBSAN checks pass, as does the registered
arithmetic CTest. Checks cover all 1024 binary basis products in each
representation, independent scalar/quartic oracles, random SIMD products,
all byte coefficients, embedded GF256 action, tails, guards, and aliases.

```sh
cmake -S spin/experiments/char2_field_spin -B build-medium -DCMAKE_BUILD_TYPE=Release
cmake --build build-medium --target spin_char2_medium_field
build-medium/spin_char2_medium_field --check
taskset -c 2 build-medium/spin_char2_medium_field 262144 101 all 1 113
```

## Proof comparison and decision

The [field-parameterized screen](MEDIUM_FIELD_PROOF.md) correctly extends
the fixed-map rank laws to genuine uniform nonzero GF32 diagonals. Seven
exact/regression tests pass. It compares the same geometry and also
evaluates the GF256 witness under GF32, avoiding an optimizer-only comparison.

At K20 and 10% target distance, GRS[64,48] with state four or eight still
has negative sampled middle-class bounds. Restricting the certified
GRS[128,87] map to state six also remains unresolved: at input density
.10 the sampled margin is -26089 bits for GF256 and -26087 for GF32.
The one-row margin can rise substantially, but that is not the class
blocking the whole proof. These are unsuccessful upper bounds, not
counterexamples to the new ensembles.

The difficult classes prefer rank tilts roughly 7--14, already permitted
by GF256. The current cover charges the rank of constraints and forgets
much of the state history; the larger coefficient set hardly helps there.
Previously proved regional-window dimension obstructions also hold over
both fields, independently of coefficient entropy.

Keep GF256 as the implementation baseline. The next useful experiment
should improve the middle-class analysis or change the source/packet
geometry before paying for full GF32 coefficients. Keep the GF32 tower
as a validated arithmetic option for a design that saves fixed-map or
routing work. Do not start a full GF32 encoder or certificate campaign
solely because its singleton bound improves.

A specific next experiment is to retain a bound on the state-form rank
between packets, using the existing positive rank-cap machinery rather
than independently maximizing Gamma over all ranks at every step. If
selected zero equations leave rank H, a packet with w active inputs and
i selected zeros on those inputs has successor rank at most
min(s,H+w-i). The present two-state cover forgets that restriction.
Generalize the checked eight-state transfer to state six, retain the joint
return cover, and first screen GRS[128,86] at densities .10 and .14.
This could tighten the proof of a cheaper GF256 inner without changing
the route or paying for GF32 arithmetic. Its continuation-cone argument
and point screens still need verification; no improvement is claimed yet.
