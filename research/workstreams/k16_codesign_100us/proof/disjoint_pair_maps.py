"""Fixed t64/s16 map with disjoint quadratic pairs and C=A^T.

Rows are the constant, six coordinate functions, and the nine quadratic
functions below. A pair means their binary sum. This is one fixed map;
the search that suggested it supplies no numerical proof premise.
"""
from collections import Counter
import hashlib
from pathlib import Path

import screen_t128 as screen

GROUPS = (((0,1),), ((0,3),), ((4,5),), ((0,2),(1,4)),
          ((0,5),(2,3)), ((1,3),(2,4)), ((0,4),(2,5)),
          ((1,5),(3,4)), ((1,2),(3,5)))


def prepare():
    rows = [(1<<64)-1]
    rows += [sum(1<<x for x in range(64) if x>>i&1) for i in range(6)]
    for group in GROUPS:
        row=0
        for i,j in group:
            row ^= sum(1<<x for x in range(64) if (x>>i&1) and (x>>j&1))
        rows.append(row)
    rows=tuple(rows)
    columns=tuple(sum(((row>>x)&1)<<i for i,row in enumerate(rows)) for x in range(64))
    maps=screen.q1.kernel_t64.s16_maps
    if (maps.binary_rank(rows)!=16 or maps.binary_rank(columns)!=16
            or any((a&b).bit_count()&1 for a in rows for b in rows)
            or any(maps.binary_rank(columns[x:x+4])!=4 for x in range(0,64,4))):
        raise ArithmeticError('rank, CA=0, or packet restriction check failed')
    images=tuple(maps.images_from_rows(rows))
    spectrum=dict(sorted(Counter(x.bit_count() for x in images).items()))
    expected={0:1,16:20,24:4640,28:13824,32:28566,36:13824,40:4640,48:20,64:1}
    if spectrum != expected or len(set(images))!=65536:
        raise ArithmeticError('fresh full state spectrum disagrees with the declared candidate')
    physical=screen.q1.kernel_t64.kernel_maps.prepare_maps(images,columns,bits=16,
        distribution='uniform_gl',birth_density='capped')
    source=Path(__file__).resolve()
    record=dict(schema='disjoint-pair-t64-s16-fixed-map-1',t=64,s=16,
        monomial_groups=[[list(pair) for pair in group] for group in GROUPS],
        expansion_rows_hex=list(map(hex,rows)),feedback_columns=list(columns),
        feedback_definition='C=A^T',feedback_times_expansion_zero=True,
        expansion_rank=16,feedback_rank=16,packet_ranks=[4]*16,
        state_count=65536,full_state_census=True,
        expansion_spectrum={str(w):n for w,n in spectrum.items()},
        feedback_kernel_spectrum=list(map(str,maps.dual_spectrum(spectrum,64,16))),
        map_sha256=physical['map_sha256'],minimum_expansion_weight=16,
        distribution='independent transitive linear family at each physical step',
        macro=dict(physical_steps=2,physical_step_bits=64,macro_step_bits=128,
            state_continuity='retained_between_halves'),
        source=dict(path=str(source),sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
    return screen.q1.kernel_t64.wrap(physical),record
