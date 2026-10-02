"""Fixed s14 map: seven affine rows and seven disjoint quadratic pairs."""
from collections import Counter
import hashlib
from pathlib import Path

import screen_t128 as screen

GROUPS=(((0,3),(1,2)),((0,4),(1,3)),((0,5),(1,4)),
        ((1,5),(2,4)),((0,1),(2,5)),((0,2),(3,4)),((2,3),(4,5)))


def prepare():
    rows=[(1<<64)-1]
    rows += [sum(1<<x for x in range(64) if x>>i&1) for i in range(6)]
    for group in GROUPS:
        row=0
        for i,j in group:
            row ^= sum(1<<x for x in range(64) if (x>>i&1) and (x>>j&1))
        rows.append(row)
    rows=tuple(rows)
    columns=tuple(sum(((r>>x)&1)<<i for i,r in enumerate(rows)) for x in range(64))
    maps=screen.q1.kernel_t64.s16_maps
    if (maps.binary_rank(rows)!=14 or maps.binary_rank(columns)!=14
            or any((a&b).bit_count()&1 for a in rows for b in rows)
            or any(maps.binary_rank(columns[x:x+4])!=4 for x in range(0,64,4))):
        raise ArithmeticError('rank, CA=0, or packet restriction failed')
    images=tuple(maps.images_from_rows(rows));spectrum=Counter(x.bit_count() for x in images)
    expected={0:1,24:1072,28:3840,32:6558,36:3840,40:1072,64:1}
    if spectrum!=expected or len(set(images))!=16384:raise ArithmeticError('fresh spectrum mismatch')
    physical=screen.q1.kernel_t64.kernel_maps.prepare_maps(images,columns,bits=14,
        distribution='uniform_gl',birth_density='capped')
    source=Path(__file__).resolve()
    record=dict(schema='seven-pair-t64-s14-fixed-map-1',physical_t=64,state_bits=14,
        monomial_groups=[[list(pair) for pair in group] for group in GROUPS],
        expansion_rows_hex=list(map(hex,rows)),feedback_columns=list(columns),
        spectrum={str(w):n for w,n in sorted(spectrum.items())},full_state_census=True,
        map_sha256=physical['map_sha256'],state_update='independent uniform GL14 per physical step',
        source_sha256={str(source):hashlib.sha256(source.read_bytes()).hexdigest()})
    return screen.q1.kernel_t64.wrap(physical),record
