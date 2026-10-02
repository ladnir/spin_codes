"""Check every selected feedback row and all compact expansion columns."""
import re
from pathlib import Path
from generate_paired import GROUPS

groups=GROUPS[:10]+GROUPS[11:]
rows=[sum((sum(p&m==m for m in group)&1)<<p for p in range(64)) for group in groups]
assert len(rows)==15 and not any((a&b).bit_count()&1 for a in rows for b in rows)
columns=[sum(((row>>p)&1)<<j for j,row in enumerate(rows)) for p in range(64)]
header=(Path(__file__).parent/'T64Paired15Map.h').read_text()
record=re.search(r'columns\[64\]=\{([^}]+)\}',header).group(1)
assert [int(x,16) for x in record.split(',')]==columns

# The wide feedback network's selected moments, each represented by its full
# 64-bit input support. Only whole 128-bit lanes are used below.
z={m:sum(1<<p for p in range(64) if p&m==m) for m in range(64)}
low=(z[0],z[1],z[2],z[3])
total=(z[4],z[8],z[16],z[32])
odd=(z[5],z[9],z[17],z[33])
upper=(z[6],z[10],z[18],z[34])
quadratic=(z[12],z[20],z[24],z[36])
both=(z[40],z[40],z[48],z[48])
both_upper=upper[:2]+(both[2],upper[3])
left2=tuple((odd+both_upper)[j] for j in (1,6,3,5))
right2=tuple(quadratic[j] for j in (0,0,0,1))
v2=left2[:2]+tuple(left2[j]^right2[j] for j in (2,3))
left3=tuple((upper+odd)[j] for j in (6,3,0))+(0,)
right3=tuple((quadratic+both)[j] for j in (3,2,4,0))
v3=tuple(left3[j]^right3[j] for j in range(3))+(0,)
v1=total[1:]+total[:1]
physical=low+v1+v2+v3
permutation=(0,1,2,7,4,5,6,3,8,9,10,11,12,13,14,15)
padded=rows+[0]
assert physical==tuple(padded[j] for j in permutation)

# Authenticate the compact ANF coefficient assignment, then the fixed zeta
# expansion, independently of the generated per-quarter C++ expressions.
fast_groups=tuple((groups+[()] if isinstance(groups,list) else groups+((),))[j] for j in permutation)
d={h:tuple(sum((sum(lane&(m&3)==(m&3) for m in group if m>>2==h)&1)<<j
                   for j,group in enumerate(fast_groups)) for lane in range(4))
   for h in range(16)}
for packet in range(16):
    for lane in range(4):
        fast=0
        for h in range(16):
            if packet&h==h:fast^=d[h][lane]
        literal=sum(((fast>>j)&1)<<permutation[j] for j in range(16))
        assert literal==columns[4*packet+lane]
print('PASS: all 64 columns, all 15 feedback rows, compact/basis indexing, dummy zero and C*A=0.')
