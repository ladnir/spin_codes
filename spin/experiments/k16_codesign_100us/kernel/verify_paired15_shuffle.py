"""Authenticate shuffle-only s15 feedback, compact rows and the zero SIMD hole."""
import re
from pathlib import Path
from generate_paired import GROUPS

permutation=(0,1,2,7,3,4,5,6,8,15,10,12,11,13,14,9)
groups=GROUPS[:10]+GROUPS[11:]
rows=[sum((sum(p&m==m for m in group)&1)<<p for p in range(64)) for group in groups]+[0]
assert not any((a&b).bit_count()&1 for a in rows for b in rows)
z={m:sum(1<<p for p in range(64) if p&m==m) for m in range(64)}
low=(z[0],z[1],z[2],z[3])
total=(z[4],z[8],z[16],z[32])
odd=(z[5],z[9],z[17],z[33])
upper=(z[6],z[10],z[18],z[34])
quadratic=(z[12],z[20],z[24],z[36])
both=(z[40],z[40],z[48],z[48])
def shuffle(a,b,imm):return (a[imm&3],a[(imm>>2)&3],b[(imm>>4)&3],b[(imm>>6)&3])
def masked(a,b,mask):return tuple(a[j]^(b[j] if mask>>(2*j)&1 else 0) for j in range(4))
left2=shuffle(odd,odd,0xb1)
left2=left2[:1]+(0,)+left2[2:]
f2=masked(left2,shuffle(quadratic,quadratic,0xc0),0xf0)
f3=masked(shuffle(upper,both,0x8d),shuffle(quadratic,upper,0x09),0x3f)
assert low+total+f2+f3==tuple(rows[j] for j in permutation)

header=(Path(__file__).parent/'T64Paired15ShuffleMap.h').read_text()
matches=re.findall(r'const auto q(\d+)=_mm512_shuffle_i32x4\(state.v\[(\d)\],state.v\[\2\],0x([0-9a-f]+)\);',header)
assert len(matches)==11
for q,word,immediate in matches:
    q=int(q)
    assert q!=10
    compact=q-(q>10)
    index=permutation.index(compact)
    assert int(word)==index//4 and int(immediate,16)==0x55*(index%4)
assert permutation[:4]==(0,1,2,7) and permutation[9]==15
print('PASS: all 15 feedback rows, state broadcasts, zero SIMD hole, compact literal order and C*A=0.')
