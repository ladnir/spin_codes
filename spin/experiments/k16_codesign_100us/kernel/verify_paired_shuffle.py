"""Authenticate the fixed shuffle-only feedback and state-coordinate mapping."""
import re
from pathlib import Path
from generate_paired import GROUPS

permutation=(0,1,2,7,3,4,5,6,8,10,11,13,12,14,15,9)
rows=[sum((sum(p&m==m for m in group)&1)<<p for p in range(64)) for group in GROUPS]
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
f2=masked(shuffle(odd,odd,0xb1),shuffle(upper,quadratic,0xc8),0xfc)
f3=masked(shuffle(upper,both,0x8d),shuffle(quadratic,upper,0x09),0x3f)
assert low+total+f2+f3==tuple(rows[j] for j in permutation)

header=(Path(__file__).parent/'T64PairedShuffleMap.h').read_text()
matches=re.findall(r'const auto q(\d+)=_mm512_shuffle_i32x4\(state.v\[(\d)\],state.v\[\2\],0x([0-9a-f]+)\);',header)
assert len(matches)==12
for q,word,immediate in matches:
    index=permutation.index(int(q))
    assert int(word)==index//4 and int(immediate,16)==0x55*(index%4)
assert permutation[:4]==(0,1,2,7)

# Simulate every bit of every matrix under the two inverse basis mappings.
for row in range(16):
    for col in range(16):
        literal=[0]*16
        literal[row]=1<<col
        fast=[sum(((literal[permutation[r]]>>permutation[c])&1)<<c for c in range(16)) for r in range(16)]
        back=[0]*16
        for r in range(16):
            for c in range(16):back[permutation[r]]|=((fast[r]>>c)&1)<<permutation[c]
        assert back==literal
print('PASS: all feedback columns, all state broadcasts, C*A=0 and both inverse matrix-basis maps.')
