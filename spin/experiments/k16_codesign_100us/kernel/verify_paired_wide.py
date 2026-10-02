"""Authenticate every physical input bit of the SIMD paired-feedback network."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
record=json.loads((ROOT/'research/workstreams/k16_codesign_100us/proof/disjoint-pairs-proposal.json').read_text())
groups=[(m,) for m in (0,1,2,4,8,16,32)]+[tuple((1<<a)|(1<<b) for a,b in g) for g in record['monomial_groups']]
rows=[sum((sum(p&m==m for m in group)&1)<<p for p in range(64)) for group in groups]
high={h:tuple(sum(1<<(4*p+lane) for p in range(16) if p&h==h) for lane in range(4)) for h in range(16)}
def xor(*vectors):
    out=[0]*4
    for v in vectors:
        out=[a^b for a,b in zip(out,v)]
    return tuple(out)
def shuffle(a,b,imm):return (a[imm&3],a[(imm>>2)&3],b[(imm>>4)&3],b[(imm>>6)&3])
def masked(a,b,mask):return tuple(a[j]^(b[j] if mask>>(2*j)&1 else 0) for j in range(4))
def transpose(a,b,c,d):
    ab0,ab1,cd0,cd1=shuffle(a,b,0x44),shuffle(a,b,0xee),shuffle(c,d,0x44),shuffle(c,d,0xee)
    return shuffle(ab0,cd0,0x88),shuffle(ab0,cd0,0xdd),shuffle(ab1,cd1,0x88),shuffle(ab1,cd1,0xdd)
def select(a,b,indices):return tuple((a+b)[j] for j in indices)
l0=masked(high[0],shuffle(high[0],high[0],0xf5),0x33)
low=masked(l0,shuffle(l0,l0,0xee),0x0f)
a,b,c,d=transpose(high[1],high[2],high[4],high[8])
odd,upper=xor(b,d),xor(c,d)
total=xor(a,c,odd)
quadratic=xor(*transpose(high[3],high[5],high[6],high[9]))
pair=xor(shuffle(high[10],high[12],0x44),shuffle(high[10],high[12],0xee))
both=xor(pair,shuffle(pair,pair,0xb1))
f0=low[:3]+(total[0],)
f1=(total[1],total[2],total[3],low[3])
f2=masked(select(odd,both,(1,6,0,3)),select(upper,quadratic,(0,0,2,4)),0xf0)
f3=xor(select(upper,odd,(1,6,3,0)),select(quadratic,both,(1,3,2,4)))
assert list(f0+f1+f2+f3)==rows
print('PASS: four-ZMM feedback equals all 16 literal paired rows on all 64 physical input bases.')
