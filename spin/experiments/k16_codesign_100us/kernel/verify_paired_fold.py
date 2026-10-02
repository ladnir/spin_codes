"""Authenticate the final totals4 circuit on all 16 independent input lanes.

This is a portable symbolic correctness check, not a benchmark. Each integer
bit names an independent 128-bit SIMD lane; XOR propagates its full support.
The source check binds the symbolic network to both frozen C++ variants.
"""
import re
from pathlib import Path

HERE=Path(__file__).resolve().parent

def shuffle(a,b,imm):
    return (a[imm&3],a[(imm>>2)&3],b[(imm>>4)&3],b[(imm>>6)&3])

def xor(a,b):
    return tuple(x^y for x,y in zip(a,b))

a,b,c,d=(tuple(1<<(4*r+j) for j in range(4)) for r in range(4))
ab=xor(shuffle(a,b,0x44),shuffle(a,b,0xee))
cd=xor(shuffle(c,d,0x44),shuffle(c,d,0xee))
got=xor(shuffle(ab,cd,0x88),shuffle(ab,cd,0xdd))
expected=tuple(v[0]^v[1]^v[2]^v[3] for v in (a,b,c,d))
assert got==expected

expected_body="""
const auto ab=_mm512_xor_si512(_mm512_shuffle_i32x4(a,b,0x44),_mm512_shuffle_i32x4(a,b,0xee));
const auto cd=_mm512_xor_si512(_mm512_shuffle_i32x4(c,d,0x44),_mm512_shuffle_i32x4(c,d,0xee));
return _mm512_xor_si512(_mm512_shuffle_i32x4(ab,cd,0x88),_mm512_shuffle_i32x4(ab,cd,0xdd));
"""
compact=lambda text:re.sub(r'\s+','',text)
for name in ('K16PairedFold.cpp','K16Paired15Fold.cpp'):
    source=(HERE/name).read_text()
    start=source.index('static SPIN_FORCEINLINE __m512i totals4(')
    begin=source.index('{',start)+1
    end=source.index('}',begin)
    body=source[begin:end]
    assert compact(body)==compact(expected_body),name
    assert body.count('_mm512_shuffle_i32x4(')==6,name
    assert body.count('_mm512_xor_si512(')==3,name

print('PASS: totals4 on all 16 independent lanes; both frozen sources match six shuffles and three XORs.')
