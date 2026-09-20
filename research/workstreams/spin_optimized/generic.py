"""Generate an opt-in XOR-element transpose without modifying SIMD hot paths."""
from pathlib import Path
import ast
import re
import sys
from k16_map import header as k16_header

p=Path(sys.argv[1]);here=Path(__file__).resolve().parent
root=here.parents[1]

def edit(s,a,b):
    assert s.count(a)==1,a
    return s.replace(a,b)

def body(text,signature):
    start=text.index(signature);begin=text.index('{',start);end=begin+1;depth=1
    while depth:
        depth+=(text[end]=='{')-(text[end]=='}');end+=1
    return text[begin+1:end-1]

h=(p/'Spin.h').read_text()
assert 'MessageLength' in h and 'forwardBits(' not in h
h=edit(h,'class Spin {','class GenericTranspose;\nclass Spin {')
h=edit(h,'    void validateSetup() const;', '''    // Owns a generic-element copy of this realized map; works after compact().
    GenericTranspose generic() const;
    void validateSetup() const;''')
(p/'Spin.h').write_text(h,newline='\n')
(p/'GenericSpin.h').write_text((here/'GenericSpin.h').read_text(),newline='\n')

# Reuse existing straight-line circuits, replacing only their element operations.
asym=(p/'AsymmetricMap.h').read_text();small=k16_header()
circuits=['#pragma once','namespace bare_spin::generic_detail {',
          'template<class E> inline E vx(const E& a,const E& b) {return E(a^b); }']
for name,t,s,text in [('Map128',128,19,asym),('Map64',64,12,small)]:
    finish=body(text,'void finish(').replace('__m128i','E')
    emit=body(text,'void emitShared(').replace('.mData','').replace('block(', 'E(')
    circuits += [f'struct {name} {{',f'static constexpr unsigned T={t},S={s};',
                 'template<class E> static inline void finish(const E* z,E* out) {'+finish+'}',
                 'template<class E,class Emit> static inline void emitShared(const E* in,E* raw,',
                 'const E* state,std::size_t base,Emit& emit) {'+emit+'}', '};']

# Check each emitted BCH linear form against the generator, not just C++ tests.
source=(p/'generated/BchCircuit.cpp').read_text()
forms={};lines=[];outputs={}
def form(expr):
    def walk(n):
        if isinstance(n,ast.Name): return forms[n.id]
        if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='vx' and len(n.args)==2:
            return walk(n.args[0])^walk(n.args[1])
        raise ValueError(ast.dump(n))
    return walk(ast.parse(expr,mode='eval').body)

for line in source.splitlines():
    m=re.fullmatch(r'const auto v(\d+)=_mm256_set_m128i\(b\[(\d+)\].mData,a\[\2\].mData\);',line)
    if m:
        index,col=map(int,m.groups());assert index==col
        forms[f'v{index}']=1<<col;lines.append(f'const auto v{index}=a[{col}];');continue
    m=re.fullmatch(r'const auto ([vo]\d+)=(vx\(.*\)|v\d+);',line)
    if m:
        name,expr=m.groups();forms[name]=form(expr);lines.append(line);continue
    m=re.fullmatch(r'x\[(\d+)\]=block\(_mm256_castsi256_si128\((o\d+)\)\);',line)
    if m:
        index=int(m[1]);outputs[index]=forms[m[2]];lines.append(f'out[{index}]={m[2]};')
words=[int(x,16) for x in re.findall(r'0x([0-9a-f]+)ULL',(p/'generated/BchCircuit.h').read_text())]
targets=[sum(words[4*i+j]<<(64*j) for j in range(4)) for i in range(128)]
assert len(outputs)==128 and [outputs[i] for i in range(128)]==targets
circuits+=['template<class E> inline void bch(const E* a,E* out) {']+lines+['}','}']
(p/'GenericCircuits.h').write_text('\n'.join(circuits)+'\n',newline='\n')

# Recover natural row-major routing from whichever schedule compact() retained.
# This one-time conversion does not enlarge Spin or alter any encoding method.
factory='''#include "GenericSpin.h"
#include <cstring>
namespace bare_spin {
static u32 read24(const u8* p) {u32 x;std::memcpy(&x,p,4);return x&0xffffff;}
static u32 natural(u32 x) {return (x&~1023U)|((x&3U)<<8)|((x>>2)&255U);}
GenericTranspose Spin::generic() const {
    std::vector<u32> route(codeBlocks());
    if(!mRoute.empty()) route=mRoute;
'''
for member in ('mSmallRoute32','mSmallRoute24','mK18Route32','mK18Route24','mK20Route32','mK20Route24','mRangeRoute32'):
    if member in h:
        read=f'{member}[i]' if member.endswith('32') else f'read24({member}.data()+3*i)'
        factory+=f'    else if(!{member}.empty()) for(std::size_t i=0;i<codeBlocks();++i) route[i]=natural({read});\n'
factory+='''    else {
        const bool packed=mRetainedLayout==Layout::Packed24;
        const auto tile=tileBlocks();
        for(std::size_t i=0;i<codeBlocks();++i) {
            const auto slot=packed?read24(mSlots24.data()+3*i):mSlots32[i];
            const auto offset=packed?read24(mOffsets24.data()+3*slot):mOffsets32[slot];
            route[i]=u32(slot&~(tile-1))+(mBchBackend==BchBackend::Avx512?natural(offset):offset);
        }
    }
    return GenericTranspose(mConfig,mK,std::move(route),mFieldRows);
}
}
'''
(p/'GenericSpin.cpp').write_text(factory,newline='\n')
print('Generated XOR-element fallback; BCH forms verified exactly.')
