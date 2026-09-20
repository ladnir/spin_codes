"""Four-row BCH forward circuit and ISA-isolated forward routing/inner."""
from pathlib import Path
import re
import sys
p=Path(sys.argv[1])
def edit(s,a,b):
    assert s.count(a)==1,a
    return s.replace(a,b)

original=(p/'generated/BchForward.cpp').read_text()
body=original[original.index('const auto o0='):original.rindex('\n}\n}')]
body,n=re.subn(r'const auto o(\d+)=_mm256_set_m128i\(b\[\1\].mData,a\[\1\].mData\);',
    r'const auto o\1=join(a[\1],a[128+\1],a[256+\1],a[384+\1]);',body)
assert n==128
body,n=re.subn(r'x\[(\d+)\]=block\(_mm256_castsi256_si128\(z\1\)\);',
    r'_mm512_storeu_si512(x+4*\1,z\1);',body)
assert n==256
body,n=re.subn(r'\ny\[(\d+)\]=block\(_mm256_extracti128_si256\(z\1,1\)\);','',body)
assert n==256
(p/'generated/BchForward512.cpp').write_text('''// Original forward XOR DAG, four independent outer rows per vector.
#include "Spin.h"
namespace bare_spin {
static inline __m512i vx(__m512i a,__m512i b) {return _mm512_xor_si512(a,b);}
static inline __m512i join(block a,block b,block c,block d) {
    auto x=_mm512_castsi128_si512(a.mData);
    x=_mm512_inserti32x4(x,b.mData,1);x=_mm512_inserti32x4(x,c.mData,2);
    return _mm512_inserti32x4(x,d.mData,3);
}
void bchForward4(const block* a,block* x) {
'''+body+'\n}\n}\n',newline='\n')
h=(p/'Spin.h').read_text()
h=edit(h,'    template<class Map, bool Packed> void runForward',
       '    template<class Map, bool Packed> void runForwardFour(const block*,block*,Workspace&) const;\n'
       '    template<class Map, bool Packed> void runForward')
(p/'Spin.h').write_text(h,newline='\n')
s=(p/'Spin.cpp').read_text()
# Direct transpose paths do not need these tables, but the forward path does.
for name in ('mBchOffsets24','mBchOffsets32'):
    s=s.replace(f'decltype({name})().swap({name});', '/* Retained for four-row forward routing. */')
start=s.index('template<class Map,bool Packed> void Spin::runForward(')
end=s.index('void Spin::forwardUnchecked',start)
kernel=s[start:end].replace('Spin::runForward(', 'Spin::runForwardFour(')
kernel=kernel.replace('mOffsets24','mBchOffsets24').replace('mOffsets32','mBchOffsets32')
kernel=edit(kernel,'''for(std::size_t j=0;j<tileSize;j+=512)
            bchForward2(in+base/2+j/2,in+base/2+j/2+128,tile+j,tile+j+256);''',
    '''for(std::size_t j=0;j<tileSize;j+=1024)
            bchForward4(in+base/2+j/2,tile+j);''')
inst=''
for name in ('Map128S19','Map64S12','Map64S12R2'):
    for packed in ('true','false'):
        inst+=f'template void Spin::runForwardFour<{name},{packed}>(const block*,block*,Workspace&) const;\n'
(p/'ForwardFast.cpp').write_text('''#include "Spin.h"
#include "Inner.h"
#include <cstring>
namespace bare_spin {
void bchForward4(const block*,block*);
static OC_FORCEINLINE u32 unpack(const u8* p) {u32 v;std::memcpy(&v,p,4);return v&0xffffff;}
'''+kernel+inst+'}\n',newline='\n')
s=edit(s,'void Spin::forwardUnchecked(const block* in,block* out,Workspace& w,Layout layout) const {',
'''void Spin::forwardUnchecked(const block* in,block* out,Workspace& w,Layout layout) const {
#if SPIN_BCH_AVX512
    if(mBchBackend==BchBackend::Avx512) {
#define FAST_FORWARD(Enum,Map) case Configuration::Enum: if(layout==Layout::Packed24) runForwardFour<Map,true>(in,out,w); else runForwardFour<Map,false>(in,out,w); return
        switch(mConfig) {
            FAST_FORWARD(T128S19,Map128S19);FAST_FORWARD(T64S12,Map64S12);FAST_FORWARD(T64S12R2,Map64S12R2);
            default:break;
        }
#undef FAST_FORWARD
    }
#endif''')
(p/'Spin.cpp').write_text(s,newline='\n')
