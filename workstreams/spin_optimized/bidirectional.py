"""Port the fast transpose to an explicit upstream snapshot; leave upstream untouched."""
from pathlib import Path
import hashlib
import json
import sys

source, common, out = map(Path,sys.argv[1:])
out.mkdir(parents=True,exist_ok=True)
(out/'generated').mkdir(exist_ok=True)
def edit(text,old,new):
    if text.count(old)!=1: raise RuntimeError(f'upstream anchor changed: {old}')
    return text.replace(old,new)
def put(name,text): (out/name).write_text(text,newline='\n')
names=('Spin.cpp','Spin.h','Block.h','Inner.h','WorkspaceRouting.h','generated/SelectedMaps.h',
       'generated/BchCircuit.cpp','generated/BchCircuit.h','generated/BchForward.cpp','generated/BchForward.h')
for name in names: put(name,(source/name).read_text())
put('generated/BchAvx512.cpp',(common/'generated/BchAvx512.cpp').read_text())
put('SOURCE_HASHES.json',json.dumps({n:hashlib.sha256((source/n).read_bytes()).hexdigest() for n in names},indent=2)+'\n')
h=(source/'Spin.h').read_text()
h=edit(h,'enum class Layout','enum class BchBackend { Auto, Avx2, Avx512 };\nbool bchAvx512Available() noexcept;\nenum class Layout')
h=edit(h,'unsigned tileRows=0);','unsigned tileRows=0, BchBackend backend=BchBackend::Auto);')
h=edit(h,'    void validateSetup() const;',
 '    void encodeInplace(block*,std::size_t,Workspace&,Layout=Layout::Packed24) const;\n'
 '    BchBackend bchBackend() const noexcept {return mBchBackend;}\n    void validateSetup() const;')
h=edit(h,'    bool mCompacted=false;',
 '    BchBackend mBchBackend=BchBackend::Avx2;\n'
 '    std::vector<u8> mBchOffsets24;\n    std::vector<u32> mBchOffsets32;\n    bool mCompacted=false;')
h=edit(h,'    template<class Map> void setupInner',
 '    template<bool Packed> void runFour(const block*,block*,Workspace&) const;\n    template<class Map> void setupInner')
put('Spin.h',h)
s=(source/'Spin.cpp').read_text()
s=edit(s,'namespace bare_spin {','''namespace bare_spin {
bool bchAvx512Available() noexcept {
#if SPIN_BCH_AVX512 && !SPIN_TEST_NO_AVX512
    return __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512vl");
#else
    return false;
#endif
}
static u32 packFour(u32 x) {return (x&~1023U)|((x&255U)<<2)|((x>>8)&3U);}
''')
s=edit(s,'unsigned tileRows)\n','unsigned tileRows,BchBackend backend)\n')
s=edit(s,'(tileRows && !std::has_single_bit(tileRows))','(tileRows && (tileRows<2 || !std::has_single_bit(tileRows)))')
s=edit(s,'    if(rows%step())','''    if(backend!=BchBackend::Auto && backend!=BchBackend::Avx2 && backend!=BchBackend::Avx512)
        throw std::invalid_argument("unknown BCH backend");
    const bool four=c==Configuration::T128S19 && mTileRows>=4 && bchAvx512Available();
    if(backend==BchBackend::Avx512 && !four) throw std::invalid_argument("unsupported AVX-512 BCH configuration");
    mBchBackend=backend!=BchBackend::Avx2 && four?BchBackend::Avx512:BchBackend::Avx2;
    if(mBchBackend==BchBackend::Avx512) {mBchOffsets24.resize(3*n+4);mBchOffsets32.resize(n);}
    if(rows%step())''')
s=edit(s,'        pack(mSlots24.data()+3*inner,slot); pack(mOffsets24.data()+3*slot,local);',
 '''        pack(mSlots24.data()+3*inner,slot); pack(mOffsets24.data()+3*slot,local);
        if(mBchBackend==BchBackend::Avx512) {
            mBchOffsets32[slot]=packFour(local);pack(mBchOffsets24.data()+3*slot,packFour(local));
        }''')
s=edit(s,'    return mSlots24.capacity()+mOffsets24.capacity()+4*(',
 '    return mBchOffsets24.capacity()+4*mBchOffsets32.capacity()+mSlots24.capacity()+mOffsets24.capacity()+4*(')
s=edit(s,'    // The packed single-row path',
 '''    if(layout==Layout::Packed24) std::vector<u32>().swap(mBchOffsets32);
    else std::vector<u8>().swap(mBchOffsets24);
    // The packed single-row path''')
s=edit(s,'void Spin::encodeUnchecked(const block* in,block* out,Workspace& w,Layout layout) const {',
 '''void Spin::encodeUnchecked(const block* in,block* out,Workspace& w,Layout layout) const {
#if SPIN_BCH_AVX512
    if(mBchBackend==BchBackend::Avx512) {
        if(layout==Layout::Packed24) runFour<true>(in,out,w); else runFour<false>(in,out,w);
        return;
    }
#endif''')
s=edit(s,'        const auto outer=mRoute[i],slot=mSlots32[i];',
 '''        const auto outer=mRoute[i],slot=mSlots32[i];
        if(mBchBackend==BchBackend::Avx512 &&
           (slot>=mBchOffsets32.size() || mBchOffsets32[slot]!=packFour(mOffsets32[slot]) ||
            unpack(mBchOffsets24.data()+3*slot)!=mBchOffsets32[slot]))
            throw std::runtime_error("packed BCH offset mismatch");''')
# This port retains the historical layout interface. Do not import the newer
# transpose target's Auto/MessageLength helpers into the upstream class.
baseline=Path(__file__).resolve().parents[2]/'workstreams/inner_design/asymmetric/bch256/weight5/implementation/Weight5Spin.cpp'
inplace=baseline.read_text().split('void Spin::encodeInplace',1)[1].split('template<class Map,bool Packed',1)[0]
s=edit(s,'void Spin::validateSetup() const {','void Spin::encodeInplace'+inplace+'void Spin::validateSetup() const {')
put('Spin.cpp',s)
start=s.index('template<class Map,bool Packed> void Spin::run(')
end=s.index('void Spin::encodeUnchecked',start)
fast=s[start:end].replace('template<class Map,bool Packed> void Spin::run(', 'template<bool Packed> void Spin::runFour(')
fast=edit(fast,'    block* values=w.buckets.data();','    using Map=Map128S19;\n    block* values=w.buckets.data();')
fast=fast.replace('mOffsets24','mBchOffsets24').replace('mOffsets32','mBchOffsets32')
fast=edit(fast,'for(std::size_t j=0;j<tileSize;j+=512)\n            bchTranspose2(tile+j,tile+j+256,out+base/2+j/2,out+base/2+j/2+128);',
 'for(std::size_t j=0;j<tileSize;j+=1024)\n            bchTranspose4(tile+j,out+base/2+j/2);')
put('Fast.cpp','''#include "Spin.h"
#include "Inner.h"
#include <cstring>
namespace bare_spin {
void bchTranspose4(const block*,block*);
static OC_FORCEINLINE u32 unpack(const u8* p) {u32 v;std::memcpy(&v,p,4);return v&0xffffff;}
'''+fast+'''
template void Spin::runFour<true>(const block*,block*,Workspace&) const;
template void Spin::runFour<false>(const block*,block*,Workspace&) const;
}
''')
print('Generated bidirectional port; original forward and WideView schedules preserved.')
