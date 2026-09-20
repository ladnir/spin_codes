"""General aligned lengths for the transpose target; preserve full-tile kernels."""
from pathlib import Path
import re
import sys
p=Path(sys.argv[1])
direct_max=int(sys.argv[2]) if len(sys.argv)>2 else 458752
assert 0<=direct_max<=((1<<32)-1)//2
def edit(s,a,b):
    assert s.count(a)==1,a
    return s.replace(a,b)
def put(n,s): (p/n).write_text(s,newline='\n')
h=(p/'Spin.h').read_text();s=(p/'Spin.cpp').read_text();f=(p/'Fast.cpp').read_text()
(p/'LengthGeometry.h').write_bytes((Path(__file__).parent/'LengthGeometry.h').read_bytes())
h='#include "LengthGeometry.h"\n'+h
assert 'forwardBits(' not in h, 'Transpose target only'
h=edit(h,'enum class Layout { Packed24, Indices32 };', '''// Explicit length in logical elements; transpose maps 2K elements to K.
struct MessageLength { std::size_t value; };
enum class Layout { Packed24, Indices32, Auto };''')
start=h.index('    Spin(Configuration');end=h.index(';',start)+1
decl=h[start:end]
h=h[:end]+'\n'+decl.replace('unsigned messageExponent','MessageLength length')+h[end:]
h=h.replace('Layout layout=Layout::Packed24','Layout layout=Layout::Auto')
h=edit(h,'    std::size_t messageBlocks() const noexcept', '''    static constexpr std::size_t maxMessageBlocks=length_geometry::maxMessageBlocks;
    bool packed24Available() const noexcept {return codeBlocks()<=(std::size_t{1}<<24);}
    Layout preferredLayout() const noexcept {
        return mCompacted?mRetainedLayout:(packed24Available()?Layout::Packed24:Layout::Indices32);
    }
    std::size_t messageBlocks() const noexcept''')
h=edit(h,'    unsigned mTileRows;', '    unsigned mTileRows;\n    bool mPartialTile=false;')
h=edit(h,'    template<class Map> void setupInner', '''    template<class Map,bool Packed,bool Quarter> void runTail(const block*,block*,Workspace&) const;
    template<class Map,bool Packed> void runFourTail(const block*,block*,Workspace&) const;
    template<class Map> void setupInner''')
old='Spin::Spin(Configuration c,unsigned exponent,u64 routeSeed,u64 coefficientSeed,unsigned tileRows,Outer outer,BchBackend backend)'
s=edit(s,old,'''static std::size_t exponentLength(unsigned exponent) {
    return length_geometry::fromExponent(exponent);
}
'''+old+'''
    :Spin(c,MessageLength{exponentLength(exponent)},routeSeed,coefficientSeed,tileRows,outer,backend) {}
Spin::Spin(Configuration c,MessageLength length,u64 routeSeed,u64 coefficientSeed,unsigned tileRows,Outer outer,BchBackend backend)''')
s=edit(s,'exponent<16 || exponent>20 || ','')
s=s.replace('supported message exponents: 16..20; tile rows must be a power of two >=2','tile rows must be a power of two >=2')
s,n=re.subn(r'    if\([^\n]*exponent!=16\)\n        throw [^\n]*;\n','',s)
assert n==1, 'certificate-only K16 restriction'
s=edit(s,'    mK=std::size_t{1}<<exponent;', '''    mK=length.value;
    length_geometry::check(mK,step());''')
s=s.replace('exponent<=16','mK<=(1U<<16)').replace('exponent<=18','mK<=(1U<<18)')
s=edit(s,'std::min<std::size_t>(selectedTile,rows)', 'std::min<std::size_t>(selectedTile,std::bit_floor(rows))')
s=edit(s,'    if(rows%step())', '    mPartialTile=(rows%mTileRows)!=0;\n    if(rows%step())')
# Configuration-specific direct kernels remain exact-size specializations.
for exponent in (18,20):
    s=s.replace(f'if(mK==(1U<<{exponent}) && mBchBackend',
                f'if(mConfig==Configuration::T128S19 && mK==(1U<<{exponent}) && mBchBackend')
    s=s.replace(f'if(mK==(1U<<{exponent})) {{runK{exponent}',
                f'if(mConfig==Configuration::T128S19 && mK==(1U<<{exponent})) {{runK{exponent}')
for suffix in ('','R2'):
    needle=f'if(mConfig==Configuration::T64S12{suffix}) {{runSmallK16{suffix}(in,out,w);return;}}'
    if needle in s:
        s=edit(s,needle,f'if(mConfig==Configuration::T64S12{suffix} && mK==(1U<<16)) {{runSmallK16{suffix}(in,out,w);return;}}')
        # Add ordinary tiled dispatch for these maps before the S19 small path.
        target=f'runFour<Map64S12{suffix}Transpose'
        loc=s.index('        if(mK==(1U<<16))',s.index('void Spin::encodeUnchecked'))
        s=s[:loc]+f'        if(mConfig==Configuration::T64S12{suffix}) {{if(layout==Layout::Packed24) {target},true>(in,out,w); else {target},false>(in,out,w);return;}}\n'+s[loc:]
s=edit(s,'    mSlots24.resize(3*n+4); mOffsets24.resize(3*n+4);',
       '    if(packed24Available()) {mSlots24.resize(3*n+4); mOffsets24.resize(3*n+4);}')
s=edit(s,'        pack(mSlots24.data()+3*inner,slot); pack(mOffsets24.data()+3*slot,offset);',
       '        if(packed24Available()) {pack(mSlots24.data()+3*inner,slot); pack(mOffsets24.data()+3*slot,offset);}')
s=edit(s,'    std::vector<u32> counts(n/tile,u32(tile));',
       '    std::vector<u32> counts((n+tile-1)/tile,u32(tile));\n    if(n%tile) counts.back()=u32(n%tile);')
s=edit(s,'if(unpack(mSlots24.data()+3*i)!=slot || unpack(mOffsets24.data()+3*slot)!=mOffsets32[slot] ||',
       'if((packed24Available() && (unpack(mSlots24.data()+3*i)!=slot || unpack(mOffsets24.data()+3*slot)!=mOffsets32[slot])) ||')
for method in ('encode','encodeInplace','encodeUnchecked','compact'):
    start=s.index(f'void Spin::{method}(');brace=s.index('{',start)+1
    extra='\n    if(layout==Layout::Auto) layout=preferredLayout();\n'
    if method!='encodeUnchecked':
        extra+='    if(layout==Layout::Packed24 && !packed24Available()) throw std::invalid_argument("Packed24 index range exceeded; use Auto or Indices32");\n'
    s=s[:brace]+extra+s[brace:]

def get_kernel(text, signature):
    start=text.index(signature);brace=text.index('{',start);end=brace+1;depth=1
    while depth:
        depth+=(text[end]=='{')-(text[end]=='}');end+=1
    return text[start:end]

# Clone only the generic-tail versions. The original full-tile bodies stay byte
# for byte unchanged; their routing loops get no min(), divide, or tail branch.
base=get_kernel(s,'template<class Map,bool Packed,bool Quarter,bool Four> void Spin::run(')
tail=base.replace('template<class Map,bool Packed,bool Quarter,bool Four> void Spin::run(',
                  'template<class Map,bool Packed,bool Quarter> void Spin::runTail(')
tail=edit(tail,'    block* values=', '    constexpr bool Four=false;\n    block* values=')
four=get_kernel(f,'template<class Map,bool Packed> void Spin::runFour(')
four_tail=four.replace('Spin::runFour(', 'Spin::runFourTail(')
def tail_loop(kernel):
    kernel=edit(kernel,'for(std::size_t base=0;base<n;base+=tileSize) {',
                'for(std::size_t base=0;base<n;base+=tileSize) {\n        const auto activeSize=std::min(tileSize,n-base);')
    return kernel.replace('j<tileSize','j<activeSize').replace('j+W5_PREFETCH<tileSize','j+W5_PREFETCH<activeSize')
tail=tail_loop(tail);four_tail=tail_loop(four_tail)
s=edit(s,'void Spin::encodeUnchecked',tail+'\nvoid Spin::encodeUnchecked')
f=f[:-2]+four_tail+'\n'+''.join(
    f'template void Spin::runFourTail<{m},{v}>(const block*,block*,Workspace&) const;\n'
    for m in ('AsymmetricMap','Map64S12Transpose','Map64S12R2Transpose') for v in ('true','false'))+'}\n'
f='#include <algorithm>\n'+f
start=s.index('void Spin::encodeUnchecked');loc=s.index('#if SPIN_BCH_AVX512',start)
dispatch='''    if(mPartialTile) {
#define TAIL(Enum,Map) case Configuration::Enum: \\
    if(mBchBackend==BchBackend::Avx512) { \\
        TAIL_FOUR(Map) \\
    } else {if(layout==Layout::Packed24) runTail<Map,true,false>(in,out,w);else runTail<Map,false,false>(in,out,w);} return
#if SPIN_BCH_AVX512
#define TAIL_FOUR(Map) if(layout==Layout::Packed24) runFourTail<Map,true>(in,out,w);else runFourTail<Map,false>(in,out,w);
#else
#define TAIL_FOUR(Map) throw std::logic_error("AVX-512 backend unavailable");
#endif
        switch(mConfig) {
            TAIL(T128S19,AsymmetricMap);TAIL(T64S12,Map64S12Transpose);TAIL(T64S12R2,Map64S12R2Transpose);
            default:throw std::logic_error("unsupported transpose map");
        }
#undef TAIL_FOUR
#undef TAIL
    }
'''
s=s[:loc]+dispatch+s[loc:]
# A runtime-length direct schedule fills the gaps between the exact-size
# specializations. Its inner map remains compile-time selected; no callback
# erasure, allocation, or per-element routing-mode branch enters encoding.
if direct_max:
    h=edit(h,'    template<class Map> void setupInner', '''    std::vector<u32> mRangeRoute32;
    template<class Map> void runRange(const block*,block*,Workspace&) const;
    template<class Map> void setupInner''')
    exclusions=[]
    if 'mSmallRoute' in h: exclusions.append('mK!=(1U<<16)')
    if 'mK18Route' in h: exclusions.append('(mConfig!=Configuration::T128S19 || mK!=(1U<<18))')
    if 'mK20Route' in h: exclusions.append('(mConfig!=Configuration::T128S19 || mK!=(1U<<20))')
    condition=' && '.join([f'mK<={direct_max}', 'mBchBackend==BchBackend::Avx512']+exclusions)
    s=edit(s,'    if(packed24Available()) {mSlots24.resize',f'''    if({condition}) {{
        mRangeRoute32.resize(n);
        for(std::size_t i=0;i<n;++i) mRangeRoute32[i]=packFour(mRoute[i]);
    }}
    if(packed24Available()) {{mSlots24.resize''')
    start=s.index('std::size_t Spin::setupBytes() const noexcept {')
    loc=s.index('    return ',start)+len('    return ')
    s=s[:loc]+'mRangeRoute32.capacity()*sizeof(u32)+'+s[loc:]
    s=edit(s,'    mRetainedLayout=layout;mCompacted=true;', '''    if(!mRangeRoute32.empty()) {
        decltype(mSlots24)().swap(mSlots24); decltype(mSlots32)().swap(mSlots32);
        decltype(mOffsets24)().swap(mOffsets24); decltype(mOffsets32)().swap(mOffsets32);
    }
    mRetainedLayout=layout;mCompacted=true;''')
    s=edit(s,'void Spin::validateSetup() const {', '''void Spin::validateSetup() const {
    if(!mCompacted && !mRangeRoute32.empty()) {
        for(std::size_t i=0;i<codeBlocks();++i)
            if(mRangeRoute32[i]!=packFour(mRoute[i])) throw std::runtime_error("range route mismatch");
    }''')
    s=edit(s,'    if(mPartialTile) {', '''#if SPIN_BCH_AVX512
    if(!mRangeRoute32.empty()) {
        switch(mConfig) {
            case Configuration::T128S19:runRange<AsymmetricMap>(in,out,w);return;
            case Configuration::T64S12:runRange<Map64S12Transpose>(in,out,w);return;
            case Configuration::T64S12R2:runRange<Map64S12R2Transpose>(in,out,w);return;
            default:throw std::logic_error("unsupported transpose map");
        }
    }
#endif
    if(mPartialTile) {''')
    kernel='''template<class Map> void Spin::runRange(const block* in,block* out,Workspace& w) const {
    block* direct=w.buckets.data();
    const auto n=codeBlocks();
    selectedReverse<Map>(in,n,mFieldRows.data(),[&](std::size_t i,block v) {
        direct[mRangeRoute32[i]]=v;
    });
    for(std::size_t j=0;j<n;j+=1024) bchTranspose4(direct+j,out+j/2);
}
'''
    f=f[:-2]+kernel+''.join(f'template void Spin::runRange<{m}>(const block*,block*,Workspace&) const;\n'
        for m in ('AsymmetricMap','Map64S12Transpose','Map64S12R2Transpose'))+'}\n'
put('Spin.h',h);put('Spin.cpp',s);put('Fast.cpp',f)
