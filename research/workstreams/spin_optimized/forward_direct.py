"""Cache-resident forward route: remove tile-to-bucket copying, preserve the map."""
from pathlib import Path
import sys

p=Path(sys.argv[1])
def edit(text,old,new):
    assert text.count(old)==1,old
    return text.replace(old,new)

h=(p/'Spin.h').read_text()
h=edit(h,'    template<class Map, bool Packed> void runForwardFour',
    '    std::vector<u32> mForwardDirect;\n'
    '    template<class Map> void runForwardDirect(const block*,block*,Workspace&) const;\n'
    '    template<class Map, bool Packed> void runForwardFour')
(p/'Spin.h').write_text(h,newline='\n')
s=(p/'Spin.cpp').read_text()
s=edit(s,'    mSlots24.resize(3*n+4);', '''    if(mK<=(1U<<18) && mBchBackend==BchBackend::Avx512) {
        mForwardDirect.resize(n);
        for(std::size_t i=0;i<n;++i) mForwardDirect[i]=packFour(mRoute[i]);
    }
    mSlots24.resize(3*n+4);''')
start=s.index('std::size_t Spin::setupBytes() const noexcept {')
index=s.index('    return ',start)+len('    return ')
s=s[:index]+'4*mForwardDirect.capacity()+'+s[index:]
s=edit(s,'void Spin::validateSetup() const {', '''void Spin::validateSetup() const {
    if(!mCompacted && !mForwardDirect.empty()) {
        for(std::size_t i=0;i<codeBlocks();++i)
            if(mForwardDirect[i]!=packFour(mRoute[i])) throw std::runtime_error("forward direct route mismatch");
    }''')
anchor='''#define FAST_FORWARD(Enum,Map) case Configuration::Enum: if(layout==Layout::Packed24) runForwardFour<Map,true>(in,out,w); else runForwardFour<Map,false>(in,out,w); return'''
s=edit(s,anchor,'''#define FAST_FORWARD(Enum,Map) case Configuration::Enum: if(!mForwardDirect.empty()) runForwardDirect<Map>(in,out,w); else if(layout==Layout::Packed24) runForwardFour<Map,true>(in,out,w); else runForwardFour<Map,false>(in,out,w); return''')
(p/'Spin.cpp').write_text(s,newline='\n')
f=(p/'ForwardFast.cpp').read_text()
assert f.endswith('}\n')
kernel='''template<class Map> void Spin::runForwardDirect(const block* in,block* out,Workspace& w) const {
    auto* values=w.buckets.data();const auto n=codeBlocks();
    for(std::size_t j=0;j<n;j+=1024) bchForward4(in+j/2,values+j);
    innerForward<Map>(n,mForwardFieldRows.data(),[&](std::size_t i) {
        return values[mForwardDirect[i]];
    },out);
}
'''
for name in ('Map128S19','Map64S12','Map64S12R2'):
    kernel+=f'template void Spin::runForwardDirect<{name}>(const block*,block*,Workspace&) const;\n'
(p/'ForwardFast.cpp').write_text(f[:-2]+kernel+'}\n',newline='\n')
