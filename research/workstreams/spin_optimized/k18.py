"""Generate an isolated K18/K20 direct-routing kernel without changing the code map."""
from pathlib import Path
import sys
p=Path(sys.argv[1]);mode=int(sys.argv[2]);assert mode in (1,2,3,4)
exponent=int(sys.argv[3]) if len(sys.argv)>3 else 18
assert exponent in (18,20)
label=f'K{exponent}';length=1 << (exponent+1)
prefetch={3:16,4:64}.get(mode,0)
packed=mode==2
def edit(text,a,b):
    assert text.count(a)==1,a
    return text.replace(a,b)
def put(name,text): (p/name).write_text(text,newline='\n')
h=(p/'Spin.h').read_text();s=(p/'Spin.cpp').read_text();f=(p/'Fast.cpp').read_text()
bidir='forwardBits(' in h
member=f'm{label}Route24' if packed else f'm{label}Route32'
typ='u8' if packed else 'u32'
index=f'unpack({member}.data()+3*i)' if packed else f'{member}[i]'
h=edit(h,'    template<class Map> void setupInner',
    f'    std::vector<{typ}> {member};\n    void run{label}(const block*,block*,Workspace&) const;\n    template<class Map> void setupInner')
setup=f'    if(mK==(1U<<{exponent}) && mBchBackend==BchBackend::Avx512) {{\n        {member}.resize('+('3*n+4' if packed else 'n')+');\n'
setup+='        for(std::size_t i=0;i<n;++i) '+(f'pack({member}.data()+3*i,packFour(mRoute[i]));' if packed else f'{member}[i]=packFour(mRoute[i]);')+'\n    }\n'
s=edit(s,'    mSlots24.resize(3*n+4);',setup+'    mSlots24.resize(3*n+4);')
start=s.index('std::size_t Spin::setupBytes() const noexcept {')
idx=s.index('    return ',start)+len('    return ')
s=s[:idx]+f'{member}.capacity()*sizeof({typ})+'+s[idx:]
tables=('mBchOffsets24','mBchOffsets32') if bidir else ('mSlots24','mSlots32','mOffsets24','mOffsets32')
clear=f'    if(!{member}.empty()) {{\n'+''.join(f'        decltype({t})().swap({t});\n' for t in tables)+'    }\n'
s=edit(s,'    mRetainedLayout=layout;mCompacted=true;',clear+'    mRetainedLayout=layout;mCompacted=true;')
s=edit(s,'void Spin::validateSetup() const {',f'''void Spin::validateSetup() const {{
    if(!mCompacted && !{member}.empty()) {{
        for(std::size_t i=0;i<codeBlocks();++i)
            if({index}!=packFour(mRoute[i])) throw std::runtime_error("{label} direct route mismatch");
    }}''')
anchor='    if(mBchBackend==BchBackend::Avx512) {\n'
idx=s.index(anchor,s.index('void Spin::encodeUnchecked'))+len(anchor)
s=s[:idx]+f'        if(mK==(1U<<{exponent})) {{run{label}(in,out,w);return;}}\n'+s[idx:]
inner='innerReverse<Map>' if bidir else 'selectedReverse<Map>'
map_type='Map128S19' if bidir else 'AsymmetricMap'
future=f'{member}[i-{prefetch}]'
hint=f'            if(i>={prefetch}) __builtin_prefetch(direct+{future},1,3);\n' if prefetch else ''
new=f'''void Spin::run{label}(const block* in,block* out,Workspace& w) const {{
    using Map={map_type};
    block* direct=w.buckets.data();
    {inner}(in,{length},mFieldRows.data(),[&](std::size_t i,block v) {{
{hint}        direct[{index}]=v;
    }});
    for(std::size_t j=0;j<{length};j+=1024) bchTranspose4(direct+j,out+j/2);
}}
'''
assert f.endswith('}\n')
f=f[:-2]+new+'}\n'
put('Spin.h',h);put('Spin.cpp',s);put('Fast.cpp',f)
print('Generated',label,'direct-routing mode',mode)
