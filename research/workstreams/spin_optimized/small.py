"""Experimental K=2^16 direct routing; operates only on generated build files."""
from pathlib import Path
import sys
p=Path(sys.argv[1]);mode=int(sys.argv[2]);assert mode in (1,2,3,4)
prefetch={3:16,4:64}.get(mode,0)
if mode in (3,4): mode=1
def edit(s,a,b):
    if s.count(a)!=1: raise RuntimeError(f'small-path anchor changed: {a}')
    return s.replace(a,b)
def put(n,t): (p/n).write_text(t,newline='\n')
h=(p/'Spin.h').read_text()
member='mSmallRoute32' if mode==1 else 'mSmallRoute24'
typ='u32' if mode==1 else 'u8'
h=edit(h,'    template<class Map> void setupInner',
    f'    std::vector<{typ}> {member};\n    void runSmall(const block*,block*,Workspace&) const;\n    template<class Map> void setupInner')
put('Spin.h',h)
s=(p/'Spin.cpp').read_text()
setup=f'    if(mK==(1U<<16) && mBchBackend==BchBackend::Avx512) {{\n        {member}.resize('+('n' if mode==1 else '3*n+4')+');\n'
setup+='        for(std::size_t i=0;i<n;++i) '+(f'{member}[i]=packFour(mRoute[i]);' if mode==1 else f'pack({member}.data()+3*i,packFour(mRoute[i]));')+'\n    }\n'
s=edit(s,'    mSlots24.resize(3*n+4);',setup+'    mSlots24.resize(3*n+4);')
key='    return mBchOffsets24.capacity()' if 'mBchOffsets24' in s else '    return mSlots24.capacity()'
s=edit(s,key,f'    return {member}.capacity()*sizeof({typ})+'+key.removeprefix('    return '))
# The direct path never reads the tiled transpose schedule after compaction.
# Bidirectional still needs the original forward/WideView tables.
tables=('mBchOffsets24','mBchOffsets32') if 'mBchOffsets24' in s else ('mSlots24','mSlots32','mOffsets24','mOffsets32')
clear=f'    if(!{member}.empty()) {{\n'
for t in tables: clear+=f'        decltype({t})().swap({t});\n'
clear+='    }\n'
s=edit(s,'    mRetainedLayout=layout;mCompacted=true;',clear+'    mRetainedLayout=layout;mCompacted=true;')
verify=f'''    if(!mCompacted && !{member}.empty()) {{
        for(std::size_t i=0;i<codeBlocks();++i)
            if('''+(f'{member}[i]' if mode==1 else f'unpack({member}.data()+3*i)')+'''!=packFour(mRoute[i]))
                throw std::runtime_error("small direct route mismatch");
    }
'''
s=edit(s,'void Spin::validateSetup() const {','void Spin::validateSetup() const {\n'+verify)
s=edit(s,'    if(mBchBackend==BchBackend::Avx512) {\n        if(layout==Layout::Packed24)',
    '    if(mBchBackend==BchBackend::Avx512) {\n        if(mK==(1U<<16)) {runSmall(in,out,w);return;}\n        if(layout==Layout::Packed24)')
put('Spin.cpp',s)
f=(p/'Fast.cpp').read_text()
inner='innerReverse<Map>' if 'using Map=Map128S19' in f else 'asymmetricReverse<Map,true,W5_MASKED>'
index=f'{member}[i]' if mode==1 else f'unpack({member}.data()+3*i)'
map_type='Map128S19' if 'using Map=Map128S19' in f else 'AsymmetricMap'
direct=f'''void Spin::runSmall(const block* in,block* out,Workspace& w) const {{
        using Map={map_type};
        block* direct=w.buckets.data();
        {inner}(in,131072,mFieldRows.data(),[&](std::size_t i,block v) {{
            direct[{index}]=v;
        }});
        for(std::size_t j=0;j<131072;j+=1024) bchTranspose4(direct+j,out+j/2);
        return;
    }}
'''
assert f.endswith('}\n')
f=f[:-2]+direct+'}\n'
put('Fast.cpp',f)
if prefetch:
    future=f'{member}[i-{prefetch}]'
    f=edit(f,f'            direct[{index}]=v;',
        f'            if(i>={prefetch}) __builtin_prefetch(direct+{future},1,3);\n            direct[{index}]=v;')
    put('Fast.cpp',f)
print(f'Generated experimental K16 direct route mode {mode}.')
