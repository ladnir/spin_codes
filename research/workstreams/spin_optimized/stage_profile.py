"""Instrument only the tiled AVX-512 kernel in a separate diagnostic build."""
from pathlib import Path
import sys
p=Path(sys.argv[1]);f=(p/'Fast.cpp').read_text()
start=f.index(' void Spin::runFour(')
brace=f.index('{',start);end=brace+1;depth=1
while depth:
    depth+=(f[end]=='{')-(f[end]=='}');end+=1
kernel=f[brace:end]
def edit(old,new):
    global kernel
    assert kernel.count(old)==1,old
    kernel=kernel.replace(old,new)
edit('    block* values=w.buckets.data();',
     '    auto& prof=spin_profile::state(); ++prof.calls; auto timer=prof.begin();\n'
     '    block* values=w.buckets.data();')
edit('    block* tile=w.tile.data();',
     '    prof.end(0,timer);\n    block* tile=w.tile.data();')
edit('    for(std::size_t base=0;base<n;base+=tileSize) {',
     '    for(std::size_t base=0;base<n;base+=tileSize) {\n        timer=prof.begin();')
if 'if constexpr(Four)' in kernel:
    edit('        if constexpr(Four)', '        prof.end(1,timer);timer=prof.begin();\n        if constexpr(Four)')
else:
    edit('        for(std::size_t j=0;j<tileSize;j+=1024)',
         '        prof.end(1,timer);timer=prof.begin();\n        for(std::size_t j=0;j<tileSize;j+=1024)')
assert kernel.endswith('    }\n}')
kernel=kernel[:-7]+'        prof.end(2,timer);\n    }\n}'
f='#include "StageProfile.h"\n'+f[:brace]+kernel+f[end:]
(p/'Fast.cpp').write_text(f,newline='\n')
(p/'StageProfile.h').write_text((Path(__file__).parent/'StageProfile.h').read_text(),newline='\n')
