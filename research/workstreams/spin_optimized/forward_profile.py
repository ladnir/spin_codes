"""Diagnostic clocks only; ordinary forward builds have no instrumentation."""
from pathlib import Path
import sys
p=Path(sys.argv[1])
for filename,method in [('Spin.cpp','runForward'),('ForwardFast.cpp','runForwardFour')]:
    path=p/filename
    if not path.exists(): continue
    s=path.read_text();start=s.index(f'template<class Map,bool Packed> void Spin::{method}(')
    brace=s.index('{',start);end=brace+1;depth=1
    while depth:
        depth+=(s[end]=='{')-(s[end]=='}');end+=1
    k=s[start:end]
    k=k.replace('    block* values=', '    auto& prof=forward_profile::state();++prof.calls;\n    block* values=',1)
    k=k.replace('for(std::size_t base=0;base<n;base+=tileSize) {',
                'for(std::size_t base=0;base<n;base+=tileSize) {\n        auto timer=forward_profile::Clock::now();')
    k=k.replace('        for(std::size_t j=0;j<tileSize;++j) {',
                '        prof.add(0,timer);timer=forward_profile::Clock::now();\n        for(std::size_t j=0;j<tileSize;++j) {')
    k=k.replace('    }\n    innerForward<Map>',
                '        prof.add(1,timer);\n    }\n    auto innerTimer=forward_profile::Clock::now();\n    innerForward<Map>')
    k=k[:-1]+'    prof.add(2,innerTimer);\n}'
    path.write_text('#include "ForwardProfile.h"\n'+s[:start]+k+s[end:],newline='\n')
(p/'ForwardProfile.h').write_bytes((Path(__file__).parent/'ForwardProfile.h').read_bytes())
