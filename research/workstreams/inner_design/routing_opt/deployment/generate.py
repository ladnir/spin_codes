"""Apply the supported routing hook to the unchanged experimental inner."""
from pathlib import Path
import hashlib
import json
HERE=Path(__file__).resolve().parent
INNER=HERE.parents[1]
ROOT=INNER.parents[1]
source=INNER/'generated/asymmetric_greedy3_2_sparse/CandidateSpin.cpp'
original=source.read_text()
record=next(r for r in json.loads((INNER/'asymmetric/IMPLEMENTATION.json').read_text())['candidates'] if r['name']=='asymmetric_greedy3_2_sparse')
assert hashlib.sha256(source.read_bytes()).hexdigest()==record['source_sha256']['CandidateSpin.cpp']
supported=(ROOT/'workstreams/bare_bch_rm2sub/Spin.cpp').read_text()
code=original.replace('#include "Spin.h"','#include "Spin.h"\n#include "WorkspaceRouting.h"',1)
start=supported.index('Spin::Workspace::Workspace(')
end=supported.index('\nstd::size_t Spin::Workspace::bytes()',start)
old='Spin::Workspace::Workspace(const Spin& code):buckets(code.codeBlocks()),tile(code.tileBlocks()) {}'
assert old in code
code=code.replace(old,supported[start:end])
start=supported.index('    if constexpr(requires {Map::tunedWorkspaceRouting;})')
end=supported.index('    for(std::size_t base=0;base<n;base+=tileSize)',supported.index('\n    }',start)+6)
hook=supported[start:end]
marker='    for(std::size_t base=0;base<n;base+=tileSize) {'
assert code.count(marker)==1
code=code.replace(marker,hook+marker)
start=supported.index('        if constexpr(workspace_routing::compiled)',supported.index('void Spin::encodeUnchecked'))
end=supported.index('        if(layout==Layout::Packed24) run<Map128S19,true,true>',start)
dispatch=supported[start:end].replace('TunedMap<Map128S19>','TunedMap<AsymmetricMap>')
marker='        if(layout==Layout::Packed24) run<AsymmetricMap,true,true>'
assert code.count(marker)==1
code=code.replace(marker,dispatch+marker)
(HERE/'AsymmetricSpin.cpp').write_text(code,encoding='utf-8',newline='\n')
(HERE/'AsymmetricMap.h').write_text((source.parent/'AsymmetricMap.h').read_text(),encoding='utf-8',newline='\n')
print('Generated deployment candidate with unchanged asymmetric inner')
