"""Generate isolated asymmetric kernels with an independent dense oracle."""
import hashlib
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
PARENT=HERE.parent
ROOT=PARENT.parents[1]
sys.path.insert(0,str(PARENT))
import generate_mixer as mixer
import generate_balanced as balanced


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    screen_path=HERE/'FEEDBACK_SCREEN.json';screen=json.loads(screen_path.read_text())
    for name,digest in screen['source_sha256'].items():assert sha(ROOT/name)==digest,name
    a_path=PARENT/'NO_CONSTANT_MAP.json';a_record=json.loads(a_path.read_text())
    a_header,a_stats=balanced.header(a_record['columns'])
    prior_path=PARENT/'BALANCED_IMPLEMENTATION.json';prior=json.loads(prior_path.read_text())
    sources=[Path(__file__),HERE/'AsymmetricInner.h',screen_path,a_path,prior_path,Path(balanced.__file__)]
    records=[]
    for name in ('mixed2_3','greedy3_2'):
        b_record=next(r for r in screen['candidates'] if r['name']==name)
        columns=b_record['b_columns'];assert all(1<=c.bit_count()<=3 for c in columns)
        code=a_header+'\nnamespace bare_spin {\nstruct AsymmetricMap : BalancedMap {\n'
        code+='static constexpr auto feedbackColumns=BalancedMap::columns;\n'
        code+='static constexpr std::array<std::uint32_t,T> columns{'+','.join(map(hex,columns))+'};\n'
        code+='static constexpr auto groupedColumns=columns;\n'
        code+='static constexpr std::array<unsigned,S> groupOrder{'+','.join(map(str,range(19)))+'};\n};\n}\n'
        for style in ('sparse','masked'):
            baseline=PARENT/'generated'/('balanced_'+style)/'CandidateSpin.cpp'
            expected=next(r['source_sha256']['CandidateSpin.cpp'] for r in prior['candidates'] if r['name']=='balanced_'+style)
            assert sha(baseline)==expected;sources.append(baseline)
            source=mixer.once(baseline.read_text(),'#include "BalancedMap.h"','#include "AsymmetricMap.h"')
            source=mixer.once(source,'#include "../../MixerInner.h"','#include "../../asymmetric/AsymmetricInner.h"')
            source=source.replace('<BalancedMap','<AsymmetricMap')
            source=mixer.once(source,'mixerReverse<Map,Quarter,','asymmetricReverse<Map,Quarter,')
            source=mixer.once(source,'"t128_s19_balanced_transvection_r1"',f'"t128_s19_asymmetric_{name}_r1"')
            source=mixer.once(source,'for(unsigned j=0;j<Map::S;++j) if((Map::columns[p]>>j)&1) syndrome[j]^=value;',
                'if constexpr(Quarter) {\n'
                '                for(unsigned j=0;j<Map::S;++j) if((Map::feedbackColumns[p]>>j)&1) syndrome[j]^=in[e*Map::T+p];\n'
                '            } else {\n'
                '                for(unsigned j=0;j<Map::S;++j) if((Map::columns[p]>>j)&1) syndrome[j]^=value;\n'
                '            }')
            directory=PARENT/'generated'/f'asymmetric_{name}_{style}';directory.mkdir(parents=True,exist_ok=True)
            (directory/'AsymmetricMap.h').write_text(code,encoding='utf-8',newline='\n')
            (directory/'CandidateSpin.cpp').write_text(source,encoding='utf-8',newline='\n')
            records.append(dict(name=directory.name,feedback_candidate=name,update_style=style,
                feedback_circuit=a_stats,source_sha256={p.name:sha(p) for p in (directory/'AsymmetricMap.h',directory/'CandidateSpin.cpp')}))
    payload=dict(status='ISOLATED_IMPLEMENTATIONS_GENERATED_NOT_A_CERTIFICATE',candidates=records,
                 source_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sources})
    (HERE/'IMPLEMENTATION.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    print(json.dumps([r['name'] for r in records]))


if __name__=='__main__':main()
