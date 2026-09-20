"""Isolated macroblock implementations; exact adjoint, no new certificate."""
import hashlib
import json
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;PARENT=HERE.parent;ROOT=PARENT.parents[1]
sys.path.insert(0,str(PARENT))
import generate_mixer as mixer
import generate_balanced as balanced


def main():
    path=HERE/'MACRO_Q1.json';screen=json.loads(path.read_text())
    for name,digest in screen['source_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    a0=json.loads((PARENT/'NO_CONSTANT_MAP.json').read_text())['columns']
    header,_=balanced.header(a0)
    template=PARENT/'generated/asymmetric_greedy3_2_sparse/CandidateSpin.cpp'
    original=template.read_text();sources=[Path(__file__),HERE/'MacroInner.h',path,template,Path(balanced.__file__)]
    records=[]
    for repeat in (2,4):
        row=next(r for r in screen['results'] if r['name']==f'repeat{repeat}_random3_seed1')
        name=f'redesign_macro{128*repeat}_r2'
        code=original.replace('#include "AsymmetricMap.h"','#include "MacroMap.h"').replace('#include "../../asymmetric/AsymmetricInner.h"','#include "../../redesign/MacroInner.h"')
        code=code.replace('<AsymmetricMap','<MacroMap').replace('asymmetricReverse<Map,Quarter,false>','macroReverse<Map,Quarter>')
        code=code.replace('"t128_s19_asymmetric_greedy3_2_r1"',f'"{name}_NOT_CERTIFIED"')
        code=mixer.once(code,'unsigned Spin::step() const noexcept {','unsigned Spin::step() const noexcept {\n    if(mOuter==Outer::Bch128x32) return MacroMap::T;')
        start=code.index('    if constexpr(Map::T==128 && Map::S==19)',code.index('void Spin::setupInner'))
        end=code.index('    mCoefficients.resize(epochs);',start)
        setup='''    if constexpr(std::is_same_v<Map,MacroMap>) {
        mCoefficients.resize(2*Map::Rounds*epochs);mFieldRows.resize(2*Map::Rounds*epochs);
        for(std::size_t e=0;e<epochs*Map::Rounds;++e) {
            u32 u;do u=u32(splitmix(seed))&((1U<<19)-1);while(!u);
            u32 v=u32(splitmix(seed))&((1U<<19)-1);
            if(std::popcount(u&v)&1) v^=u&-u;
            mCoefficients[2*e]=mFieldRows[2*e]=u;
            mCoefficients[2*e+1]=mFieldRows[2*e+1]=v;
        }
        return;
    }
'''
        code=code[:start]+setup+code[end:]
        code=mixer.once(code,'e<codeBlocks()/step();++e','e<codeBlocks()/step()*MacroMap::Rounds;++e')
        start=code.index('        for(unsigned j=0;j<Map::S;++j) {\n            block v=syndrome[j];',code.index('void Spin::oracle'))
        end=code.index('        state=next;',start)+len('        state=next;')
        update='''        if constexpr(Quarter) {
            // Independent dense transpose of each rank-one factor, reverse order.
            for(unsigned round=Map::Rounds;round-->0;) {
                const auto u=mCoefficients[2*(e*Map::Rounds+round)];
                const auto v=mCoefficients[2*(e*Map::Rounds+round)+1];
                for(unsigned j=0;j<Map::S;++j) {
                    auto value=state[j];
                    if((v>>j)&1) for(unsigned b=0;b<Map::S;++b) if((u>>b)&1) value^=state[b];
                    next[j]=value;
                }
                state=next;
            }
            for(unsigned j=0;j<Map::S;++j) state[j]^=syndrome[j];
        } else {
            for(unsigned j=0;j<Map::S;++j) {
                auto value=syndrome[j];const auto mask=referenceMultiply<Map>(1U<<j,mCoefficients[e]);
                for(unsigned b=0;b<Map::S;++b) if((mask>>b)&1) value^=state[b];
                next[j]=value;
            }
            state=next;
        }'''
        code=code[:start]+update+code[end:]
        map_header=header+'\nnamespace bare_spin {\nstruct MacroMap : BalancedMap {\n'
        map_header+=f'static constexpr unsigned T={128*repeat},Rounds=2;\n'
        map_header+='static constexpr std::array<std::uint32_t,T> columns{'+','.join(map(hex,row['b_columns']))+'};\n'
        map_header+='static constexpr std::array<std::uint32_t,T> feedbackColumns{'+','.join(map(hex,row['a_columns']))+'};\n};\n}\n'
        directory=PARENT/'generated'/name;directory.mkdir(parents=True,exist_ok=True)
        for filename,text in (('CandidateSpin.cpp',code),('MacroMap.h',map_header)):
            (directory/filename).write_text(text,encoding='utf-8',newline='\n')
        records.append(dict(name=name,map_record=row['name'],rounds=2,
            generated_sha256={f:hashlib.sha256((directory/f).read_bytes()).hexdigest() for f in ('CandidateSpin.cpp','MacroMap.h')}))
    payload=dict(status='ISOLATED_MACRO_IMPLEMENTATIONS_NOT_CERTIFIED',candidates=records,
        source_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (HERE/'MACRO_IMPLEMENTATION.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    print(json.dumps([r['name'] for r in records]))


if __name__=='__main__':main()
