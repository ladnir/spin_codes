"""Isolated one-transvection variants; certified baseline is never edited.

The dense-row, sparse-update and masked-update variants sample identical maps.
The oracle uses raw u/v and forms matrix rows independently of the hot masks.
"""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
BARE=ROOT/'workstreams/bare_bch_rm2sub'


def once(text,old,new):
    assert text.count(old)==1,old
    return text.replace(old,new)


def main():
    baseline=(BARE/'Spin.cpp').read_text()
    records=[]
    for style in ('rows','sparse','masked'):
        name='mixer_'+style;directory=HERE/'generated'/name;directory.mkdir(parents=True,exist_ok=True)
        source=once(baseline,'#include "Inner.h"','#include "Inner.h"\n#include "../../MixerInner.h"')
        source=once(source,'#include "../rate_quarter_bch/implementation/generated/QuarterCircuit.h"','#include "QuarterCircuit.h"')
        source=once(source,'const char* Spin::name() const noexcept {',
                    'const char* Spin::name() const noexcept {\n    if(mOuter==Outer::Bch128x32) return "t128_s19_transvection_r1";')
        setup='''
    if constexpr(Map::T==128 && Map::S==19) {
        if(mOuter==Outer::Bch128x32) {
            mCoefficients.resize(2*epochs);mFieldRows.resize(FIELD_SIZE*epochs);
            for(std::size_t e=0;e<epochs;++e) {
                u32 u;
                do u=u32(splitmix(seed))&((1U<<Map::S)-1); while(!u);
                u32 v=u32(splitmix(seed))&((1U<<Map::S)-1);
                if(std::popcount(u&v)&1) v^=u&-u;
                mCoefficients[2*e]=u;mCoefficients[2*e+1]=v;
                FIELD_ASSIGN
            }
            return;
        }
    }
'''
        if style=='rows':
            setup=setup.replace('FIELD_SIZE','Map::S').replace('FIELD_ASSIGN',
                'for(unsigned j=0;j<Map::S;++j) mFieldRows[e*Map::S+j]=(1U<<j)^(((v>>j)&1)?u:0);')
        else:
            setup=setup.replace('FIELD_SIZE','2').replace('FIELD_ASSIGN','''u32 grouped=0;
                if constexpr(SPIN_GROUPED_A) {
                    for(unsigned j=0;j<Map::S;++j) grouped|=((u>>Map::groupOrder[j])&1)<<j;
                } else grouped=u;
                mFieldRows[2*e]=grouped;mFieldRows[2*e+1]=v;''')
            source=once(source,'innerReverse<Map>(in,codeBlocks(),mFieldRows.data(),',
                        f'mixerReverse<Map,Quarter,{str(style=="masked").lower()}>(in,codeBlocks(),mFieldRows.data(),')
        source=once(source,'const auto epochs=codeBlocks()/Map::T;',
                    'const auto epochs=codeBlocks()/Map::T;'+setup)
        source=once(source,'const u32 mask=referenceMultiply<Map>(1U<<j,mCoefficients[e]);','''u32 mask;
            if constexpr(Quarter) {
                const u32 u=mCoefficients[2*e],v=mCoefficients[2*e+1];
                mask=1U<<j;
                if((v>>j)&1) mask^=u;
            } else mask=referenceMultiply<Map>(1U<<j,mCoefficients[e]);''')
        source=once(source,'void Spin::validateSetup() const {','''void Spin::validateSetup() const {
    if(!mCompacted && mOuter==Outer::Bch128x32) {
        for(std::size_t e=0;e<codeBlocks()/step();++e) {
            const auto u=mCoefficients[2*e],v=mCoefficients[2*e+1];
            if(!u || (u>>19) || (v>>19) || (std::popcount(u&v)&1))
                throw std::runtime_error("invalid transvection sample");
        }
    }''')
        path=directory/'CandidateSpin.cpp';path.write_text(source,encoding='utf-8',newline='\n')
        records.append(dict(name=name,source_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    inputs=[Path(__file__),HERE/'MixerInner.h',BARE/'Spin.cpp',BARE/'Spin.h',BARE/'Inner.h',BARE/'generated/SelectedMaps.h']
    payload=dict(status='UNCERTIFIED_NEW_MIXER_IMPLEMENTATION',rounds=1,
                 source_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                 candidates=records)
    (HERE/'MIXER_IMPLEMENTATION.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(records,indent=2))


if __name__=='__main__':main()
