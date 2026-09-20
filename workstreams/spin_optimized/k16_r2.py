"""Add the explicit two-round K16 option after k16.py, preserving old kernels."""
from pathlib import Path
import sys
p=Path(sys.argv[1]);here=Path(__file__).resolve().parent
def edit(text,a,b):
    assert text.count(a)==1,a
    return text.replace(a,b)
def put(name,text): (p/name).write_text(text,newline='\n')
for name in ('ImtRounds.h','K16R2Inner.h'): put(name,(here/name).read_text())
h=(p/'Spin.h').read_text();s=(p/'Spin.cpp').read_text();f=(p/'Fast.cpp').read_text()
bidir='forwardBits(' in h
h=edit(h,'T64S12 };','T64S12, T64S12R2 };')
s=edit(s,'static_cast<unsigned>(Configuration::T64S12)','static_cast<unsigned>(Configuration::T64S12R2)')
# Eligibility is selected outside the kernels. Preserve wideView's restriction.
s=s.replace('c==Configuration::T64S12','(c==Configuration::T64S12 || c==Configuration::T64S12R2)')
s=s.replace('c!=Configuration::T64S12','(c!=Configuration::T64S12 && c!=Configuration::T64S12R2)')
s=s.replace('mConfig!=Configuration::T64S12','(mConfig!=Configuration::T64S12 && mConfig!=Configuration::T64S12R2)')
s=edit(s,'case Configuration::T64S12:return 64;','case Configuration::T64S12:case Configuration::T64S12R2:return 64;')
s=edit(s,'case Configuration::T64S12:return 12;','case Configuration::T64S12:case Configuration::T64S12R2:return 12;')
s=edit(s,'    if(mConfig==Configuration::T64S12) return Map64S12::name;',
    '    if(mConfig==Configuration::T64S12R2) return "imt_t64_s12_subspace_v1_r2";\n'
    '    if(mConfig==Configuration::T64S12) return Map64S12::name;')
s=edit(s,'#include "K16Inner.h"','#include "K16Inner.h"\n#include "K16R2Inner.h"')
f=edit(f,'#include "K16Inner.h"','#include "K16Inner.h"\n#include "K16R2Inner.h"')

# Both setups use a flat array of (u,v) pairs, grouped by epoch then round.
begin=s.index('template<class Map> void Spin::setupInner')
end=s.index('Spin::Workspace::Workspace',begin)
setup=s[begin:end]
setup=edit(setup,'const auto epochs=codeBlocks()/Map::T;',
    'const auto epochs=codeBlocks()/Map::T;\n    constexpr unsigned rounds=isTwoRoundMap<Map>?2:1;')
setup=setup.replace('resize(2*epochs)','resize(2*rounds*epochs)')
assert setup.count('for(std::size_t e=0;e<epochs;++e)')==2
setup=setup.replace('for(std::size_t e=0;e<epochs;++e)','for(std::size_t e=0;e<rounds*epochs;++e)',1)
s=s[:begin]+setup+s[end:]

if bidir:
    s=edit(s,'void Spin::validateSetup() const {', '''void Spin::validateSetup() const {
    if(!mCompacted && (mConfig==Configuration::T64S12 || mConfig==Configuration::T64S12R2)) {
        const auto count=2*(mConfig==Configuration::T64S12R2?2:1)*codeBlocks()/64;
        if(mCoefficients.size()!=count || mFieldRows.size()!=count || mForwardFieldRows.size()!=count)
            throw std::runtime_error("invalid IMT schedule length");
        for(std::size_t i=0;i<count;i+=2) {
            const auto u=mCoefficients[i],v=mCoefficients[i+1];
            if(!u || (u>>12) || (v>>12) || (std::popcount(u&v)&1) ||
               mFieldRows[i]!=u || mFieldRows[i+1]!=v || mForwardFieldRows[i]!=u || mForwardFieldRows[i+1]!=v)
                throw std::runtime_error("invalid IMT schedule sample");
        }
    }''')
    s=edit(s,'case Configuration::T64S12: setupInner<Map64S12>(coefficientSeed); break;',
        'case Configuration::T64S12R2: setupInner<Map64S12R2>(coefficientSeed); break;\n'
        '        case Configuration::T64S12: setupInner<Map64S12>(coefficientSeed); break;')
    for macro in ('RUN_CASE','FORWARD_CASE'):
        s=edit(s,f'{macro}(T64S12,Map64S12);',f'{macro}(T64S12R2,Map64S12R2); {macro}(T64S12,Map64S12);')
    for method in ('oracle','forwardOracle'):
        anchor=f'case Configuration::T64S12:{method}<Map64S12>(in,out);break;'
        s=edit(s,anchor,f'case Configuration::T64S12R2:{method}<Map64S12R2>(in,out);break;\n        '+anchor)
    # Both dense oracles use rows of the transpose, but apply them differently.
    anchor='if constexpr(isImtMap<Map>) mask=(1U<<j)^(((mCoefficients[2*e+1]>>j)&1)?mCoefficients[2*e]:0);'
    assert s.count(anchor)==2
    s=s.replace(anchor,'if constexpr(isTwoRoundMap<Map>) mask=imtR2ReferenceRow(j,mCoefficients.data()+4*e);\n            else '+anchor)
    s=edit(s,'    if(mConfig==Configuration::T64S12) forwardBitsMap<Map64S12>(in,out,scratch);',
        '    if(mConfig==Configuration::T64S12R2) forwardBitsMap<Map64S12R2>(in,out,scratch);\n'
        '    else if(mConfig==Configuration::T64S12) forwardBitsMap<Map64S12>(in,out,scratch);')
    anchor='''        u32 next=0;
        for(unsigned j=0;j<S;++j) {'''
    s=edit(s,anchor,'''        u32 mixed=state;
        if constexpr(isTwoRoundMap<Map>) {
            const auto* pair=mForwardFieldRows.data()+4*(base/T);
            mixed^=(0U-u32(std::popcount(mixed&pair[1])&1))&pair[0];
            mixed^=(0U-u32(std::popcount(mixed&pair[3])&1))&pair[2];
        }
        u32 next=0;
        for(unsigned j=0;j<S;++j) {''')
    anchor='''int((state>>j)&1)^int(((mForwardFieldRows[2*(base/T)]>>j)&1) &
                (std::popcount(state&mForwardFieldRows[2*(base/T)+1])&1))'''
    s=edit(s,anchor,'''([&] {if constexpr(isTwoRoundMap<Map>) return int((mixed>>j)&1);
                else return int((state>>j)&1)^int(((mForwardFieldRows[2*(base/T)]>>j)&1) &
                (std::popcount(state&mForwardFieldRows[2*(base/T)+1])&1));}())''')
    inner=(p/'Inner.h').read_text().replace('#include "Map64S12.h"','#include "ImtRounds.h"')
    inner=edit(inner,'std::is_same_v<Map,Map64S12>;','std::is_same_v<Map,Map64S12> || isTwoRoundMap<Map>;')
    inner=inner.replace('if constexpr(std::is_same_v<Map,Map64S12>)','if constexpr(std::is_same_v<Map,Map64S12> || isTwoRoundMap<Map>)')
    for first,second in [('v','u'),('u','v')]:
        anchor=('            auto v=fieldRows[2*epoch+1];auto dot' if first=='v'
                else '            auto u=fieldRows[2*epoch];auto dot')
        start=inner.index(anchor)
        end=inner.index('            for(unsigned j=0;j<Map::S;++j) state[j]',start)
        original=inner[start:end]
        calls=('imtStep(state,fieldRows[4*epoch+1],fieldRows[4*epoch]);\n'
               '                imtStep(state,fieldRows[4*epoch+3],fieldRows[4*epoch+2]);') if first=='v' else (
               'imtStep(state,fieldRows[4*epoch+2],fieldRows[4*epoch+3]);\n'
               '                imtStep(state,fieldRows[4*epoch],fieldRows[4*epoch+1]);')
        inner=inner[:start]+'            if constexpr(isTwoRoundMap<Map>) {\n                '+calls+'\n            } else {\n'+original+'            }\n'+inner[end:]
    put('Inner.h',inner)
    # Duplicate only the dispatch wrapper; shared map and kernel circuits stay shared.
    start=f.index('template<bool Packed> void Spin::runK16Four(')
    end=f.index('template void Spin::runK16Four',start)
    new=f[start:end].replace('runK16Four','runK16R2Four').replace('inner64Reverse','inner64R2Reverse')
    new+='template void Spin::runK16R2Four<true>(const block*,block*,Workspace&) const;\n'
    new+='template void Spin::runK16R2Four<false>(const block*,block*,Workspace&) const;\n'
    f=f[:-2]+new+'}\n'
    h=edit(h,'    template<bool Packed> void runK16Four',
        '    template<bool Packed> void runK16R2Four(const block*,block*,Workspace&) const;\n    template<bool Packed> void runK16Four')
    dispatch='if(layout==Layout::Packed24) runK16R2Four<true>(in,out,w); else runK16R2Four<false>(in,out,w);return;'
else:
    s=edit(s,'    if((c==Configuration::T64S12 || c==Configuration::T64S12R2)) setupInner<Map64S12Transpose>(coefficientSeed);',
        '    if(c==Configuration::T64S12R2) setupInner<Map64S12R2Transpose>(coefficientSeed);\n'
        '    else if(c==Configuration::T64S12) setupInner<Map64S12Transpose>(coefficientSeed);')
    for name in ('s','f'):
        text=s if name=='s' else f
        text=edit(text,'    if constexpr(Map::T==64 && Map::S==12) inner64Reverse(in,n,masks,emit);',
            '    if constexpr(isTwoRoundMap<Map>) inner64R2Reverse(in,n,masks,emit);\n'
            '    else if constexpr(Map::T==64 && Map::S==12) inner64Reverse(in,n,masks,emit);')
        if name=='s': s=text
        else: f=text
    s=edit(s,'    if(mConfig==Configuration::T64S12) oracle<Map64S12Transpose,false>(in,out);',
        '    if(mConfig==Configuration::T64S12R2) oracle<Map64S12R2Transpose,false>(in,out);\n'
        '    else if(mConfig==Configuration::T64S12) oracle<Map64S12Transpose,false>(in,out);')
    s=edit(s,'                if((v>>j)&1) mask^=u;',
        '                if((v>>j)&1) mask^=u;\n'
        '                if constexpr(isTwoRoundMap<Map>) mask=imtR2ReferenceRow(j,mCoefficients.data()+4*e);')
    s=edit(s,'for(std::size_t e=0;e<codeBlocks()/step();++e) {',
        'for(std::size_t e=0;e<(mConfig==Configuration::T64S12R2?2:1)*codeBlocks()/step();++e) {')
    s=s.replace('(u>>19) || (v>>19)','(u>>state()) || (v>>state())')
    start=s.index('    if(mConfig==Configuration::T64S12) {\n',s.index('void Spin::encodeUnchecked'))
    end=s.index('    if(layout==Layout::Packed24) run<AsymmetricMap',start)
    s=s[:start]+s[start:end].replace('T64S12','T64S12R2').replace('Map64S12Transpose','Map64S12R2Transpose')+s[start:]
    f=f[:-2]+'''template void Spin::runFour<Map64S12R2Transpose,true>(const block*,block*,Workspace&) const;
template void Spin::runFour<Map64S12R2Transpose,false>(const block*,block*,Workspace&) const;
}
'''
    dispatch='if(layout==Layout::Packed24) runFour<Map64S12R2Transpose,true>(in,out,w); else runFour<Map64S12R2Transpose,false>(in,out,w);return;'

if 'void runSmallK16(' in h:
    h=edit(h,'    void runSmallK16(', '    void runSmallK16R2(const block*,block*,Workspace&) const;\n    void runSmallK16(')
    start=f.index('void Spin::runSmallK16(')
    end=f.index('        return;\n    }',start)+len('        return;\n    }')
    f=f[:-2]+f[start:end].replace('runSmallK16','runSmallK16R2').replace('inner64Reverse','inner64R2Reverse')+'\n}\n'
    dispatch='runSmallK16R2(in,out,w);return;'
anchor='    if(mBchBackend==BchBackend::Avx512) {\n'
idx=s.index(anchor,s.index('void Spin::encodeUnchecked'))+len(anchor)
s=s[:idx]+'        if(mConfig==Configuration::T64S12R2) {'+dispatch+'}\n'+s[idx:]
put('Spin.h',h);put('Spin.cpp',s);put('Fast.cpp',f)
print('Integrated explicit K16 T64S12R2 option; defaults unchanged.')
