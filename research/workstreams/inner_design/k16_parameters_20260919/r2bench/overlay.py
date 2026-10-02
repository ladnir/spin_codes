"""Mechanical overlay of a fresh standalone build; never changes source defaults."""
from pathlib import Path
import sys

p=Path(sys.argv[1])
def edit(text,a,b):
    assert text.count(a)==1,a
    return text.replace(a,b)

s=(p/'Spin.cpp').read_text()
assert 'forwardBits(' not in s
s=edit(s,'mCoefficients.resize(2*epochs);mFieldRows.resize(2*epochs);',
    'constexpr unsigned rounds=Map::S==12?2:1;\n'
    '            mCoefficients.resize(2*rounds*epochs);mFieldRows.resize(2*rounds*epochs);')
s=edit(s,'for(std::size_t e=0;e<epochs;++e) {\n                u32 u;',
    'for(std::size_t e=0;e<rounds*epochs;++e) {\n                u32 u;')
s=edit(s,'if(mConfig==Configuration::T64S12) return Map64S12::name;',
    'if(mConfig==Configuration::T64S12) return "imt_t64_s12_subspace_v1_r2";')
s=edit(s,'for(std::size_t e=0;e<codeBlocks()/step();++e) {',
    'for(std::size_t e=0;e<(mConfig==Configuration::T64S12?2:1)*codeBlocks()/step();++e) {')
s=edit(s,'(u>>19) || (v>>19)','(u>>state()) || (v>>state())')
anchor='''                  const u32 u=mCoefficients[2*e],v=mCoefficients[2*e+1];
                  mask=1U<<j;
                  if((v>>j)&1) mask^=u;'''
# Match indentation independently of historical formatting.
anchor='\n'.join(line[2:] for line in anchor.splitlines())
s=edit(s,anchor,'''                mask=1U<<j;
                if constexpr(Map::S==12) {
                    // Row j of T1^T T2^T. Independent dense oracle: form
                    // the composite row before taking its dot with state.
                    for(unsigned r=0;r<2;++r) {
                        const u32 u=mCoefficients[4*e+2*r],v=mCoefficients[4*e+2*r+1];
                        if(std::popcount(mask&v)&1) mask^=u;
                    }
                } else {
                    const u32 u=mCoefficients[2*e],v=mCoefficients[2*e+1];
                    if((v>>j)&1) mask^=u;
                }''')
(p/'Spin.cpp').write_text(s)
h=(p/'K16Inner.h').read_text()
start=h.index('            auto u=masks[2*epoch]')
end=h.index('            for(unsigned j=0;j<Map::S;++j) state[j]',start)
h=h[:start]+'''            // Forward update is T2 T1 state + Bx. Reverse applies
            // T2^T first, then T1^T, and adds the A^T input only once.
            auto u=masks[4*epoch+2];auto dot=_mm_setzero_si128();
            while(u) {const auto j=std::countr_zero(u);u&=u-1;dot=_mm_xor_si128(dot,state[j]);}
            auto v=masks[4*epoch+3];
            while(v) {const auto j=std::countr_zero(v);v&=v-1;state[j]=_mm_xor_si128(state[j],dot);}
            u=masks[4*epoch];dot=_mm_setzero_si128();
            while(u) {const auto j=std::countr_zero(u);u&=u-1;dot=_mm_xor_si128(dot,state[j]);}
            v=masks[4*epoch+1];
            while(v) {const auto j=std::countr_zero(v);v&=v-1;state[j]=_mm_xor_si128(state[j],dot);}
'''+h[end:]
(p/'K16Inner.h').write_text(h)
print('Isolated t64s12r2 setup, dense oracle, and reverse kernel generated.')
