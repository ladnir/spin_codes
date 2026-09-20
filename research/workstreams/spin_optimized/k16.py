"""Add an explicit certified K16 option to generated encoders, preserving defaults."""
from pathlib import Path
import sys
from k16_map import header

p = Path(sys.argv[1])
here = Path(__file__).resolve().parent
def edit(s, a, b):
    if s.count(a) != 1:
        raise RuntimeError(f'K16 integration anchor changed: {a}')
    return s.replace(a, b)
def put(n, t):
    (p/n).write_text(t, newline='\n')

put('Map64S12.h', header())
put('K16Inner.h', (here/'K16Inner.h').read_text())
h = (p/'Spin.h').read_text()
bidir = 'forwardBits(' in h
h = edit(h, 'T128S19, T256S14 };', 'T128S19, T256S14, T64S12 };')
s = (p/'Spin.cpp').read_text()
f = (p/'Fast.cpp').read_text()
for text_name in ('s', 'f'):
    text = s if text_name == 's' else f
    text = edit(text, '#include "Inner.h"', '#include "Inner.h"\n#include "K16Inner.h"')
    if text_name == 's': s = text
    else: f = text
s = edit(s, '    mK=std::size_t{1}<<exponent;',
         '    if(c==Configuration::T64S12 && exponent!=16)\n'
         '        throw std::invalid_argument("T64S12 is certified and supported only at K=2^16");\n'
         '    mK=std::size_t{1}<<exponent;')
s = edit(s, 'static_cast<unsigned>(Configuration::T256S14)', 'static_cast<unsigned>(Configuration::T64S12)')
s = edit(s, 'case Configuration::T64S20:return 64;', 'case Configuration::T64S20:case Configuration::T64S12:return 64;')
s = edit(s, 'case Configuration::T64S16:return 16;', 'case Configuration::T64S12:return 12; case Configuration::T64S16:return 16;')
s = edit(s, 'const char* Spin::name() const noexcept {',
         'const char* Spin::name() const noexcept {\n'
         '    if(mConfig==Configuration::T64S12) return Map64S12::name;')

if bidir:
    h = edit(h, '    template<bool Packed> void runFour',
             '    template<bool Packed> void runK16Four(const block*,block*,Workspace&) const;\n'
             '    template<bool Packed> void runFour')
    h = edit(h, '    template<class Map> void setupInner',
             '    template<class Map> void forwardBitsMap(const u64*,u64*,u64*) const;\n'
             '    template<class Map> void setupInner')
    s = edit(s, 'const bool four=c==Configuration::T128S19 &&',
             'const bool four=(c==Configuration::T128S19 || c==Configuration::T64S12) &&')
    s = edit(s, 'case Configuration::T128S19: setupInner<Map128S19>(coefficientSeed); break;',
             'case Configuration::T64S12: setupInner<Map64S12>(coefficientSeed); break;\n'
             '        case Configuration::T128S19: setupInner<Map128S19>(coefficientSeed); break;')
    # Explicit map trait; never infer the recurrence from state dimension.
    inner = (p/'Inner.h').read_text()
    inner = edit(inner, '#include <utility>', '#include <utility>\n#include <type_traits>\n#include "Map64S12.h"')
    inner = edit(inner, 'namespace bare_spin {', 'namespace bare_spin {\n'
                 'template<class Map> inline constexpr bool isImtMap=\n'
                 '    std::is_same_v<Map,Map128S19> || std::is_same_v<Map,Map64S12>;')
    inner = inner.replace('Map::S==19', 'isImtMap<Map>')
    inner = edit(inner, '            imtReversePoints<Map>(in+base,values,state,base,emit,std::make_index_sequence<Map::T>{});',
                 '            if constexpr(std::is_same_v<Map,Map64S12>) Map::emitShared(in+base,values,state,base,emit);\n'
                 '            else imtReversePoints<Map>(in+base,values,state,base,emit,std::make_index_sequence<Map::T>{});')
    put('Inner.h', inner)
    s = s.replace('Map::S==19', 'isImtMap<Map>').replace('Map::S!=19', '!isImtMap<Map>')
    s = s.replace('((1U<<19)-1)', '((1U<<Map::S)-1)')
    s = edit(s, 'if(mConfig!=Configuration::T128S19) std::vector<u32>().swap(mRoute);',
             'if(mConfig!=Configuration::T128S19 && mConfig!=Configuration::T64S12) std::vector<u32>().swap(mRoute);')
    for macro in ('RUN_CASE', 'FORWARD_CASE'):
        s = edit(s, f'{macro}(T128S19,Map128S19);', f'{macro}(T64S12,Map64S12); {macro}(T128S19,Map128S19);')
    for method in ('oracle', 'forwardOracle'):
        anchor = f'case Configuration::T128S19:{method}<Map128S19>(in,out);break;'
        s = edit(s, anchor, f'case Configuration::T64S12:{method}<Map64S12>(in,out);break;\n        '+anchor)
    start = f.index('template<bool Packed> void Spin::runFour(')
    end = f.index('template void Spin::runFour', start)
    new = f[start:end].replace('Spin::runFour(', 'Spin::runK16Four(').replace('using Map=Map128S19;', 'using Map=Map64S12;')
    new = new.replace('innerReverse<Map>', 'inner64Reverse')
    new += '\ntemplate void Spin::runK16Four<true>(const block*,block*,Workspace&) const;\n'
    new += 'template void Spin::runK16Four<false>(const block*,block*,Workspace&) const;\n'
    f = f[:-2]+new+'}\n'
    fast_dispatch = 'if(layout==Layout::Packed24) runK16Four<true>(in,out,w); else runK16Four<false>(in,out,w);'

    # Generalize the fixed-width bit-packed recurrence at compile time. Old T128
    # code remains a separate instantiation; T64 has one word and three nibbles.
    s = edit(s, 'if(mConfig!=Configuration::T128S19 || !in || !out || !scratch ||',
             'if((mConfig!=Configuration::T128S19 && mConfig!=Configuration::T64S12) || !in || !out || !scratch ||')
    marker = '    // 16 KiB of four-bit BCH tables.'
    idx = s.index(marker)
    s = s[:idx]+'''    if(mConfig==Configuration::T64S12) forwardBitsMap<Map64S12>(in,out,scratch);
    else forwardBitsMap<Map128S19>(in,out,scratch);
}
template<class Map> void Spin::forwardBitsMap(const u64* in,u64* out,u64* scratch) const {
    constexpr unsigned S=Map::S,T=Map::T,W=T/64,G=(S+3)/4;
'''+s[idx:]
    idx = s.index('template<class Map> void Spin::forwardBitsMap')
    bits = s[idx:]
    bits = bits.replace('std::array<u64,2>,19', 'std::array<u64,W>,S').replace('std::array<u64,2>,16>,5', 'std::array<u64,W>,16>,G')
    bits = bits.replace('j<19', 'j<S').replace('p<128', 'p<T').replace('g<5', 'g<G').replace('4*g+b<19', '4*g+b<S').replace('w<2', 'w<W')
    bits = bits.replace('Map128S19::', 'Map::').replace('base+=128', 'base+=T').replace('base+128', 'base+T').replace('base/128', 'base/T')
    bits = edit(bits, 'const auto a=mRoute[base+p],b=mRoute[base+64+p];', 'const auto a=mRoute[base+p];')
    bits = edit(bits, '            x1|=((scratch[b/64]>>(b%64))&1)<<p;',
                '            if constexpr(W==2) {const auto b=mRoute[base+64+p];x1|=((scratch[b/64]>>(b%64))&1)<<p;}')
    bits = edit(bits, 'y0^=t[0]; y1^=t[1];', 'y0^=t[0]; if constexpr(W==2) y1^=t[1];')
    bits = edit(bits, 'out[base/64]=y0; out[base/64+1]=y1;', 'out[base/64]=y0; if constexpr(W==2) out[base/64+1]=y1;')
    bits = edit(bits, 'std::popcount(x1&feedbackMasks[j][1])', '([&] {if constexpr(W==2) return std::popcount(x1&feedbackMasks[j][1]); else return 0;}())')
    s = s[:idx]+bits
else:
    s = edit(s, 'if(c!=Configuration::T128S19 || outer!=Outer::Bch256x128)',
             'if((c!=Configuration::T128S19 && c!=Configuration::T64S12) || outer!=Outer::Bch256x128)')
    s = edit(s, '    setupInner<AsymmetricMap>(coefficientSeed);',
             '    if(c==Configuration::T64S12) setupInner<Map64S12Transpose>(coefficientSeed);\n'
             '    else setupInner<AsymmetricMap>(coefficientSeed);')
    s = edit(s, 'if constexpr(Map::T==128 && Map::S==19)',
             'if constexpr((Map::T==128 && Map::S==19) || (Map::T==64 && Map::S==12))')
    # Both variants call exactly their original reverse kernel; selection folds
    # away at compile time, including in the ISA-specific translation unit.
    wrapper = '''
template<class Map,class Emit> OC_FORCEINLINE void selectedReverse(
    const block* in,std::size_t n,const u32* masks,Emit&& emit) {
    if constexpr(Map::T==64 && Map::S==12) inner64Reverse(in,n,masks,emit);
    else asymmetricReverse<Map,true,W5_MASKED>(in,n,masks,emit);
}
'''
    s = s.replace('asymmetricReverse<Map,true,W5_MASKED>', 'selectedReverse<Map>')
    f = f.replace('asymmetricReverse<Map,true,W5_MASKED>', 'selectedReverse<Map>')
    s = edit(s, 'namespace bare_spin {', 'namespace bare_spin {\n'+wrapper)
    f = edit(f, 'namespace bare_spin {', 'namespace bare_spin {\n'+wrapper)
    s = edit(s, '    oracle<AsymmetricMap,false>(in,out);',
             '    if(mConfig==Configuration::T64S12) oracle<Map64S12Transpose,false>(in,out);\n'
             '    else oracle<AsymmetricMap,false>(in,out);')
    s = edit(s, '    if(layout==Layout::Packed24) run<AsymmetricMap,true,false>(in,out,w);',
             '    if(mConfig==Configuration::T64S12) {\n'
             '        if(layout==Layout::Packed24) run<Map64S12Transpose,true,false>(in,out,w);\n'
             '        else run<Map64S12Transpose,false,false>(in,out,w); return;\n    }\n'
             '    if(layout==Layout::Packed24) run<AsymmetricMap,true,false>(in,out,w);')
    f = f[:-2]+'''template void Spin::runFour<Map64S12Transpose,true>(const block*,block*,Workspace&) const;
template void Spin::runFour<Map64S12Transpose,false>(const block*,block*,Workspace&) const;
}
'''
    fast_dispatch = 'if(layout==Layout::Packed24) runFour<Map64S12Transpose,true>(in,out,w); else runFour<Map64S12Transpose,false>(in,out,w);'

# The new configuration is dispatched once, not inside an epoch or route loop.
dispatch = '    if(mBchBackend==BchBackend::Avx512) {\n'
idx = s.index(dispatch, s.index('void Spin::encodeUnchecked'))+len(dispatch)
small_dispatch = ''
if 'void runSmall(' in h:
    h = edit(h, '    void runSmall(', '    void runSmallK16(const block*,block*,Workspace&) const;\n    void runSmall(')
    start = f.index('void Spin::runSmall(')
    end = f.index('        return;\n    }', start)+len('        return;\n    }')
    new = f[start:end].replace('Spin::runSmall(', 'Spin::runSmallK16(')
    new = new.replace('using Map=Map128S19;', 'using Map=Map64S12;').replace('using Map=AsymmetricMap;', 'using Map=Map64S12Transpose;')
    new = new.replace('innerReverse<Map>', 'inner64Reverse').replace('selectedReverse<Map>', 'inner64Reverse')
    f = f[:-2]+new+'\n}\n'
    small_dispatch = 'runSmallK16(in,out,w);return;'
s = s[:idx]+'        if(mConfig==Configuration::T64S12) {'+(small_dispatch or fast_dispatch+'return;')+'}\n'+s[idx:]
put('Spin.h', h)
put('Spin.cpp', s)
put('Fast.cpp', f)
print('Integrated explicit K16 T64S12 option; existing configurations unchanged.')
