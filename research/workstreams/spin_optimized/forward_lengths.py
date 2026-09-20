"""Natural lengths and 32-bit wide routing, preserving complete-tile kernels."""
from pathlib import Path
import re
import sys

p = Path(sys.argv[1])


def edit(text, old, new):
    if text.count(old) != 1:
        raise RuntimeError(f"Forward length anchor changed: {old}")
    return text.replace(old, new)


def function(text, signature):
    start = text.index(signature)
    end = text.index('{', start) + 1
    depth = 1
    while depth:
        depth += (text[end] == '{') - (text[end] == '}')
        end += 1
    return text[start:end]


def tail(kernel, name):
    kernel = edit(kernel, f'Spin::{name}(', f'Spin::{name}Tail(')
    kernel = edit(kernel, 'for(std::size_t base=0;base<n;base+=tileSize) {',
                  'for(std::size_t base=0;base<n;base+=tileSize) {\n'
                  '        const auto activeSize=std::min(tileSize,n-base);')
    return kernel.replace('<tileSize', '<activeSize')


def put(name, text):
    (p/name).write_text(text, newline='\n')


h = (p/'Spin.h').read_text()
s = (p/'Spin.cpp').read_text()
(p/'LengthGeometry.h').write_bytes((Path(__file__).parent/'LengthGeometry.h').read_bytes())
h = '#include "LengthGeometry.h"\n' + h
h = edit(h, 'enum class Layout { Packed24, Indices32 };',
         'struct MessageLength {std::size_t value;};\nenum class Layout { Packed24, Indices32, Auto };')
start = h.index('    Spin(Configuration')
end = h.index(';', start) + 1
h = h[:end] + '\n' + h[start:end].replace('unsigned messageExponent', 'MessageLength length') + h[end:]
h = h.replace('Layout layout=Layout::Packed24', 'Layout layout=Layout::Auto')
h = edit(h, '    std::size_t messageBlocks() const noexcept', '''    static constexpr std::size_t maxMessageBlocks=length_geometry::maxMessageBlocks;
    bool packed24Available() const noexcept {return codeBlocks()<=(std::size_t{1}<<24);}
    Layout preferredLayout() const noexcept {
        return mCompacted?mRetainedLayout:(packed24Available()?Layout::Packed24:Layout::Indices32);
    }
    std::size_t messageBlocks() const noexcept''')
h = edit(h, '    unsigned mTileRows;', '    unsigned mTileRows;\n    bool mPartialTile=false;')
signature = 'Spin::Spin(Configuration c,unsigned exponent,u64 routeSeed,u64 coefficientSeed,unsigned tileRows,BchBackend backend)'
s = edit(s, signature, signature + '''
    :Spin(c,MessageLength{length_geometry::fromExponent(exponent)},routeSeed,coefficientSeed,tileRows,backend) {}
Spin::Spin(Configuration c,MessageLength length,u64 routeSeed,u64 coefficientSeed,unsigned tileRows,BchBackend backend)''')
s = edit(s, '    const unsigned minimum=c==Configuration::T128S19?14:16;\n', '')
s = edit(s, 'exponent<minimum || exponent>20 || ', '')
s = s.replace('supported message exponents: S19 14..20, others 16..20; ', '')
s, count = re.subn(r'    if\([^\n]*exponent!=16\)\n        throw [^\n]*;\n', '', s)
if count != 1:
    raise RuntimeError('K16-only restriction changed')
s = edit(s, '    mK=std::size_t{1}<<exponent;', '    mK=length.value;\n    length_geometry::check(mK,step());')
s = s.replace('exponent<=18', 'mK<=(1U<<18)')
s = edit(s, 'std::min<std::size_t>(selectedTile,rows)', 'std::min<std::size_t>(selectedTile,std::bit_floor(rows))')
s = edit(s, '    if(rows%step())', '    mPartialTile=(rows%mTileRows)!=0;\n    if(rows%step())')
s = edit(s, 'mBchOffsets24.resize(3*n+4);mBchOffsets32.resize(n);',
         'if(packed24Available()) mBchOffsets24.resize(3*n+4);mBchOffsets32.resize(n);')
s = edit(s, '    mSlots24.resize(3*n+4); mOffsets24.resize(3*n+4);',
         '    if(packed24Available()) {mSlots24.resize(3*n+4); mOffsets24.resize(3*n+4);}')
s = edit(s, '    std::vector<u32> counts(n/tile,u32(tile));',
         '    std::vector<u32> counts((n+tile-1)/tile,u32(tile));\n    if(n%tile) counts.back()=u32(n%tile);')
s = edit(s, '        pack(mSlots24.data()+3*inner,slot); pack(mOffsets24.data()+3*slot,local);',
         '        if(packed24Available()) {pack(mSlots24.data()+3*inner,slot); pack(mOffsets24.data()+3*slot,local);}')
s = edit(s, 'pack(mBchOffsets24.data()+3*slot,packFour(local));',
         'if(packed24Available()) pack(mBchOffsets24.data()+3*slot,packFour(local));')
s = edit(s, 'unpack(mBchOffsets24.data()+3*slot)!=mBchOffsets32[slot]',
         '(packed24Available() && unpack(mBchOffsets24.data()+3*slot)!=mBchOffsets32[slot])')
s = edit(s, 'if(unpack(mSlots24.data()+3*i)!=slot || unpack(mOffsets24.data()+3*slot)!=mOffsets32[slot] ||',
         'if((packed24Available() && (unpack(mSlots24.data()+3*i)!=slot || unpack(mOffsets24.data()+3*slot)!=mOffsets32[slot])) ||')
for method in ('encode', 'encodeInplace', 'encodeUnchecked', 'forward', 'forwardUnchecked', 'compact'):
    start = s.index(f'void Spin::{method}(')
    brace = s.index('{', start) + 1
    extra = '\n    if(layout==Layout::Auto) layout=preferredLayout();\n'
    if 'Unchecked' not in method:
        extra += '    if(layout==Layout::Packed24 && !packed24Available()) throw std::invalid_argument("Packed24 index range exceeded; use Auto or Indices32");\n'
    s = s[:brace] + extra + s[brace:]

# Generate separate tail versions: no per-element min or tail branch in full tiles.
for name in ('run', 'runForward'):
    kernel = function(s, f'template<class Map,bool Packed> void Spin::{name}(')
    s = edit(s, kernel, kernel + '\n' + tail(kernel, name))
    h = edit(h, '    template<class Map> void setupInner',
             f'    template<class Map,bool Packed> void {name}Tail(const block*,block*,Workspace&) const;\n'
             '    template<class Map> void setupInner')

fast_names = []
for filename in ('Fast.cpp', 'ForwardFast.cpp'):
    if not (p/filename).exists():
        continue
    if filename == 'ForwardFast.cpp' and 'void runForwardFour(' not in h:
        continue
    text = (p/filename).read_text()
    names = ('runFour', 'runK16Four', 'runK16R2Four') if filename == 'Fast.cpp' else ('runForwardFour',)
    for name in names:
        match = re.search(r'template<[^\n]+> void Spin::' + name + r'\(', text)
        if not match:
            continue
        kernel = function(text, match[0])
        text = edit(text, kernel, kernel + '\n' + tail(kernel, name))
        prototype = re.search(r'    template<[^\n]+> void ' + name + r'\([^\n]+;', h)[0]
        h = edit(h, prototype, prototype + '\n' + prototype.replace(name+'(', name+'Tail('))
        inst = re.findall(r'template void Spin::' + name + r'<[^\n]+;', text)
        if not inst:
            raise RuntimeError(f'Missing instantiations for {name}')
        text = text[:-2] + '\n'.join(x.replace(name+'<', name+'Tail<') for x in inst) + '\n}\n'
        fast_names.append(name)
    put(filename, '#include <algorithm>\n' + text)

# Only the exact-size transpose direct kernels may assume fixed n.
for suffix in ('', 'R2'):
    old = f'if(mConfig==Configuration::T64S12{suffix}) {{runSmallK16{suffix}(in,out,w);return;}}'
    if old in s:
        s = edit(s, old, old.replace(')', ' && mK==(1U<<16))', 1))
for exp in (18, 20):
    s = s.replace(f'if(mK==(1U<<{exp})) {{runK{exp}',
                  f'if(mConfig==Configuration::T128S19 && mK==(1U<<{exp})) {{runK{exp}')
for cfg, name in (('T64S12', 'runK16Four'), ('T64S12R2', 'runK16R2Four')):
    if name in fast_names:
        anchor = '        if(layout==Layout::Packed24) runFour<true>(in,out,w); else runFour<false>(in,out,w);'
        s = edit(s, anchor, f'        if(mConfig==Configuration::{cfg}) {{if(layout==Layout::Packed24) {name}<true>(in,out,w);else {name}<false>(in,out,w);return;}}\n'+anchor)

maps = [('T64S16','Map64S16'), ('T64S20','Map64S20'), ('T128S19','Map128S19'),
        ('T256S14','Map256S14'), ('T64S12','Map64S12'), ('T64S12R2','Map64S12R2')]
for method, name in (('encodeUnchecked', 'run'), ('forwardUnchecked', 'runForward')):
    partial = 'mPartialTile'
    if method == 'forwardUnchecked' and 'mForwardDirect' in h:
        partial += ' && mForwardDirect.empty()'  # Direct routing has no tile tail.
    dispatch = f'    if({partial}) {{\n#if SPIN_BCH_AVX512\n        if(mBchBackend==BchBackend::Avx512) {{\n'
    for cfg, mapname in maps:
        fast = {'T128S19':'runFour','T64S12':'runK16Four','T64S12R2':'runK16R2Four'}.get(cfg)
        if method == 'forwardUnchecked':
            fast = 'runForwardFour' if fast else None
        if fast not in fast_names:
            continue
        args = f'{mapname},' if fast == 'runForwardFour' else ''
        dispatch += f'            if(mConfig==Configuration::{cfg}) {{if(layout==Layout::Packed24) {fast}Tail<{args}true>(in,out,w);else {fast}Tail<{args}false>(in,out,w);return;}}\n'
    dispatch += '        }\n#endif\n        switch(mConfig) {\n'
    for cfg, mapname in maps:
        dispatch += f'            case Configuration::{cfg}:if(layout==Layout::Packed24) {name}Tail<{mapname},true>(in,out,w);else {name}Tail<{mapname},false>(in,out,w);return;\n'
    dispatch += '        }\n    }\n'
    start = s.index(f'void Spin::{method}(')
    loc = s.index('    if(layout==Layout::Auto) layout=preferredLayout();', start)
    loc = s.index('\n', loc) + 1
    s = s[:loc] + dispatch + s[loc:]

if 'wideForwardView()' in h:
    h = edit(h, '        Configuration configuration=Configuration::T128S19;',
             '        Configuration configuration=Configuration::T128S19;\n'
             '        const u32* slots32=nullptr;\n        const u32* offsets32=nullptr;')
    s = edit(s, ' || mSlots24.empty())\n        throw std::invalid_argument("wide forward requires S19 or K16 R2 and Packed24");',
             ' || (mSlots24.empty() && mSlots32.empty()))\n        throw std::invalid_argument("wide forward requires S19 or K16 R2 routing");')
    s = edit(s, 'return {codeBlocks(),tileBlocks(),mSlots24.data(),mOffsets24.data(),mForwardFieldRows.data(),mConfig};',
             'return {codeBlocks(),tileBlocks(),mSlots24.empty()?nullptr:mSlots24.data(),mOffsets24.data(),mForwardFieldRows.data(),mConfig,mSlots32.data(),mOffsets32.data()};')
    text = (p/'wide/WideKernel.h').read_text()
    text = edit(text, 'template<class Map> static void encode(', 'template<class Map,bool Packed> static void encode(')
    text = edit(text, 'tile[unpack(view.offsets+3*(base+j))]',
                'tile[Packed?unpack(view.offsets+3*(base+j)):view.offsets32[base+j]]')
    text = edit(text, 'values[unpack(view.slots+3*i)]', 'values[Packed?unpack(view.slots+3*i):view.slots32[i]]')
    kernel = function(text, 'template<class Map,bool Packed> static void encode(')
    extra = edit(kernel, 'static void encode(', 'static void encodeTail(')
    extra = edit(extra, 'for(std::size_t base=0;base<view.n;base+=view.tile) {',
                 'for(std::size_t base=0;base<view.n;base+=view.tile) {\n        const auto activeSize=std::min(view.tile,view.n-base);')
    extra = extra.replace('j<view.tile', 'j<activeSize')
    wrapper = '''template<class Map> static void dispatch(Spin::WideView v,const Vec* in,Vec* out,Workspace& w) {
    if(v.n%v.tile) {
        if(v.slots) encodeTail<Map,true>(v,in,out,w);else encodeTail<Map,false>(v,in,out,w);
    } else {
        if(v.slots) encode<Map,true>(v,in,out,w);else encode<Map,false>(v,in,out,w);
    }
}
'''
    text = edit(text, kernel, kernel+'\n'+extra+'\n'+wrapper)
    text = text.replace('HC_NS::encode<', 'HC_NS::dispatch<')
    put('wide/WideKernel.h', '#include <algorithm>\n'+text)
put('Spin.h', h)
put('Spin.cpp', s)
print('Forward and wide natural lengths; Auto/32-bit routing; separate tail kernels.')
