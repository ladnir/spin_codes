"""Import the exact frozen mode-52 circuits into the independent core module.

This development-only generator reads the pinned research sources. The library
sources it emits have no research include or build dependency.
"""
from pathlib import Path
import hashlib
import re
import sys

here = Path(__file__).resolve().parent
repo = here.parents[2]
research = repo / 'spin/experiments/k16_codesign_100us'
check = '--check' in sys.argv
def source_hash(name):
    # Hash the LF-normalized source so Windows Git checkouts remain reproducible.
    return hashlib.sha256((research/name).read_text(encoding='utf-8').encode('utf-8')).hexdigest()
def emit(name, text):
    path = here/name
    expected = text.encode('utf-8')
    if check:
        if not path.exists() or path.read_text(encoding='utf-8') != text:
            raise SystemExit(f'{name} is stale; run src/paired15/import_frozen.py')
    else:
        path.write_bytes(expected)
fold = (research / 'kernel/K16Paired15Fold.cpp').read_text()
common = fold[fold.index('struct Four {'):fold.index('struct ByteRouteWide {')]
shuffle = (research / 'kernel/T64Paired15ShuffleMap.h').read_text()
shuffle = shuffle[shuffle.index('template<bool RawFeedback'):]
shuffle = shuffle[:shuffle.rfind('}')]
shuffle = shuffle.replace('template<bool RawFeedback=false,class Emit>', 'template<bool Gather,bool RawFeedback=false,class Emit>')
shuffle = shuffle.replace('const Block* raw,std::size_t packetBase', 'const Block* raw,const std::uint32_t* route,std::size_t packetBase')
shuffle = re.sub(r'raw\+(\d+)', lambda m: 'raw+(Gather?route[%d]:%s)' % (int(m[1])//4,m[1]), shuffle)
emit('Paired15FastCommon.h', '''// Imported exact mode-52 paired-s15 state/feedback circuits.
// Regenerate with src/paired15/import_frozen.py.
''' + '// K16Paired15Fold.cpp SHA256: ' + source_hash('kernel/K16Paired15Fold.cpp') + '\n'
    + '// T64Paired15ShuffleMap.h SHA256: ' + source_hash('kernel/T64Paired15ShuffleMap.h') + '\n' + '''
#pragma once
#include "Paired15.h"
#include "../packet/PacketInnerFast.h"
namespace spin::detail::paired15::fast {
using namespace detail::packet::fast;
''' + common + shuffle + '\n}\n')

native = (research / 'kernel/NativeOuterFinish.h').read_text()
native = native.replace('#include "K16CodeSign.h"', '#include "Paired15.h"')
native = native.replace('spin::research::k16codesign::nativeouter', 'spin::detail::paired15::nativeouter')
emit('NativeOuterFinish.h', '// Frozen NativeOuterFinish.h SHA256: ' + source_hash('kernel/NativeOuterFinish.h') + '\n' + native)

outer = (research / 'outer_variants/OuterVariants.cpp').read_text()
outer = outer[outer.index('using Coeff ='):outer.index('template<std::size_t... Symbol>')]
outer = outer.replace('rs::tower32byte::quadraticMultiply', 'quadraticMultiply')
field = '''struct Pair { __m512i lo,hi; };
static SPIN_FORCEINLINE Pair quadraticMultiply(__m512i x0,__m512i x1,const std::uint8_t* c) {
    const auto p0=_mm512_gf2p8mul_epi8(x0,_mm512_set1_epi8(char(c[0])));
    const auto p1=_mm512_gf2p8mul_epi8(x1,_mm512_set1_epi8(char(c[1])));
    const auto p2=_mm512_gf2p8mul_epi8(_mm512_xor_si512(x0,x1),_mm512_set1_epi8(char(c[2])));
    return {_mm512_xor_si512(p0,p1),_mm512_xor_si512(p0,p2)};
}
'''
emit('Paired15TransposeOuter.h', '''// Imported selected fieldLoopSharedParity circuit (mode 52).
''' + '// OuterVariants.cpp SHA256: ' + source_hash('outer_variants/OuterVariants.cpp') + '\n' + '''
#pragma once
#include "NativeOuterFinish.h"
#include "../packet/PacketLargeInner.h"
namespace spin::detail::paired15::outerfast {
''' + field + outer + '''
static SPIN_NOINLINE void group(const Block* __restrict in,Block* __restrict out,
                               const Coeff* __restrict coeff) {
    alignas(64) __m512i packed[64];
    // Preserve the measured compact symbol loop; full unrolling lost.
    for(unsigned symbol=0;symbol<16;++symbol) mix(in+16*symbol,packed+4*symbol,coeff[symbol]);
    finish<0,true>(packed,out);finish<2,true>(packed,out);
}
inline void run(const Block* in,Block* out,const Plan& p) {
    for(std::size_t g=0;g<p.groups;++g)
        group(in+groupStride*g,out+128*g,p.outerField.data()+16*g);
}
}
''')
print(('Verified' if check else 'Imported') + ' selected mode-52 state and outer circuits.')
