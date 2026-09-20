"""Build the integrated half-rate encoder without changing submission source pins."""
from pathlib import Path
import importlib.util
import re
import sys

ROOT = Path(__file__).resolve().parents[2]
BARE = ROOT / 'workstreams/bare_bch_rm2sub'
BASE = ROOT / 'workstreams/inner_design/asymmetric/bch256/weight5/implementation'
OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
(OUT / 'generated').mkdir(exist_ok=True)

def put(name, text):
    (OUT / name).write_text(text, newline='\n')

def replace(text, old, new):
    if text.count(old) != 1:
        raise RuntimeError(f'upstream integration anchor changed: {old}')
    return text.replace(old, new)

header = (BARE / 'Spin.h').read_text()
header = replace(header, 'enum class Outer',
    'enum class BchBackend { Auto, Avx2, Avx512 };\n'
    '// Auto falls back safely; explicit Avx512 rejects unsupported CPUs/builds/tiles.\n'
    'bool bchAvx512Available() noexcept;\nenum class Outer')
header = replace(header, 'Outer outer=Outer::Bch256x128);',
    'Outer outer=Outer::Bch256x128, BchBackend backend=BchBackend::Auto);')
header = replace(header, '    const char* name() const noexcept;',
    '    BchBackend bchBackend() const noexcept { return mBchBackend; }\n'
    '    const char* name() const noexcept;')
header = replace(header, '    bool mCompacted=false;',
    '    BchBackend mBchBackend=BchBackend::Avx2;\n    bool mCompacted=false;')
header = replace(header, 'bool Packed, bool Quarter> void run',
    'bool Packed, bool Quarter, bool Four=false> void run')
header = replace(header, '    template<class Map> void setupInner',
    '    template<class Map,bool Packed> void runFour(const block*,block*,Workspace&) const;\n'
    '    template<class Map> void setupInner')
put('Spin.h', header)

spin = (BASE / 'Weight5Spin.cpp').read_text()
spin = replace(spin, 'namespace bare_spin {', '''namespace bare_spin {
void bchTranspose4(const block*,block*);
bool bchAvx512Available() noexcept {
#if SPIN_BCH_AVX512 && !SPIN_TEST_NO_AVX512
    // GCC/Clang also check OS support for the extended register state.
    return __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512vl");
#else
    return false;
#endif
}
static u32 packFour(u32 x) noexcept {
    return (x&~1023U)|((x&255U)<<2)|((x>>8)&3U);
}
static u32 unpackFour(u32 x) noexcept {
    return (x&~1023U)|((x&3U)<<8)|((x>>2)&255U);
}''')
spin = replace(spin, 'unsigned tileRows,Outer outer)', 'unsigned tileRows,Outer outer,BchBackend backend)')
spin = replace(spin, '    if(rows%step())', '''    if(backend!=BchBackend::Auto && backend!=BchBackend::Avx2 && backend!=BchBackend::Avx512)
        throw std::invalid_argument("unknown BCH backend");
    const bool four=mTileRows>=4 && bchAvx512Available();
    if(backend==BchBackend::Avx512 && !four)
        throw std::invalid_argument("AVX-512 BCH requires an enabled build, supported CPU, and >=4 tile rows");
    mBchBackend=(backend!=BchBackend::Avx2 && four)?BchBackend::Avx512:BchBackend::Avx2;
    if(rows%step())''')
spin = replace(spin, 'mSlots32[inner]=slot; mOffsets32[slot]=local;\n        pack(mSlots24.data()+3*inner,slot); pack(mOffsets24.data()+3*slot,local);',
    'const u32 offset=mBchBackend==BchBackend::Avx512?packFour(local):local;\n'
    '        mSlots32[inner]=slot; mOffsets32[slot]=offset;\n'
    '        pack(mSlots24.data()+3*inner,slot); pack(mOffsets24.data()+3*slot,offset);')
spin = replace(spin, 'template<class Map,bool Packed,bool Quarter> void Spin::run',
    'template<class Map,bool Packed,bool Quarter,bool Four> void Spin::run')
spin = replace(spin, '        if constexpr(Quarter) {\n            for(std::size_t j=0;j<tileSize;j+=256)',
    '        if constexpr(Four) {\n'
    '            for(std::size_t j=0;j<tileSize;j+=1024)\n'
    '                bchTranspose4(tile+j,out+base/2+j/2);\n'
    '        } else if constexpr(Quarter) {\n            for(std::size_t j=0;j<tileSize;j+=256)')
spin = replace(spin, '    if(layout==Layout::Packed24) run<AsymmetricMap,true,false>(in,out,w);',
    '''#if SPIN_BCH_AVX512
    if(mBchBackend==BchBackend::Avx512) {
        if(layout==Layout::Packed24) runFour<AsymmetricMap,true>(in,out,w);
        else runFour<AsymmetricMap,false>(in,out,w);
        return;
    }
#endif
    if(layout==Layout::Packed24) run<AsymmetricMap,true,false>(in,out,w);''')
spin = replace(spin, 'mOffsets32[slot]!=outer',
    '(mBchBackend==BchBackend::Avx512?unpackFour(mOffsets32[slot]):mOffsets32[slot])!=outer')
put('Spin.cpp', spin)
# Instantiate the same complete hot loop in its own ISA-specific translation unit.
# Setup and AVX2 instantiations remain in Spin.cpp; the constructor gates entry.
start = spin.index('template<class Map,bool Packed,bool Quarter,bool Four> void Spin::run')
end = spin.index('void Spin::encodeUnchecked',start)
fast = spin[start:end].replace(
    'template<class Map,bool Packed,bool Quarter,bool Four> void Spin::run',
    'template<class Map,bool Packed> void Spin::runFour')
fast = replace(fast, '    block* values=w.buckets.data();',
    '    constexpr bool Four=true,Quarter=false;\n    block* values=w.buckets.data();')
put('Fast.cpp', '''#include "Spin.h"
#include "Inner.h"
#include "AsymmetricMap.h"
#include "Weight5Inner.h"
#include "generated/BchCircuit.h"
#include "QuarterCircuit.h"
#include <cstring>
namespace bare_spin {
void bchTranspose4(const block*,block*);
static OC_FORCEINLINE u32 unpack(const u8* p) {u32 v;std::memcpy(&v,p,4);return v&0xffffff;}
''' + fast + '''
template void Spin::runFour<AsymmetricMap,true>(const block*,block*,Workspace&) const;
template void Spin::runFour<AsymmetricMap,false>(const block*,block*,Workspace&) const;
}
''')
# Copy headers so quoted includes cannot accidentally select the legacy Spin class.
for name in ('Inner.h','WorkspaceRouting.h','generated/SelectedMaps.h',
             'generated/BchCircuit.h','generated/BchCircuit.cpp'):
    put(name, (BARE / name).read_text())
for name in ('AsymmetricMap.h','Weight5Inner.h'):
    put(name, (BASE / name).read_text())
put('QuarterCircuit.h', (ROOT / 'workstreams/rate_quarter_bch/implementation/generated/QuarterCircuit.h').read_text())

# Reuse the exact symbolic verifier and emitter, but synthesize only the winner.
spec = importlib.util.spec_from_file_location('bch_emit', ROOT / 'workstreams/inner_design/bch_20260919/generate.py')
emit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(emit)
emit.HERE = OUT
words = [int(x,16) for x in re.findall(r'0x([0-9a-f]+)ULL', (BARE / 'generated/BchCircuit.h').read_text())]
targets = [sum(words[4*i+j] << (64*j) for j in range(4)) for i in range(128)]
emit.paar.DIMENSION = 256
circuit = emit.paar.synthesize(0,4,targets)
emit.emit('BchAvx512',circuit.gates,circuit.output_signals,targets,128,width=512)
path = OUT / 'generated/BchAvx512.cpp'
text, count = re.subn(r'const auto v(\d+)=_mm512_inserti32x4\(.*;',
    lambda m: f'const auto v{m[1]}=_mm512_loadu_si512(a+4*{m[1]});', path.read_text())
if count != 256:
    raise RuntimeError('four-row load generation failed')
put('generated/BchAvx512.cpp',text)
print('Generated integrated encoder; BCH output forms verified exactly.')
