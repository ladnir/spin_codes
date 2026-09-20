"""Import width-lifted forward circuits without editing the upstream snapshot."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

source, out = map(Path, sys.argv[1:])
if not (out / 'Spin.h').exists():
    raise RuntimeError('Generate the bidirectional integration first')
names = ('Wide256.cpp', 'Wide512.cpp', 'WideKernel.h', 'generate_wide.py', 'libote-LICENSE')
wide = out / 'wide'
(wide / 'generated').mkdir(parents=True, exist_ok=True)
for name in names:
    (wide / name).write_bytes((source / name).read_bytes())
# Generate from the unchanged S19 source, not our subsequently patched Inner.h
# (which also contains K16-specific traits that upstream's extractor does not lift).
inputs = ('generated/BchForward.cpp', 'generated/SelectedMaps.h', 'Inner.h')
for name in inputs:
    (wide / name).write_bytes((source / name).read_bytes())
subprocess.run([sys.executable, str(wide / 'generate_wide.py')], check=True)

def edit(text, old, new):
    if text.count(old) != 1:
        raise RuntimeError(f'wide integration anchor changed: {old}')
    return text.replace(old, new)

def lift(text):
    for a, b in (('__m128i', 'Vec'), ('_mm_setzero_si128()', 'Ops::zero()'),
                 ('_mm_xor_si128', 'Ops::vx'), ('alignas(32)', 'alignas(64)')):
        text = text.replace(a, b)
    return text

# Lift the same forward recurrence used by the integrated 128-bit path. In
# particular retain T2(T1(state)) followed by feedback, not two feedback updates.
inner = (out / 'Inner.h').read_text()
start = inner.index('template<class Map,std::size_t P,class Gather>')
end = inner.index('template<u32 Mask> OC_FORCEINLINE __m128i imtSparseSum')
forward = lift(inner[start:end])
forward = forward.replace('isImtMap<Map>', '(Map::Rounds>0)').replace('isTwoRoundMap<Map>', '(Map::Rounds==2)')
forward = forward.replace('block* out', 'Vec* out').replace('gather(base+P).mData', 'gather(base+P)')
forward = forward.replace('out[base+P]=block(v);', 'Ops::store(out+base+P,v);')
forward = forward.replace('values[p]=v.mData; out[base+p]=v;', 'values[p]=v; Ops::store(out+base+p,v);')
forward = forward.replace('OC_FORCEINLINE', 'static OC_FORCEINLINE')
rounds = (out / 'ImtRounds.h').read_text()
step = lift(rounds[rounds.index('OC_FORCEINLINE void imtStep'):rounds.index('// Independent dense-reference')])
step = step.replace('OC_FORCEINLINE', 'static OC_FORCEINLINE')
small = (out / 'Map64S12.h').read_text()
small = small[small.index('struct Map64S12 {'):small.index('template<class Emit>')]
small = lift(small).replace('struct Map64S12 {',
    'template<class Ops> struct WideMap64S12R2 {\nusing Vec=typename Ops::Vec;\nstatic constexpr unsigned Rounds=2;') + '};\n'
header_path = wide / 'generated/WideCircuit.h'
header = header_path.read_text()
header = edit(header, 'template<class Ops> struct WideMap128S19 {',
    'template<class Ops> struct WideMap128S19 {\nstatic constexpr unsigned Rounds=1;')
header = edit(header, 'template<class Ops> struct WideCircuit {', small + '\ntemplate<class Ops> struct WideCircuit {')
start = header.index('template<class Map,std::size_t P,class Gather>')
assert header.endswith('\n};\n}\n')
header = header[:start] + step + forward + '\n};\n}\n'
assert '__m128i' not in header and '.mData' not in header
header_path.write_text(header, newline='\n')
manifest_path = wide / 'generated/WIDE_MANIFEST.json'
manifest = json.loads(manifest_path.read_text())
manifest['integrated_inputs'] = {n: hashlib.sha256((out/n).read_bytes()).hexdigest()
                                 for n in ('Map64S12.h', 'Inner.h', 'ImtRounds.h')}
manifest['output_sha256'] = hashlib.sha256(header_path.read_bytes()).hexdigest()
manifest['description'] = 'Width lifting of S19 and certified K16 two-round forward maps; no new XOR synthesis.'
manifest_path.write_text(json.dumps(manifest, indent=2)+'\n', newline='\n')

# Preserve wideView()'s original S19-only contract for upstream consumers.
h = (out/'Spin.h').read_text()
h = edit(h, 'const u32* fieldRows; // IMT: consecutive (u,v) masks per epoch.',
         'const u32* fieldRows; // IMT: (u,v) pairs, epoch then round.\n'
         '        Configuration configuration=Configuration::T128S19;')
h = edit(h, '    WideView wideView() const;', '    WideView wideView() const;\n    WideView wideForwardView() const;')
(out/'Spin.h').write_text(h, newline='\n')
s = (out/'Spin.cpp').read_text()
s = edit(s, 'Spin::WideView Spin::wideView() const {', '''Spin::WideView Spin::wideForwardView() const {
    if((mConfig!=Configuration::T128S19 && mConfig!=Configuration::T64S12R2) || mSlots24.empty())
        throw std::invalid_argument("wide forward requires S19 or K16 R2 and Packed24");
    return {codeBlocks(),tileBlocks(),mSlots24.data(),mOffsets24.data(),mForwardFieldRows.data(),mConfig};
}
Spin::WideView Spin::wideView() const {''')
(out/'Spin.cpp').write_text(s, newline='\n')
kernel = (wide/'WideKernel.h').read_text().replace('->wideView()', '->wideForwardView()')
kernel = edit(kernel, 'using Map=WideMap128S19<Ops>;', '')
kernel = edit(kernel, 'static void encode(Spin::WideView', 'template<class Map> static void encode(Spin::WideView')
kernel = edit(kernel, '''        bare_spin::HC_NS::encode(v,static_cast<const bare_spin::HC_NS::Vec*>(in),
                                  static_cast<bare_spin::HC_NS::Vec*>(out),w);''', '''        using namespace bare_spin;
        using Ops=HC_NS::Ops;
        const auto* input=static_cast<const HC_NS::Vec*>(in);auto* output=static_cast<HC_NS::Vec*>(out);
        if(v.configuration==Configuration::T64S12R2) HC_NS::encode<WideMap64S12R2<Ops>>(v,input,output,w);
        else HC_NS::encode<WideMap128S19<Ops>>(v,input,output,w);''')
(wide/'WideKernel.h').write_text(kernel, newline='\n')
# Retain Hypercat's actual streaming 4x4 assembly kernel, not a scalar control.
# Only its unsupported-shape fallback changes: this benchmark requires 16 planes.
text = (source / 'pcs_kernels.cpp').read_text()
start = text.index('extern "C" HC_AVX512 void hc_pcs_assemble512_stream(')
brace = text.index('{', start)
depth = 1
end = brace + 1
while depth:
    depth += (text[end] == '{') - (text[end] == '}')
    end += 1
body = text[start:end].replace('HC_AVX512 ', '')
old = 'hc_pcs_assemble(input,output,n,stride,count,columns); return;'
assert body.count(old) == 1
body = body.replace(old, 'throw std::invalid_argument("assembly requires aligned four-word groups");')
(out / 'Column128.cpp').write_text(
    '// From Hypercat pcs_kernels.cpp; see wide/libote-LICENSE and WIDE_IMPORT.json.\n'
    '#include "Block.h"\n#include <cstddef>\n#include <cstdint>\n#include <stdexcept>\n'
    'using osuCrypto::block;\n' + body + '\n', newline='\n')
(out / 'WIDE_IMPORT.json').write_text(json.dumps({
    'inputs': {n: hashlib.sha256((source/n).read_bytes()).hexdigest()
               for n in (*names, *inputs, 'pcs_kernels.cpp')},
    'column128': 'Streaming kernel extracted verbatim except unsupported-shape fallback throws.',
    'wide': 'Upstream outer/routing schedules plus width-lifted integrated forward recurrence; dispatch once for S19 or K16 R2.',
    'integration_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
}, indent=2) + '\n', newline='\n')
