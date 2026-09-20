"""Maintainer-only snapshot import; ordinary package builds never run Python.

Input: configured spin_optimized build using ForwardRecommended.cmake.
Writes SPIN's private kernel snapshot, generic circuit headers, and generation
record. The library's MIT license is maintained separately, not copied from
consumer checkouts.
"""
from pathlib import Path
import hashlib
import sys

if len(sys.argv) != 2:
    raise SystemExit('usage: import_kernels.py CONFIGURED_SPIN_BUILD')
build = Path(sys.argv[1])
package = Path(__file__).resolve().parents[1]
dest = package / 'src/kernels'
dest.mkdir(parents=True, exist_ok=True)


def rewrite(s):
    return (s.replace('namespace bare_spin', 'namespace spin::detail::kernel')
            .replace('bare_spin::', 'spin::detail::kernel::')
            .replace('using namespace bare_spin;', 'using namespace spin::detail::kernel;')
            .replace('osuCrypto', 'spin::detail::storage')
            .replace('hc_spin_wide_', 'spin_internal_wide_')
            .replace('OC_FORCEINLINE', 'SPIN_FORCEINLINE')
            .replace('__attribute__((noinline))', 'SPIN_NOINLINE')
            .replace('__builtin_cpu_supports("avx512f")', 'cpu_avx512f()')
            .replace('__builtin_cpu_supports("avx512vl")', 'cpu_avx512vl()'))


files = ('Spin.h', 'Spin.cpp', 'Inner.h', 'ImtRounds.h', 'K16Inner.h',
         'K16R2Inner.h', 'Map64S12.h', 'LengthGeometry.h', 'WorkspaceRouting.h',
         'Fast.cpp', 'ForwardFast.cpp', 'generated/BchCircuit.h',
         'generated/BchCircuit.cpp', 'generated/BchForward.h',
         'generated/BchForward.cpp', 'generated/BchAvx512.cpp',
         'generated/BchForward512.cpp', 'generated/SelectedMaps.h',
         'wide/Wide256.cpp', 'wide/Wide512.cpp', 'wide/WideKernel.h',
         'wide/generated/WideCircuit.h')
hashes = []
for name in files:
    source = build / 'bidirectional' / name
    raw = source.read_bytes()
    text = rewrite(raw.decode())
    if name == 'Spin.h':
        text = text.replace('private:\n', 'private:\n    friend struct Access;\n', 1)
        text = text.replace('    std::vector<u32> mForwardDirect;',
            '    std::vector<u32> mRangeRoute32;\n'
            '    template<class Map> void runRange(const block*,block*,Workspace&) const;\n'
            '    std::vector<u32> mForwardDirect;')
    if name == 'Spin.cpp':
        anchor = '    if(packed24Available()) {mSlots24.resize'
        assert text.count(anchor) == 1
        text = text.replace(anchor, '''    // Standalone transpose's measured range-direct schedule. Reuse the
    // forward direct table where present, without extending forward dispatch.
    if(mK<=458752 && mBchBackend==BchBackend::Avx512 && mForwardDirect.empty()) {
        mRangeRoute32.resize(n);
        for(std::size_t i=0;i<n;++i) mRangeRoute32[i]=packFour(mRoute[i]);
    }
''' + anchor)
        text = text.replace('    return 4*mForwardDirect.capacity()',
                            '    return 4*mRangeRoute32.capacity()+4*mForwardDirect.capacity()')
        anchor = 'void Spin::encodeUnchecked(const block* in,block* out,Workspace& w,Layout layout) const {'
        text = text.replace(anchor, anchor + '''
#if SPIN_BCH_AVX512
    if(mK<=458752 && mBchBackend==BchBackend::Avx512 && mK!=(1U<<16) &&
       (mConfig!=Configuration::T128S19 || mK!=(1U<<18))) {
        switch(mConfig) {
        case Configuration::T128S19:runRange<Map128S19>(in,out,w);return;
        case Configuration::T64S12:runRange<Map64S12>(in,out,w);return;
        case Configuration::T64S12R2:runRange<Map64S12R2>(in,out,w);return;
        default:break;
        }
    }
#endif
''')
    if name == 'Fast.cpp':
        assert text.endswith('}\n')
        text = text[:-2] + '''// Range-direct transpose from the standalone target; compile-time map selection.
template<class Map> void Spin::runRange(const block* in,block* out,Workspace& w) const {
    block* direct=w.buckets.data();
    const auto n=codeBlocks();
    const auto* route=mForwardDirect.empty()?mRangeRoute32.data():mForwardDirect.data();
    auto emit=[&](std::size_t i,block v) {direct[route[i]]=v;};
    if constexpr(std::is_same_v<Map,Map64S12>) inner64Reverse(in,n,mFieldRows.data(),emit);
    else if constexpr(std::is_same_v<Map,Map64S12R2>) inner64R2Reverse(in,n,mFieldRows.data(),emit);
    else innerReverse<Map>(in,n,mFieldRows.data(),emit);
    for(std::size_t j=0;j<n;j+=1024) bchTranspose4(direct+j,out+j/2);
}
template void Spin::runRange<Map128S19>(const block*,block*,Workspace&) const;
template void Spin::runRange<Map64S12>(const block*,block*,Workspace&) const;
template void Spin::runRange<Map64S12R2>(const block*,block*,Workspace&) const;
}
'''
    out = dest / name
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, newline='\n')
    hashes.append((name, hashlib.sha256(raw).hexdigest(),
                   hashlib.sha256(out.read_bytes()).hexdigest()))

# These are public template implementation details, with no SIMD or consumer types.
circuits = (build / 'generated/GenericCircuits.h').read_text().replace(
    'namespace bare_spin::generic_detail', 'namespace spin::detail::generic')
inc = package / 'include/spin/detail'
inc.mkdir(parents=True, exist_ok=True)
circuits = circuits.replace('template<class E> static inline void finish(const E* z,E* out) {',
    'template<class E,class Op> static inline void finish(const E* z,E* out,const Op& vx) {')
circuits = circuits.replace('template<class E,class Emit>', 'template<class E,class Emit,class Op>')
circuits = circuits.replace('std::size_t base,Emit& emit) {', 'std::size_t base,Emit& emit,const Op& vx) {')
circuits = circuits.replace('template<class E> inline void bch(const E* a,E* out) {',
    'template<class E,class Op> inline void bch(const E* a,E* out,const Op& vx) {')
(inc / 'GenericCircuits.h').write_text(circuits, newline='\n')

original = (package.parent / 'workstreams/spin_optimized/GenericSpin.h').read_text()
start = original.index('    template<class E> static inline void imt(')
end = original.rindex('\n};')
body = original[start:end].replace('generic_detail::', 'detail::generic::')
body = body.replace('mRoute', 'data_->route').replace('mMasks', 'data_->masks')
body = body.replace('mK', 'data_->k').replace('codeBlocks()', 'code_size()')
body = body.replace('w.values', 'w.values_')
body = body.replace('template<class E>', 'template<class E,class Op>')
body = body.replace('class E> static', 'class E,class Op> static')
body = body.replace('XorElement E>', 'ValueElement E,class Op>')
body = body.replace('u32 u,u32 v)', 'u32 u,u32 v,const Op& op)')
body = body.replace('(E* x)', '(E* x,const Op& op)')
body = body.replace('Workspace<E>& w)', 'Workspace<E>& w,const Op& op)')
body = body.replace('detail::generic::vx(', 'op(')
body = body.replace('zetaStage<T,D,Base+2*D>(x)', 'zetaStage<T,D,Base+2*D>(x,op)')
body = body.replace('zetaStage<T,D>(x)', 'zetaStage<T,D>(x,op)')
body = body.replace('zeta<T,D/2>(x)', 'zeta<T,D/2>(x,op)')
body = body.replace('state.data(),base,emit)', 'state.data(),base,emit,op)')
body = body.replace('zeta<Map::T>(raw.data())', 'zeta<Map::T>(raw.data(),op)')
body = body.replace('syndrome.data())', 'syndrome.data(),op)')
body = body.replace('masks[3])', 'masks[3],op)').replace('masks[1])', 'masks[1],op)')
body = body.replace('out+128*row)', 'out+128*row,op)')
# Included inside GenericTranspose. Kernel bodies keep the existing compile-time shape.
(inc / 'GenericMethods.h').write_text(body + '\n', newline='\n')

lines = ['# Kernel generation record', '',
         'SPIN is an independently licensed MIT library; see [LICENSE](LICENSE).',
         'The packaged kernels were generated from the SPIN research workspace at commit',
         '`3ced6dc0`, using `ForwardRecommended.cmake`. The wide-kernel development snapshot',
         'was staged at commit `2582a9fb366d28750ef92afdd9b2dbe9538552c6`.',
         "These references record SPIN's development and packaging history, not attribution",
         'of the library to its consumer projects.', '',
         'The import changes namespaces, internal symbol prefixes, and compiler',
         'portability spellings, and adds a private setup accessor. It also ports',
         'the standalone transpose range-direct dispatch through K=458752, reusing',
         'the forward direct route table where present. It does not',
         'resynthesize circuits or change the full/partial-tile kernel schedules.',
         'Block storage and capability detection are package-owned files.',
         'The generic transpose circuits come from the same configured build.', '',
         'The package has no runtime or build dependency on its consumer projects.', '',
         '| Source file | Import input SHA-256 | Imported SHA-256 |',
         '|---|---|---|']
lines += [f'| `{n}` | `{a}` | `{b}` |' for n, a, b in hashes]
(package / 'PROVENANCE.md').write_text('\n'.join(lines) + '\n', newline='\n')
print(f'Imported {len(files)} private kernel files and generic circuits.')
