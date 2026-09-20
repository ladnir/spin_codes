"""Isolated half-rate implementations of the certified weight-five map."""
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import candidate
model = candidate.model
ROOT = model.ROOT
INNER = model.ASYMMETRIC.parent
sys.path.insert(0,str(INNER))
import generate_balanced as balanced


def once(text,old,new):
    assert text.count(old) == 1, (old[:80],text.count(old))
    return text.replace(old,new,1)


def between(text,start,end,new):
    a=text.index(start);b=text.index(end,a)
    return text[:a]+new+text[b:]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    certificate_path=HERE.parent/'FULL_M20_VERIFIED.json'
    cert=model.base.read(certificate_path)
    model.authenticate(cert)
    assert cert['full_distance_proved'] and cert['target_bits']==40
    engine=candidate.Engine(20)
    assert engine.identity()==cert['instance']
    code,stats=balanced.header(engine.a_columns)
    code+='\nnamespace bare_spin {\nstruct AsymmetricMap : BalancedMap {\n'
    code+='static constexpr auto feedbackColumns=BalancedMap::columns;\n'
    code+='static constexpr std::array<std::uint32_t,T> columns{'+','.join(map(hex,engine.columns))+'};\n'
    code+='static constexpr auto groupedColumns=columns;\n'
    code+='static constexpr std::array<unsigned,S> groupOrder{'+','.join(map(str,range(19)))+'};\n'
    # A second exact-map emission strategy shares XOR expressions globally.
    balanced.shared.paar.DIMENSION=19
    circuit=balanced.shared.paar.synthesize(0,2,engine.columns)
    symbolic=[1<<j for j in range(19)]
    for a,b in circuit.gates:
        symbolic.append(symbolic[a]^symbolic[b])
    code+='template<class Emit> static OC_FORCEINLINE void emitShared(const block* in,__m128i* raw,const __m128i* state,std::size_t base,Emit& emit) {\n'
    code+='\n'.join(f'const auto v{j}=state[{j}];' for j in range(19))+'\n'
    code+='\n'.join(f'const auto v{i+19}=vx(v{a},v{b});' for i,(a,b) in enumerate(circuit.gates))+'\n'
    for p in reversed(range(128)):
        signals=circuit.output_signals[p]
        rebuilt=0
        for i in signals:rebuilt^=symbolic[i]
        assert rebuilt==engine.columns[p]
        code+=f'raw[{p}]=in[{p}].mData; emit(base+{p},block(vx(raw[{p}],{balanced.shared.expression(signals)})));\n'
    code+='}\n'
    order,before,after=balanced.shared.grouping(engine.columns,19)
    grouped=[sum(((c>>bit)&1)<<j for j,bit in enumerate(order)) for c in engine.columns]
    assert all(sum(((c>>j)&1)<<bit for j,bit in enumerate(order))==original for c,original in zip(grouped,engine.columns))
    code+='template<class Emit> static OC_FORCEINLINE void emitGrouped(const block* in,__m128i* raw,const __m128i* state,std::size_t base,Emit& emit) {\n'
    code+='alignas(32) __m128i grouped[19],table[5][16];\n'
    code+='\n'.join(f'grouped[{j}]=state[{bit}];' for j,bit in enumerate(order))+'\ntables<19>(grouped,table);\n'
    for p in reversed(range(128)):
        code+=f'raw[{p}]=in[{p}].mData; emit(base+{p},block(vx(raw[{p}],fixedSum<{hex(grouped[p])}>(table))));\n'
    code+='}\n};\n}\n'
    (HERE/'AsymmetricMap.h').write_text(code,encoding='utf-8',newline='\n')
    kernel=(model.ASYMMETRIC/'AsymmetricInner.h').read_text()
    kernel=once(kernel,'} else asymmetricPoints<Map>(in+base,raw,state,base,emit,std::make_index_sequence<Map::T>{});',
        '} else {\n#if W5_SHARED_EMISSION==2\n'
        '                Map::emitGrouped(in+base,raw,state,base,emit);\n#elif W5_SHARED_EMISSION==1\n'
        '                Map::emitShared(in+base,raw,state,base,emit);\n#else\n'
        '                asymmetricPoints<Map>(in+base,raw,state,base,emit,std::make_index_sequence<Map::T>{});\n#endif\n            }')
    (HERE/'Weight5Inner.h').write_text(kernel,encoding='utf-8',newline='\n')
    base_path=INNER/'generated/asymmetric_greedy3_2_sparse/CandidateSpin.cpp'
    code=base_path.read_text()
    code=once(code,'#include "../../asymmetric/AsymmetricInner.h"','#include "Weight5Inner.h"\n#include "WorkspaceRouting.h"')
    marker='    mK=std::size_t{1}<<exponent;'
    code=once(code,marker,'    if(c!=Configuration::T128S19 || outer!=Outer::Bch256x128)\n'
        '        throw std::invalid_argument("isolated weight-five half-rate instance only");\n'+marker)
    code=between(code,'    switch(c) {','\ntemplate<class Map> void Spin::setupInner',
        '    setupInner<AsymmetricMap>(coefficientSeed);\n}\n')
    # All uses of this runtime check in the inherited candidate select its new
    # mixer, setup validator, and name. Remove dispatch bodies separately below.
    code=code.replace('mOuter==Outer::Bch128x32','mOuter==Outer::Bch256x128')
    old='Spin::Workspace::Workspace(const Spin& code):buckets(code.codeBlocks()),tile(code.tileBlocks()) {}'
    workspace='''Spin::Workspace::Workspace(const Spin& code):buckets(code.codeBlocks()),tile(code.tileBlocks()) {
    if(workspace_routing::eligible(true,code.messageBlocks())) {
        workspace_routing::adviseOwned(tile.data(),tile.size()*sizeof(block));
        workspace_routing::adviseOwned(buckets.data(),buckets.size()*sizeof(block));
    }
}'''
    code=once(code,old,workspace)
    code=once(code,'asymmetricReverse<Map,Quarter,false>','asymmetricReverse<Map,true,W5_MASKED>')
    code=between(code,'void Spin::encodeUnchecked(', 'void Spin::validateSetup()',
        '''void Spin::encodeUnchecked(const block* in,block* out,Workspace& w,Layout layout) const {
    if(layout==Layout::Packed24) run<AsymmetricMap,true,false>(in,out,w);
    else run<AsymmetricMap,false,false>(in,out,w);
}
''')
    a=code.index('template<class Map,bool Quarter> void Spin::oracle')
    b=code.index('void Spin::reference(',a)
    code=code[:a]+code[a:b].replace('if constexpr(Quarter)','if constexpr(!Quarter)')+code[b:]
    code=code[:code.index('void Spin::reference(')]+'''void Spin::reference(const block* in,block* out) const {
    if(mCompacted) throw std::logic_error("oracle schedules were discarded");
    oracle<AsymmetricMap,false>(in,out);
}
}
'''
    code=once(code,'"t128_s19_asymmetric_greedy3_2_r1"','"t128_s19_weight5_seed0_r1"')
    code=once(code,'if(j+32<tileSize) {','if constexpr(W5_PREFETCH>0) if(j+W5_PREFETCH<tileSize) {')
    code=code.replace('base+j+32','base+j+W5_PREFETCH')
    code=once(code,'_mm_prefetch(reinterpret_cast<const char*>(tile+future),_MM_HINT_T0);',
        '__builtin_prefetch(tile+future,W5_WRITE_PREFETCH,3);')
    (HERE/'Weight5Spin.cpp').write_text(code,encoding='utf-8',newline='\n')
    # Give the old inner the same half-rate workspace/prefetch option as a control.
    bare=ROOT/'workstreams/bare_bch_rm2sub'
    baseline=(bare/'Spin.cpp').read_text()
    baseline=once(baseline,'#include "../rate_quarter_bch/implementation/generated/QuarterCircuit.h"','#include "QuarterCircuit.h"')
    baseline=once(baseline,'workspace_routing::eligible(code.outerLength()==128,code.messageBlocks())','workspace_routing::eligible(true,code.messageBlocks())')
    baseline=once(baseline,'_mm_prefetch(reinterpret_cast<const char*>(tile+future),_MM_HINT_T0);','__builtin_prefetch(tile+future,1,3);')
    (HERE/'BaselineTuned.cpp').write_text(baseline,encoding='utf-8',newline='\n')
    test=(bare/'correctness.cpp').read_text()
    test=once(test,'for(unsigned c=0;c<4;++c)','for(unsigned c=2;c<3;++c)')
    (HERE/'correctness.cpp').write_text(test,encoding='utf-8',newline='\n')
    sources=[Path(__file__),base_path,model.ASYMMETRIC/'AsymmetricInner.h',
        Path(balanced.__file__),Path(balanced.shared.__file__),Path(balanced.shared.paar.__file__),certificate_path,
        bare/'Spin.cpp',bare/'Spin.h',bare/'Inner.h',bare/'WorkspaceRouting.h',bare/'benchmark.cpp',bare/'correctness.cpp',
        bare/'generated/BchCircuit.cpp',bare/'generated/BchCircuit.h']
    payload=dict(status='SYMBOLICALLY_CHECKED_WEIGHT5_IMPLEMENTATION_NOT_BENCHMARKED',
        certificate_sha256=digest(certificate_path),instance=engine.identity(),feedback_circuit=stats,
        direct_emission_xors=512,shared_emission_xors=circuit.xor_count,
        grouped_emission_lookups=after,original_grouping_lookups=before,
        generated_sha256={name:digest(HERE/name) for name in ('AsymmetricMap.h','Weight5Inner.h','Weight5Spin.cpp','BaselineTuned.cpp','correctness.cpp')},
        source_sha256={p.relative_to(ROOT).as_posix():digest(p) for p in sources})
    (HERE/'IMPLEMENTATION.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    print('Generated certified map; emission XORs: direct',512,'shared',circuit.xor_count)


if __name__=='__main__':main()
