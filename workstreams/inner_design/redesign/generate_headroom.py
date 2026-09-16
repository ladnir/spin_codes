"""Generate an explicitly invalid-distance no-inner performance diagnostic."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
PARENT=HERE.parent
ROOT=PARENT.parents[1]


def main():
    source=PARENT/'generated/asymmetric_greedy3_2_sparse/CandidateSpin.cpp'
    manifest=json.loads((PARENT/'asymmetric/IMPLEMENTATION.json').read_text())
    record=next(c for c in manifest['candidates'] if c['name']=='asymmetric_greedy3_2_sparse')
    assert hashlib.sha256(source.read_bytes()).hexdigest()==record['source_sha256']['CandidateSpin.cpp']
    code=source.read_text()
    old='asymmetricReverse<Map,Quarter,false>(in,codeBlocks(),mFieldRows.data(),[&](std::size_t i,block v) {\n        if constexpr(Packed) values[unpack(mSlots24.data()+3*i)]=v;\n        else values[mSlots32[i]]=v;\n    });'
    new='// INVALID-DISTANCE DIAGNOSTIC: remove the inner, retain routing and outer.\n    for(std::size_t i=codeBlocks();i-->0;) {\n        if constexpr(Packed) values[unpack(mSlots24.data()+3*i)]=in[i];\n        else values[mSlots32[i]]=in[i];\n    }'
    assert code.count(old)==1;code=code.replace(old,new)
    code=code.replace('"t128_s19_asymmetric_greedy3_2_r1"','"DIAGNOSTIC_NO_INNER_NOT_A_SPIN_CODE"')
    variants={'redesign_no_inner':code}
    unrolled='''// INVALID-DISTANCE DIAGNOSTIC: retain the 128-position unrolled emission.
    auto emit=[&](std::size_t i,block v) {
        if constexpr(Packed) values[unpack(mSlots24.data()+3*i)]=v;
        else values[mSlots32[i]]=v;
    };
    for(std::size_t epoch=codeBlocks()/128;epoch-->0;)
        [&]<std::size_t... P>(std::index_sequence<P...>) {
            (emit(epoch*128+127-P,in[epoch*128+127-P]),...);
        }(std::make_index_sequence<128>{});'''
    variants['redesign_no_inner_unrolled']=code.replace(new,unrolled).replace('DIAGNOSTIC_NO_INNER_NOT_A_SPIN_CODE','DIAGNOSTIC_UNROLLED_NO_INNER_NOT_A_SPIN_CODE')
    # Hypothetical contiguous outer inputs: removes inner AND routing.
    # It is not a distance-preserving rewrite or a speedup prediction.
    start=code.index('    block* values=w.buckets.data();',code.index('void Spin::run('))
    end=code.index('\nvoid Spin::encodeUnchecked',start)
    outer='''    if constexpr(Quarter) {
        for(std::size_t j=0;j<codeBlocks();j+=256)
            quarterTranspose2(in+j,in+j+128,out+j/4,out+j/4+32);
    } else throw std::runtime_error("outer-only diagnostic supports quarter rate only");
}
'''
    variants['redesign_outer_only']=(code[:start]+outer+code[end:]).replace('DIAGNOSTIC_NO_INNER_NOT_A_SPIN_CODE','DIAGNOSTIC_OUTER_ONLY_NOT_A_SPIN_CODE')
    header=source.parent/'AsymmetricMap.h'
    generated={}
    for name,text in variants.items():
        old_reference='    if(mOuter==Outer::Bch128x32) {oracle<AsymmetricMap,true>(in,out);return;}'
        assert text.count(old_reference)==1
        data='in' if name=='redesign_outer_only' else 'routed.data()'
        route='' if name=='redesign_outer_only' else '''        std::vector<block> routed(codeBlocks());
        for(std::size_t i=0;i<codeBlocks();++i) routed[mRoute[i]]=in[i];
'''
        reference='''    if(mOuter==Outer::Bch128x32) {
'''+route+f'''        const block* source={data};
        for(std::size_t row=0;row<mK/32;++row) for(unsigned j=0;j<32;++j) {{
            block value(0,0);
            for(unsigned c=0;c<128;++c) if((QuarterRows[j][c/64]>>(c%64))&1) value^=source[row*128+c];
            out[row*32+j]=value;
        }}
        return;
    }}'''
        text=text.replace(old_reference,reference)
        directory=PARENT/'generated'/name;directory.mkdir(parents=True,exist_ok=True)
        (directory/'CandidateSpin.cpp').write_text(text,encoding='utf-8',newline='\n')
        (directory/'AsymmetricMap.h').write_text(header.read_text(),encoding='utf-8',newline='\n')
        generated[name]={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in directory.iterdir() if p.is_file()}
    sources=[Path(__file__),source,header,PARENT/'asymmetric/IMPLEMENTATION.json']
    result=dict(status='NO_INNER_TIMING_DIAGNOSTIC_NOT_VALID_DISTANCE',
        source_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        generated_sha256=generated)
    (HERE/'HEADROOM_IMPLEMENTATION.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(list(variants)))


if __name__=='__main__':main()
