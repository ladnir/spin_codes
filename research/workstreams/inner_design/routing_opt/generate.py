"""Exact-map routing variants; isolate all changes from supported sources."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
SOURCE = PARENT / 'generated/asymmetric_greedy3_2_sparse/CandidateSpin.cpp'


def main():
    original = SOURCE.read_text()
    prior = json.loads((PARENT/'asymmetric/IMPLEMENTATION.json').read_text())
    record = next(x for x in prior['candidates'] if x['name']=='asymmetric_greedy3_2_sparse')
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==record['source_sha256']['CandidateSpin.cpp']
    start = original.index('        for(std::size_t j=0;j<tileSize;++j) {')
    end = original.index('        if constexpr(Quarter) {', start)
    variants = []
    for mode, distance in [('scatter',0),('scatter',32),('scatter',128),('write',32),('write',128),('gather',0),('gather',64),('scalar',0),('scalar',128),('scalar',256),('scalarwrite',32),('scalarwrite',128),('innerwrite',128),('innerwrite',512)]:
        name=f'routeopt_{mode}{distance}'
        code=original
        if mode=='innerwrite':
            old='''        if constexpr(Packed) values[unpack(mSlots24.data()+3*i)]=v;'''
            new=f'''        if(i>={distance}) {{
            const auto next=Packed?unpack(mSlots24.data()+3*(i-{distance})):mSlots32[i-{distance}];
            __builtin_prefetch(values+next,1,3);
        }}
'''+old
            assert code.count(old)==1
            code=code.replace(old,new)
        elif mode.startswith('scalar'):
            body=original[start:end]
            if not distance:
                begin=body.index('            if(j+32<tileSize) {')
                finish=body.index('            const auto offset=',begin)
                body=body[:begin]+body[finish:]
            else:
                body=body.replace('j+32',f'j+{distance}')
                if mode=='scalarwrite':
                    body=body.replace('_mm_prefetch(reinterpret_cast<const char*>(tile+future),_MM_HINT_T0);','__builtin_prefetch(tile+future,1,3);')
            code=code[:start]+body+code[end:]
        elif mode=='gather':
            body='''        if constexpr(Quarter) {
            for(std::size_t j=0;j<tileSize;j+=256) {
                for(unsigned p=0;p<256;p+=8) {
                    auto point=[&]<unsigned K>() {
                        const auto pos=base+j+p+K;
                        PREFETCH
                        const auto off=Packed?unpack(mOffsets24.data()+3*pos):mOffsets32[pos];
                        tile[p+K]=values[base+off];
                    };
                    CALLS
                }
                quarterTranspose2(tile,tile+128,out+base/4+j/4,out+base/4+j/4+32);
            }
            continue;
        }
'''
            pf=f'''if(j+p+K+{distance}<tileSize) {{
                            const auto off=Packed?unpack(mOffsets24.data()+3*(pos+{distance})):mOffsets32[pos+{distance}];
                            _mm_prefetch(reinterpret_cast<const char*>(values+base+off),_MM_HINT_T0);
                        }}''' if distance else ''
            body=body.replace('PREFETCH',pf).replace('CALLS','\n                    '.join(f'point.template operator()<{k}>();' for k in range(8)))
            code=code[:start]+body+code[start:]
            old='mSlots32[inner]=slot; mOffsets32[slot]=local;\n        pack(mSlots24.data()+3*inner,slot); pack(mOffsets24.data()+3*slot,local);'
            new='''mSlots32[inner]=slot;
        const auto offsetIndex=mOuter==Outer::Bch128x32?outer:slot;
        const auto offsetValue=mOuter==Outer::Bch128x32?slot-bucket*u32(tile):local;
        mOffsets32[offsetIndex]=offsetValue;
        pack(mSlots24.data()+3*inner,slot);pack(mOffsets24.data()+3*offsetIndex,offsetValue);'''
            assert old in code
            code=code.replace(old,new)
            old='''        if(unpack(mSlots24.data()+3*i)!=slot || unpack(mOffsets24.data()+3*slot)!=mOffsets32[slot] ||
           (slot/tileBlocks())*tileBlocks()+mOffsets32[slot]!=outer)'''
            new='''        const auto index=mOuter==Outer::Bch128x32?outer:slot;
        const auto target=mOuter==Outer::Bch128x32?slot:outer;
        if(unpack(mSlots24.data()+3*i)!=slot || unpack(mOffsets24.data()+3*index)!=mOffsets32[index] ||
           (slot/tileBlocks())*tileBlocks()+mOffsets32[index]!=target)'''
            assert old in code
            code=code.replace(old,new)
        else:
            body='''        for(std::size_t j=0;j<tileSize;j+=8) {
            auto point=[&]<unsigned K>() {
                const auto pos=base+j+K;
                PREFETCH
                const auto offset=Packed?unpack(mOffsets24.data()+3*pos):mOffsets32[pos];
                tile[offset]=values[pos];
            };
            CALLS
        }
'''
            pf=f'''if(j+K+{distance}<tileSize) {{
                    const auto future=Packed?unpack(mOffsets24.data()+3*(pos+{distance})):mOffsets32[pos+{distance}];
                    __builtin_prefetch(tile+future,{1 if mode=='write' else 0},3);
                }}''' if distance else ''
            body=body.replace('PREFETCH',pf).replace('CALLS','\n            '.join(f'point.template operator()<{k}>();' for k in range(8)))
            code=code[:start]+body+code[end:]
        code=code.replace('"t128_s19_asymmetric_greedy3_2_r1"',f'"{name}"')
        directory=PARENT/'generated'/name
        directory.mkdir(parents=True,exist_ok=True)
        (directory/'CandidateSpin.cpp').write_text(code,encoding='utf-8',newline='\n')
        (directory/'AsymmetricMap.h').write_text((SOURCE.parent/'AsymmetricMap.h').read_text(),encoding='utf-8',newline='\n')
        variants.append(dict(name=name,mode=mode,prefetch=distance,sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in directory.iterdir() if p.is_file()}))
    (HERE/'IMPLEMENTATION.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),candidates=variants),indent=2)+'\n')
    print(' '.join(x['name'] for x in variants))


if __name__=='__main__':
    main()
