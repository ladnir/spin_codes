"""Try exact-map direct gathering; no change to permutations or code ensemble."""
import hashlib
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent;PARENT=HERE.parent;ROOT=PARENT.parents[1]


def main():
    path=PARENT/'generated/asymmetric_greedy3_2_sparse/CandidateSpin.cpp';original=path.read_text()
    manifest=json.loads((PARENT/'asymmetric/IMPLEMENTATION.json').read_text())
    previous=next(r for r in manifest['candidates'] if r['name']=='asymmetric_greedy3_2_sparse')
    assert hashlib.sha256(path.read_bytes()).hexdigest()==previous['source_sha256']['CandidateSpin.cpp']
    variants=[]
    for distance in (0,64):
        name=f'redesign_gather_pf{distance}';code=original
        marker='    switch(c) {\n        case Configuration::T64S16:'
        assert code.count(marker)==1
        code=code.replace(marker,'''    if(mOuter==Outer::Bch128x32) {
        for(std::size_t i=0;i<n;++i) {
            const auto outer=mRoute[i];mSlots32[outer]=u32(i);
            pack(mSlots24.data()+3*outer,u32(i));
        }
    }
'''+marker)
        marker='template<class Map,bool Packed,bool Quarter> void Spin::run(const block* in,block* out,Workspace& w) const {\n'
        body='''    if constexpr(Quarter) {
        block* values=w.buckets.data();
        asymmetricReverse<Map,true,false>(in,codeBlocks(),mFieldRows.data(),[&](std::size_t i,block v){values[i]=v;});
        block* tile=w.tile.data();const auto n=codeBlocks();
        for(std::size_t j=0;j<n;j+=256) {
            for(unsigned p=0;p<256;p+=8) {
                auto gather=[&]<unsigned K>() {
                    const auto position=j+p+K;
                    PREFETCH
                    const auto offset=Packed?unpack(mSlots24.data()+3*position):mSlots32[position];
                    tile[p+K]=values[offset];
                };
                gather.template operator()<0>();gather.template operator()<1>();
                gather.template operator()<2>();gather.template operator()<3>();
                gather.template operator()<4>();gather.template operator()<5>();
                gather.template operator()<6>();gather.template operator()<7>();
            }
            quarterTranspose2(tile,tile+128,out+j/4,out+j/4+32);
        }
        return;
    }
'''
        prefetch=f'''if(position+{distance}<n) {{
                        const auto future=Packed?unpack(mSlots24.data()+3*(position+{distance})):mSlots32[position+{distance}];
                        _mm_prefetch(reinterpret_cast<const char*>(values+future),_MM_HINT_T0);
                    }}''' if distance else ''
        assert code.count(marker)==1;code=code.replace(marker,marker+body.replace('PREFETCH',prefetch))
        marker='    std::vector<bool> seen(codeBlocks()),slots(codeBlocks());'
        validate='''    if(mOuter==Outer::Bch128x32) {
        std::vector<bool> seen(codeBlocks()),rowRegion(codeBlocks());
        const auto rows=mK/outerDimension();
        for(std::size_t i=0;i<codeBlocks();++i) {
            const auto outer=mRoute[i];
            if(outer>=codeBlocks() || seen[outer] || mSlots32[outer]!=i || unpack(mSlots24.data()+3*outer)!=i)
                throw std::runtime_error("invalid direct inverse route");
            seen[outer]=true;
            const auto key=(outer/128)*128+i/rows;
            if(rowRegion[key]) throw std::runtime_error("duplicate row/region");
            rowRegion[key]=true;
        }
        return;
    }
'''
        assert code.count(marker)==1;code=code.replace(marker,validate+marker)
        code=code.replace('"t128_s19_asymmetric_greedy3_2_r1"',f'"asymmetric_greedy3_2_exact_gather_pf{distance}"')
        directory=PARENT/'generated'/name;directory.mkdir(parents=True,exist_ok=True)
        (directory/'CandidateSpin.cpp').write_text(code,encoding='utf-8',newline='\n')
        (directory/'AsymmetricMap.h').write_text((path.parent/'AsymmetricMap.h').read_text(),encoding='utf-8',newline='\n')
        variants.append(dict(name=name,prefetch_distance=distance,
            generated_sha256={f:hashlib.sha256((directory/f).read_bytes()).hexdigest() for f in ('CandidateSpin.cpp','AsymmetricMap.h')}))
    sources=[Path(__file__),path,path.parent/'AsymmetricMap.h',PARENT/'asymmetric/IMPLEMENTATION.json']
    payload=dict(status='EXACT_SCHEDULE_CHANGE_PENDING_TESTS',candidates=variants,
        source_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (HERE/'ROUTING_IMPLEMENTATION.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    print(json.dumps([r['name'] for r in variants]))


if __name__=='__main__':main()
