"""Two-level bucket schedule; same route and inner, extra cache-local pass."""
import hashlib
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
PARENT=HERE.parent
SOURCE=PARENT/'generated/asymmetric_greedy3_2_sparse/CandidateSpin.cpp'


def main():
    original=SOURCE.read_text()
    variants=[]
    for small in (4096,16384,65536):
        name=f'routeopt_hier{small}'
        code=original
        code=code.replace('mOffsets24.resize(3*n+4);','mOffsets24.resize(5*n+4);')
        code=code.replace('mOffsets32.resize(n);','mOffsets32.resize(2*n);')
        # Leave the original schedule construction intact, then replace its offsets.
        marker='    switch(c) {\n        case Configuration::T64S16:'
        setup=f'''    if(mOuter==Outer::Bch128x32) {{
        const auto sub=std::min<std::size_t>({small},tile);
        std::vector<u32> old(mOffsets32.begin(),mOffsets32.begin()+n);
        for(std::size_t base=0;base<n;base+=tile) {{
            std::vector<u32> count(tile/sub,0);
            for(std::size_t j=0;j<tile;++j) {{
                const auto destination=old[base+j];
                const auto bucket=destination/sub;
                const u32 slot=u32(bucket*sub+count[bucket]++);
                mOffsets32[2*(base+j)]=slot;
                mOffsets32[2*(base+slot)+1]=destination%sub;
                pack(mOffsets24.data()+5*(base+j),slot);
                const u32 local=destination%sub;
                auto* p=mOffsets24.data()+5*(base+slot)+3;
                p[0]=u8(local);p[1]=u8(local>>8);
            }}
        }}
    }}
'''
        assert code.count(marker)==1
        code=code.replace(marker,setup+marker)
        marker='    for(std::size_t base=0;base<n;base+=tileSize) {\n'
        body=f'''    if constexpr(Quarter) {{
        const auto sub=std::min<std::size_t>({small},tileSize);
        for(std::size_t base=0;base<n;base+=tileSize) {{
            // Partition a coarse bucket into sequential sub-bucket streams.
            for(std::size_t j=0;j<tileSize;++j) {{
                const auto off=Packed?unpack(mOffsets24.data()+5*(base+j)):mOffsets32[2*(base+j)];
                tile[off]=values[base+j];
            }}
            for(std::size_t b=0;b<tileSize;b+=sub) {{
                for(std::size_t j=0;j<sub;++j) {{
                    const auto pos=base+b+j;
                    u32 off;
                    if constexpr(Packed) {{
                        const auto* p=mOffsets24.data()+5*pos+3;
                        off=u32(p[0])|(u32(p[1])<<8);
                    }} else off=mOffsets32[2*pos+1];
                    values[base+b+off]=tile[b+j];
                }}
                for(std::size_t j=0;j<sub;j+=256)
                    quarterTranspose2(values+base+b+j,values+base+b+j+128,
                        out+(base+b+j)/4,out+(base+b+j)/4+32);
            }}
        }}
        return;
    }}
'''
        assert code.count(marker)==1
        code=code.replace(marker,body+marker)
        old='''        if(unpack(mSlots24.data()+3*i)!=slot || unpack(mOffsets24.data()+3*slot)!=mOffsets32[slot] ||
           (slot/tileBlocks())*tileBlocks()+mOffsets32[slot]!=outer)
            throw std::runtime_error("packed route mismatch");'''
        new=f'''        if(mOuter==Outer::Bch128x32) {{
            const auto base=(slot/tileBlocks())*tileBlocks();
            const auto sub=std::min<std::size_t>({small},tileBlocks());
            const auto intermediate=mOffsets32[2*slot];
            const auto* p=mOffsets24.data()+5*(base+intermediate)+3;
            const auto local=u32(p[0])|(u32(p[1])<<8);
            if(intermediate>=tileBlocks() || local>=sub ||
               unpack(mSlots24.data()+3*i)!=slot || unpack(mOffsets24.data()+5*slot)!=intermediate ||
               mOffsets32[2*(base+intermediate)+1]!=local ||
               base+(intermediate/sub)*sub+local!=outer)
                throw std::runtime_error("hierarchical route mismatch");
        }} else {{
'''+old+'''\n        }'''
        assert old in code
        code=code.replace(old,new)
        code=code.replace('"t128_s19_asymmetric_greedy3_2_r1"',f'"{name}"')
        directory=PARENT/'generated'/name
        directory.mkdir(exist_ok=True)
        (directory/'CandidateSpin.cpp').write_text(code,encoding='utf-8',newline='\n')
        (directory/'AsymmetricMap.h').write_text((SOURCE.parent/'AsymmetricMap.h').read_text(),encoding='utf-8',newline='\n')
        variants.append(dict(name=name,sub_blocks=small,sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in directory.iterdir() if p.is_file()}))
    (HERE/'HIERARCHICAL.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),candidates=variants),indent=2)+'\n')


if __name__=='__main__':main()
