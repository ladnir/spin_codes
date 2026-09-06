"""Compare normalized-pair counts with explicit affine closure, including stabilizers."""
import json
import sys
from pathlib import Path
import numpy as np
from certify_bch_canonical_shell import field_tables,canonical
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]


def build():
    tables=field_tables()
    coords,bitlut,mult,inv=tables
    pos={int(x):i for i,x in enumerate(coords)}
    # Include subfields and affine flats with nontrivial stabilizers.
    supports=[[0,1],[0,1,2,3],list(range(8)),list(range(16)),
        [x for x in range(256) if int(mult[x,x])==x],
        [x for x in range(256) if int(mult[int(mult[x,x]),int(mult[x,x])])==x]]
    old=json.loads((ROOT/'generated/bch256_hull_dual_weight18_complete.json').read_text())
    supports.extend(r['field_support'] for r in old['orbits'][:3])
    results=[]
    for support in supports:
        word=sum(1<<pos[x] for x in support)
        best,stab=canonical(word,len(support),tables)
        images=np.zeros((255,256,4),dtype='<u8')
        shifts=np.arange(256,dtype=np.uint8)[None,:]
        for x in support:
            images|=bitlut[mult[1:,x,None]^shifts]
        unique=np.unique(images.reshape(-1,4),axis=0)
        assert len(unique)==65280//stab
        normalized=unique[((unique[:,0]&1)==1)&((unique[:,3]>>np.uint64(63))==1)]
        fullbest=normalized[np.lexsort(normalized.T)[0]]
        assert best==sum(int(fullbest[i])<<(64*i) for i in range(4))
        for a,b in ((1,1),(2,37),(191,255)):
            transformed=sum(1<<pos[int(mult[a,x])^b] for x in support)
            assert canonical(transformed,len(support),tables)==(best,stab)
        results.append(dict(weight=len(support),field_support=support,stabilizer=stab,orbit_size=len(unique)))
    # Every orbit in the explicit complete shell has a distinct canonical key.
    keys=set()
    for row in old['orbits']:
        key,stab=canonical(int(row['seed_hex'],16),18,tables)
        assert key not in keys and 65280//stab==row['orbit_size']
        keys.add(key)
    data=(ROOT/'generated/hull_weight18_canonical/canonical.bin').read_bytes()
    assert keys=={int.from_bytes(data[i:i+32],'little') for i in range(0,len(data),36)}
    return dict(classification='Exact regression of canonical normalization against full affine enumeration',
        cases=results,weight18_orbits_cross_checked=len(keys),all_passed=True,
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),ROOT/'code/certify_bch_canonical_shell.py',
             ROOT/'generated/bch256_hull_dual_weight18_complete.json',
             ROOT/'generated/hull_weight18_canonical/canonical.bin')})


if __name__=='__main__':
    value=build()
    output=ROOT/'generated/bch256_affine_canonical_regression.json'
    if '--verify' in sys.argv:
        assert value==json.loads(output.read_text())
    else:
        write_new(output,value)
    print(json.dumps({k:v for k,v in value.items() if k!='source_sha256'},indent=2))
