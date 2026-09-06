"""Independently verify canonical affine orbits and count their full union."""
import argparse
import json
import sys
from pathlib import Path
import numpy as np
sys.dont_write_bytecode=True
from affine_wambach import gf_mul,gf_pow
from certify_bch_hull_dual_weight14 import build as invariance_build
from prepare_bch_hull_anchor_probe import build as anchor_build
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]


def field_tables():
    coords=np.asarray([gf_pow(2,i) for i in range(255)]+[0],dtype=np.uint8)
    assert len(set(map(int,coords)))==256
    pos={int(x):i for i,x in enumerate(coords)}
    bitlut=np.zeros((256,4),dtype='<u8')
    for x in range(256):
        i=pos[x]; bitlut[x,i//64]=np.uint64(1)<<np.uint64(i%64)
    mult=np.asarray([[gf_mul(a,x) for x in range(256)] for a in range(256)],dtype=np.uint8)
    inv=np.zeros(256,dtype=np.uint8)
    for a in range(1,256):
        hit=np.flatnonzero(mult[a]==1)
        assert len(hit)==1
        inv[a]=hit[0]
    return coords,bitlut,mult,inv


def canonical(word,weight,tables):
    coords,bitlut,mult,inv=tables
    support=np.asarray([coords[i] for i in range(256) if word>>i&1],dtype=np.uint8)
    assert len(support)==weight
    first=np.repeat(support,weight)
    second=np.tile(support,weight)
    keep=first!=second
    first=first[keep]; second=second[keep]
    scalars=inv[first^second]
    mapped=mult[scalars[:,None],support[None,:]^first[:,None]]
    images=np.bitwise_or.reduce(bitlut[mapped],axis=1)
    # np.lexsort's final key is primary: high limb first for integer order.
    order=np.lexsort(images.T)
    best=images[order[0]]
    stabilizer=int(np.count_nonzero(np.all(images==best,axis=1)))
    value=sum(int(best[i])<<(64*i) for i in range(4))
    assert stabilizer>0 and 65280%stabilizer==0
    return value,stabilizer


def build(weight,folder_name):
    assert 0<weight<128 and weight%2==0 and folder_name.isidentifier()
    invariance=invariance_build()
    assert invariance==json.loads((ROOT/'generated/bch256_hull_dual_weight14.json').read_text())
    assert invariance['affine_invariance_generators_verified']
    model,_=anchor_build()
    anchor=model['metadata']['Lhull_anchor']
    basis_a=[int(x,16) for x in anchor['hull_basis_hex']]
    structure=json.loads((ROOT/'generated/bch256_hull_structure.json').read_text())
    basis_q=[int(x,16) for x in structure['codes']['HQ']['basis_hex']]
    expected=int(anchor['dual_half_spectrum'][str(weight)])
    source=ROOT/'generated'/folder_name/'canonical.bin'
    data=source.read_bytes()
    assert len(data)>0 and len(data)%36==0
    tables=field_tables()
    total=inside=0; previous=-1; records=[]
    for i in range(0,len(data),36):
        word=int.from_bytes(data[i:i+32],'little')
        stab=int.from_bytes(data[i+32:i+36],'little')
        assert word>previous and word.bit_count()==weight
        previous=word
        assert all((word&r).bit_count()%2==0 for r in basis_a)
        assert canonical(word,weight,tables)==(word,stab)
        size=65280//stab
        in_q=all((word&r).bit_count()%2==0 for r in basis_q)
        total+=size; inside+=size if in_q else 0
        assert total<=expected
        records.append(dict(canonical_hex=hex(word),stabilizer=stab,orbit_size=size,contained_in_HQdual=in_q))
    return dict(classification='Exact canonical affine-orbit verification, completeness checked against known shell size',
        weight=weight,folder_name=folder_name,affine_orbits=len(records),orbits=records,
        known_Adual_shell_count=expected,distinct_Adual_words_verified=total,
        full_Adual_shell_enumerated=total==expected,HQdual_shell_count_exact=inside if total==expected else None,
        HQdual_shell_count_lower_bound=inside,HQdual_shell_count_upper_bound=inside+expected-total,
        candidate_search_completeness_assumed=False,original_M22_target_closed=False,
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),ROOT/'code/affine_wambach.py',ROOT/'code/certify_bch_hull_dual_weight14.py',
             ROOT/'code/prepare_bch_hull_anchor_probe.py',ROOT/'generated/bch256_hull_dual_weight14.json',
             ROOT/'generated/bch256_hull_structure.json',source)})


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('weight',type=int)
    parser.add_argument('folder_name')
    parser.add_argument('receipt_name')
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    assert args.receipt_name.isidentifier()
    value=build(args.weight,args.folder_name)
    output=ROOT/'generated'/(args.receipt_name+'.json')
    if args.verify:
        assert value==json.loads(output.read_text())
    else:
        write_new(output,value)
    print(json.dumps({k:v for k,v in value.items() if k not in ('orbits','source_sha256')},indent=2))
