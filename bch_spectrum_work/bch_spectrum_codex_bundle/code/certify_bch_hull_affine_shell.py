"""General exact orbit certificate for candidates in the known hull supercode."""
import argparse
import hashlib
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


def build(weight,seed_folder):
    assert weight%2==0 and 0<weight<128 and seed_folder.isidentifier()
    invariance=invariance_build()
    assert invariance==json.loads((ROOT/'generated/bch256_hull_dual_weight14.json').read_text())
    assert invariance['affine_invariance_generators_verified']
    model,_=anchor_build()
    anchor=model['metadata']['Lhull_anchor']
    basis_a=[int(x,16) for x in anchor['hull_basis_hex']]
    structure=json.loads((ROOT/'generated/bch256_hull_structure.json').read_text())
    basis_q=[int(x,16) for x in structure['codes']['HQ']['basis_hex']]
    expected=int(anchor['dual_half_spectrum'][str(weight)])
    seed_path=ROOT/'generated'/seed_folder/'seeds.bin'
    data=seed_path.read_bytes()
    assert len(data)%32==0
    seeds={int.from_bytes(data[i:i+32],'little') for i in range(0,len(data),32)}
    assert len(seeds)*32==len(data) and seeds
    assert all(w.bit_count()==weight and all((w&r).bit_count()%2==0 for r in basis_a) for w in seeds)
    coords=[gf_pow(2,i) for i in range(255)]+[0]
    pos={x:i for i,x in enumerate(coords)}
    bitlut=np.zeros((256,4),dtype='<u8')
    for x in range(256):
        i=pos[x]
        bitlut[x,i//64]=np.uint64(1)<<np.uint64(i%64)
    mult=np.asarray([[gf_mul(a,x) for x in range(256)] for a in range(256)],dtype=np.uint8)
    assert all(np.unique(mult[a]).size==256 for a in range(1,256))
    shifts=np.arange(256,dtype=np.uint8)[None,:]
    remaining=seeds.copy()
    images=np.zeros((255*256,4),dtype='<u8')
    records=[]
    total=inside=0
    while remaining:
        seed=min(remaining)
        support=[coords[i] for i in range(256) if (seed>>i)&1]
        in_q=all((seed&r).bit_count()%2==0 for r in basis_q)
        images.fill(0)
        for x in support:
            y=np.bitwise_xor(mult[1:,x,None],shifts).reshape(-1)
            images|=bitlut[y]
        unique=np.unique(images,axis=0)
        words={int(a)|(int(b)<<64)|(int(c)<<128)|(int(d)<<192) for a,b,c,d in unique}
        assert len(words)==len(unique) and seed in words
        assert all(w.bit_count()==weight for w in words)
        hit=len(remaining&words)
        remaining-=words
        total+=len(words)
        inside+=len(words) if in_q else 0
        assert total<=expected
        records.append(dict(seed_hex=hex(seed),field_support=support,orbit_size=len(words),
            input_seeds_in_orbit=hit,contained_in_HQdual=in_q,
            word_array_sha256=hashlib.sha256(unique.tobytes()).hexdigest()))
    return dict(classification='Exact affine-orbit verification, with completeness checked against the known shell count',
        weight=weight,seed_folder=seed_folder,candidate_seeds=len(seeds),affine_orbits=len(records),orbits=records,
        known_Adual_shell_count=expected,distinct_Adual_words_verified=total,
        full_Adual_shell_enumerated=total==expected,HQdual_shell_count_exact=inside if total==expected else None,
        HQdual_shell_count_lower_bound=inside,seed_search_completeness_assumed=False,
        orbit_digest_convention='unique 4-column little-endian uint64 rows sorted lexicographically low limb first',
        original_M22_target_closed=False,source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),ROOT/'code/affine_wambach.py',ROOT/'code/certify_bch_hull_dual_weight14.py',
             ROOT/'code/prepare_bch_hull_anchor_probe.py',ROOT/'generated/bch256_hull_dual_weight14.json',
             ROOT/'generated/bch256_hull_structure.json',seed_path)})


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('weight',type=int)
    parser.add_argument('seed_folder')
    parser.add_argument('receipt_name')
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    assert args.receipt_name.isidentifier()
    output=ROOT/'generated'/(args.receipt_name+'.json')
    value=build(args.weight,args.seed_folder)
    if args.verify:
        assert json.loads(output.read_text())==value
    else:
        write_new(output,value)
    print(json.dumps({k:v for k,v in value.items() if k not in ('orbits','source_sha256')},indent=2))
