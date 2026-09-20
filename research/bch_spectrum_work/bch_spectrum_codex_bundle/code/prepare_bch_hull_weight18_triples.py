"""Prepare exact subspace indicators, syndromes, and scalar-orbit representatives."""
import json
import struct
from pathlib import Path
from certify_bch_hull_dual_weight16 import subspaces
from prepare_bch_hull_anchor_probe import build as anchor_build
from affine_wambach import gf_pow
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/hull_weight18_triples'


def build():
    model,_=anchor_build()
    basis=[int(x,16) for x in model['metadata']['Lhull_anchor']['hull_basis_hex']]
    coords=[gf_pow(2,i) for i in range(255)]+[0]
    pos={x:i for i,x in enumerate(coords)}
    bits=[1<<pos[x] for x in range(256)]
    columns=[sum(((word>>pos[x])&1)<<j for j,word in enumerate(basis)) for x in range(256)]
    words=[]
    syndromes=[]
    for members in subspaces(3):
        w=s=0
        for x in members:
            w|=bits[x]
            s^=columns[x]
        assert w.bit_count()==8 and s!=0
        words.append(w)
        syndromes.append(s)
    index={w:i for i,w in enumerate(words)}
    assert len(index)==len(words)==97155
    unseen=set(range(len(words)))
    reps=[]
    mask=(1<<255)-1
    while unseen:
        i=min(unseen)
        reps.append(i)
        w=words[i]
        orbit=set()
        for _ in range(255):
            orbit.add(index[w])
            w=(w&(1<<255))|(((w&mask)<<1)&mask)|((w>>254)&1)
        assert w==words[i] and len(orbit)==255 and orbit<=unseen
        unseen-=orbit
    assert len(reps)==381
    binary=struct.pack('<II',len(words),len(reps))
    binary+=b''.join(w.to_bytes(32,'little')+s.to_bytes(16,'little') for w,s in zip(words,syndromes))
    binary+=struct.pack('<'+'I'*len(reps),*reps)
    receipt=dict(classification='Exact indicators and scalar orbit partition for a candidate search',
        subspaces=len(words),scalar_orbits=len(reps),orbit_size=255,
        candidate_search_not_claimed_exhaustive_for_weight18=True,
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),ROOT/'code/certify_bch_hull_dual_weight16.py',
             ROOT/'code/prepare_bch_hull_anchor_probe.py',ROOT/'code/affine_wambach.py')})
    return binary,receipt


if __name__=='__main__':
    binary,receipt=build()
    FOLDER.mkdir(exist_ok=False)
    with (FOLDER/'input.bin').open('xb') as f:
        f.write(binary)
    receipt['input_sha256']=sha(FOLDER/'input.bin')
    write_new(FOLDER/'inputs.json',receipt)
    print(json.dumps(receipt,indent=2))
