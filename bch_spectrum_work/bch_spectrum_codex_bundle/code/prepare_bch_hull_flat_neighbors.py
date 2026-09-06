"""Inputs for candidate neighbors of verified low-weight orbit seeds."""
import json
import struct
from pathlib import Path
from affine_wambach import gf_pow
from certify_bch_hull_dual_weight16 import subspaces
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/hull_flat_neighbors20'
if __name__=='__main__':
    coords=[gf_pow(2,i) for i in range(255)]+[0]
    pos={x:i for i,x in enumerate(coords)}
    bits=[1<<pos[x] for x in range(256)]
    flats=[sum(bits[x] for x in members) for members in subspaces(4)]
    assert len(set(flats))==len(flats)==200787
    p14=ROOT/'generated/bch256_hull_dual_weight14.json'
    p18=ROOT/'generated/bch256_hull_dual_weight18_complete.json'
    seeds=[int(r['seed_hex'],16) for p in (p14,p18) for r in json.loads(p.read_text())['orbits']]
    assert len(seeds)==98
    binary=struct.pack('<II',len(flats),len(seeds))
    binary+=b''.join(w.to_bytes(32,'little') for w in flats+seeds)
    FOLDER.mkdir(exist_ok=False)
    with (FOLDER/'input.bin').open('xb') as f:
        f.write(binary)
    value=dict(linear_four_spaces=len(flats),seed_orbits=len(seeds),candidate_search_not_claimed_exhaustive=True,
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),ROOT/'code/certify_bch_hull_dual_weight16.py',p14,p18,FOLDER/'input.bin')})
    write_new(FOLDER/'inputs.json',value)
    print(json.dumps(value,indent=2))
