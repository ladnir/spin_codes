"""Enumerate a known weight-14 shell via pairs of binary three-dimensional flats.

No classification of minimum words is assumed: completeness is checked against
the exact known shell cardinality of hull(L71)^perp.
"""
import hashlib
import itertools
import json
import sys
from pathlib import Path
sys.dont_write_bytecode=True
from audit_bch_hulls import contains,echelon
from affine_wambach import gf_mul,gf_pow
from prepare_bch_hull_anchor_probe import build as anchor_build
from run_higher_endpoint_preflight import sha,write_new
ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'generated/bch256_hull_dual_weight14.json'


def build():
    model,_=anchor_build()
    anchor=model['metadata']['Lhull_anchor']
    basis_a=[int(x,16) for x in anchor['hull_basis_hex']]
    structure=json.loads((ROOT/'generated/bch256_hull_structure.json').read_text())
    basis_q=[int(x,16) for x in structure['codes']['HQ']['basis_hex']]
    basis_p=[int(x,16) for x in structure['codes']['HP']['basis_hex']]
    expected=int(anchor['dual_half_spectrum']['14'])
    assert expected==65280
    assert all(int(anchor['dual_half_spectrum'][str(w)])==0 for w in range(2,14,2))
    coordinates=[gf_pow(2,i) for i in range(255)]+[0]
    position={x:i for i,x in enumerate(coordinates)}
    bits=[1<<position[x] for x in range(256)]

    def columns(basis):
        return [sum(((word>>position[x])&1)<<j for j,word in enumerate(basis)) for x in range(256)]
    ca,cq,cp=(columns(b) for b in (basis_a,basis_q,basis_p))

    def permute(word,images):
        out=0
        while word:
            low=word&-word
            out|=bits[images[coordinates[low.bit_length()-1]]]
            word^=low
        return out
    for basis in (basis_a,basis_q,basis_p):
        for images in ([x^1 for x in range(256)],[gf_mul(2,x) for x in range(256)]):
            assert all(contains(basis,permute(word,images)) for word in basis)

    # Reduced row echelon bases with their leading bits ordered increasingly.
    buckets={}
    candidates=set()
    spaces=0
    collisions=0
    for pivots in itertools.combinations(range(8),3):
        options=[]
        for pivot in pivots:
            free=[j for j in range(pivot) if j not in pivots]
            options.append([(1<<pivot)|sum(((mask>>k)&1)<<j for k,j in enumerate(free))
                            for mask in range(1<<len(free))])
        for a,b,c in itertools.product(*options):
            members=(0,a,b,a^b,c,a^c,b^c,a^b^c)
            syndrome=0
            word=0
            for x in members:
                syndrome^=ca[x]
                word|=bits[x]
            assert word.bit_count()==8
            for other in buckets.get(syndrome,[]):
                candidate=word^other
                assert candidate.bit_count()==14
                candidates.add(candidate)
                collisions+=1
            buckets.setdefault(syndrome,[]).append(word)
            spaces+=1
    gaussian=(255*254*252)//(7*6*4)
    assert spaces==gaussian==97155
    assert candidates

    covered=set()
    inside_q=set()
    inside_p=set()
    orbit_records=[]
    for seed in sorted(candidates):
        if seed in covered:
            continue
        support=[coordinates[i] for i in range(256) if (seed>>i)&1]
        orbit=set()
        seed_q=all((seed&b).bit_count()%2==0 for b in basis_q)
        seed_p=all((seed&b).bit_count()%2==0 for b in basis_p)
        for a in range(1,256):
            scaled=[gf_mul(a,x) for x in support]
            for b in range(256):
                word=sa=sq=sp=0
                for x in scaled:
                    y=x^b
                    word|=bits[y]
                    sa^=ca[y]
                    sq^=cq[y]
                    sp^=cp[y]
                assert word.bit_count()==14 and sa==0
                assert (sq==0)==seed_q and (sp==0)==seed_p
                orbit.add(word)
                if sq==0:
                    inside_q.add(word)
                if sp==0:
                    inside_p.add(word)
        assert not (covered&orbit)
        covered.update(orbit)
        orbit_records.append(dict(seed_hex=hex(seed),field_support=support,orbit_size=len(orbit),
                                  contained_in_HQdual=seed_q,contained_in_HPdual=seed_p))
        assert len(covered)<=expected
        if len(covered)==expected:
            break
    digest=hashlib.sha256(b''.join(word.to_bytes(32,'little') for word in sorted(covered))).hexdigest()
    complete=len(covered)==expected
    # HQdual subset Adual, so a complete larger shell determines its subcode shell.
    assert all(contains(basis_q,x) for x in basis_a)
    return dict(classification='Exact affine-orbit shell enumeration relative to a published spectrum anchor',
                three_subspaces_enumerated=spaces,matching_syndrome_pairs=collisions,
                distinct_seed_words=len(candidates),orbits=orbit_records,
                known_Adual14=expected,distinct_Adual14_words_verified=len(covered),
                full_Adual14_shell_enumerated=complete,word_set_sha256=digest,
                HQdual_A14_exact=len(inside_q) if complete else None,
                HQdual_distance_at_least_16=complete and len(inside_q)==0,
                affine_invariance_generators_verified=True,original_M22_target_closed=False,
                source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
                    (Path(__file__),ROOT/'code/prepare_bch_hull_anchor_probe.py',
                     ROOT/'code/audit_bch_hulls.py',ROOT/'code/affine_wambach.py',
                     ROOT/'generated/bch256_hull_structure.json')})


if __name__=='__main__':
    result=build()
    if '--verify' in sys.argv:
        assert json.loads(OUTPUT.read_text())==result
    else:
        write_new(OUTPUT,result)
    print(json.dumps(result,indent=2))
