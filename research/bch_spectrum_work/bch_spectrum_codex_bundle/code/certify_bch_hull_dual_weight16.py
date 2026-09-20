"""Complete the known supercode's weight-16 shell by explicit affine flats."""
import itertools
import json
import sys
from pathlib import Path
sys.dont_write_bytecode=True
from affine_wambach import gf_pow
from certify_bch_hull_dual_weight14 import build as weight14_build
from prepare_bch_hull_anchor_probe import build as anchor_build
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'generated/bch256_hull_dual_weight16.json'


def subspaces(dimension):
    for pivots in itertools.combinations(range(8),dimension):
        options=[]
        for pivot in pivots:
            free=[j for j in range(pivot) if j not in pivots]
            options.append([(1<<pivot)|sum(((mask>>k)&1)<<j for k,j in enumerate(free))
                            for mask in range(1<<len(free))])
        for basis in itertools.product(*options):
            members=[0]
            for x in basis:
                members += [y^x for y in members]
            yield tuple(members)


def build():
    prior=weight14_build()
    assert prior==json.loads((ROOT/'generated/bch256_hull_dual_weight14.json').read_text())
    assert prior['HQdual_A14_exact']==0 and prior['affine_invariance_generators_verified']
    model,_=anchor_build()
    anchor=model['metadata']['Lhull_anchor']
    basis_a=[int(x,16) for x in anchor['hull_basis_hex']]
    structure=json.loads((ROOT/'generated/bch256_hull_structure.json').read_text())
    basis_q=[int(x,16) for x in structure['codes']['HQ']['basis_hex']]
    coords=[gf_pow(2,i) for i in range(255)]+[0]
    pos={x:i for i,x in enumerate(coords)}
    bits=[1<<pos[x] for x in range(256)]
    ca=[sum(((word>>pos[x])&1)<<j for j,word in enumerate(basis_a)) for x in range(256)]
    cq=[sum(((word>>pos[x])&1)<<j for j,word in enumerate(basis_q)) for x in range(256)]

    def describe(members):
        word=sa=sq=0
        for x in members:
            word|=bits[x]
            sa^=ca[x]
            sq^=cq[x]
        return word,sa,sq

    count4=inside4=0
    first4=None
    for members in subspaces(4):
        word,sa,sq=describe(members)
        assert word.bit_count()==16 and sa==0
        count4+=1
        if sq==0:
            inside4+=1
            if first4 is None:
                first4=hex(word)
    assert count4==200787
    # The already checked affine invariance gives all 16 cosets of each
    # four-space, with constant membership in each dual code.
    all_affine4=16*count4
    inside_affine4=16*inside4

    buckets={}
    pairs=[]
    count3=0
    for members in subspaces(3):
        word,sa,sq=describe(members)
        assert word.bit_count()==8 and sa!=0
        for other in buckets.get(sa,[]):
            assert set(members)&set(other)=={0}
            pairs.append((other,members,sa))
        buckets.setdefault(sa,[]).append(members)
        count3+=1
    assert count3==97155 and len(pairs)==255

    def cosets(members,expected):
        unseen=(1<<256)-1
        out=[]
        while unseen:
            rep=(unseen&-unseen).bit_length()-1
            translated=[x^rep for x in members]
            natural=sum(1<<x for x in translated)
            assert unseen&natural==natural
            unseen^=natural
            word,sa,sq=describe(translated)
            assert sa==expected
            out.append((rep,word,sq))
        assert len(out)==32
        return out

    special=set()
    inside_special=set()
    weight14=set()
    for u,v,syndrome in pairs:
        sumspace={x^y for x in u for y in v}
        assert len(sumspace)==64
        for a,wu,squ in cosets(u,syndrome):
            for b,wv,sqv in cosets(v,syndrome):
                word=wu^wv
                if a^b in sumspace:
                    assert word.bit_count()==14 and squ^sqv!=0
                    weight14.add(word)
                else:
                    assert word.bit_count()==16
                    # Its affine span has dimension seven: U,V span six,
                    # and the coset offset lies outside that sum. It is not
                    # an affine four-flat, so the two counted classes are disjoint.
                    special.add(word)
                    if squ^sqv==0:
                        inside_special.add(word)
    assert len(weight14)==65280
    assert len(special)==195840
    expected=int(anchor['dual_half_spectrum']['16'])
    assert all_affine4+len(special)==expected==3408432
    exact=inside_affine4+len(inside_special)
    return dict(classification='Complete exact supercode weight-16 enumeration; published shell cardinality is the completeness check',
                linear_four_spaces=count4,affine_four_flats=all_affine4,
                linear_three_spaces=count3,matching_three_space_pairs=len(pairs),
                skew_affine_three_flat_unions=len(special),two_classes_disjoint_by_affine_span=True,
                known_Adual16=expected,full_Adual16_shell_enumerated=True,
                HQdual16_affine_four_flat_count=inside_affine4,
                HQdual16_skew_union_count=len(inside_special),HQdual_A16_exact=exact,
                first_HQdual_four_flat_hex=first4,HQdual_distance_at_least_18=exact==0,
                original_M22_target_closed=False,
                source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
                    (Path(__file__),ROOT/'code/certify_bch_hull_dual_weight14.py',
                     ROOT/'generated/bch256_hull_dual_weight14.json',ROOT/'code/prepare_bch_hull_anchor_probe.py',
                     ROOT/'code/affine_wambach.py',ROOT/'generated/bch256_hull_structure.json')})


if __name__=='__main__':
    result=build()
    if '--verify' in sys.argv:
        assert json.loads(OUTPUT.read_text())==result
    else:
        write_new(OUTPUT,result)
    print(json.dumps(result,indent=2))
