"""Exact translation automorphisms preserve P and its p37 syndrome fibers."""
import argparse
from pathlib import Path
import sys
import bridge as base
sys.path.insert(0,str(base.BCH/'code'))
from bch_quotient import generator_polynomial,binary_poly_divmod
from affine_wambach import gf_pow,ALPHA,MODULUS


def run(verify=False):
    g=generator_polynomial(37);assert 255-(g.bit_length()-1)==131
    basis=[g<<i for i in range(131)]
    basis=[w|((w.bit_count()%2)<<255) for w in basis]
    coordinates=[gf_pow(ALPHA,i) for i in range(255)]+[0]
    assert len(set(coordinates))==256
    inverse={x:i for i,x in enumerate(coordinates)}
    powers=[gf_pow(x,37) for x in coordinates]
    def syndrome(word):
        value=0
        while word:
            bit=word&-word;value^=powers[bit.bit_length()-1];word^=bit
        return value
    checked=0
    for b in (1<<i for i in range(8)):
        permutation=[inverse[x^b] for x in coordinates]
        assert sorted(permutation)==list(range(256))
        for word in basis:
            moved=sum(((word>>i)&1)<<permutation[i] for i in range(256))
            assert moved.bit_count()%2==0
            assert binary_poly_divmod(moved&((1<<255)-1),g)[1]==0
            assert syndrome(moved)==syndrome(word)
            checked+=1
    result=dict(status='EXACT_BCH_TRANSLATION_SYMMETRY',field_modulus=MODULUS,
        P_dimension=131,translation_generators=[1<<i for i in range(8)],basis_checks=checked,
        p37_preserved=True,consequence='Every p37 fiber in P is preserved. C and every weight-defined subset of C are coordinate-transitive.',
        local_sha256={str(Path(__file__).relative_to(base.HERE)):base.sha(Path(__file__))},
        outer_sha256={str(p.relative_to(base.BCH)):base.sha(p) for p in [base.BCH/'code/bch_quotient.py',base.BCH/'code/affine_wambach.py']})
    path=base.HERE/'generated/bch_translation_symmetry.json'
    if verify:assert result==base.read(path)
    else:base.write_new(path,result)
    print('Verified',checked,'basis images: all 256 translations preserve the fixed C and its tail T.',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
