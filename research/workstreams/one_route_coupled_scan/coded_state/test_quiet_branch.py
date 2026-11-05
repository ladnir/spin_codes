"""Small-field checks for column-scaled B and all-quiet field descent."""
import itertools
from check_maps import Field, Maps, weight
from quiet_certificate import KAPPA, cover, N


def main():
    e = Field(2,0x7)
    maps = Maps(e,1)
    b = maps.b_matrix[0]
    fixtures = 0
    for x in itertools.product(range(4),repeat=2):
        if not any(x):
            continue
        quiet = 0
        for scales in itertools.product(range(1,4),repeat=2):
            scaled = [e.table[a][z] for a,z in zip(scales,x)]
            quiet += maps.bmap(scaled) == [0]
        assert quiet <=3  # 1/(|E|-1) with s=1
        if weight(x)==1:
            assert quiet==0
        fixtures += 1
    # Uniform state refresh survives every nonzero column scaling.
    refresh = 0
    for scales in itertools.product(range(1,4),repeat=2):
        for q in range(1,4):
            aq = maps.a([q])
            for x in itertools.product(range(4),repeat=2):
                counts = [0]*4
                for d in itertools.product(range(4),repeat=2):
                    y = [a ^ e.table[c][z] for a,c,z in zip(x,d,aq)]
                    next_q = q ^ maps.bmap([e.table[r][z] for r,z in zip(scales,y)])[0]
                    counts[next_q] +=1
                assert counts==[4]*4
                refresh +=1
    print(f'Quiet-probability fixtures={fixtures}; scaled refresh fixtures={refresh}: PASS',flush=True)

    # GF4 is embedded explicitly into GF16; no reuse of numerical byte labels.
    f = Field(4,0x13)
    root = next(a for a in range(2,16) if f.table[a][a]^a^1==0)
    embedding = [0,1,root,root^1]
    for a,c in itertools.product(range(4),repeat=2):
        assert embedding[e.table[a][c]] == f.table[embedding[a]][embedding[c]]
    # Two RS[4,2] outer blocks; region j contains the two j-th evaluations.
    # Identity regional orders and unit scales deliberately give a nonzero
    # quiet subcode, so the field-descent test is not vacuous.
    generator = []
    for point in range(4):
        generator.extend(([1,point,0,0],[0,0,1,point]))
    def minimum(field,embed):
        matrix = [[embed[z] for z in row] for row in generator]
        bb = [embed[z] for z in b]
        least,words = 9,0
        for message in itertools.product(range(field.size),repeat=4):
            if not any(message):
                continue
            x = field.apply(matrix,message)
            if any(field.table[bb[0]][x[2*j]] ^ field.table[bb[1]][x[2*j+1]] for j in range(3)):
                continue
            least = min(least,weight(x))
            words +=1
        return least,words
    base = minimum(e,list(range(4)))
    extension = minimum(f,embedding)
    assert base==(6,15) and extension==(6,255),(base,extension)
    print(f'Nonvacuous quiet-code field descent: GF4 {base}, GF16 {extension}: PASS',flush=True)
    assert KAPPA.numerator==1279 and KAPPA.denominator==20
    assert 255**160 >= 2**1279
    boxes,cells = cover(491*N//1000)
    assert cells==15714
    assert sum(hi-lo+1 for lo,hi,_,_ in boxes)==cells
    print('Integer coefficient bound and all-quiet coverage count PASS')


if __name__=='__main__':
    main()
