"""Outward output-weight bound for one rank-one four-row group.

Uniform shared 256-coordinate shuffle, four-column macroregions, IMT(128,19).
Includes all 2048 possible groups and their nonempty row subsets. Does not
cover higher-rank tuples or messages spanning groups. Writes no data files.
"""
import argparse
from collections import Counter,defaultdict
from itertools import product
from math import comb
from flint import arb,arb_mat,ctx,fmpq_mat

import bch_joint_support as joint
from group_moment import maps
from group_rank_one_verify import TILTS,up


def shapes_and_masks():
    shapes=[]
    masks=defaultdict(list)
    for mask in range(65536):
        shape=tuple(sorted(((mask>>(4*j))&15).bit_count() for j in range(4)))
        shapes.append(shape)
        active=[v for v in shape if v]
        if active and len(set(active))==1:
            masks[active[0],len(active)].append(mask)
    assert len(masks)==16
    for (a,b),row in masks.items():
        assert len(row)==comb(4,b)*comb(4,a)**b
    return shapes,masks


def intersections(shape,a):
    """Count masks by occupied-column count and intersection with shape."""
    counts={(0,0):1}
    for weight in shape:
        updated=Counter(counts)  # do not activate this column
        for (b,overlap),count in counts.items():
            for r in range(max(0,a+weight-4),min(a,weight)+1):
                updated[b+1,overlap+r]+=count*comb(weight,r)*comb(4-weight,a-r)
        counts=updated
    return counts


def census():
    images,columns,spectrum=maps()
    shapes,masks=shapes_and_masks()
    windows=defaultdict(Counter)
    for image in images[1:]:
        weight=image.bit_count()
        for start in range(0,128,16):
            windows[weight][shapes[(image>>start)&65535]]+=1
    for v,count in spectrum.items():
        assert sum(windows[v].values())==8*count
    moments={key:defaultdict(Counter) for key in masks}
    for a in range(1,5):
        for shape in set(shapes):
            counts=intersections(shape,a)
            for v in spectrum:
                multiplicity=windows[v][shape]
                if not multiplicity:
                    continue
                for (b,overlap),count in counts.items():
                    if b:
                        moments[a,b][v][v+a*b-2*overlap]+=multiplicity*count
    maximum={}
    cancellation={key:defaultdict(Counter) for key in masks}
    for key,allowed in sorted(masks.items()):
        atoms=Counter()
        for start in range(0,128,16):
            for mask in allowed:
                remaining,value=mask,0
                while remaining:
                    bit=remaining&-remaining
                    value^=columns[start+bit.bit_length()-1]
                    remaining^=bit
                assert value  # every nonempty rank-one block activates zero
                atoms[value]+=1
                image=images[value]
                cancellation[key][image.bit_count()][(image^(mask<<start)).bit_count()]+=1
        maximum[key]=max(atoms.values())
        choices=8*len(allowed)
        assert sum(atoms.values())==choices
        assert sum(sum(row.values()) for row in cancellation[key].values())==choices
        for v,count in spectrum.items():
            assert sum(moments[key][v].values())==choices*count
        print('Exact orbit',key,'choices',choices,'maximum syndrome multiplicity',maximum[key],flush=True)
    return spectrum,moments,maximum,cancellation


def transfers(data,tilt,a):
    spectrum,moments,maximum,cancellation=data
    levels=sorted(spectrum)
    n,m=len(levels)+2,(1 << 19)-1
    powers=[up((-arb(tilt)*w).exp()) for w in range(145)]
    empty=[[arb(0) for _ in range(n)] for _ in range(n)]
    empty[0][0]=arb(1)
    empty[1][1]=powers[min(levels)]/2
    for j,w in enumerate(levels):
        empty[1][j+2]=powers[min(levels)]*spectrum[w]/(2*m)
    for i,v in enumerate(levels):
        empty[i+2][i+2]+=powers[v]/2
        for j,w in enumerate(levels):
            empty[i+2][j+2]+=powers[v]*spectrum[w]/(2*m)
    result=[arb_mat([[up(x) for x in row] for row in empty])]
    for b in range(1,5):
        key=a,b
        choices=8*comb(4,b)*comb(4,a)**b
        active=[[arb(0) for _ in range(n)] for _ in range(n)]
        active[0][1]=powers[a*b]
        d=powers[min(levels)-a*b]
        min_cancel=min(w for row in cancellation[key].values() for w in row)
        cd=min(up(d),up(powers[min_cancel]*maximum[key]/choices))
        active[1][0]=cd/2+d/(2*m)
        active[1][1]=d/2
        for j,w in enumerate(levels):
            active[1][j+2]=d*spectrum[w]/(2*m)
        for i,v in enumerate(levels):
            moment=up(sum((count*powers[w] for w,count in moments[key][v].items()),arb(0))/(choices*spectrum[v]))
            cancel=up(sum((count*powers[w] for w,count in cancellation[key][v].items()),arb(0))/(choices*spectrum[v]))
            active[i+2][0]=cancel/2+moment/(2*m)
            active[i+2][1]=moment/2
            for j,w in enumerate(levels):
                active[i+2][j+2]=moment*spectrum[w]/(2*m)
        result.append(arb_mat([[up(x) for x in row] for row in active]))
    return result


def regions(epoch,epochs=256):
    n=epoch[0].nrows()
    empty=arb_mat([[int(i==j) for j in range(n)] for i in range(n)])
    active=[arb_mat(n,n) for _ in epoch[1:]]
    for _ in range(epochs):
        active=[old*epoch[0]+empty*one for old,one in zip(active,epoch[1:])]
        empty=empty*epoch[0]
    result=[empty]+[row/epochs for row in active]
    return [arb_mat([[up(row[i,j]) for j in range(n)] for i in range(n)]) for row in result]


def conditional_numerators(region,length=64,width=4,matrix=arb_mat,round_up=up):
    """Yield [future support weight, first active occupancy b-1] matrices."""
    assert len(region)==width+1
    n=region[0].nrows()
    current=matrix([[1]*n])
    first=matrix([[region[b][0,j] for b in range(1,width+1)] for j in range(n)])
    transposed=[row.transpose()*comb(width,b) for b,row in enumerate(region)]
    for remaining in range(length):
        yield remaining,current*first
        if remaining+1==length:
            break
        shifted=[current*row for row in transposed]
        count=current.nrows()
        current=matrix([[round_up(sum((row[w-b,j] for b,row in enumerate(shifted) if 0<=w-b<count),0))
                         for j in range(n)] for w in range(count+width)])


def self_test():
    shapes,masks=shapes_and_masks()
    for pattern in (0,65535,0xaaaa,0x5555,0x1234,0x00ff,0xff00,0x8001,0x369c,0x1357):
        for a in range(1,5):
            actual=intersections(shapes[pattern],a)
            for b in range(1,5):
                expected=Counter((pattern&mask).bit_count() for mask in masks[a,b])
                assert expected=={r:count for (bb,r),count in actual.items() if bb==b and count}
    for width,length in ((2,4),(4,2),(4,64)):
        for w in range(1,width*length+1):
            assert sum(comb(width,b)*comb(width*l,w-b) for l in range(length) for b in range(1,width+1)
                       if 0<=w-b<=width*l)==comb(width*length,w)
    region=[fmpq_mat([[1,0],[0,2]]),fmpq_mat([[1,2],[3,4]]),fmpq_mat([[2,3],[1,4]])]
    # Independently enumerate the active epoch's two possible positions.
    epoch=[arb_mat([[int(row[i,j]) for j in range(2)] for i in range(2)]) for row in region]
    aggregated=regions(epoch,epochs=2)
    assert aggregated[0]==epoch[0]*epoch[0]
    for b in (1,2):
        assert aggregated[b]==(epoch[0]*epoch[b]+epoch[b]*epoch[0])/2
    for remaining,values in conditional_numerators(region,4,2,fmpq_mat,lambda x:x):
        direct=[[0]*2 for _ in range(2*remaining+1)]
        for support in product((0,1),repeat=2*remaining):
            tail=fmpq_mat([[1],[1]])
            for i in reversed(range(remaining)):
                tail=region[sum(support[2*i:2*i+2])]*tail
            for b in (1,2):
                direct[sum(support)][b-1]+=(region[b]*tail)[0,0]
        assert values==fmpq_mat(direct)
    print('Exact orbit intersections, region averaging, first-macroregion counts, and small transfer enumeration passed',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--rows',type=int,nargs='+',default=[1,2,3,4])
    parser.add_argument('--test-only',action='store_true')
    args=parser.parse_args()
    assert args.precision>=128 and args.rows and all(1<=a<=4 for a in args.rows)
    assert len(set(args.rows))==len(args.rows)
    self_test()
    if args.test_only:
        return
    spectrum=joint.authenticated_caps()
    data=census()
    # Authentication imports historical modules that change FLINT precision.
    ctx.prec=args.precision
    assert ctx.prec==args.precision
    weights=[w for w in range(1,257) if spectrum[w]]
    total=[[[arb(0) for _ in weights] for _ in range(4)] for _ in range(64)]
    for a in args.rows:
        best=[[[arb(1) for _ in weights] for _ in range(4)] for _ in range(64)]
        for tilt in TILTS:
            factor=(arb(tilt)*209715).exp()
            for remaining,values in conditional_numerators(regions(transfers(data,tilt,a))):
                for b in range(1,5):
                    for index,w in enumerate(weights):
                        if 0<=w-b<=4*remaining:
                            value=up(values[w-b,b-1]*factor/comb(4*remaining,w-b))
                            best[remaining][b-1][index]=min(best[remaining][b-1][index],value)
            print('Outward replay',ctx.prec,'bits: row copies',a,'tilt',tilt,flush=True)
        for remaining in range(64):
            for b in range(4):
                for index in range(len(weights)):
                    total[remaining][b][index]=up(total[remaining][b][index]+comb(4,a)*best[remaining][b][index])
    upper=arb(0)
    shells=[]
    for index,w in enumerate(weights):
        probability=sum((comb(4,b)*comb(4*l,w-b)*min(arb(1),total[l][b-1][index])
                         for l in range(64) for b in range(1,5) if 0<=w-b<=4*l),arb(0))/comb(256,w)
        shell=up(2048*spectrum[w]*probability)
        shells.append((w,shell))
        upper+=shell
    upper=up(upper)
    print('Outward union upper:',upper,flush=True)
    print('Margin (display only):',-upper.log()/arb(2).log(),flush=True)
    print('Dominant shells:',sorted(shells,key=lambda p:float(p[1]),reverse=True)[:5],flush=True)
    if set(args.rows)=={1,2,3,4} and 0<upper<arb(2)**-40:
        print('VERIFIED: rank-one, one-active-group output weight <=209715 union < 2^-40',flush=True)
    else:
        print('The complete rank-one 40-bit target has NOT been verified by this run.',flush=True)
    print('Full SPIN certificate: NO (higher ranks and multiple active groups remain).')


if __name__=='__main__':
    main()
