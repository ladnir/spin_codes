"""Binary64 rank-four screen retaining total outer weight as well as support.

Integer outer counts and exact invariant tests; numerical moments are not
outward-certified. No benchmark and no production change.
"""
import argparse
from itertools import combinations,product
from math import comb, log, log2
import numpy as np
from flint import fmpz_poly

from bch_joint_support import authenticated_caps, support_caps
from basis_lattice import improve_caps
from shortened_bound import dimension_caps
from joint_support import span, gf2_rank
from rank_two_moment import census, region_transfers
from rank_three_flags import TILTS


def support_weight_counts(spectrum,caps,step):
    n=len(spectrum)-1
    upper=np.arange(step,4*n+step,step)
    result=np.zeros((len(upper),n+1),dtype=object)
    for u in range(n+1):
        coefficients=(fmpz_poly([0]+spectrum[1:u+1])**4).coeffs()
        for j,hi in enumerate(upper):
            lo=max(0,int(hi)-step+1)
            count=sum(int(x) for x in coefficients[lo:int(hi)+1])
            result[j,u]=min(caps[u],count)
    return upper,result


def refine_flag_counts(spectrum,caps3,caps4,dimensions,upper,counts,step,h_shell=None):
    n=len(spectrum)-1
    outer_cdf=np.cumsum(np.array(spectrum,dtype=object))
    result=counts.copy()
    for u in range(n+1):
        by_v=[]
        for v in range(u+1):
            spaces=caps3[v]//2520
            if h_shell is not None:
                spaces=min(spaces,h_shell[v])
            if not spaces:
                by_v.append(0)
                continue
            assert dimensions[v]>=3
            extensions=(comb(n-v,u-v)*(1<<dimensions[v])-(8 if u==v else 0))//8
            extensions=min(extensions,int(outer_cdf[(2*u-v)//2]))
            by_v.append(min(caps4[u],1344*spaces*extensions))
        prefix=[0]
        for count in by_v:
            prefix.append(prefix[-1]+count)
        for j,hi in enumerate(upper):
            lo=int(hi)-step+1
            # Exactly u-v columns are all ones; the other v occupied columns
            # have weights 1,2,3. Thus 4u-3v <= W <= 4u-v.
            first=max(0,-(-(4*u-int(hi))//3))
            last=min(u,4*u-lo)
            flag=prefix[last+1]-prefix[first] if first<=last else 0
            result[j,u]=min(result[j,u],flag)
    return result


def self_test():
    for rows in ([1,2,4,8],[0x97,0x4b,0x2d,0x1e]):
        words=span(rows)
        n=max(words).bit_length()
        spectrum=[sum(w.bit_count()==i for w in words) for i in range(n+1)]
        exact=np.zeros((4*n+1,n+1),dtype=object)
        for tup in product(words[1:],repeat=4):
            union=tup[0]|tup[1]|tup[2]|tup[3]
            total=sum(w.bit_count() for w in tup)
            assert total==sum(sum((w>>i)&1 for w in tup) for i in range(n))
            if gf2_rank(tup)==4:
                exact[total,union.bit_count()]+=1
        caps=np.cumsum(np.sum(exact,axis=0))
        hspaces={tuple(sorted(span(b))) for b in combinations(words[1:],3) if gf2_rank(b)==3}
        hcounts=[0]*(n+1)
        for h in hspaces:
            union=0
            for w in h:
                union|=w
            hcounts[union.bit_count()]+=1
        caps3=np.cumsum(np.array(hcounts,dtype=object))*2520
        dimensions=[0]*(n+1)
        for support in range(1<<n):
            size=sum(not(w&~support) for w in words)
            dimensions[support.bit_count()]=max(dimensions[support.bit_count()],size.bit_length()-1)
        for step in (1,2,4):
            upper,counts=support_weight_counts(spectrum,caps,step)
            counts=refine_flag_counts(spectrum,caps3,caps,dimensions,upper,counts,step)
            for j,hi in enumerate(upper):
                for u in range(n+1):
                    assert sum(exact[max(0,hi-step+1):hi+1,u])<=counts[j,u]
    print('Joint support/total-weight counts passed exhaustive small-code checks',flush=True)


def probabilities(data,upper,tilts,penalties,memory=False,columns=4):
    regions=256//columns
    if memory:
        import memory_moment
        extra=memory_moment.memory_data(data,windows=32//columns)
        memory_moment.self_test(data,extra)
    best=np.ones((len(upper),regions,columns,257))
    for tilt in map(float,tilts):
        if memory:
            transfers=memory_moment.regions(data,extra,tilt,steps=64*columns)
        else:
            shapes,empty,active=region_transfers(data,tilt,windows=32//columns,epochs=64*columns)
        for penalty in penalties:
            size=7+len(extra[0]) if memory else 7
            if not memory:
                weighted=active*np.array([penalty**sum(s) for s in shapes])[:,None,None]
                region=[empty]+[np.max(weighted[[i for i,s in enumerate(shapes) if sum(w!=0 for w in s)==b]],axis=0)
                                for b in range(1,columns+1)]
                first=np.stack([r[0] for r in region[1:]],axis=1)
                transposed=[r.T*comb(columns,b) for b,r in enumerate(region)]
            current=np.ones((1,size))
            for l in range(regions):
                if memory:
                    applied=memory_moment.actions(current,transfers,penalty)
                    moments=np.stack([a[:,0] for a in applied[1:]],axis=1)
                else:
                    moments=current@first
                for b in range(1,columns+1):
                    last=min(256,columns*l+b)
                    if last<72:
                        continue
                    u=np.arange(72,last+1)
                    base=np.log(np.maximum(moments[u-b,b-1],1e-300))+tilt*209715-np.array([log(comb(columns*l,int(w-b))) for w in u])
                    bounds=np.exp(np.minimum(0,base[None,:]-upper[:,None]*log(penalty)))
                    best[:,l,b-1,72:last+1]=np.minimum(best[:,l,b-1,72:last+1],bounds)
                if l==regions-1:
                    break
                shifted=[a*comb(columns,b) for b,a in enumerate(applied)] if memory else [current@r for r in transposed]
                current=np.zeros((len(current)+columns,size))
                for b,row in enumerate(shifted):
                    current[b:b+len(row)]+=row
        print('Total-weight diagnostic tilt',tilt,flush=True)
    result=np.ones((len(upper),257))
    for u in range(72,257):
        result[:,u]=sum(comb(columns,b)*comb(columns*l,u-b)/comb(256,u)*best[:,l,b-1,u]
                        for l in range(regions) for b in range(1,columns+1) if 0<=u-b<=columns*l)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--step',type=int,default=16,choices=(1,2,4,8,16,32))
    parser.add_argument('--fine',action='store_true')
    parser.add_argument('--memory',action='store_true')
    parser.add_argument('--flags',action='store_true')
    parser.add_argument('--oa',action='store_true')
    parser.add_argument('--full-window',action='store_true')
    parser.add_argument('--columns',type=int,default=4,choices=(2,4))
    args=parser.parse_args()
    if args.full_window and args.columns!=4:
        parser.error('full-window screen currently supports four columns only')
    print('Screen settings:',vars(args),flush=True)
    self_test()
    spectrum=authenticated_caps()
    dimensions=dimension_caps()
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimensions),spectrum)
    upper,counts=support_weight_counts(spectrum,caps[3],args.step)
    h_shell=None
    if args.oa:
        import joint_oa_caps
        joint_oa_caps.self_test()
        oa4=joint_oa_caps.shell_caps(256,128,4)
        oa3=joint_oa_caps.shell_caps(256,128,3)
        h_shell=[x//168 for x in oa3]
        for u in range(257):
            for j in range(len(upper)):
                counts[j,u]=min(counts[j,u],oa4[u])
    if args.flags:
        counts=refine_flag_counts(spectrum,caps[2],caps[3],dimensions,upper,counts,args.step,h_shell)
    tilts=TILTS if args.fine else TILTS[::2]
    penalties=(1,.98,.95,.9,.875,.8,.75,.625,.5,.375,.25) if args.fine else (1,.95,.875,.75,.5,.25)
    if args.columns==2:
        import two_column_moment
        data=two_column_moment.census()
    else:
        data=census()
    if args.full_window:
        from full_window_moment import lift
        data=lift(data)
    p=probabilities(data,upper,tilts,penalties,args.memory,args.columns)
    terms=[(2048*int(counts[j,u])*p[j,u],u,int(w),log2(int(counts[j,u])),log2(p[j,u]))
           for j,w in enumerate(upper) for u in range(72,257) if counts[j,u] and p[j,u]]
    total=sum(x for x,*_ in terms)
    print('Total-weight rank-four log2 union',log2(total),flush=True)
    print('Dominant (u,weight upper,log2 term,log2 count,log2 p)',[(u,w,log2(x),c,p) for x,u,w,c,p in sorted(terms,reverse=True)[:12]],flush=True)
    print('Binary64 diagnostic only; no rank-four certificate or multiple-group coverage.')


if __name__=='__main__':
    main()
