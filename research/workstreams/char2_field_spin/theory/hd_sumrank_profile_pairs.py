"""FLOAT paired-packet witnesses retaining low-occupancy tile profiles.

High-occupancy profiles are not discarded. For each such count, replace the
conditional matrix by its entrywise maximum over ALL profiles at that count.
This enlarges every product. Fractional powers are applied after multiplying
two packet matrices, preserving the intermediate-state sum.
"""

from math import comb,log2
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp

from hd_sumrank_fast_screen import packet_logs,profile_data
from hd_randomized_outer_state_profile import power,total


def coefficients(z,alphas,*,state=8,mixer='feedback_refresh',cutoff=16,packet_source=None):
    a,j,probability,ranks,floor,choose=profile_data(state)
    raw=packet_logs(z,state,mixer) if packet_source is None else packet_source(z,state)
    keep=j<=cutoff
    degrees=np.concatenate((j[keep],np.arange(cutoff+1,65)))
    weight=np.concatenate((probability[keep]+choose[j[keep]],choose[cutoff+1:]))
    matrices=np.concatenate((raw[keep],np.stack([raw[j==k].max(axis=0) for k in range(cutoff+1,65)])))
    result=[np.zeros((129,2,2)) for alpha in alphas]
    # Bound peak memory by batching only the left profile index.
    for start in range(0,len(degrees),32):
        left=matrices[start:start+32]
        degree=(degrees[start:start+32,None]+degrees[None,:]).reshape(-1)
        probability=weight[start:start+32,None]+weight[None,:]
        for r in range(2):
            for s in range(2):
                entry=np.logaddexp2(left[:,r,0,None]+matrices[None,:,0,s],left[:,r,1,None]+matrices[None,:,1,s])
                for alpha,coeff in zip(alphas,result):
                    coeff[:,r,s]+=np.bincount(degree,weights=np.exp2(probability+alpha*entry).reshape(-1),minlength=129)
    with np.errstate(divide='ignore'):
        return [np.log2(coeff) for coeff in result]


def screen(rows,z,alpha,coeff,*,regional=False,width=128,outer_length=256):
    if outer_length not in (64,128,256):raise ValueError('Unsupported rate-half RS outer')
    region_length=2**21//outer_length;dimension=outer_length//2
    lc=log2(comb(region_length,rows))
    prefix=lc+alpha*(rows*(16*outer_length-dimension*log2(65535))-log2(65535)-838860*log2(z))
    def local(odds):
        return logsumexp((coeff+(np.arange(width+1)*odds)[:,None,None])*np.log(2),axis=0)/np.log(2)
    if regional:
        entries=np.zeros((2,2)); markers=np.zeros((2,2))
        for r in range(2):
            for s in range(2):
                def entry(odds):
                    value=power(tuple(tuple(row) for row in local(odds)),region_length//width)[r][s]
                    return value-rows*odds-lc
                chosen=minimize_scalar(entry,bounds=(-24,20),method='bounded')
                entries[r,s]=chosen.fun; markers[r,s]=chosen.x
        value=prefix+total(power(tuple(tuple(row) for row in entries),outer_length)[0])
        marker=markers.tolist()
    else:
        def objective(odds):
            return prefix+total(power(tuple(tuple(row) for row in local(odds)),2**21//width)[0])-outer_length*(lc+rows*odds)
        chosen=minimize_scalar(objective,bounds=(-24,20),method='bounded')
        value=float(chosen.fun); marker=float(chosen.x)
    return dict(active_rows=rows,mark=z,alpha=alpha,marker_log_odds=marker,regional=regional,
                log2_upper_screen=value,directed=False,whole_code_certificate=False)


def exact_region_screen(rows,z,alpha,coeff,*,width=128,outer_length=256):
    from hd_sumrank_regions import product
    proposal=screen(rows,z,alpha,coeff,width=width,outer_length=outer_length)
    region_length=2**21//outer_length;dimension=outer_length//2
    odds=proposal['marker_log_odds']
    p=1/(1+2**(-odds))
    degree=min(rows,width)
    shifted=coeff[:degree+1]+(np.arange(degree+1)*log2(p)+(width-np.arange(degree+1))*log2(1-p))[:,None,None]
    peak=np.max(shifted)
    base=(np.exp2(shifted-peak),peak)
    size=coeff.shape[1]
    result=(np.eye(size).reshape((1,size,size)),0.0)
    exponent=region_length//width
    while exponent:
        if exponent&1: result=product(result,base,rows)
        exponent>>=1
        if exponent: base=product(base,base,rows)
    value,scale=result
    if not np.all(value[rows]>0):
        raise ArithmeticError('FLOAT paired regional coefficient underflow')
    lc=log2(comb(region_length,rows))
    mass=lc+rows*log2(p)+(region_length-rows)*log2(1-p)
    region=np.log2(value[rows])+scale-mass
    bound=(lc+alpha*(rows*(16*outer_length-dimension*log2(65535))-log2(65535)-838860*log2(z))
           +total(power(tuple(tuple(row) for row in region),outer_length)[0]))
    return dict(active_rows=rows,mark=z,alpha=alpha,marker_log_odds=odds,
                log2_upper_screen=bound,exact_regional_coefficient=True,
                directed=False,whole_code_certificate=False)


if __name__=='__main__':
    import argparse,json
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows',type=int,nargs='+',default=(128,512,1024))
    parser.add_argument('--state',type=int,default=8)
    parser.add_argument('--mixer',choices=('feedback_refresh','double_refresh','uniform_pre'),default='feedback_refresh')
    parser.add_argument('--cutoff',type=int,default=20)
    parser.add_argument('--alphas',type=float,nargs='+',default=(.04,.05,.06,.07,.08,.09,.1,.125))
    parser.add_argument('--marks',type=float,nargs='+')
    parser.add_argument('--exact-regions',action='store_true')
    args=parser.parse_args()
    for r in args.rows:
        best=None
        zs=args.marks or (2**(-r/512),2**(-r/480),2**(-r/544))
        for z in zs:
            cs=coefficients(z,args.alphas,state=args.state,mixer=args.mixer,cutoff=args.cutoff)
            for alpha,coeff in zip(args.alphas,cs):
                for regional in (False,True):
                    result=screen(r,z,alpha,coeff,regional=regional)
                    result.update(state=args.state,mixer=args.mixer,cutoff=args.cutoff)
                    if best is None or result['log2_upper_screen']<best['log2_upper_screen']: best=result
                if args.exact_regions:
                    result=exact_region_screen(r,z,alpha,coeff)
                    result.update(state=args.state,mixer=args.mixer,cutoff=args.cutoff)
                    if best is None or result['log2_upper_screen']<best['log2_upper_screen']: best=result
        print(json.dumps(best),flush=True)
