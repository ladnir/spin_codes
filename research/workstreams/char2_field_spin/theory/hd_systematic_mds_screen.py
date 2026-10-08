"""FLOAT screen: Y=X+A Q, Q'=M Q+T A^T X with A an RS[w,s] map.

For 2s<=w on an additive w-point GF256 subspace, A^T A=0 and every
coordinate restriction of A has its MDS rank. No tile profiles are needed.
This is a design screen, never a distance certificate.
"""
from math import comb,log2
import numpy as np
from scipy.optimize import minimize,minimize_scalar
from hd_sumrank_fast_screen import bound
from hd_sumrank_profile_pairs import exact_region_screen,screen as paired_marker_screen


def packet_logs(z,state,mixer=True,feedback=True,residual=False,width=64):
    q=65536.;dim=q**state-1
    if not feedback and not mixer:raise ValueError('The residual operator retains the rank-one mixer')
    lite=feedback=='lite'
    rank1=feedback in ('rank1','rank1_cone')
    size=2 if feedback and not residual and not lite and not rank1 else 3
    values=np.zeros((width+1,size,size))
    for j in range(width+1):
        k=sum(comb(j,c)*z**(j-c)*(1-z)**c*q**(-c-min(state,j-c)) for c in range(j+1))
        t=(z+(1-z)/q)**j
        v=t*sum(comb(width-j,d)*z**(width-j-d)*(1-z)**d*(q**(state-d)-1)/dim
                for d in range(min(state,width+1-j)))
        floor=max(0,width+1-j-state)
        if feedback and not lite and not rank1:
            gamma=min(k,q**(-min(state,j))*z**floor,dim*v)
            live=(v+gamma/q if j else v*(q-1)/q+gamma/q) if mixer else (v+gamma if j else gamma)
            if not residual:values[j]=((k,max(0,t-k)),(v/dim if j else 0,live))
            else:
                cap=min(t*z**floor,dim*v)
                values[j]=((k,max(0,t-k),0),(v/dim if j else 0,live,0),
                           (cap/dim if j else 0,cap if j else cap*(q-1)/q,gamma/q))
        else:
            # The established zero/flat/residual operator applies because
            # C A=0 makes I+lambda A C invertible for every scalar lambda.
            # Unlike the two-component proof, births are arbitrary residual.
            rank=min(state,j);eps=1/q;c=1-eps
            live=min(t*z**floor,dim*v)
            gamma=min(t,z**floor,dim*v)
            gamma0=min(max(0,t-q**(-j)),z**floor*(1-q**(-rank)),dim*v)
            # Under a uniform j-subset, a fixed nonzero state belongs to
            # span(C_U) with probability <=j/64 for j<state (MDS span lemma).
            # K, the syndrome rank, and the output floor are independent of
            # the particular subset, so the return bound can retain this.
            span_cap=j/width if j<state else 1
            return_factor=eps*span_cap+c*(q**rank-1)/dim
            zero=min(k*return_factor,q**(-rank)*z**floor*return_factor,q**(state-1)/dim*live) if j else 0
            values[j]=((k,0,max(0,t-q**(-j))),
                       ((c*v+eps*gamma0)/dim if j else 0,c*v+eps*gamma,0),
                       (zero,c*live,eps*live))
            if lite:
                # Candidate T=D2 (I+cJ) D1, c not 0 or 1, even state.
                # Dense H and independent nonzero diagonals give a good
                # full-support part <=b U_* and residual mass <=delta.
                # The old flat-column bound applies to every invertible T.
                if state%2:raise ValueError('The rank-one H candidate needs even state')
                b=dim/(q-1)**state;delta=min(1,state/(q-1))
                quiet_cap=min(k,q**(-rank)*z**floor,dim*v)
                birth=max(0,t-k)
                good=c+eps*b
                values[j,0]=(k,b*birth,delta*birth)
                values[j,2]=(min(zero,good*live/dim+eps*delta*live) if j else 0,
                             good*live,eps*min(live,quiet_cap+delta*live))
            if rank1:
                # Independent rank-one M and T. Each one-vector law is a
                # (1-1/q) uniform part plus 1/q scalar-line-preserving part.
                # Two simultaneous line branches cost 1/q^2, retaining
                # I+lambda A C invertibility in the flat-column bound.
                quiet_cap=min(k,q**(-rank)*z**floor,dim*v)
                birth=max(0,t-k);good=1-eps*eps
                factor=eps*eps*span_cap+good*(q**rank-1)/dim
                ret=min(k*factor,q**(-rank)*z**floor*factor,
                        (good/dim+eps*eps)*live) if j else 0
                values[j]=((k,c*birth,eps*birth),
                           ((good*v+eps*eps*gamma0)/dim if j else 0,
                            good*v+eps*(quiet_cap+eps*gamma) if j else c*v+eps*quiet_cap,0),
                           (ret,good*live if j else c*live,eps*min(live,quiet_cap+eps*live)))
    if feedback=='rank1_cone':
        # Every row creates residual <=B/(q-1). Preserve that invariant
        # before fractional powers; normalized live mass is (1+eta) B.
        eta=1/(q-1)
        if not np.all(values[:,:,2]<=eta*values[:,:,1]*(1+1e-12)):
            raise ArithmeticError('Invalid almost-flat cone invariant')
        values=np.stack((values[:,0,0],(1+eta)*values[:,0,1],
                         (values[:,1,0]+eta*values[:,2,0])/(1+eta),
                         values[:,1,1]+eta*values[:,2,1]),axis=1).reshape((width+1,2,2))
    with np.errstate(divide='ignore'):return np.log2(values)


def coefficients(logs,alpha,right=None):
    width=len(logs)-1
    if right is None:right=logs
    choose=np.array([log2(comb(width,j)) for j in range(width+1)])
    degree=(np.arange(width+1)[:,None]+np.arange(width+1)[None,:]).reshape(-1)
    weight=choose[:,None]+choose[None,:]
    size=logs.shape[1];coeff=np.zeros((2*width+1,size,size))
    for r in range(size):
        for s in range(size):
            entry=np.logaddexp2(logs[:,r,0,None]+right[None,:,0,s],logs[:,r,1,None]+right[None,:,1,s])
            for h in range(2,size):entry=np.logaddexp2(entry,logs[:,r,h,None]+right[None,:,h,s])
            coeff[:,r,s]=np.bincount(degree,weights=np.exp2(weight+alpha*entry).reshape(-1),minlength=2*width+1)
    with np.errstate(divide='ignore'):return np.log2(coeff)


def optimize(rows,state=4,iterations=80,regional=False,mixer=True,alternate=False,feedback=True,refresh_alternate=False,width=64,packet_source=None,pair_source=None,outer_length=256):
    if width not in (8,16,32,64,128,256) or not 1<=state<=min(32,width//2):
        raise ValueError('Unsupported self-orthogonal RS screening dimensions')
    if (alternate or refresh_alternate) and not regional:
        raise ValueError('Alternating ensembles require paired composition')
    if pair_source and (alternate or refresh_alternate):
        raise ValueError('Custom pair source requires ordinary paired composition')
    if refresh_alternate and (alternate or not mixer or not feedback):
        raise ValueError('Select just one alternating ensemble')
    best=None;cached={}
    def objective(theta):
        nonlocal best
        ell,alpha=theta
        if not .0001<=ell<=16 or not .005<=alpha<=1:return 1e12
        z=2**(-ell)
        if z not in cached:
            cached.clear();cached[z]=(packet_source(z,state,width) if packet_source else packet_logs(z,state,False if alternate else mixer,False if refresh_alternate else feedback,width=width))
            if alternate:cached['right']=packet_logs(z,state,True,width=width)
            if refresh_alternate:cached['right']=packet_logs(z,state,True,True,True,width)
        if regional or pair_source or outer_length!=256:
            evaluator=exact_region_screen if regional else paired_marker_screen
            try:result=evaluator(rows,z,alpha,pair_source(z,alpha,state,width) if pair_source else coefficients(cached[z],alpha,cached.get('right')),width=2*width,outer_length=outer_length)
            except ArithmeticError:return 1e12
            value=result['log2_upper_screen'];odds=result['marker_log_odds']
        else:
            result=minimize_scalar(lambda p:bound(rows,z,alpha,p,alpha*cached[z],tiles=width//8),bounds=(-24,20),method='bounded')
            value=float(result.fun);odds=float(result.x)
        if not np.isfinite(value):return 1e12
        if best is None or value<best['log2_upper_screen']:
            best=dict(active_rows=rows,state=state,packet_width=width,outer_length=outer_length,mark=z,alpha=alpha,
                      marker_log_odds=odds,log2_upper_screen=value,regional=regional,
                      construction='systematic-global-MDS',state_mixer=mixer,mix_alternate=alternate,
                      feedback_refresh=feedback,refresh_alternate=refresh_alternate,directed=False,whole_code_certificate=False)
        return value
    seed=max(.002,rows*outer_length/131072)
    minimize(objective,(seed,.0625),method='Nelder-Mead',options={'maxiter':iterations})
    if best is None:
        minimize(objective,(min(2.,seed),.0625),method='Nelder-Mead',options={'maxiter':iterations})
    if best is None:raise ArithmeticError('No numerically usable FLOAT proposal')
    return best


if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rows',type=int,nargs='+',default=(16,128,512,1024))
    p.add_argument('--width',type=int,default=64)
    p.add_argument('--outer-length',type=int,choices=(64,128,256),default=256)
    p.add_argument('--states',type=int,nargs='+',default=(4,8,16))
    p.add_argument('--iterations',type=int,default=80)
    p.add_argument('--regional',action='store_true')
    p.add_argument('--no-mixer',action='store_true')
    p.add_argument('--alternate',action='store_true')
    p.add_argument('--no-feedback',action='store_true')
    p.add_argument('--refresh-alternate',action='store_true')
    p.add_argument('--lite-feedback',action='store_true')
    p.add_argument('--rank-one-feedback',action='store_true')
    p.add_argument('--rank-one-cone',action='store_true')
    args=p.parse_args()
    for s in args.states:
        for r in args.rows:print(json.dumps(optimize(r,s,args.iterations,args.regional,not args.no_mixer,args.alternate,
                                                   'rank1_cone' if args.rank_one_cone else 'rank1' if args.rank_one_feedback else 'lite' if args.lite_feedback else not args.no_feedback,args.refresh_alternate,args.width,outer_length=args.outer_length)),flush=True)
