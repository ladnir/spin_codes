"""Experimental inner: independent state refresh after adding feedback.

Keep y=x+E(q), but use q'=R_post(R_pre(q)+B(x)). The current inner
omits R_post. This changes the encoder/setup distribution. The additional
updates fix zero and have lazy/uniform action on every nonzero state.
The outgoing envelope is therefore obtained by right-composing the
existing epoch operator with the positive map below. No performance or
whole-code certificate is claimed by this binary64 screen.
"""
import argparse

from flint import arb, arb_mat

import attack_cache
import mass_density_screen as screen
from mixing_attack import retarget
from occupancy_memory import Z, F, M, C, U
from group_rank_one_verify import up


def refresh(spectrum, rounds):
    """Eleven-coordinate envelope for an independent output-state refresh.

F, M, and U entries bound masses of separate submeasures. C/L entries
are auxiliary bounds for the M submeasure, not additional source mass.
"""
    if (rounds is not None and (type(rounds) is not int or rounds<0)
            or len(spectrum)!=5 or any(type(n) is not int or n<=0 for n in spectrum.values())):
        raise ValueError('nonnegative update count and five nonempty state classes required')
    alpha=arb(0) if rounds is None else arb(2)**(-rounds)
    beta=1-alpha
    state_count=sum(spectrum.values())
    matrix=arb_mat(11,11)
    matrix[Z,Z]=1
    for source in range(1,11):matrix[source,source]=alpha
    for source in (F,M,*range(U,U+5)):
        for target,n in enumerate((spectrum[v] for v in sorted(spectrum)),U):
            matrix[source,target]=up(matrix[source,target]+beta*n/state_count)
    return matrix


def compose(base,spectrum,rounds):
    post=refresh(spectrum,rounds)
    return [screen.baseline.rounded(t*post) for t in base]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tilts',nargs='+',default=['.072','.088','.096','.104','.12'])
    parser.add_argument('--pre',type=int,nargs='+',default=[2,3])
    parser.add_argument('--post',nargs='+',default=['1','2','3','ideal'])
    parser.add_argument('--groups',type=int,nargs='+',default=[96,128])
    parser.add_argument('--supports',type=int,nargs='+',default=[200])
    parser.add_argument('--directory',default='tmp/four-bit-attack-cache')
    args=parser.parse_args()
    try:post_counts=[None if r=='ideal' else int(r) for r in args.post]
    except ValueError:parser.error('post count must be an integer or ideal')
    if (any(not 2<=r<=32 for r in args.pre)
            or any(r is not None and not 0<=r<=32 for r in post_counts)
            or any(not 1<=q<=2048 for q in args.groups)
            or any(not 38<=u<=256 for u in args.supports)):
        parser.error('invalid update count or selected event')
    print('NEW PRE/POST INNER; SELECTED BINARY64 SCREEN, NOT A CERTIFICATE',flush=True)
    data=attack_cache.build(args.tilts,directory=args.directory)
    caps=screen.baseline.authenticated_caps()
    cdf=screen.baseline.integer_cdf(screen.baseline.weighted_cdf_upper(caps,1<<128,full_weight=screen.Q(10,9)))
    shells=screen.baseline.weighted_union_shells(caps,full_weight=screen.Q(10,9))
    counts={u:min(screen.Q(cdf[u]),shells[u]) for u in args.supports}
    best={}
    for tilt,(base,_,details) in data.items():
        spectrum=details['spectrum']
        for before in args.pre:
            pre=retarget(base,spectrum,tilt,2,before)
            for after in post_counts:
                operators=compose(pre,spectrum,after)
                region=screen.float_placement(screen.as_array(operators),max(args.groups))
                for q in args.groups:
                    for u,count in counts.items():
                        value,p=screen.score(region,q,u,count,tilt,cutoff=193986)
                        print('PREPOST q/u/tilt/pre/post',q,u,tilt,before,after,'log2',value,'p',p,flush=True)
                        key=q,u,before,after
                        if key not in best or value<best[key][0]:best[key]=(value,tilt,p)
    for key,row in best.items():print('PREPOST BEST',key,row,flush=True)


if __name__=='__main__':main()
