"""Refine the state envelope using exact full feedback distributions.

The CLI is a binary64 selected-point diagnostic. Coefficient construction
uses outward Arb arithmetic; no full-code certificate is claimed here.
"""
import argparse
from math import log
import numpy as np
from flint import arb,ctx
from scipy.optimize import minimize_scalar

from full_feedback_census import census,check_prior_zeros
from group_moment import maps
from group_rank_one_verify import up
from occupancy_memory import Z,F,M,C,U
from occupancy_sensitivity import float_placement
from occupancy_screen import matrix_for_probabilities
from two_group_screen import log_power_moment,log_binomial_mass
from mature_tail import TAIL_TERMINAL,L48,L56


def coefficients(data,spectrum,tilt,penalty,rounds=2,input_penalty='1',odd_penalty=1):
    """Return shape-maximized lazy coefficients and separate refresh slots.

    rounds=None is the uniform-refresh limit, used only by diagnostics.
    """
    powers=[(-arb(tilt)*w).exp() for w in range(129)]
    lazy=arb(0) if rounds is None else arb(2)**(-rounds)
    candidates={}
    for weights,(zero,peak,denominator,classes) in data.items():
        j=len(weights)
        if j<2:continue
        weight=sum(weights)
        scale=arb(penalty)**weights.count(4)*arb(input_penalty)**weight*arb(odd_penalty)**sum(w%2 for w in weights)
        all_peak=arb(max(zero,peak))/denominator
        nz_peak=arb(peak)/denominator
        zero_probability=arb(zero)/denominator
        row={(Z,Z):powers[weight]*zero_probability,
             (Z,M):powers[weight]*(1-zero_probability),
             (Z,C):powers[weight]*nz_peak,
             (F,C):lazy*powers[max(0,48-weight)]*all_peak,
             (F,Z):lazy*powers[max(0,48-weight)]*nz_peak}
        # If q+Bx=0, output is E(Bx)+x. Its weight is at least |v-weight|
        # when E(Bx) has weight v. Source zero is excluded from C and U.
        row[C,Z]=lazy*sum((arb(count)*powers[abs(v-weight)] for v,count in classes.items()),arb(0))/denominator
        for i,v in enumerate(sorted(spectrum)):
            row[U+i,Z]=lazy*powers[abs(v-weight)]*classes[v]/(denominator*spectrum[v])
            row[U+i,C]=lazy*powers[abs(v-weight)]*min(all_peak,arb(1)/spectrum[v])
        for target,cut in ((L48,48),(L56,56)):
            row[Z,target]=powers[weight]*sum(count for v,count in classes.items() if v<=cut)/denominator
        for (source,target),value in row.items():
            key=j,source,target
            value=up(value*scale)
            candidates[key]=max(candidates.get(key,arb(0)),value)
    return candidates


def refine(base,data,spectrum,tilt,penalty,rounds=2,input_penalty='1',odd_penalty=1):
    """Preserve existing coordinates and take smaller valid scalar bounds.

    base must already use the stated number of updates. Fresh/uniform
    zero bounds below concern the lazy contribution; add a refresh upper
    from the unchanged uniform-class columns before taking a minimum.
    """
    result=[t*1 for t in base]
    for (j,source,target),lazy in coefficients(data,spectrum,tilt,penalty,rounds,input_penalty,odd_penalty).items():
        if j>=len(result) or target>=result[j].ncols():continue
        value=lazy
        if target==Z and (source==F or U<=source<U+5):
            refresh=min(up(base[j][source,U+i]/spectrum[v]) for i,v in enumerate(sorted(spectrum)))
            value=up(value+refresh)
        result[j][source,target]=min(result[j][source,target],value)
    return result


def check_pairs(data):
    from occupancy_model import local_data
    single,pairs,_=local_data(2)
    for weights,(zero,peak,den,classes) in data.items():
        if len(weights)!=2:continue
        a,b=weights
        choices,z,maximum,histograms=pairs[1][(a,),(b,)]
        assert zero*choices==z*den and peak*choices==maximum*den
        assert all(classes[v]*choices==sum(histograms[v].values())*den for v in classes)
    print('Full distribution peaks and expansion classes match independent pair census',flush=True)


def main():
    from flint import arb_mat
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot',default='tmp/r2-q64-sensitivity.npz')
    parser.add_argument('--maximum',type=int,choices=range(2,11),default=6)
    args=parser.parse_args();ctx.prec=192
    with np.load(args.snapshot,allow_pickle=False) as saved:
        base=saved['base'];outer=float(saved['outer'])
        groups,support,tilt,penalty=map(str,saved['parameters'])
    groups=int(groups);support=int(support)
    assert base.shape==(min(groups,32)+1,11,11) and groups>=4
    data=census(args.maximum)
    if args.maximum>=4:check_prior_zeros(data)
    check_pairs(data)
    spectrum=maps()[2]
    # The cache is binary64 and not a certificate input. Reconstructing Arb
    # matrices here cannot recover its rounding history or certify the result.
    matrices=[arb_mat([[arb(float(x)) for x in row] for row in t]) for t in base]
    changed=refine(matrices,data,spectrum,tilt,penalty)
    revised=np.array([[[float(t[i,j]) for j in range(11)] for i in range(11)] for t in changed])
    assert (revised<=base+1e-15).all()
    for limit in [0]+list(range(2,args.maximum+1)):
        operators=base.copy();operators[:limit+1]=revised[:limit+1]
        region=float_placement(operators,groups)
        def objective(z):
            p=1/(1+np.exp(-z))
            return (log_power_moment(matrix_for_probabilities(region,[p]*groups),256,TAIL_TERMINAL)
                    +float(tilt)*209715-groups*log_binomial_mass(256,support,p))
        fit=minimize_scalar(objective,bounds=(-8.,16.),method='bounded')
        print('DIAGNOSTIC NOT CERTIFICATE: full feedback through',limit,
              'log2 score',(fit.fun+outer)/log(2),'p',1/(1+np.exp(-fit.x)),flush=True)


if __name__=='__main__':main()
