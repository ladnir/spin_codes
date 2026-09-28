"""Counterfactual sensitivity to binomial-shaped dual shells, NOT a proof.

This standalone diagnostic substitutes unproved shell caps only to decide
whether further work on dual shells can plausibly repair the current bound.
Its results must not be used by a certificate or the production encoder.
"""
from math import comb,log2
from dual_sandwich import refined_spectrum
from dual_moments import refine_bch
from shortened_bound import dimension_caps
from shortening_polynomial import improve_dimensions
from dual_shortening import improve_dimensions as dual_dimensions
from bch_joint_support import authenticated_caps,support_caps
from basis_lattice import improve_caps
from shortening_moments import improve
from occupancy_allones import weighted_cdf


def main():
    print('COUNTERFACTUAL DIAGNOSTIC. Binomial shell caps are NOT proved for this BCH dual.',flush=True)
    spectrum=authenticated_caps()
    dims=dual_dimensions(improve_dimensions(dimension_caps()))
    primal=improve_caps(support_caps(spectrum,g=4,dimensions=dims),spectrum)
    primal,_=improve(primal,dims)
    checked=refined_spectrum()
    before=refine_bch(primal,dims,dual_spectrum=checked)
    guessed=checked[:]
    for w in range(30,227,2):
        # Round up the even-binomial reference, intersect with proven caps.
        guessed[w]=min(checked[w],-(-comb(256,w)//(1<<127)))
    assert guessed==guessed[::-1]
    after=refine_bch(primal,dims,dual_spectrum=guessed)
    old=weighted_cdf(spectrum,before,dims,'.75',prefix_flags=True)
    new=weighted_cdf(spectrum,after,dims,'.75',prefix_flags=True)
    # Previously measured selected points: identical u in each of 64 groups,
    # lambda=.032, rho=.75, all refinements, full_feedback=6. For a singleton
    # support interval, changing the CDF changes the log bound by exactly
    # 64*log2(new[u]/old[u]), independent of the auxiliary Bernoulli parameter.
    prior={128:2839.5478020247538,144:2408.779213684673,160:2259.1081345521698}
    for u in (96,112,128,144,160,176,192):
        shift=64*(log2(new[u])-log2(old[u]))
        print('COUNTERFACTUAL support',u,'weighted log2 count old/new',log2(old[u]),log2(new[u]),
              'q64 shift',shift,'projected selected point',prior[u]+shift if u in prior else None,flush=True)
    print('No certificate: neither the guessed shells nor projected scores are BCH guarantees.',flush=True)


if __name__=='__main__':main()
