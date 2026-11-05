"""FLOAT pre-alpha grouping for rank-one refresh after feedback.

The 17-state packet cap has arbitrary masses on masks0..15 and a separate
uniform-live density component16. Pair placements preserve the matching
{(0,1),(2,3)}, whose order-eight symmetry gives six mask orbits below.
This module changes neither the construction nor existing proof defaults.

Validity interface. Condition on the ENTIRE packet-type sequence (s,f).
Given those types, local placements and fair swaps are independent across
packets, including independence of future placements from incoming history.
Average the unpowered labelled caps conditional on each type. For source
orbit I and target orbit J, set H[I,J]=max_i sum_j E[M[i,j]|type], with
i in I and j in J. This bounds every incoming mixture in I; equality of
source rows is unnecessary. The singleton flat orbit retains its density
meaning, and the zero orbit contains only the actual zero state.

For the complete marked moment Z, Jensen and positive path expansion give
E[Z**alpha|types] <= (E[Z|types])**alpha
 <= (e0 product(H) 1)**alpha <= e0 product(H entrywise**alpha) 1.
Therefore the packet polynomial coefficient is N[s,f]*H[s,f]**alpha,
where N[s,f]=multinomial(32;32-s-f,s,f). N is NOT fractionalized.
The outer comparison uses uniform-F active coordinates, not the XOR
forced-nonzero comparison. All arithmetic here is FLOAT, not directed.

The single frozen (r,a,b)=(1024,400,400) gate below gives
+101453.78397569829 bits, versus +99844.04545436724 for the original
17-state cap at the same z, alpha, and markers. The seven-state witness
therefore fails and is worse by1609.73852133 bits. No fitting was done.
This is a proof-screen failure, not a counterexample to the construction.
Reproduce with: python -B hd_rankone_orbit_screen.py
"""
from math import comb, log2

import numpy as np
from scipy.special import logsumexp

from hd_grouped_packet_screen import LIMIT, log_multinomial
from hd_rankone_mask_screen import pair_coefficients, log_power
from hd_shifted_outer_density import TOEPLITZ_OUTER_MEASURE, log_density_float
from hd_shifted_pairs_count import one_pairing_count


ORBITS=((0,), (1,2,4,8), (3,12), (5,6,9,10),
        (7,11,13,14), (15,), (16,))
ORBIT_NAMES=('zero','one','two-within-paired-half','two-across-halves',
             'three','four','flat')
FROZEN_MEDIUM=dict(weight=1024,full=400,z=.4828132651072611,
    alpha=.10118724647150915,
    marker_log2=(-3.710160740811969,-2.660319071541368))


def lump_log_matrix(log_matrix):
    """Sum target components, then maximize source rows, without alpha."""
    values=np.asarray(log_matrix,dtype=float)
    if (values.shape!=(17,17) or np.any(np.isnan(values))
            or np.any(np.isposinf(values))):
        raise ValueError('Finite or negative-infinite 17-state log cap required')
    result=np.empty((7,7))
    for j,targets in enumerate(ORBITS):
        sums=logsumexp(values[:,targets]*np.log(2),axis=1)/np.log(2)
        for i,sources in enumerate(ORBITS):
            result[i,j]=np.max(sums[list(sources)])
    return result


def coefficients_from_unpowered(unpowered,alpha):
    """Input logs include actual placement multiplicity; output logs do too."""
    if not np.isfinite(alpha) or not 0<alpha<=1:
        raise ValueError('Fractional exponent must be in (0,1]')
    result={}
    for (singleton,full),matrix in unpowered.items():
        if (type(singleton) is not int or type(full) is not int
                or min(singleton,full)<0 or singleton+full>32):
            raise ValueError('Invalid packet type')
        log_count=log2(comb(32,singleton)*comb(32-singleton,full))
        mean=lump_log_matrix(matrix)-log_count
        result[singleton,full]=log_count+alpha*mean
    return result


def grouped_coefficients(z,alpha,maximum=64):
    """One bounded labelled build; no marker or witness optimization."""
    if (not np.isfinite(z) or not 0<z<=1 or not np.isfinite(alpha)
            or not 0<alpha<=1 or type(maximum) is not int or not 0<=maximum<=64):
        raise ValueError('Invalid mark, exponent, or packet cutoff')
    return coefficients_from_unpowered(pair_coefficients(z,1.,maximum),alpha)


def fixed_symmetric_gate(weight,full,z,alpha,marker_log2):
    """Single positive-marker bound for a=b=full; never a whole-code proof."""
    if type(weight) is not int or type(full) is not int:
        raise ValueError('Integer source profile required')
    histogram=(4096-weight+full,weight-2*full,full)
    if min(histogram)<=0:
        raise ValueError('Interior symmetric histogram required')
    odds=np.asarray(marker_log2,dtype=float)
    if odds.shape!=(2,) or not np.all(np.isfinite(odds)):
        raise ValueError('Two finite fixed coefficient markers required')
    coefficients=grouped_coefficients(z,alpha)
    keys=np.array(tuple(coefficients))
    logs=np.array(tuple(coefficients.values()))
    matrix=logsumexp((logs+np.dot(keys,odds)[:,None,None])*np.log(2),axis=0)/np.log(2)
    end=log_power(matrix,32768)[0]
    prefix=(log2(one_pairing_count(4096,weight,full))
        +alpha*(weight*log_density_float(TOEPLITZ_OUTER_MEASURE)-log2(255)-LIMIT*log2(z))
        -256*log_multinomial(histogram)-256*np.dot(histogram[1:],odds))
    value=prefix+float(logsumexp(end*np.log(2))/np.log(2))
    return dict(active_groups=weight,a=full,b=full,z=z,alpha=alpha,
        marker_log2=odds.tolist(),log2_bound_FLOAT=value,evaluations=1,
        count_method='single-matching-upper',directed=False,
        whole_code_certificate=False,components=7,orbit_names=list(ORBIT_NAMES),
        conditional_average_before_alpha=True,source_orbit_rule='entrywise-row-maximum',
        refresh='rank-one-after-feedback',outer_measure=TOEPLITZ_OUTER_MEASURE)


if __name__=='__main__':
    import json
    print(json.dumps(fixed_symmetric_gate(**FROZEN_MEDIUM),indent=2,allow_nan=False))
