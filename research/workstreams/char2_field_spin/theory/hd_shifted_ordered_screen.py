"""FLOAT diagnosis of potential loss using ordered regional matrix powers.

For fixed scalar support weight r, evaluate every feasible histogram pair
(a,b) with its actual positive two-state matrix product (R_a R_b)**128.
The regional matrices are coefficient upper bounds supplied by the caller.
This removes common-potential slack but not regional or support-count slack.
All arithmetic in this module is FLOAT; no result is a certificate.
"""
import numpy as np

from hd_shifted_run_potential import log_run_counts,log_joint_counts
from hd_shifted_swap_potential import full_counts


def product(left,right):
    a,b,c,d = left
    e,f,g,h = right
    return (np.logaddexp2(a+e,b+g),np.logaddexp2(a+f,b+h),
            np.logaddexp2(c+e,d+g),np.logaddexp2(c+f,d+h))


def ordered_log_screen(pairs,weight,log_regions,*,repetitions=256,joint_counts=False):
    """Broadcast over all ordered matching counts, keeping state chronology."""
    if not isinstance(repetitions,int) or repetitions<2 or repetitions%2:
        raise ValueError('An even positive number of regional repetitions is required')
    counts = tuple(full_counts(pairs,weight))
    entries = np.asarray([np.asarray(log_regions[f]).reshape(4) for f in counts])
    if np.any(np.isnan(entries)) or np.any(np.isposinf(entries)):
        raise ValueError('Finite or negative-infinite matrix logs required')
    left = tuple(entries[:,i,None] for i in range(4))
    right = tuple(entries[None,:,i] for i in range(4))
    base = product(left,right)
    exponent,result = repetitions//2,None
    while exponent:
        if exponent&1:
            result = base if result is None else product(result,base)
        exponent >>= 1
        if exponent:
            base = product(base,base)
    multiplicities = (log_joint_counts if joint_counts else log_run_counts)(pairs,weight)
    terms = multiplicities+np.logaddexp2(result[0],result[1])
    indices = np.unravel_index(np.argmax(terms),terms.shape)
    return dict(log2_support_sum_screen=float(np.logaddexp2.reduce(terms.ravel())),
                largest_profile=(weight,counts[indices[0]],counts[indices[1]]),
                largest_profile_term_screen=float(terms[indices]),
                matrix_pairs_evaluated=terms.size,regional_histograms=len(counts),
                support_count_method='joint-upper' if joint_counts else 'cyclic-runs',
                common_potential=False,directed=False,whole_code_certificate=False)
