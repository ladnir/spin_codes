"""Composable quiet/singleton rank-cap refinement for full-Cauchy packets.

States 1..s are upper bounds on the surviving state-form rank, together
with actual witness-state liveness. They are not disjoint rank events.
Domination is on nonnegative continuations nondecreasing on these caps;
state0 is unconstrained relative to them. No complete distance certificate.
"""

from fractions import Fraction as F
import argparse
import json
import math

from state_rank_cap_transfer import rank_cap_weights


def singleton_return(cap, y, tau, q=65536, state=4, width=4):
    """Return-first affine-torus bound, uniform over incoming ranks<=cap."""
    if not (1 <= cap <= state and 0 < y <= 1 and 1 <= tau <= q-1):
        raise ValueError("invalid cap or marks")
    denominator = F(q-1) if isinstance(tau, (F, int)) else q-1
    return y**width * denominator**(cap-state) / tau**cap


def local_matrix(weight, y, tau, q=65536, state=4, width=4):
    """A positive cap operator; multi-input packets retain the old cover.

    The new proof applies to any square full-Cauchy width/state geometry.
    The multi-input fallback is scoped to the existing width/state4 helper.
    Fresh independent nonzero field diagonals are required at every packet.
    """
    if not (width == state and 1 <= state <= 4 and 0 <= weight <= width
            and 0 < y <= 1 and 1 <= tau <= q-1):
        raise ValueError("invalid square packet or marks")
    if isinstance(y, int):
        y = F(y)
    if isinstance(tau, int):
        tau = F(tau)
    zero = y*0 + tau*0
    result = [[zero for _ in range(state+1)] for _ in range(state+1)]
    if weight >= 2:
        from high_distance_fields import local_matrices
        bits = q.bit_length()-1
        if (state, width) != (4,4) or 1 << bits != q:
            raise ValueError("multi-input fallback requires the checked 4/4 geometry")
        old = local_matrices(y,tau,width,width,state,bits)[weight]
        result[0][0], result[0][state] = old[0]
        for cap in range(1,state+1):
            result[cap][0], result[cap][state] = old[1]
        return tuple(map(tuple,result))

    if weight == 0:
        result[0][0] = y**0
    else:
        # A full-Cauchy column and its systematic state-injection column
        # are both nonzero in every coordinate. A singleton birth is rank1.
        result[0][1] = y**width

    for cap in range(1,state+1):
        if weight:
            result[cap][0] = singleton_return(cap,y,tau,q,state,width)
        for selected in range(width+1):
            factor = math.comb(width,selected)*y**(width-selected)*(1-y)**selected
            if not factor:
                continue
            if weight and selected == 0:
                result[cap][min(state,cap+1)] += factor
                continue
            # Quiet: all selected rows constrain old state.
            # Singleton: eliminate its input column; expose the remaining
            # superregular Schur rows before the one fresh input scalar.
            roots = selected-weight
            if roots >= state:
                continue
            for h,mass in enumerate(rank_cap_weights(cap,roots,tau,q,state),1):
                result[cap][h] += factor*mass/tau**weight
    return tuple(map(tuple,result))


def cone_generators(state=4):
    """Extreme rays of v>=0, v1<=...<=vs; v0 has no ordering condition."""
    yield (1,)+(0,)*state
    for first in range(1,state+1):
        yield (0,)+(0,)*(first-1)+(1,)*(state-first+1)


def log_polynomial_product(left,right,limit):
    """Small research-only positive matrix polynomial multiplication."""
    import numpy as np
    from scipy.special import logsumexp
    count = min(limit+1,len(left)+len(right)-1)
    states = left.shape[1]
    out = np.full((count,states,states),-math.inf)
    for degree,matrix in enumerate(left):
        length = min(len(right),count-degree)
        if length <= 0:
            break
        for row in range(states):
            contributions = logsumexp(matrix[row][None,:,None]+right[:length],axis=1)
            out[degree:degree+length,row,:] = np.logaddexp(
                out[degree:degree+length,row,:],contributions)
    return out


def log_terminal(matrix,count):
    import numpy as np
    from scipy.special import logsumexp
    value = np.full(len(matrix),-math.inf)
    value[0] = 0.
    for _ in range(count):
        value = logsumexp(value[:,None]+matrix,axis=0)
    return float(logsumexp(value))


def fixed_bound(point,value,sparse=True):
    """One floating K20 [256,128], 4/4 keep4, GF65536, 40% screen."""
    import numpy as np
    from cauchy_firstmoment_screen import log_choose,log_terminal_power
    from gf8_rank_history_screen import modal_log_penalty,source_row_log
    from high_distance_expansion_screen import profile
    from sparse_active_interval_screen import regional_from_coefficients
    g = profile(256,128,(4,),4,4,F(2,5),20)
    states = 5
    identity = np.full((1,states,states),-math.inf)
    np.fill_diagonal(identity[0],0.)
    if sparse:
        active = int(value)
        logc,logtau = map(float,point)
        logy = -math.exp(logc)*active/g['output_n']
        tau = math.exp(logtau)
        matrices = np.asarray([local_matrix(w,math.exp(logy),tau)
                               for w in range(min(4,active)+1)],dtype=float)
        with np.errstate(divide='ignore'):
            factor = np.log(matrices)
        for w in range(len(factor)):
            factor[w] += math.log(math.comb(4,w))
        result = identity
        steps = g['rows']//4
        while steps:
            if steps & 1:
                result = log_polynomial_product(result,factor,active)
            steps >>= 1
            if steps:
                factor = log_polynomial_product(factor,factor,active)
        regional = regional_from_coefficients(result,active,g['rows'],logtau)
        return (log_choose(g['rows'],active)-(128*active+1)*logtau
                -g['output_threshold']*logy+log_terminal(regional,256))/math.log(2)
    logx,logy,logtau,logtilt = map(float,point)
    matrices = np.asarray([local_matrix(w,math.exp(logy),math.exp(logtau))
                           for w in range(5)],dtype=float)
    packet = np.einsum('w,wij->ij',[math.comb(4,w)*math.exp(w*logx)
                                   for w in range(5)],matrices)
    return (256*modal_log_penalty(g['rows'])
            +g['rows']*source_row_log(value/g['input_n'],logtau,logtilt,256,128,True)
            +log_terminal_power(packet,g['input_n']//4)
            -value*(logx+logtilt)-g['output_threshold']*logy-logtau)/math.log(2)


def replay():
    from hd_fast_screens import CASES
    from high_distance_expansion_screen import bound
    for case in CASES:
        name,b,k,width,state,pattern,sparse,target,marks,expected = case
        if not name.startswith('2K-state4-'):
            continue
        old = -bound(marks,target,sparse,bits=16,length=b,dimension=k,
                     pattern=pattern,width=width,state=state,distance=F(2,5),k_log=20)
        new = -fixed_bound(marks,target,sparse)
        yield dict(name=name,log_parameters=marks,target=target,
                   old_margin_bits=old,new_margin_bits=new,improvement_bits=new-old,
                   rigorous=False,whole_code=False,optimized=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay',action='store_true')
    args = parser.parse_args()
    if args.replay:
        for receipt in replay():
            print(json.dumps(receipt),flush=True)
