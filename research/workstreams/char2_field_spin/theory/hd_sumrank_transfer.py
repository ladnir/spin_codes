"""Exact local mask moments for the sum-rank candidate; FLOAT global screens.

These routines do not certify the original FFT64 maps. Universal local
ranks here require the distinct-conjugacy LRS blocks in hd_sumrank_maps.
"""

from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial, log2, prod

from hd_randomized_outer_state_screen import NEG, add, total
from hd_randomized_outer_state_profile import power


def _advance(dp, local):
    state=len(dp)-1
    following = [0]*(state+1)
    for r, value in enumerate(dp):
        for s, coefficient in enumerate(local):
            following[min(state, r+s)] += value*coefficient
    return following


@lru_cache(maxsize=20)
def profiles(mark, field=65536, input_mark=F(1), state=4, maximum=64, with_live=False,tiles=8,with_gamma=False,feedthrough='cauchy'):
    """(occupancies,j,multiplicity,K,T,V), all three moments exact rational.

    K: zero previous state and zero next state.
    T: zero previous state, unrestricted next state.
    V: uniform NONZERO previous state, unrestricted next state.
    """
    z, input_mark = F(mark), F(input_mark)
    if not 0 < z <= 1 or input_mark < 1 or tiles not in (1,2,4,8) or not 1<=state<=4*tiles:
        raise ValueError('Invalid weight mark')
    if with_gamma and (not with_live or input_mark!=1):
        raise ValueError('The quiet live-state mask cap requires unmarked live profiles')
    if feedthrough not in ('cauchy','identity'):
        raise ValueError('Unknown feedthrough')
    quiet_local, flat_local, joint_local = [], [], []
    for a in range(9):
        quiet, flat, joint = [0]*5, [0]*5, []
        if feedthrough=='cauchy':
            terms=((zeros,min(a,zeros),max(0,zeros-a),comb(8,zeros)) for zeros in range(9))
        else:
            terms=((c+d,c,d,comb(a,c)*comb(8-a,d)) for c in range(a+1) for d in range(9-a))
        for zeros,overlap,outside,multiplicity in terms:
            residual=a-overlap
            coefficient=(multiplicity*z.numerator**(8-zeros)
                         *(z.denominator-z.numerator)**zeros*field**residual)
            quiet[min(4,residual)]+=coefficient
            flat[min(4,outside)]+=coefficient
            joint.append((min(4,residual),min(4,outside),coefficient))
        quiet_local.append(quiet)
        flat_local.append(flat)
        joint_local.append(joint)
    if input_mark != 1:
        def transform(local):
            return [[sum(comb(a, b)*input_mark.numerator**b
                         *(input_mark.denominator-input_mark.numerator)**(a-b)*local[b][r]
                         for b in range(a+1)) for r in range(5)] for a in range(9)]
        quiet_local, flat_local = transform(quiet_local), transform(flat_local)
    common = z.denominator**(8*tiles)
    result, counts = [], [0]*9

    def visit(profile, start, j, multiplicity, quiet, flat, free, joint):
        if len(profile) == tiles:
            multiplicity *= factorial(tiles)//prod(factorial(n) for n in counts)
            denominator = common*input_mark.denominator**j
            k = F(sum(value*field**(state-r) for r, value in enumerate(quiet)), denominator*field**(j+state))
            t = F(free, denominator*field**j)
            v = F(sum(value*(field**(state-r)-1) for r, value in enumerate(flat)),
                  denominator*field**j*(field**state-1))
            if not F(1, field**j) <= k <= t <= input_mark**j or not 0 < v <= input_mark**j:
                raise ArithmeticError('Invalid exact profile moments')
            record=(profile, j, multiplicity, k, t, v)
            if with_live:
                if input_mark!=1:
                    raise ValueError('The maskwise live cap requires unmarked inputs')
                record+=(F(sum(flat[:-1]),denominator*field**j),)
                if with_gamma:
                    record+=(F(sum(value*field**(state-r) for (r,e),value in joint.items()),
                               denominator*field**(j+state)),)
            result.append(record)
            return
        for a in range(start, 9):
            if j+a>maximum:
                break
            counts[a] += 1
            following={}
            if with_gamma:
                for (r,e),value in joint.items():
                    for dr,de,coefficient in joint_local[a]:
                        if e+de<state:
                            key=(min(state,r+dr),e+de)
                            following[key]=following.get(key,0)+value*coefficient
            visit(profile+(a,), a, j+a, multiplicity*comb(8, a),
                  _advance(quiet, quiet_local[a]), _advance(flat, flat_local[a]),
                  free*sum(quiet_local[a]),following)
            counts[a] -= 1

    visit((), 0, 0, 1, (1,)+(0,)*state, (1,)+(0,)*state, 1,{(0,0):1})
    return tuple(result)


@lru_cache(None)
def minimum_live_weight(profile,state=4):
    """Exact minimum over arbitrary potential values and a NONZERO state."""
    best=[0]+[-1000]*(state-1)
    for a in profile:
        following=[-1000]*state
        for spent,value in enumerate(best):
            for zeros in range(9):
                cost=min(4,max(0,zeros-a))
                if spent+cost<state:
                    following[spent+cost]=max(following[spent+cost],value+zeros)
        best=following
    return 8*len(profile)-max(best)


def operator(profile, j, k, t, v, mark, field=65536, state=4,live_cap=None):
    """One positive (zero, flat, residual) transfer, rational upper bounds."""
    z, q = F(mark), field
    d, eps, c = q**state-1, F(1, q), F(q-1, q)
    syndrome_rank=min(state,sum(min(4,a) for a in profile))
    zeros_allowed=8*((state-1)//4)+(state-1)%4
    untouched = z**max(0, 8*profile.count(0)-zeros_allowed)
    support_floor=z**minimum_live_weight(profile,state)
    live = min(t*untouched, d*v,support_floor)
    if live_cap is not None:
        live=min(live,live_cap)
    gamma = min((F(1, q)+c*z)**j, support_floor, d*v)
    gamma_zero = min((F(1, q)+c*z)**j-F(1, q**j),
                     support_floor*(1-F(1, q**syndrome_rank)), d*v)
    return_factor = eps+c*F(q**syndrome_rank-1, d)
    return_zero = (min(k*return_factor, support_floor*F(1, q**syndrome_rank)*return_factor,
                       F(q**(state-1), d)*live) if j else F(0))
    return ((k, F(0), t-F(1, q**j)),
            ((c*v+eps*gamma_zero)/d if j else F(0), c*v+eps*gamma, F(0)),
            (return_zero, c*live, eps*live))


@lru_cache(maxsize=24)
def operator_logs(mark, mixer='rank_one', input_mark=F(1), state=4, maximum=64,tiles=8):
    if mixer not in ('rank_one', 'uniform_pre', 'feedback_refresh'):
        raise ValueError('Unknown mixer')
    result = []
    for weighted, unweighted in zip(profiles(F(mark), input_mark=F(input_mark),state=state,maximum=maximum,tiles=tiles), profiles(F(mark),state=state,maximum=maximum,with_live=True,tiles=tiles,with_gamma=mixer=='feedback_refresh')):
        profile, j, multiplicity, k, t, v = weighted
        _, _, _, old_k, old_t, old_v,live_cap = unweighted[:7]
        quiet_live_cap = unweighted[7] if mixer=='feedback_refresh' else old_k
        factor = F(input_mark)**j
        if mixer == 'rank_one':
            matrix = [list(row) for row in operator(profile, j, old_k, old_t, old_v, mark,state=state,live_cap=live_cap)]
            matrix = [[value*factor for value in row] for row in matrix]
            matrix[0] = [k, F(0), t-k]
            # The uniform refresh acts on the exact marked flat moment.
            q, d = 65536, 65536**state-1
            eps, c = F(1, q), F(q-1, q)
            rank=min(state,sum(min(4,a) for a in profile))
            untouched = F(mark)**minimum_live_weight(profile,state)
            gamma = min(factor*(F(1,q)+c*F(mark))**j, factor*untouched, d*v)
            gamma0 = min(factor*((F(1,q)+c*F(mark))**j-F(1,q**j)),
                         factor*untouched*(1-F(1,q**rank)), d*v)
            matrix[1] = [(c*v+eps*gamma0)/d if j else F(0), c*v+eps*gamma, F(0)]
        elif mixer=='uniform_pre':
            q, d = 65536, 65536**state-1
            rank=min(state,sum(min(4,a) for a in profile))
            untouched = F(mark)**minimum_live_weight(profile,state)
            zero = (min(v, factor*old_k*F(q**rank-1, d),
                        factor*((F(1, q)+F(q-1, q)*F(mark))**j-F(1, q**j))/d,
                        factor*untouched*F(1, q**rank)*F(q**rank-1, d)) if j else F(0))
            matrix = ((k, t-k), (zero, v))
        else:
            q,d=65536,65536**state-1
            rank=min(state,sum(min(4,a) for a in profile))
            untouched=F(mark)**minimum_live_weight(profile,state)
            gamma=min(factor*old_k,factor*quiet_live_cap,factor*untouched/F(q**rank),d*v)
            # Feedback-active branches become flat before translation by M Q.
            # Feedback-zero branches use the old rank-one flat-density law.
            bulk=(F(q-1,q)*v if j==0 else v)+gamma/q
            matrix=((k,t-k),(v/d if j else F(0),bulk))
        logarithms = tuple(tuple(log2(value.numerator)-log2(value.denominator) if value else NEG
                                 for value in row) for row in matrix)
        result.append((j, log2(multiplicity)-log2(comb(8*tiles, j)), logarithms))
    return tuple(result)


@lru_cache(maxsize=120)
def fractional_packets(mark, alpha, *, quiet_only=False, mixer='rank_one', input_mark=F(1),state=4,tiles=8):
    z, a = F(mark), F(alpha)
    if not 0 < a <= 1:
        raise ValueError('Invalid exponent')
    if quiet_only:
        quiet, free = [NEG]*(8*tiles+1), [NEG]*(8*tiles+1)
        for profile, j, multiplicity, k, t, v in profiles(z, input_mark=F(input_mark),state=state,tiles=tiles):
            weight = log2(multiplicity)
            quiet[j] = add(quiet[j], weight+float(a)*(log2(k.numerator)-log2(k.denominator)))
            free[j] = add(free[j], weight+float(a)*(log2(t.numerator)-log2(t.denominator)))
        return (tuple(x-log2(comb(8*tiles, j)) for j, x in enumerate(quiet)),
                tuple(x-log2(comb(8*tiles, j)) for j, x in enumerate(free)))
    size = 3 if mixer == 'rank_one' else 2
    result = [[[NEG]*size for _ in range(size)] for j in range(8*tiles+1)]
    for j, weight, matrix in operator_logs(z, mixer, F(input_mark),state,tiles=tiles):
        for r in range(size):
            for s in range(size):
                value = matrix[r][s]
                result[j][r][s] = add(result[j][r][s], weight+float(a)*value)
    return tuple(tuple(tuple(row) for row in matrix) for matrix in result)


def screen(rows, mark, alpha, *, quiet_only=False, mixer='rank_one', input_mark=F(1),state=4,tiles=8):
    width=8*tiles
    steps=2**21//width
    if quiet_only:
        from hd_randomized_outer_quiet_exact import weighted_logsum_mean
        quiet,free=fractional_packets(mark,alpha,quiet_only=True,input_mark=input_mark,state=state,tiles=tiles)
        ordinary=tuple(log2(comb(width,j))+value for j,value in enumerate(quiet))
        terminal=tuple(log2(comb(width,j))+value for j,value in enumerate(free))
        lower,upper=-128.0,128.0
        for _ in range(80):
            mid=(lower+upper)/2
            _,mean=weighted_logsum_mean(ordinary,mid)
            _,last_mean=weighted_logsum_mean(terminal,mid)
            if (steps-1)*mean+last_mean<256*rows: lower=mid
            else: upper=mid
        odds=(lower+upper)/2
        first,_=weighted_logsum_mean(ordinary,odds)
        last,_=weighted_logsum_mean(terminal,odds)
        value=(float(alpha)*(rows*(4096-128*log2(65535))-log2(65535)-838860*log2(mark)-129*rows*log2(input_mark))
               +(steps-1)*first+last-256*rows*odds-255*log2(comb(8192,rows)))
        return dict(active_rows=rows,mark=str(mark),alpha=str(alpha),input_mark=str(input_mark),
                    state=state,packet_width=width,marker_log_odds=odds,log2_upper_screen=value,
                    directed=False,whole_code_certificate=False)
    from scipy.optimize import minimize_scalar
    packets = fractional_packets(F(mark), F(alpha), mixer=mixer, input_mark=F(input_mark),state=state,tiles=tiles)
    size = len(packets[0])
    choose = log2(comb(8192, rows))
    prefix = choose+float(alpha)*(rows*(4096-128*log2(65535))-log2(65535)-838860*log2(mark)-129*rows*log2(input_mark))

    def objective(odds):
        matrix = tuple(tuple(total(log2(comb(width, j))+j*odds+packets[j][r][s]
                                   for j in range(width+1)) for s in range(size)) for r in range(size))
        return prefix+total(power(matrix, steps)[0])-256*(choose+rows*odds)

    selected = minimize_scalar(objective, bounds=(-32, 32), method='bounded', options={'xatol':1e-7})
    return dict(active_rows=rows, mark=str(mark), alpha=str(alpha), marker_log_odds=float(selected.x),
                log2_upper_screen=float(selected.fun), quiet_only=False, mixer=mixer, input_mark=str(input_mark),state=state,packet_width=width,
                directed=False, whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quiet-only', action='store_true')
    parser.add_argument('--mixer', choices=('rank_one', 'uniform_pre','feedback_refresh'), default='rank_one')
    parser.add_argument('--input-marks', type=F, nargs='+', default=(F(1),))
    parser.add_argument('--state', type=int, default=4)
    parser.add_argument('--tiles',type=int,default=8)
    parser.add_argument('--rows', type=int, nargs='+', default=(16, 128, 512, 1024, 1536, 2048, 4096, 6144))
    parser.add_argument('--marks',type=F,nargs='+',default=(F(63,64),F(7,8),F(1,2),F(1,4),F(1,8),F(1,16),F(1,64),F(1,256),F(1,4096)))
    parser.add_argument('--alphas',type=F,nargs='+',default=(F(1,64),F(1,32),F(1,16),F(1,8),F(1,4),F(1,2),F(1)))
    args = parser.parse_args()
    best = {r:None for r in args.rows}
    for z in args.marks:
        for v in args.input_marks:
            for a in args.alphas:
                for r in args.rows:
                    result = screen(r, z, a, quiet_only=args.quiet_only, mixer=args.mixer, input_mark=v,state=args.state,tiles=args.tiles)
                    if best[r] is None or result['log2_upper_screen'] < best[r]['log2_upper_screen']:
                        best[r] = result
            print(json.dumps(dict(completed_mark=str(z), input_mark=str(v), best=list(best.values()))), flush=True)
