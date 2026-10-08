"""Directed middle replay for the higher-rate fixed-puncture family.

Default profile: GF256 outer[128,80], full-Cauchy width8/state8, keep6.
The frozen outer[128,84]/keep7 certificate is neither edited nor relabeled.
Rational witnesses are proposals; directed positive arithmetic alone
accepts every integer cell in the stated input-weight interval.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb

from gf8_outer_rate_middle import conditioning
from gf8_punctured_screen import local_matrices
from gf8_rate_climb_profile import metadata
from return_cover_certificate import Witness,interval,terminal_power,row_upper
from state8_density_certificate import integer_span,refined_cells,retarget_witness
from subfield_certificate import HIGH,LOW,ZERO,ONE,integer_margin


LOW.Emax=HIGH.Emax=1_000_000_000
LOW.Emin=HIGH.Emin=-1_000_000_000
WITNESSES_80_6=tuple((F(density),Witness(*(F(value) for value in values))) for density,values in (
    ('3/100', ('886/2877', '330/511', '3331/500', '21965/4557')),
    ('1/20', ('4702/7901', '2737/5906', '73449/9998', '8567/3079')),
    ('1/10', ('4421/2550', '859/4682', '47465/4542', '5072/4263')),
    ('7/50', ('20762/6449', '401/4681', '131890/8799', '1669/2350')),
    ('1/5', ('49937/7357', '39/1330', '251247/9178', '1601/4501')),
    ('3/10', ('147606/7649', '31/5359', '794386/9809', '757/6294')),
    ('2/5', ('256751/5068', '9/7027', '255', '11/268')),
    ('1/2', ('113917/2143', '8/8311', '255', '48/1747')),
    ('3/5', ('130987/2390', '4/5191', '255', '138/7501')),
    ('7/10', ('107518/1921', '3/4667', '255', '10/841')),
    ('4/5', ('83949/1477', '2/3627', '255', '53/7576'))))
WITNESSES_67_5=tuple((F(density),Witness(*(F(value) for value in values))) for density,values in (
    ('3/100', ('2946/9385', '5367/9059', '73004/9593', '427/64')),
    ('1/20', ('4988/8275', '1989/4970', '34341/4000', '21237/5537')),
    ('1/10', ('7258/4323', '557/3969', '108039/8335', '9370/5759')),
    ('7/50', ('29114/9923', '386/6145', '38779/2025', '6578/6753')),
    ('1/5', ('20239/3637', '109/5070', '350529/9943', '4789/9562')),
    ('3/10', ('38698/2999', '5/1079', '877378/8983', '1544/8253')),
    ('2/5', ('202374/7823', '8/6401', '255', '585/7654')),
    ('1/2', ('33029/1227', '9/9560', '255', '65/1271')),
    ('3/5', ('230108/8335', '6/7949', '255', '248/7247')),
    ('7/10', ('205828/7327', '6/9517', '255', '58/2625')),
    ('4/5', ('283092/9949', '1/1849', '255', '31/2391')),
    ('9/10', ('246943/8594', '2/4223', '255', '29/4967'))))


def powers(base,length):
    result=[ONE]
    for _ in range(length):
        result.append(result[-1]*base)
    return result


def source_upper(lower,upper,witness,length=128,parity=48):
    """Same positive termwise envelope, with linearly tabulated powers."""
    if not (isinstance(lower,F) and isinstance(upper,F) and 0<lower<=upper<1):
        raise ValueError('exact interior density interval required')
    beta=witness.tau*witness.tilt
    if beta<2:
        return row_upper(lower,upper,witness,length,parity)
    bp=powers(interval(beta),length)
    endpoints={a:(powers(interval(a),length),powers(interval(1-a),length)) for a in (lower,upper)}
    tail=ZERO
    for h in range(parity+1,length+1):
        a=min(upper,max(lower,F(h,length)))
        if a in endpoints:
            yes,no=endpoints[a]
            term=yes[h]*no[length-h]
        else:
            term=interval(a)**h*interval(1-a)**(length-h)
        tail+=comb(length,h)*bp[h]*term
    return endpoints[lower][1][length]+tail/interval(witness.tau)**parity


def exact_transfer(witness,outputs=6):
    if any(not isinstance(getattr(witness,name),F) for name in ('x','y','tau','tilt')):
        raise ValueError('exact rational packet witnesses required')
    if not (witness.x>0 and 0<witness.y<=1 and 1<=witness.tau<=255 and witness.tilt>0):
        raise ValueError('invalid rational marks')
    local=local_matrices(witness.y,witness.tau,8,outputs,8)
    result=[[F(0),F(0)],[F(0),F(0)]]
    for weight,matrix in enumerate(local):
        factor=comb(8,weight)*witness.x**weight
        for i in range(2):
            for j in range(2):
                if not isinstance(matrix[i][j],F):
                    raise TypeError('packet arithmetic is not exact')
                result[i][j]+=factor*matrix[i][j]
    return result


def exact_period_transfer(witness,pattern):
    if not pattern or any(not isinstance(p,int) or not 1<=p<=8 for p in pattern):
        raise ValueError('invalid chronological packet period')
    result=[[F(1),F(0)],[F(0),F(1)]]
    for outputs in pattern:
        following=exact_transfer(witness,outputs)
        result=[[sum(result[i][k]*following[k][j] for k in range(2)) for j in range(2)]
                for i in range(2)]
    return result


def prefactor(witness,g):
    pattern=tuple(g['output_pattern']) if 'output_pattern' in g else (g['outputs'],)
    if g['input_n']%(8*len(pattern)):
        raise ValueError('packet period does not divide the raw input length')
    matrix=[[interval(value) for value in row] for row in exact_period_transfer(witness,pattern)]
    return (interval(conditioning(g['rows']))**128
            *terminal_power(matrix,g['input_n']//(8*len(pattern)))
            /(interval(witness.tau)*interval(witness.y)**g['output_threshold']))


def cell_upper(lo,hi,g,witness,factor):
    if not 0<lo<=hi<g['input_n']:
        raise ValueError('invalid integer middle cell')
    source=source_upper(F(lo,g['input_n']),F(hi,g['input_n']),witness,
                        length=128,parity=128-g['outer'][1])
    base=witness.x*witness.tilt
    endpoint=lo if base>=1 else hi
    return (hi-lo+1)*factor*source**g['rows']/interval(base)**endpoint


def certify(dimension=80,outputs=6,first=F(3,100),last=F(4,5),k_log=20,
            initial_grid=1000,cell_bits=128,witnesses=None):
    if witnesses is None:
        defaults={(80,6):WITNESSES_80_6,(67,5):WITNESSES_67_5}
        if (dimension,outputs) not in defaults:
            raise ValueError('supply rational witnesses for this nondefault profile')
        witnesses=defaults[dimension,outputs]
    return certify_profile(metadata(dimension,outputs,k_log),first,last,witnesses,initial_grid,cell_bits)


def certify_profile(g,first,last,witnesses,initial_grid=1000,cell_bits=128):
    if (not isinstance(first,F) or not isinstance(last,F) or initial_grid<1
            or cell_bits<40 or not witnesses):
        raise ValueError('exact interval and valid replay parameters required')
    if any(not isinstance(anchor,F) or not 0<anchor<1 for anchor,_ in witnesses):
        raise ValueError('exact witness anchors required')
    dimension=g['outer'][1]
    n=g['input_n']
    first_weight,last_weight=integer_span(n,first,last)
    factors=tuple(prefactor(w,g) for _,w in witnesses)
    def evaluate(lo,hi):
        middle=F(lo+hi,2*n)
        candidates=sorted(range(len(witnesses)),key=lambda i:abs(witnesses[i][0]-middle))[:3]
        return min(((i,cell_upper(lo,hi,g,retarget_witness(*witnesses[i],middle),factors[i]))
                    for i in candidates),key=lambda item:item[1].hi)
    maximum=interval(F(1,2**cell_bits)).lo
    span=(n+initial_grid-1)//initial_grid
    next_weight,cells,total,largest=first_weight,0,ZERO,ZERO.hi
    digest,uses=sha256(),[0]*len(witnesses)
    for lo,hi,index,value in refined_cells(first_weight,last_weight,span,evaluate,maximum):
        if lo!=next_weight or hi<lo:
            raise ArithmeticError('middle coverage is not contiguous')
        next_weight=hi+1
        cells+=1
        uses[index]+=1
        total+=value
        largest=max(largest,value.hi)
        digest.update(f'{lo}:{hi}:{index}\n'.encode('ascii'))
    if next_weight!=last_weight+1:
        raise ArithmeticError('middle coverage ends early')
    bits=integer_margin(total.hi)
    if bits<40:
        raise ArithmeticError(f'middle component below40bits: {bits}')
    return dict(**g,weights=[first_weight,last_weight],covered_weights=last_weight-first_weight+1,
        checked_cells=cells,density_min=str(first),density_max_exclusive=str(last),
        failure_upper=str(total.hi),margin_bits_lower=bits,
        minimum_cell_margin_bits_lower=integer_margin(largest),partition_sha256=digest.hexdigest(),
        witness_cell_counts=uses,conditioning='central-mode-reciprocal^128',precision_digits=80,
        source_polynomial=f'empty plus feasible support weights{129-dimension}..128; beta>=2, otherwise relaxed',
        complete=True,rigorous=True,whole_code=False,
        scope='complete middle inverse-input-weight interval only')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first',type=F,default=F(3,100))
    parser.add_argument('--last',type=F,default=F(4,5))
    parser.add_argument('--cell-bits',type=int,default=128)
    args=parser.parse_args()
    print(json.dumps(certify(first=args.first,last=args.last,cell_bits=args.cell_bits)),flush=True)
