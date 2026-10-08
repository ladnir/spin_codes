"""Directed middle replay for GRS64/40, Cauchy4/4, fixed keep3.

Every local coefficient is rational and every accepted integer interval
uses outward rounding. The earlier length128 certificates are unchanged.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import comb

from gf8_outer_rate_middle import conditioning
from gf8_punctured_screen import local_matrices
from gf8_rate_climb_middle_certificate import source_upper
from gf8_speed_profile import metadata
from return_cover_certificate import Witness,interval,terminal_power
from state8_density_certificate import integer_span,refined_cells,retarget_witness
from subfield_certificate import HIGH,LOW,ZERO,integer_margin


LOW.Emax=HIGH.Emax=1_000_000_000
LOW.Emin=HIGH.Emin=-1_000_000_000
WITNESSES=tuple((F(anchor),Witness(*(F(v) for v in values))) for anchor,values in (
    ('3/100',('3097/9748','3655/5679','72067/8789','37768/8651')),
    ('1/20',('3556/5693','1048/2299','60896/7097','23787/9253')),
    ('1/10',('14180/7761','1092/6191','104654/9225','10897/9636')),
    ('7/50',('31286/9339','118/1431','50627/3204','3538/5171')),
    ('1/5',('67323/9698','261/9134','198916/7065','3364/9683')),
    ('3/10',('160727/8222','43/7540','186019/2257','62/523')),
    ('2/5',('394584/7787','3/2345','255','197/4807')),
    ('1/2',('517261/9730','8/8315','255','89/3231')),
    ('3/5',('299015/5456','2/2595','255','57/3077')),
    ('7/10',('250517/4476','2/3111','255','47/3905')),
    ('4/5',('519556/9141','5/9068','255','29/4062'))))


def exact_transfer(witness,g):
    if any(not isinstance(getattr(witness,name),F) for name in ('x','y','tau','tilt')):
        raise ValueError('exact rational marks required')
    if not(witness.x>0 and witness.tilt>0 and 0<witness.y<=1 and 1<=witness.tau<=255):
        raise ValueError('invalid exact marks')
    local=local_matrices(witness.y,witness.tau,g['width'],g['outputs'],g['state'])
    matrix=[[F(0),F(0)],[F(0),F(0)]]
    for w,term in enumerate(local):
        factor=comb(g['width'],w)*witness.x**w
        for i in range(2):
            for j in range(2):
                if not isinstance(term[i][j],F):
                    raise TypeError('local transfer escaped rational arithmetic')
                matrix[i][j]+=factor*term[i][j]
    return matrix


def prefactor(witness,g):
    matrix=[[interval(v) for v in row] for row in exact_transfer(witness,g)]
    return (interval(conditioning(g['rows']))**g['outer'][0]
            *terminal_power(matrix,g['input_n']//g['width'])
            /(interval(witness.tau)*interval(witness.y)**g['output_threshold']))


def cell_upper(lo,hi,g,witness,factor):
    if not 0<lo<=hi<g['input_n']:
        raise ValueError('invalid middle weight interval')
    length,dimension=g['outer']
    source=source_upper(F(lo,g['input_n']),F(hi,g['input_n']),witness,length,length-dimension)
    base=witness.x*witness.tilt
    endpoint=lo if base>=1 else hi
    return (hi-lo+1)*factor*source**g['rows']/interval(base)**endpoint


def certify(dimension=40,k_log=20,first=F(3,100),last=F(4,5),
            initial_grid=1000,cell_bits=128):
    if not isinstance(first,F) or not isinstance(last,F) or cell_bits<40 or initial_grid<1:
        raise ValueError('invalid exact middle replay parameters')
    g=metadata(dimension,k_log)
    n=g['input_n']
    lo0,hi0=integer_span(n,first,last)
    factors=tuple(prefactor(w,g) for _,w in WITNESSES)
    def evaluate(lo,hi):
        middle=F(lo+hi,2*n)
        choices=sorted(range(len(WITNESSES)),key=lambda i:abs(WITNESSES[i][0]-middle))[:3]
        return min(((i,cell_upper(lo,hi,g,retarget_witness(*WITNESSES[i],middle),factors[i]))
                    for i in choices),key=lambda item:item[1].hi)
    maximum=interval(F(1,1<<cell_bits)).lo
    span=(n+initial_grid-1)//initial_grid
    expected,total,cells=lo0,ZERO,0
    digest=sha256()
    uses=[0]*len(WITNESSES)
    largest=ZERO.hi
    for lo,hi,index,value in refined_cells(lo0,hi0,span,evaluate,maximum):
        if lo!=expected or hi<lo:
            raise ArithmeticError('middle interval coverage gap')
        expected=hi+1
        total+=value
        cells+=1
        uses[index]+=1
        largest=max(largest,value.hi)
        digest.update(f'{lo}:{hi}:{index}\n'.encode('ascii'))
    if expected!=hi0+1:
        raise ArithmeticError('middle interval ended early')
    bits=integer_margin(total.hi)
    if bits<40:
        raise ArithmeticError('middle failure bound exceeds40bit target')
    return dict(**g,weights=[lo0,hi0],covered_weights=hi0-lo0+1,checked_cells=cells,
        density_min=str(first),density_max_exclusive=str(last),failure_upper=str(total.hi),
        margin_bits_lower=bits,minimum_cell_margin_bits_lower=integer_margin(largest),
        partition_sha256=digest.hexdigest(),witness_cell_counts=uses,
        conditioning='central-mode-reciprocal^64',precision_digits=80,
        source_polynomial=f'empty plus feasible support weights{65-dimension}..64; beta>=2, otherwise relaxed',
        complete=True,rigorous=True,whole_code=False,scope='complete middle inverse-input-weight interval only')


if __name__=='__main__':
    import json
    print(json.dumps(certify()),flush=True)
