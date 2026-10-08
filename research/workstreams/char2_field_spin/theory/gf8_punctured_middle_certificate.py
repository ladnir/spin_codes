"""Directed middle certificate for the frozen seven-of-eight GF256 code.

The packet cover is evaluated with exact rational witnesses. Positive
source and matrix powers are rounded outward at 80 decimal digits.
All integer inverse-input weights in the reported interval are covered.
This component is not a whole-code certificate.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb

from gf8_outer_rate_middle import conditioning
from gf8_punctured_profile import metadata
from gf8_punctured_screen import local_matrices
from return_cover_certificate import Witness, interval, row_upper, terminal_power
from state8_density_certificate import integer_span, refined_cells, retarget_witness
from subfield_certificate import HIGH, LOW, ZERO, integer_margin


LOW.Emax = HIGH.Emax = 1_000_000_000
LOW.Emin = HIGH.Emin = -1_000_000_000
WITNESSES = tuple((F(density), Witness(*(F(v) for v in values))) for density,values in (
    ('1/50', ('1179/6175', '6235/7954', '61627/9963', '9042/1381')),
    ('3/100', ('2039/6600', '6699/9767', '56159/8782', '25990/6023')),
    ('1/20', ('268/443', '4689/9151', '28420/4071', '18865/7552')),
    ('1/10', ('13441/7073', '1295/5993', '78757/8041', '3545/3309')),
    ('7/50', ('2999/758', '882/8971', '22514/1553', '61/97')),
    ('1/5', ('37107/3433', '29/1029', '50249/1620', '2521/8752')),
    ('3/10', ('525285/6566', '2/919', '2101945/9297', '64/1145')),
    ('2/5', ('305713/3105', '12/9265', '255', '101/3006')),
    ('9/20', ('339799/3344', '2/1799', '255', '105/3827'))))


def exact_transfer(witness):
    if any(not isinstance(getattr(witness,name),F) for name in ('x','y','tau','tilt')):
        raise ValueError('exact rational packet witnesses required')
    if not (witness.x>0 and 0<witness.y<=1 and 1<=witness.tau<=255 and witness.tilt>0):
        raise ValueError('invalid packet witnesses')
    matrices = local_matrices(witness.y,witness.tau,8,7,8)
    result = [[F(0),F(0)],[F(0),F(0)]]
    for weight,matrix in enumerate(matrices):
        factor = comb(8,weight)*witness.x**weight
        for i in range(2):
            for j in range(2):
                if not isinstance(matrix[i][j],F):
                    raise TypeError('local transfer escaped exact arithmetic')
                result[i][j] += factor*matrix[i][j]
    return result


def source_upper(lower,upper,witness,length=128,parity=44):
    """Uniform positive source bound, retaining the MDS support cutoff.

    beta>=2 permits the audited Bernoulli homogenization. Each nonempty
    term alpha^h(1-alpha)^(length-h) is maximized at h/length, clipped
    into the cell. The empty term is maximized at its lower endpoint.
    """
    if not (isinstance(lower,F) and isinstance(upper,F) and 0<lower<=upper<1):
        raise ValueError('exact interior density cell required')
    beta = witness.tau*witness.tilt
    if beta < 2:
        return row_upper(lower,upper,witness,length,parity)
    empty = interval(1-lower)**length
    bp = [interval(beta)**h for h in range(length+1)]
    endpoint_powers = {}
    for point in (lower,upper):
        endpoint_powers[point] = ([interval(point)**h for h in range(length+1)],
                                  [interval(1-point)**h for h in range(length+1)])
    tail = ZERO
    for h in range(parity+1,length+1):
        point = min(upper,max(lower,F(h,length)))
        if point in endpoint_powers:
            yes,no = endpoint_powers[point]
            term = yes[h]*no[length-h]
        else:
            term = interval(point)**h*interval(1-point)**(length-h)
        tail += comb(length,h)*bp[h]*term
    return empty+tail/interval(witness.tau)**parity


def prefactor(witness,g):
    matrix = [[interval(value) for value in row] for row in exact_transfer(witness)]
    return (interval(conditioning(g['rows']))**128
            *terminal_power(matrix,g['input_n']//8)
            /(interval(witness.tau)*interval(witness.y)**g['output_threshold']))


def cell_upper(lo,hi,g,witness,factor):
    if not 0<lo<=hi<g['input_n']:
        raise ValueError('invalid middle integer cell')
    source = source_upper(F(lo,g['input_n']),F(hi,g['input_n']),witness)
    base = witness.x*witness.tilt
    endpoint = lo if base>=1 else hi
    return (hi-lo+1)*factor*source**g['rows']/interval(base)**endpoint


def certify(first=F(1,50),last=F(9,20),initial_grid=1000,cell_bits=128,
            witnesses=WITNESSES,k_log=20):
    if (not isinstance(first,F) or not isinstance(last,F) or initial_grid<1
            or cell_bits<40 or not witnesses):
        raise ValueError('exact interval and valid replay settings required')
    if any(not isinstance(anchor,F) or not 0<anchor<1 for anchor,_ in witnesses):
        raise ValueError('exact witness anchors required')
    g = metadata(k_log)
    n = g['input_n']
    first_weight,last_weight = integer_span(n,first,last)
    factors = tuple(prefactor(w,g) for _,w in witnesses)

    def evaluate(lo,hi):
        middle = F(lo+hi,2*n)
        candidates = sorted(range(len(witnesses)),key=lambda i:abs(witnesses[i][0]-middle))[:3]
        return min(((i,cell_upper(lo,hi,g,retarget_witness(*witnesses[i],middle),factors[i]))
                    for i in candidates),key=lambda item:item[1].hi)

    maximum = interval(F(1,2**cell_bits)).lo
    span = (n+initial_grid-1)//initial_grid
    next_weight,cells,total,largest = first_weight,0,ZERO,ZERO.hi
    digest,uses = sha256(),[0]*len(witnesses)
    for lo,hi,index,value in refined_cells(first_weight,last_weight,span,evaluate,maximum):
        if lo != next_weight or hi<lo:
            raise ArithmeticError('middle cells are not contiguous')
        next_weight = hi+1
        cells += 1
        uses[index] += 1
        total += value
        largest = max(largest,value.hi)
        digest.update(f'{lo}:{hi}:{index}\n'.encode('ascii'))
    if next_weight != last_weight+1:
        raise ArithmeticError('middle coverage ends early')
    bits = integer_margin(total.hi)
    if bits<40:
        raise ArithmeticError(f'middle component below40bits: {bits}')
    return dict(**g,weights=[first_weight,last_weight],covered_weights=last_weight-first_weight+1,
        checked_cells=cells,density_min=str(first),density_max_exclusive=str(last),
        failure_upper=str(total.hi),margin_bits_lower=bits,
        minimum_cell_margin_bits_lower=integer_margin(largest),partition_sha256=digest.hexdigest(),
        witness_cell_counts=uses,conditioning='central-mode-reciprocal^128',precision_digits=80,
        source_polynomial='empty plus feasible support weights45..128; beta>=2, otherwise relaxed',
        complete=True,rigorous=True,whole_code=False,
        scope='complete middle inverse-input-weight interval only')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first',type=F,default=F(1,50))
    parser.add_argument('--last',type=F,default=F(9,20))
    parser.add_argument('--cell-bits',type=int,default=128)
    args=parser.parse_args()
    print(json.dumps(certify(args.first,args.last,cell_bits=args.cell_bits)),flush=True)
