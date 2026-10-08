"""Exact K20 dense-class certificate for GRS[128,87], width16/state8.

Frozen rational witnesses cover every inverse-input density in [.24,1].
Only integers and Fractions participate in acceptance; no optimizer or native
encoder is invoked. This certifies a base-field class, not the whole code.
"""

import argparse
from fractions import Fraction as F
import hashlib
import json


GRID = 200
FIRST = 48
# Previously proposed for GRS[32,22], now frozen and independently rechecked.
_TILT_TEXT = '''
209/15 163/12 53/4 155/12 63/5 123/10 12 82/7 103/9 123/11
153/14 171/16 167/16 143/14 10 127/13 86/9 131/14 119/13 143/16
114/13 103/12 101/12 107/13 113/14 87/11 31/4 91/12 119/16 51/7
50/7 7 103/15 101/15 33/5 97/15 89/14 81/13 55/9 6
47/8 75/13 17/3 61/11 49/9 16/3 21/4 36/7 81/16 79/16
34/7 62/13 14/3 55/12 9/2 53/12 13/3 17/4 25/6 41/10
4 63/16 58/15 19/5 41/11 51/14 43/12 7/2 31/9 44/13
53/16 13/4 51/16 25/8 46/15 3 47/16 26/9 17/6 36/13
19/7 8/3 34/13 23/9 5/2 27/11 12/5 33/14 23/10 9/4
11/5 28/13 19/9 33/16 2 2 27/14 17/9 11/6 9/5
7/4 12/7 5/3 18/11 19/12 14/9 3/2 22/15 23/16 7/5
15/11 4/3 13/10 5/4 11/9 19/16 15/13 9/8 11/10 16/15 17/16
'''
TILTS = tuple(map(F,_TILT_TEXT.split()))+(F(1),)*31


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def geometry():
    rows = 64*((2**20+87*64-1)//(87*64))
    n = 128*rows
    return dict(rows=rows,k=87*rows,n=n,threshold=(n+9)//10-1,
                dense_first=(6*n+24)//25)


def interval_base(length, dimension, index, tilt):
    """Endpoint upper bound on the unnormalized MDS row bracket."""
    if not (0 <= index < GRID and tilt >= 1):
        raise ValueError('invalid interval or tilt')
    lower,upper = F(index,GRID),F(index+1,GRID)
    return ((1-lower)**length + F(1,256**(length-dimension))
            *(1+(F(256,255)*tilt-1)*upper)**length)


def certify_interval(index, tilt):
    """Check the small witness, fourth-power lift, and target inequality."""
    small = interval_base(32,22,index,tilt)
    target = interval_base(128,87,index,tilt)
    return (small**GRID * 2**(256*32) <= tilt**(index*32)
            and target <= small**4
            and target**GRID * 2**(256*128) <= tilt**(index*128))


def weight_intervals(n):
    """Disjoint integer cells inside the closed density cover."""
    return tuple(((index*n+GRID-1)//GRID,
                  min(n,((index+1)*n+GRID-1)//GRID-1) if index<GRID-1 else n)
                 for index in range(FIRST,GRID))


def certify(show_tilts=False):
    g = geometry()
    require(g == dict(rows=12096,k=1052352,n=1548288,threshold=154828,
                      dense_first=371590),'unexpected finite geometry')
    require(len(TILTS)==GRID-FIRST,'incomplete rational witness list')
    intervals = weight_intervals(g['n'])
    require(intervals[0][0]==g['dense_first'] and intervals[-1][1]==g['n'],
            'incorrect dense endpoints')
    require(all(lo<=hi for lo,hi in intervals)
            and all(left[1]+1==right[0] for left,right in zip(intervals,intervals[1:])),
            'gap or overlap in integer coverage')
    for index,tilt in enumerate(TILTS,FIRST):
        require(certify_interval(index,tilt),f'failed exact density cell {index}')
    require(F(10,9)**100*2295**10 < 2**127,'failed radius-.1 ball bound')
    require(g['rows']+1 <= 2**g['rows'].bit_length(),'failed conditioning bound')
    margin = g['n']//100-128*g['rows'].bit_length()
    require(margin==13690,'unexpected dense margin')
    witness_text = ';'.join(f'{index}:{tilt}' for index,tilt in enumerate(TILTS,FIRST))
    result = dict(outer=[128,87],q=256,width=16,state=8,k_log=20,
        requested_k=2**20,rows=g['rows'],padded_source=g['k'],output_length=g['n'],
        forbidden_output_threshold=g['threshold'],dense_weights=[g['dense_first'],g['n']],
        covered_weights=sum(hi-lo+1 for lo,hi in intervals),dense_intervals_checked=len(TILTS),
        membership_exponent='32/25',ball_exponent='127/100',
        dense_margin_bits_lower=margin,failure_upper_power_of_two=-margin,
        witness_sha256=hashlib.sha256(witness_text.encode('ascii')).hexdigest(),
        dense_rigorous=True,cutoff_cover_complete=True,whole_code_certified=False,
        scope='GF256 dense inverse-input class; independent invertible inner',
        construction='fixed GRS; independent regional permutations and nonzero input diagonals')
    if show_tilts:
        result['rational_tilts'] = [(index,tilt.numerator,tilt.denominator)
                                   for index,tilt in enumerate(TILTS,FIRST)]
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--show-tilts',action='store_true')
    args = parser.parse_args()
    print(json.dumps(certify(args.show_tilts)),flush=True)
