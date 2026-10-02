"""Transfer shape-weighted operators to independent uniform GF16 packets.

The change of measure is on the *new* packet values, not on BCH words.
No weighted outer spectrum or lane-permutation certificate is imported.
"""
from copy import copy
from fractions import Fraction as Q
from math import comb
import numpy as np
from flint import arb,arb_mat
import sparse_cover


def normalizer(full_penalty,weight_tilt=1):
    rho,a=Q(full_penalty),Q(weight_tilt)
    if not 0<rho<=1 or a<=0:raise ValueError('positive shape weights required')
    return sum((Q(comb(4,w),15)/(a**w*rho**int(w==4)) for w in range(1,5)),Q(0))


def build_operators(args,penalties,weight_tilt='1'):
    """Each region with j active groups contains exactly j active packets."""
    weighted=copy(args);weighted.penalties=list(penalties);weighted.weight_tilt=weight_tilt
    raw=sparse_cover.sparse.build_operators(weighted)
    result={}
    for (tilt,rho),(matrices,_) in raw.items():
        c=normalizer(rho,weight_tilt);factor=arb(c.numerator)/c.denominator
        exact=[sparse_cover.sparse.up(factor**j)*matrix for j,matrix in enumerate(matrices)]
        exact=[arb_mat([[sparse_cover.sparse.up(matrix[i,j]) for j in range(matrix.ncols())]
                       for i in range(matrix.nrows())]) for matrix in exact]
        result[tilt,rho]=exact,[np.array([[float(matrix[i,j]) for j in range(matrix.ncols())]
                                       for i in range(matrix.nrows())]) for matrix in exact]
        print('GF16 shape change of measure',tilt,'rho',rho,'a',weight_tilt,'normalizer',str(c),flush=True)
    return result
