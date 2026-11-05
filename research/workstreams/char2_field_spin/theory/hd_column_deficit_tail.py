"""Positive truncated deficit tails for up to eight original input rows.

Each supplied column matrix polynomial has activity degrees0..a. Reverse
these degrees to deficits a-w. Raise each matrix entry to r, then compose
sixteen columns in chronological order. Truncation at total deficit7ar
retains exactly the relaxed RS activity condition >=9ar. All coefficients
are nonnegative; no high-degree term can return to the retained prefix.

The input polynomials may be exact norms or a proved upper-tail/activity
majorant. Establish that premise before calling this numerical helper.
FLOAT only: no RS normalization, source count, row-set sum, or certificate.
"""
import numpy as np

from hd_pair_state_tail_general import (entrywise_power,matrix_polynomial_product,
                                        tail_from_polynomial)


def column_deficit_tail(columns,arity,groups,*,check_time=None):
    if (type(arity) is not int or not 1<=arity<=8 or type(groups) is not int
            or not 1<=groups<=18 or len(columns)!=16):
        raise ValueError('Sixteen columns, arity1..8, and groups1..18 required')
    limit=7*arity*groups
    shape=np.asarray(columns[0]).shape
    if len(shape)<3 or shape[-3:]!=(arity+1,2,2):
        raise ValueError('Activity polynomial shape (...,arity+1,2,2) required')
    word=np.full(shape[:-3]+(1,2,2),-np.inf)
    word[...,0,0,0]=word[...,0,1,1]=0.
    for column in columns:
        if check_time: check_time()
        column=np.asarray(column,dtype=float)
        if column.shape!=shape or np.any(np.isnan(column)) or np.any(np.isposinf(column)):
            raise ValueError('Matching finite log-polynomials, with -inf for zero, required')
        factor=entrywise_power(column[...,::-1,:,:],groups,check_time=check_time)
        word=matrix_polynomial_product(word,factor,check_time=check_time)
        word=word[...,:limit+1,:,:]
    return tail_from_polynomial(word,0)


def metadata():
    return dict(arithmetic='positive-FLOAT-log2-not-directed',column_count=16,
        deficit_degree='a minus original input activity or a proved upper surrogate',
        retained_maximum_deficit='7*a*r',minimum_activity='9*a*r',
        chronological_matrix_product=True,coefficient_subtraction=False,
        original_per_row_floors_relaxed=True,source_count_included=False,
        rs_normalization_included=False,whole_code_certificate=False)
