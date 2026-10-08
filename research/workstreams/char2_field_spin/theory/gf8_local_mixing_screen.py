"""Exact local-mixer laws and obstruction checks, not a SPIN certificate.

The candidate is independent GL(2,q) mixing of packet-output pairs. A
uniform local matching admits a positive fixed-word selected-zero kernel.
That kernel is not an automatic replacement for the source-rank transfer.
No native code, setup timing, or performance measurement is included.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
import json
from math import comb


def gl_pair_nonzero_moment(y, q=256):
    """Output-weight moment of a nonzero pair under independent uniform GL2."""
    y=F(y)
    if q<2 or not 0<y<=1:
        raise ValueError('valid field size and output mark required')
    return (2*y+(q-1)*y*y)/(q+1)


def fixed_mds_pair_moment(weight, y, q=256):
    """All-nonzero invertible 2x2 matrix after independent nonzero diagonals."""
    y=F(y)
    if weight not in (0,1,2) or q<3 or not 0<y<=1:
        raise ValueError('valid pair weight, field size and mark required')
    if weight==0:
        return F(1)
    if weight==1:
        return y*y
    return y*y+F(2,q-1)*y*(1-y)


def matching_zero_coefficients(y, q=256, width=16):
    """Positive coefficients c_z in h_t=sum_z C(width-t,z)c_z.

    Each pair receives an independent uniform GL2. The pairing is uniform
    over all perfect matchings. The input vector is fixed before sampling.
    """
    if width<=0 or width%2:
        raise ValueError('positive even packet width required')
    g=gl_pair_nonzero_moment(y,q)
    pairs=width//2
    coefficients=[F(0)]*(width+1)
    for j in range(pairs+1):
        coefficients[2*j]=F(comb(pairs,j),comb(width,2*j))*g**(pairs-j)*(1-g)**j
    return tuple(coefficients)


def matching_moment(weight,y,q=256,width=16):
    if not 0<=weight<=width:
        raise ValueError('weight outside packet')
    coefficients=matching_zero_coefficients(y,q,width)
    return sum(comb(width-weight,z)*coefficients[z] for z in range(width-weight+1))


def gl_pair_rank_moment(form_rank, selected_rows, tau, q=256):
    """Exact E[tau^-rank(G_selected L)] for rank(L) in {0,1,2}.

    G is fresh uniform GL(2,q). This scalar moment does not encode the
    quotient's orientation or a following packet's state-cap continuation.
    """
    tau=F(tau)
    if form_rank not in (0,1,2) or selected_rows not in (0,1,2) or tau<1 or q<2:
        raise ValueError('invalid rank-moment parameters')
    if not selected_rows or not form_rank:
        return F(1)
    if selected_rows==2:
        return tau**(-form_rank)
    if form_rank==2:
        return 1/tau
    return F(1,q+1)+F(q,q+1)/tau


def matchings(labels):
    labels=tuple(labels)
    if not labels:
        yield ()
        return
    if len(labels)%2:
        raise ValueError('even label count required')
    first=labels[0]
    for index in range(1,len(labels)):
        partner=labels[index]
        remaining=labels[1:index]+labels[index+1:]
        for rest in matchings(remaining):
            yield ((first,partner),)+rest


def archived_flat_witness():
    """Adjacent GL2/GL4 output blocks preserve the existing four-zero flat."""
    from cheaper_inner_map import coefficient_matrix
    from state_capacity_reference import matvec,rank
    matrix=coefficient_matrix(0)
    state=tuple(1 if j==0 else 80 if j==4 else 0 for j in range(8))
    output=matvec(matrix,state)
    zeros=tuple(i for i,x in enumerate(output) if not x)
    assert zeros==(0,1,2,3)
    restricted=tuple(tuple(matrix[i][j] for j in (0,4)) for i in zeros)
    assert rank(restricted)==1
    return dict(state_support=[0,4],state_values=[1,80],zero_outputs=list(zeros),
                restricted_rank=1,zero_event_probability='1/255',
                survives_adjacent_pair_output_mixing=True,
                survives_adjacent_width4_output_mixing=True,
                whole_code_bad_word=False)


def accounting():
    """Four 128-bit records per ZMM; static pair swap plus two GF products."""
    return dict(packet_records=16,record_bits=128,baseline_gfni=70,
        extra_gfni=8,extra_lane_shuffles=4,extra_xors=4,
        fresh_field_coefficients=32,postmix_gfni=78,
        excludes=['local matching','coefficient materialization','loads','stores'],
        native_measurement=False,full_certificate=False)


if __name__=='__main__':
    print(json.dumps(dict(cost=accounting(),preserved_obstruction=archived_flat_witness())),flush=True)
    for weight in (1,2,4,8,12,16):
        y=F(9,10)
        value=matching_moment(weight,y)
        print(json.dumps(dict(input_weight=weight,old_moment=float(y**weight),
            ideal_matching_gl2_moment=float(value),fixed_word_only=True,
            whole_code_certificate=False)),flush=True)
