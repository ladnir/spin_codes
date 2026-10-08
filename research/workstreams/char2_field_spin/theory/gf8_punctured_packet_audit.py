"""Independent fixed-row puncturing algebra/layout audit; no bound optimizer.

The full8/8 packet is reused unchanged. Only its firstm emitted coordinates
are retained. Checks include global rank, independence from input diagonals,
literal adjoints, and a four-packet coordinate-major payload layout.
"""
from random import Random
import json

from gf8_packet_maps import MUL,determinant,matvec,rank
from gf8_mds_packet import fullstate_matrices,fullstate_factored_packet,fullstate_factored_four_packet_adjoint


def forward(values,input_diagonals,state_diagonals,outputs=7):
    state=(0,)*8
    result=[]
    for start,dx,ds in zip(range(0,len(values),8),input_diagonals,state_diagonals):
        emitted,state=fullstate_factored_packet(values[start:start+8],state,dx,ds)
        result.extend(emitted[:outputs])
    return tuple(result)


def adjoint(values,input_diagonals,state_diagonals,outputs=7):
    state=(0,)*8
    packets=len(input_diagonals)
    result=[]
    for packet in range(packets-1,-1,-1):
        padded=tuple(values[packet*outputs:(packet+1)*outputs])+(0,)*(8-outputs)
        emitted,state=fullstate_factored_packet(padded,state,input_diagonals[packet],state_diagonals[packet],True)
        result.append(emitted)
    return tuple(v for row in reversed(result) for v in row)


def matrix(input_diagonals,state_diagonals,outputs=7):
    size=8*len(input_diagonals)
    columns=[]
    for index in range(size):
        unit=tuple(int(i==index) for i in range(size))
        columns.append(forward(unit,input_diagonals,state_diagonals,outputs))
    return tuple(zip(*columns))


def pack_four(emitted,outputs=7):
    """7 or6 ZMMs; each has four records in forward packet chronology."""
    return tuple(emitted[packet][coordinate] for coordinate in range(outputs) for packet in range(4))


def unpack_four(payload,outputs=7):
    """Return reverse-chronological8-coordinate adjoint packets, zero-filled."""
    return tuple(tuple(payload[4*coordinate+3-packet] for coordinate in range(outputs))+(0,)*(8-outputs)
                 for packet in range(4))


def check():
    rng=Random(7084105)
    global_cases=packing_cases=0
    f,_,_,_=fullstate_matrices()
    for outputs in (4,6,7,8):
        # A fixed square pivot block works at every packet. The Cauchy
        # property proves the same for any equally sized column set.
        assert determinant(tuple(row[:outputs] for row in f[:outputs]))
        for packets in (1,2,4,5):
            for _ in range(4):
                dx=tuple(tuple(rng.randrange(1,256) for _ in range(8)) for _ in range(packets))
                ds=tuple(tuple(rng.randrange(1,256) for _ in range(8)) for _ in range(packets))
                direct=matrix(dx,ds,outputs)
                h=matrix(((1,)*8,)*packets,ds,outputs)
                diagonal=tuple(v for row in dx for v in row)
                assert direct==tuple(tuple(MUL[v][d] for v,d in zip(row,diagonal)) for row in h)
                assert rank(direct)==outputs*packets
                pivots=tuple(8*p+j for p in range(packets) for j in range(outputs))
                assert determinant(tuple(tuple(row[j] for j in pivots) for row in direct))
                output=tuple(rng.randrange(256) for _ in range(outputs*packets))
                assert adjoint(output,dx,ds,outputs)==matvec(tuple(zip(*direct)),output)
                global_cases+=1
    for outputs in (6,7):
        for _ in range(128):
            emitted=tuple(tuple(rng.randrange(256) for _ in range(8)) for _ in range(4))
            payload=pack_four(emitted,outputs)
            reverse=unpack_four(payload,outputs)
            assert reverse==tuple(tuple(row[:outputs])+(0,)*(8-outputs) for row in reversed(emitted))
            dx=tuple(tuple(rng.randrange(1,256) for _ in range(8)) for _ in range(4))
            ds=tuple(tuple(rng.randrange(1,256) for _ in range(8)) for _ in range(4))
            terminal=tuple(rng.randrange(256) for _ in range(8))
            got,carry=fullstate_factored_four_packet_adjoint(reverse,terminal,dx,ds)
            expected=[]
            for values,d1,d2 in zip(reverse,dx,ds):
                values,terminal=fullstate_factored_packet(values,terminal,d1,d2,True)
                expected.append(values)
            assert got==tuple(expected) and carry==terminal
            packing_cases+=1
    return dict(global_rank_input_diagonal_factorization_and_transpose=global_cases,
                four_packet_puncture_zero_fill_and_chronology=packing_cases,
                setup_lane_count=16,same_scalar_code_in_every_lane=True,native_work=False)


if __name__=='__main__':print(json.dumps(check()))
