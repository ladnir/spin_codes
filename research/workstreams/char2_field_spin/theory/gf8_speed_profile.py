"""Canonical geometry for the smaller-outer/smaller-packet speed candidate."""
from fractions import Fraction as F

from gf8_punctured_screen import profile
from gf8_speed_maps import matrix_digest,FACTORED_LAMBDA


EXPECTED_MATRIX_SHA256='ecd9399dee2d87cde4c12a66704a7884d7a590d38603a3de18cfeb07292fabd3'


def metadata(dimension=40,k_log=20):
    if dimension not in (39,40) or k_log not in (16,18,20):
        raise ValueError('unsupported speed profile')
    g=profile(4,3,64,dimension,k_log)
    if matrix_digest()!=EXPECTED_MATRIX_SHA256:
        raise ArithmeticError('speed profile matrix identity changed')
    return dict(construction=f'gf256-full8-cauchy-w4-s4-keep3-grs64-{dimension}-regional-nonzero-diagonals',
        outer=[64,dimension],q=256,width=4,state=4,outputs=3,kept_outputs=[0,1,2],eta=2,
        k_log=k_log,requested_k=2**k_log,padded_k=g['k'],rows=g['rows'],
        input_n=g['input_n'],output_n=g['output_n'],output_threshold=g['threshold'],
        distance='1/10',nominal_rate=str(F(dimension,48)),requested_rate=str(F(2**k_log,g['output_n'])),
        matrix_sha256=matrix_digest(),matrix_digest_encoding='row-major F,A,C bytes followed by eta byte',
        field_polynomial='0x11b',initial_state='zero',terminal_state='discarded',
        input_layout='outer-column-major, then independent regional permutations',
        input_diagonals='independent uniform nonzero GF256 scalars',
        state_diagonals='independent uniform nonzero GF256 scalars',
        routing='one independent uniform permutation of the rows per outer coordinate',
        lane_interpretation='sixteen parallel scalar GF256 codewords; shared sampled code',
        factored_lambda=FACTORED_LAMBDA)


if __name__=='__main__':
    import json
    print(json.dumps(metadata()))
