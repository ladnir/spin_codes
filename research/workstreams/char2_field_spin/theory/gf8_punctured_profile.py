"""Shared exact identity for the candidate three-quarter-rate punctured code.

This metadata identifies the construction, not whether a complete
certificate or measured implementation has been obtained.
"""
from fractions import Fraction as F
from hashlib import sha256
import json

from gf8_mds_packet import fullstate_matrices
from gf8_punctured_screen import profile


CONSTRUCTION='gf256-full16-cauchy-w8-s8-keep7-grs128-84-regional-nonzero-diagonals'
EXPECTED_MATRIX_SHA256='6e3cf3e5d6b5876838afc5b2d7230f931974f363f9f67b402843f677521450a1'


def matrix_digest():
    # Ordered row-major bytes F(64), A(64), C(64), followed by eta(1).
    f,a,c,eta=fullstate_matrices(normalized=True)
    encoded=bytes(v for matrix in (f,a,c) for row in matrix for v in row)+bytes((eta,))
    digest=sha256(encoded).hexdigest()
    if digest!=EXPECTED_MATRIX_SHA256:
        raise ValueError('frozen punctured map identity changed; requires a new profile audit')
    return digest


def metadata(k_log=20):
    if k_log not in (16,18,20):
        raise ValueError('research profile covers requested K16,K18,K20 geometry only')
    g=profile(k_log=k_log)
    return dict(construction=CONSTRUCTION,outer=[128,84],q=256,width=8,state=8,
        outputs=7,kept_outputs=list(range(7)),eta=2,k_log=k_log,requested_k=2**k_log,
        padded_k=g['k'],rows=g['rows'],input_n=g['input_n'],output_n=g['output_n'],
        output_threshold=g['threshold'],distance='1/10',nominal_rate='3/4',
        requested_rate=str(F(2**k_log,g['output_n'])),matrix_sha256=matrix_digest(),
        matrix_digest_encoding='row-major F,A,C bytes followed by eta byte',
        field_polynomial='0x11b',initial_state='zero',terminal_state='discarded',
        input_layout='outer-column-major, then independent regional permutations',
        input_diagonals='independent uniform nonzero GF256 scalars',
        state_diagonals='independent uniform nonzero GF256 scalars',
        routing='one independent uniform permutation of the rows per outer coordinate',
        lane_interpretation='sixteen parallel scalar GF256 codewords; shared sampled code')


if __name__=='__main__':
    print(json.dumps(metadata()))
