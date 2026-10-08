"""Identities for the higher-rate full-Cauchy punctured research family.

This module does not certify a parameter choice. Frozen [128,84]/keep7
receipts remain in gf8_punctured_profile and are not relabeled here.
"""
from fractions import Fraction as F

from gf8_punctured_profile import metadata as frozen_metadata, matrix_digest
from gf8_punctured_screen import profile


def metadata(dimension=84, outputs=7, k_log=20):
    if (not isinstance(dimension,int) or not isinstance(outputs,int)
            or not 1 <= outputs <= 8 or not 1 <= dimension < min(128,16*outputs)
            or k_log not in (16,18,20)):
        raise ValueError('invalid rate-climb geometry')
    g=profile(width=8,outputs=outputs,length=128,dimension=dimension,k_log=k_log)
    result=frozen_metadata(k_log)
    result.update(
        construction=f'gf256-full16-cauchy-w8-s8-keep{outputs}-grs128-{dimension}-regional-nonzero-diagonals',
        outer=[128,dimension],outputs=outputs,kept_outputs=list(range(outputs)),
        requested_k=2**k_log,padded_k=g['k'],rows=g['rows'],
        input_n=g['input_n'],output_n=g['output_n'],output_threshold=g['threshold'],
        nominal_rate=str(F(dimension,16*outputs)),
        requested_rate=str(F(2**k_log,g['output_n'])),matrix_sha256=matrix_digest())
    return result


if __name__=='__main__':
    import argparse,json
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dimension',type=int,default=84)
    parser.add_argument('--outputs',type=int,default=7)
    args=parser.parse_args()
    print(json.dumps(metadata(args.dimension,args.outputs)))
