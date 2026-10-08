"""Directed dense wrapper for a fixed chronological puncturing period."""
from fractions import Fraction as F
import json

from gf8_periodic_climb_profile import metadata
from gf8_rate_climb_dense import certify_profile


def certify(dimension=84,pattern=(6,6,6,7),first_density=F(17,20),k_log=20,
            cell_bits=128,initial_grid=200,precision=192,include_witnesses=False):
    return certify_profile(metadata(dimension,pattern,k_log),first_density,cell_bits,
                           initial_grid,precision,include_witnesses)


if __name__=='__main__':
    result=certify()
    result.pop('failure_upper_dyadic_hex')
    print(json.dumps(result),flush=True)
