"""Search or replay the three-type cover with no-constant expansion maps."""
import argparse
from pathlib import Path

import parameter_no_constant as maps
import verify_parameter_dense_v2 as replay


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--input',type=Path)
    p.add_argument('--minimum',type=int,default=65)
    p.add_argument('--nodes',type=int,default=1023)
    a = p.parse_args()
    with maps.use():
        if a.input:
            replay.verify(a.input.resolve(),a.output.resolve())
        else:
            maps.dense.run(a.output.resolve(),128,64,20,20,a.minimum,a.nodes)
