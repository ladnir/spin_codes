"""Recompute or replay all Q1 parameter slices with no-constant IMT maps."""
import argparse
from pathlib import Path

import parameter_no_constant as maps
import verify_parameter_q1 as replay


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--input',type=Path)
    a = p.parse_args()
    with maps.use():
        if a.input:
            replay.verify(a.input.resolve(),a.output.resolve())
        else:
            maps.grid.run(a.output.resolve())
