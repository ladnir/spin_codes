"""Search or replay sparse IMT bounds with no-constant expansion maps."""
import argparse
from pathlib import Path

import parameter_no_constant as maps
import parameter_sparse_balanced as compositions
import parameter_sparse_range as ranges
import verify_parameter_sparse as composition_replay
import verify_parameter_range as range_replay


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--component',choices=('compositions','range'),required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--input',type=Path)
    a = p.parse_args()
    with maps.use():
        if a.component == 'compositions':
            if a.input:
                composition_replay.verify(a.input.resolve(),a.output.resolve())
            else:
                compositions.run(a.output.resolve(),4)
        elif a.input:
            range_replay.verify(a.input.resolve(),a.output.resolve())
        else:
            ranges.run(a.output.resolve(),128,64,20,20,5,64)
