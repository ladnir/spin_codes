"""Extend or replay the fixed no-constant IMT pilot's sparse range."""
import argparse
from pathlib import Path

import parameter_no_constant as maps
import parameter_sparse_range as producer
import verify_parameter_range as replay


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--input',type=Path)
    parser.add_argument('--first',type=int,default=65)
    parser.add_argument('--last',type=int,default=256)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(output)
    with maps.use():
        if args.input:
            replay.verify(args.input.resolve(),output)
        else:
            producer.run(output,128,64,20,20,args.first,args.last)
