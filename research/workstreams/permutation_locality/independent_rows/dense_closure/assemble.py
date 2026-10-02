"""Combine a replayed complete dense atlas with the earlier sparse theorem."""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
from flint import arb,ctx

import kernel
from mixture import actual_components,N
from atlas import Model,SCHEMA,replay

# Conservative published upper from ../LOW_OCCUPANCIES.md, all q=1..58,
# same maps, two updates, and cutoff 209715. This driver relies on that
# earlier certificate; it does not rerun its 58 individual support covers.
SPARSE_UPPER=Q('6.786362e-14')
SPARSE_EXTENSION=Q('5.789821e-40')  # q=59..64 at cutoff 104857; SPARSE_HANDOFF.md.


def sparse_upper(minimum_groups,threshold):
    if type(minimum_groups) is not int or type(threshold) is not int or threshold<0:
        raise ValueError('integer occupancy threshold and nonnegative integer cutoff required')
    if minimum_groups==59 and threshold<=209715:return SPARSE_UPPER
    if minimum_groups==65 and threshold<=104857:return SPARSE_UPPER+SPARSE_EXTENSION
    raise ValueError('no matching complete sparse theorem for this occupancy handoff and cutoff')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('dense',type=Path)
    parser.add_argument('--precision',type=int,default=256)
    args=parser.parse_args();ctx.prec=args.precision
    record=json.loads(args.dense.read_text())
    if record.get('schema')!=SCHEMA:parser.error('matching four-bit atlas schema required')
    try:sparse=sparse_upper(record.get('minimum_groups'),record.get('threshold'))
    except ValueError as error:parser.error(str(error))
    if record.get('unresolved') or record.get('screen_only'):
        parser.error('a complete outward atlas is required; this record is incomplete or screening-only')
    components=actual_components(record['central_bits'],Q(record['theta']))
    data=kernel.actual();ctx.prec=args.precision
    model=Model(components,data,record['minimum_groups'],record['threshold'],[1,*map(Q,record['input_tilt'])])
    print('OUTWARD PRECISION',ctx.prec,flush=True)
    dense=replay(model,record)
    total=kernel.up(dense+kernel.aq(sparse))
    if not 0<total<arb(2)**-40:raise ArithmeticError('complete sum misses the 40-bit target')
    print('COMPLETE FOUR-BIT CERTIFICATE; ideal independent setup; two updates',flush=True)
    print('minimum distance >=',record['threshold']+1,'of',N,flush=True)
    print('bad-setup probability upper',total,flush=True)
    print('margin lower enclosure',-total.log()/arb(2).log(),flush=True)
    print('Uses the documented sparse theorem for q=1..',record['minimum_groups']-1,
          '; dense q>=',record['minimum_groups'],'replayed here.',flush=True)


if __name__=='__main__':main()
