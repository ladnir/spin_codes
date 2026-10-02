"""Source-bound memo of rebuilt local bounds and exact census details.

Arb values are stored as exact upper endpoints. The memo is not a global
certificate. Omitting the directory always rebuilds independently.
"""
import hashlib
import json
import os
from pathlib import Path
import tempfile
from fractions import Fraction

from flint import arb,arb_mat,ctx
import mass_density_screen as screen
from local_family import pack,unpack,operator_cache


def encode(value):
    if isinstance(value,arb):return {'arb':pack(value)}
    if isinstance(value,arb_mat):
        return {'matrix':[[encode(value[i,j]) for j in range(value.ncols())] for i in range(value.nrows())]}
    if isinstance(value,Fraction):return {'fraction':[value.numerator,value.denominator]}
    if isinstance(value,dict):return {'dict':[[encode(k),encode(v)] for k,v in value.items()]}
    if isinstance(value,tuple):return {'tuple':[encode(x) for x in value]}
    if isinstance(value,list):return [encode(x) for x in value]
    if value is None or type(value) in (int,str,bool):return value
    raise TypeError(type(value))


def decode(value):
    if isinstance(value,list):return [decode(x) for x in value]
    if not isinstance(value,dict):return value
    if len(value)!=1:raise ValueError('invalid memo node')
    tag,data=next(iter(value.items()))
    if tag=='arb':return unpack(data)
    if tag=='matrix':return arb_mat([[decode(x) for x in row] for row in data])
    if tag=='fraction':return Fraction(*data)
    if tag=='tuple':return tuple(decode(x) for x in data)
    if tag=='dict':return {decode(k):decode(v) for k,v in data}
    raise ValueError('unknown memo type')


def sources():
    digest=hashlib.sha256(operator_cache.source_digest().encode())
    for path in (Path(__file__),Path(screen.__file__)):
        digest.update(path.read_bytes())
    return digest.hexdigest()


def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':')).encode()


def key(tilt,penalty,rounds,precision):
    return dict(schema=1,source=sources(),tilt=str(Fraction(tilt)),penalty=str(Fraction(penalty)),
                rounds=rounds,precision=precision,windows=32,dimension=11,fresh_classes=True)


def build(tilts,penalty='.9',rounds=2,precision=192,directory=None):
    ctx.prec=precision
    tilts=tuple(dict.fromkeys(map(str,tilts)));result={};keys={};paths={}
    for tilt in tilts:
        keys[tilt]=key(tilt,penalty,rounds,precision)
        if directory is None:continue
        paths[tilt]=Path(directory)/(hashlib.sha256(canonical(keys[tilt])).hexdigest()+'.json')
        if not paths[tilt].exists():continue
        record=json.loads(paths[tilt].read_bytes());body=record['body']
        if body['key']!=keys[tilt] or record['checksum']!=hashlib.sha256(canonical(body)).hexdigest():
            raise ValueError('memo source/parameter/checksum mismatch')
        result[tilt]=decode(body['data'])
        print('ATTACK LOCAL MEMO hit',tilt,flush=True)
    missing=[t for t in tilts if t not in result]
    if missing:
        grid,feedback,details=screen.epoch_grid(missing,penalty,precision,rounds=rounds,details=True)
        for tilt in missing:
            result[tilt]=(grid[tilt],feedback,details[tilt])
            if directory is None:continue
            body={'key':keys[tilt],'data':encode(result[tilt])}
            record={'body':body,'checksum':hashlib.sha256(canonical(body)).hexdigest()}
            paths[tilt].parent.mkdir(parents=True,exist_ok=True)
            with tempfile.NamedTemporaryFile(dir=paths[tilt].parent,delete=False,suffix='.tmp') as stream:
                temporary=Path(stream.name);stream.write(canonical(record))
            os.replace(temporary,paths[tilt])
            print('ATTACK LOCAL MEMO saved exact endpoints',tilt,flush=True)
    return result


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tilts',nargs='+',default=['.056','.072','.088'])
    parser.add_argument('--directory',default='tmp/four-bit-attack-cache')
    args=parser.parse_args()
    build(args.tilts,directory=args.directory)
