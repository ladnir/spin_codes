"""Reusable outward local bounds for the experimental witness searches.

The optional memo binds every generator source, map, precision and option.
It stores exact upper endpoints, not binary64 proposals or certificates.
Omit the cache directory to rebuild all censuses independently.
"""
import hashlib
import json
import os
from pathlib import Path
import tempfile

from flint import arb, arb_mat, ctx

from mass_density_screen import epoch_grid, baseline
from mass_density import conditioned_coefficients, rational
import operator_cache


def source_digest():
    digest=hashlib.sha256(operator_cache.source_digest().encode())
    for path in sorted(Path(__file__).resolve().parent.glob('*.py')):
        digest.update(path.name.encode()+b'\0'+path.read_bytes()+b'\0')
    return digest.hexdigest()


def parameters(tilt,penalty,precision,maximum,exact_feedback,spectral_feedback,
               joint_four,density_through,sources):
    tilt,penalty=map(rational,(tilt,penalty))
    if (tilt<=0 or not 0<penalty<=1 or precision<128 or not 6<=maximum<=32
            or exact_feedback not in (0,7,8,9,10) or spectral_feedback not in (0,7,8,9,10)
            or (exact_feedback and spectral_feedback) or density_through not in (0,7,8,9,10)):
        raise ValueError('invalid local family parameters')
    return dict(schema=1,stage='independent-row-local-family',sources=sources,
        tilt=str(tilt),penalty=str(penalty),precision=precision,maximum=maximum,
        exact_feedback=exact_feedback,spectral_feedback=spectral_feedback,
        joint_four=bool(joint_four),density_through=density_through,
        windows=32,dimension=11,rounds=2,input_weight='1',full_feedback=6,
        window_histogram=8,joint_cancellation=True,column_density=True,
        feedback_density=6,mature_tail=64)


def cache_path(directory,key):
    return Path(directory)/('local-family-'+hashlib.sha256(operator_cache.encoded(key)).hexdigest()+'.json')


def pack(value):
    if not value.is_finite() or not value>=0:
        raise ValueError('finite nonnegative local bound required')
    return [int(x) for x in value.upper().man_exp()]


def unpack(pair):
    if (not isinstance(pair,list) or len(pair)!=2
            or any(type(x) is not int for x in pair) or pair[0]<0):
        raise ValueError('invalid local dyadic endpoint')
    value=arb(pair[0])*arb(2)**pair[1]
    if not value.is_exact() or not value.is_finite():
        raise ValueError('precision does not preserve the cached endpoint')
    return value


def save(directory,key,family):
    base,coefficients=family
    if (len(base)!=33 or any(t.nrows()!=11 or t.ncols()!=11 for t in base)
            or set(coefficients)!=set(range(1,key['maximum']+1))
            or any(len(row)!=3 for row in coefficients.values())):
        raise ValueError('invalid local family dimensions')
    body=dict(key=key,base=[[[pack(t[i,j]) for j in range(11)] for i in range(11)] for t in base],
              coefficients=[[pack(x) for x in coefficients[j]] for j in range(1,key['maximum']+1)])
    record=dict(body=body,checksum=hashlib.sha256(operator_cache.encoded(body)).hexdigest())
    path=cache_path(directory,key)
    path.parent.mkdir(parents=True,exist_ok=True)
    # Only publish a complete record. A killed calculation leaves no partial hit.
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=path.stem+'-',suffix='.tmp',delete=False) as stream:
        temporary=Path(stream.name)
        stream.write(operator_cache.encoded(record))
    try:
        os.replace(temporary,path)
    finally:
        if temporary.exists():
            temporary.unlink()


def load(directory,key):
    path=cache_path(directory,key)
    if not path.exists():
        return None
    if ctx.prec<key['precision']:
        raise ValueError('cache precision exceeds active precision')
    record=json.loads(path.read_bytes()); body=record['body']
    if body['key']!=key or hashlib.sha256(operator_cache.encoded(body)).hexdigest()!=record['checksum']:
        raise ValueError('local family key/checksum mismatch')
    rows=body['base']; coefficients=body['coefficients']
    if (len(rows)!=33 or any(len(t)!=11 or any(len(row)!=11 for row in t) for t in rows)
            or len(coefficients)!=key['maximum'] or any(len(row)!=3 for row in coefficients)):
        raise ValueError('cached local family dimensions mismatch')
    return ([arb_mat([[unpack(x) for x in row] for row in t]) for t in rows],
            {j+1:tuple(unpack(x) for x in row) for j,row in enumerate(coefficients)})


def build(tilts,penalty,precision=192,maximum=16,*,exact_feedback=0,
          spectral_feedback=0,joint_four=False,density_through=0,directory=None):
    """Return tilt -> (unsplit base, coupled mass coefficients).

All refinement steps precede any coupled-column replacement. A caller can
therefore keep a new scalar density bound or choose a fixed convex mixture.
"""
    tilts=tuple(dict.fromkeys(map(str,tilts)))
    if not tilts:
        raise ValueError('at least one tilt required')
    maximum=max(maximum,exact_feedback,spectral_feedback)
    sources=source_digest() if directory is not None else None
    keys={tilt:parameters(tilt,penalty,precision,maximum,exact_feedback,spectral_feedback,
                         joint_four,density_through,sources) for tilt in tilts}
    ctx.prec=precision
    result={}
    if directory is not None:
        for tilt,key in keys.items():
            cached=load(directory,key)
            if cached is not None:
                result[tilt]=cached
                print('LOCAL FAMILY exact cache hit',tilt,flush=True)
    missing=[tilt for tilt in tilts if tilt not in result]
    if not missing:
        return result
    grid,census=epoch_grid(missing,penalty,precision)
    through=6
    if exact_feedback or spectral_feedback:
        from full_feedback_refinement import refine
        if exact_feedback:
            from exact_feedback import build as feedback_build
            census=feedback_build(exact_feedback,census)
        else:
            from spectral_feedback import build as feedback_build,combined_records
            census=combined_records(feedback_build(spectral_feedback,census),census)
        spectrum=baseline.prepare_inputs()[3]
        grid={tilt:refine(base,census,spectrum,tilt,penalty,2,'1') for tilt,base in grid.items()}
        through=exact_feedback or spectral_feedback
    if joint_four:
        from joint_four import census as joint_census
        from cancellation_joint import check_feedback,refine
        joint=joint_census(); check_feedback(joint,census)
        grid={tilt:refine(base,joint,tilt,penalty,2,'1') for tilt,base in grid.items()}
    if density_through:
        from density_extend import build as density_build,refine
        records=density_build(density_through,missing,precision=precision)
        spectrum=baseline.prepare_inputs()[3]
        grid={tilt:refine(base,records[tilt],spectrum,penalty,maximum=density_through)
              for tilt,base in grid.items()}
    for tilt,base in grid.items():
        result[tilt]=base,conditioned_coefficients(census,through,maximum,tilt,penalty)
    if directory is not None:
        if source_digest()!=sources:
            print('LOCAL FAMILY sources changed during build; memo not written',flush=True)
        else:
            for tilt in missing:
                save(directory,keys[tilt],result[tilt])
                print('LOCAL FAMILY saved exact endpoints',tilt,flush=True)
    return result
