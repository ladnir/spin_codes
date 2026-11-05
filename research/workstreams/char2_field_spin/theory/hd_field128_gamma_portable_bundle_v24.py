"""Package exact proposals and frozen helpers; no numerical evaluation."""
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import shutil
import zipfile

import hd_field128_gamma_portable_bundle_v23 as previous
import hd_field128_gamma_midpoint_weight_band_typed_v2_audit as audit

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
LIMIT = 5 << 20
RESERVE = 1 << 30
ORIGINS = (
    (152, 'edge', '4f2d236edf57f68119381141207b4390dc8516404487c95607b0828c5fea53bb'),
    (192, 'geometric', '30e4f828deb33c41c2da2be7597c93589731fb4770926c8fa0f4204ea98c90ba'),
    (224, 'edge', '671f1563baad138785fb8175307edac9200a6c5a89d7b07ef8fce684e9f99955'),
    (288, 'edge', '5f42ecffde8e1cf34594b7770fc504b674415884ec3e97f43347f5a1dc8fef5f'),
    (376, 'edge', '0637f97c3af6df508386c3e99ba87187ff87042124729c91e56f27887ed0cc5b'),
)


def main():
    pins = audit.source_pins()
    for name in (Path(__file__).name, Path(previous.__file__).name,
                 'hd_field128_gamma_portable_band_runner_v24.py',
                 'hd_field128_gamma_portable_band_runner_v24_tests.py'):
        pins[name] = sha256((HERE/name).read_bytes()).hexdigest()
    payload, inputs, seeds, pending = {}, {}, [], []
    for name, digest in sorted(pins.items()):
        path = (HERE/name).resolve()
        raw = path.read_bytes()
        if sha256(raw).hexdigest() != digest:
            raise ValueError('Helper changed: '+name)
        payload[path.relative_to(ROOT).as_posix()] = raw
        inputs[path] = raw
    for weight, kind, digest in ORIGINS:
        path = HERE/f'hd_field128_gamma30_d40_r{weight}_{kind}_pruned_v22_20261010.json'
        raw, seed_raw = previous.portable_seed(path, digest, weight)
        inputs[path] = raw
        name = f'hd_field128_gamma30_d40_r{weight}_portable_v24_20261010.json'
        payload['theory/'+name] = seed_raw
        pending.append((HERE/name, seed_raw))
        seeds.append(dict(active_groups=weight,
            origin=dict(path=str(path), sha256=digest),
            seed=dict(path=name, sha256=sha256(seed_raw).hexdigest()),
            exact_literals_equal=True, numerical_bounds_transferred=False))
    manifest = dict(schema='gamma-portable-transport-v24',
        construction=audit.CONSTRUCTION, exact_literal_origins=seeds,
        files=[dict(member=name, bytes=len(raw), sha256=sha256(raw).hexdigest())
               for name, raw in sorted(payload.items())],
        source_sha256=dict(sorted(pins.items())),
        maximum_total_local_write_bytes=LIMIT,
        numerical_evaluation=False, whole_code_certificate=False)
    manifest_raw = previous.canonical(manifest)
    payload['transport-v24.json'] = manifest_raw
    if sum(len(raw)+1024 for raw in payload.values()) > LIMIT:
        raise OSError('Transport payload exceeds small packaging budget')
    buffer = BytesIO()
    with zipfile.ZipFile(buffer, 'w', compression=zipfile.ZIP_DEFLATED) as stream:
        for name, raw in sorted(payload.items()):
            stream.writestr(name, raw)
    archive_raw = buffer.getvalue()
    pending.extend((
        (HERE/'hd_field128_gamma_portable_transport_v24_20261010.json', manifest_raw),
        (HERE/'hd_field128_gamma_portable_transport_v24_20261010.zip', archive_raw),
    ))
    total = sum(len(raw) for path, raw in pending)
    if total > LIMIT or shutil.disk_usage(HERE).free-total < RESERVE:
        raise OSError('Bounded housekeeping write preflight failed')
    if any(path.exists() for path, raw in pending):
        raise FileExistsError('Refusing to replace transport or seed artifacts')
    if any(path.read_bytes() != raw for path, raw in inputs.items()):
        raise RuntimeError('A transport input changed')
    for path, raw in pending:
        if shutil.disk_usage(HERE).free-len(raw) < RESERVE:
            raise OSError('Housekeeping write would breach reserve')
        with path.open('xb') as stream:
            stream.write(raw)
    if any(path.read_bytes() != raw for path, raw in inputs.items()):
        raise RuntimeError('A transport input changed while writing')
    with zipfile.ZipFile(BytesIO(archive_raw)) as stream:
        if set(stream.namelist()) != set(payload) or any(
                stream.read(name) != raw for name, raw in payload.items()):
            raise RuntimeError('Archive membership or bytes changed')
    print(json.dumps(dict(archive=str(pending[-1][0]),
        archive_sha256=sha256(archive_raw).hexdigest(),
        archive_bytes=len(archive_raw), total_local_write_bytes=total,
        manifest_sha256=sha256(manifest_raw).hexdigest(), seeds=seeds)), flush=True)


if __name__ == '__main__':
    main()
