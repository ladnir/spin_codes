"""Create a small byte-pinned transport bundle; never generate numerical caps."""
from fractions import Fraction as F
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import shutil
import zipfile

import hd_field128_gamma_midpoint_weight_band_typed_v2_audit as audit

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESERVE = 600 << 20
MAX_ARCHIVE = 5 << 20
ORIGINS = {
    68: 'b8bc16786c138030739ce024db46a35d937fcd01af124052d36de117028ee01a',
    96: '101691c609bc6f8c2dbe64c91059ae618c180dad3e7b65e7bc4077507116fcba',
}


def canonical(value):
    return (json.dumps(value, indent=2, allow_nan=False) + '\n').encode()


def write_new(path, raw):
    if len(raw) > MAX_ARCHIVE or shutil.disk_usage(HERE).free - len(raw) < RESERVE:
        raise OSError('Small transport output would breach600MiB local reserve')
    with path.open('xb') as stream:
        stream.write(raw)


def portable_seed(origin, digest, weight):
    raw = origin.read_bytes()
    if sha256(raw).hexdigest() != digest:
        raise ValueError('Original exact witness byte pin differs')
    row = json.loads(raw)
    if (row['active_groups'] != weight or row['state_dimension'] != 30
            or row['target_distance'] != '2/5' or row['construction'] != audit.CONSTRUCTION):
        raise ValueError('Wrong construction/weight in original witness')
    witness = {key: row['witness'][key] for key in ('z', 'alpha', 'markers')}
    z, alpha = F(witness['z']), F(witness['alpha'])
    if (not 0 < z < 1 or not 0 < alpha <= 1
            or not 1 <= len(witness['markers']) <= 128
            or any(len(pair) != 2 or any(F(x) <= 0 for x in pair)
                   for pair in witness['markers'])):
        raise ValueError('Invalid positive literal witness')
    witness.update(grid=128, state_dimension=30)
    seed = dict(schema='gamma-portable-exact-witness-v1', active_groups=weight,
        construction=audit.CONSTRUCTION, field_order=1 << 128, state_dimension=30, target_distance='2/5',
        witness=witness, source=dict(path=str(origin.resolve()), sha256=digest),
        origin_role='Exact literal proposal only; rebuild all numerical coefficients and bounds',
        directed=False, selected_weight_certificate=False, whole_code_certificate=False)
    return raw, canonical(seed)


def main():
    pins = audit.source_pins()
    for name in ('hd_field128_gamma_portable_band_runner.py',
                 'hd_field128_gamma_portable_band_runner_tests.py', Path(__file__).name):
        pins[name] = sha256((HERE / name).read_bytes()).hexdigest()
    payload, inputs, origins = {}, {}, []
    for name, digest in sorted(pins.items()):
        path = (HERE / name).resolve()
        member = path.relative_to(ROOT).as_posix()
        raw = path.read_bytes()
        if sha256(raw).hexdigest() != digest:
            raise ValueError('Helper bytes changed: ' + name)
        payload[member], inputs[path] = raw, raw
    pending_seeds = []
    for weight, digest in ORIGINS.items():
        origin = HERE / f'hd_field128_gamma30_d40_r{weight}_edge_pruned_v22_20261010.json'
        raw, seed_raw = portable_seed(origin, digest, weight)
        inputs[origin] = raw
        seed_name = f'hd_field128_gamma30_d40_r{weight}_portable_v23_20261010.json'
        payload['theory/' + seed_name] = seed_raw
        pending_seeds.append((HERE / seed_name, seed_raw))
        origins.append(dict(active_groups=weight, origin=dict(path=str(origin), sha256=digest),
            seed=dict(path=seed_name, sha256=sha256(seed_raw).hexdigest()),
            exact_literals_equal=True, numerical_bounds_transferred=False))
    manifest = dict(schema='gamma-portable-transport-v1', construction=audit.CONSTRUCTION,
        files=[dict(member=name, sha256=sha256(raw).hexdigest(), bytes=len(raw))
               for name, raw in sorted(payload.items())],
        exact_literal_origins=origins,
        source_sha256=dict(sorted(pins.items())),
        local_free_reserve_bytes=RESERVE, maximum_archive_bytes=MAX_ARCHIVE,
        numerical_free_reserve_bytes=750 << 20,
        audit_location='Full typed-v2 provenance audit after byte-identical return to Windows',
        numerical_cap_reuse=False, whole_code_certificate=False)
    manifest_raw = canonical(manifest)
    manifest_path = HERE / 'hd_field128_gamma_portable_transport_v23_20261010.json'
    payload['transport.json'] = manifest_raw
    archive = HERE / 'hd_field128_gamma_portable_transport_v23_20261010.zip'
    outputs = [path for path, unused in pending_seeds] + [manifest_path, archive]
    if any(path.exists() for path in outputs):
        raise FileExistsError('Refusing to replace a transport artifact')
    if any(path.read_bytes() != raw for path, raw in inputs.items()):
        raise RuntimeError('Transport source changed')
    # Bound packaging memory first, then preflight all exact output sizes.
    if sum(len(raw) + 1024 for raw in payload.values()) > MAX_ARCHIVE:
        raise OSError('Uncompressed transport payload exceeds bounded packaging limit')
    packed = BytesIO()
    with zipfile.ZipFile(packed, 'w', compression=zipfile.ZIP_DEFLATED) as stream:
        for name, raw in sorted(payload.items()):
            stream.writestr(name, raw)
    archive_raw = packed.getvalue()
    total = len(archive_raw) + len(manifest_raw) + sum(len(raw) for path, raw in pending_seeds)
    if len(archive_raw) > MAX_ARCHIVE or shutil.disk_usage(HERE).free - total < RESERVE:
        raise OSError('Bounded bundle would breach local reserve')
    for path, raw in pending_seeds:
        write_new(path, raw)
    write_new(manifest_path, manifest_raw)
    write_new(archive, archive_raw)
    if any(path.read_bytes() != raw for path, raw in inputs.items()):
        raise RuntimeError('Transport source changed during packaging')
    with zipfile.ZipFile(archive) as stream:
        if set(stream.namelist()) != set(payload):
            raise RuntimeError('Archive membership differs')
        if any(stream.read(name) != raw for name, raw in payload.items()):
            raise RuntimeError('Archive changed file contents')
    print(json.dumps(dict(archive=str(archive), bytes=archive.stat().st_size,
        sha256=sha256(archive.read_bytes()).hexdigest(), file_count=len(payload),
        manifest_sha256=sha256(manifest_raw).hexdigest(), seeds=origins)), flush=True)


if __name__ == '__main__':
    main()
