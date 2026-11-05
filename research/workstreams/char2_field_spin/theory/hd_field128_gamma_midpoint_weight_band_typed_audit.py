"""Typed saved-band provenance adapter; frozen count/tree/union arithmetic.

JSON provenance objects are parsed recursively. Explicit Python and Markdown
links are checked as source bytes, never executed. This is not authentication
or a coefficient, regional-matrix, or cell-transfer replay.
"""
import ast
from hashlib import sha256
import json
import os
from pathlib import Path

import hd_field128_gamma_midpoint_weight_band_audit as previous

HERE = Path(__file__).resolve().parent
PREVIOUS_SHA256 = 'b3099b35281ad95597de024912a7e52356b1df96b16c9af08c3c11a52d6cf94f'
producer, base = previous.producer, previous.base
require, canonical, CONSTRUCTION = previous.require, previous.canonical, previous.CONSTRUCTION


def source_pins():
    require(sha256(Path(previous.__file__).read_bytes()).hexdigest() == PREVIOUS_SHA256,
            'Frozen numerical auditor changed')
    pins = previous.source_pins()
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def check_source_pins(pins):
    previous.check_source_pins(pins)


class Inputs(previous.Inputs):
    def __init__(self):
        super().__init__()
        self.non_json_pins = {}

    def load(self, path, expected=None):
        path = Path(path).resolve()
        require(path.suffix.lower() == '.json', 'JSON input must have .json extension')
        value, pin = super().load(path, expected)
        require(isinstance(value, (dict, list)), 'JSON provenance must be an object or list')
        return value, pin

    def load_source(self, path, expected):
        path = Path(path).resolve()
        suffix = path.suffix.lower()
        require(suffix in ('.py', '.md'), 'Unsupported non-JSON provenance extension')
        raw = path.read_bytes()
        digest = sha256(raw).hexdigest()
        require(digest == expected, 'Non-JSON source byte pin mismatch')
        text = raw.decode('utf-8-sig')
        try:
            json.loads(text)
        except json.JSONDecodeError:
            pass
        else:
            raise ValueError('JSON content with a non-JSON source extension')
        if suffix == '.py':
            ast.parse(text, filename=str(path))
        pin = dict(path=str(path), sha256=digest, kind='python-source' if suffix == '.py' else 'markdown-source')
        prior = self.non_json_pins.get(str(path))
        require(prior is None or prior == pin, 'Non-JSON input changed during audit')
        self.non_json_pins[str(path)] = pin
        return pin

    def check(self):
        super().check()
        for path, pin in self.non_json_pins.items():
            require(sha256(Path(path).read_bytes()).hexdigest() == pin['sha256'],
                    'Pinned non-JSON input changed: '+path)


def _source(sources, name, digest):
    require(name not in sources or sources[name] == digest, 'Conflicting helper source pin: '+name)
    require(sha256((HERE/name).read_bytes()).hexdigest() == digest, 'Helper source changed: '+name)
    sources[name] = digest


def _provenance(value, inputs, sources):
    if isinstance(value, dict):
        if isinstance(value.get('path'), str) and isinstance(value.get('sha256'), str):
            path = Path(value['path'])
            path = (path if path.is_absolute() else HERE/path).resolve()
            if path.suffix.lower() == '.json':
                if str(path) not in inputs.pins:
                    linked, unused = inputs.load(path, value['sha256'])
                    _provenance(linked, inputs, sources)
                else:
                    require(inputs.pins[str(path)]['sha256'] == value['sha256'],
                            'Conflicting JSON provenance pins')
            else:
                pin = inputs.load_source(path, value['sha256'])
                _source(sources, Path(os.path.relpath(path, HERE)).as_posix(), pin['sha256'])
        for name, child in value.items():
            if (name == 'source_sha256' or name.endswith('_source_sha256')) and isinstance(child, dict):
                for source, digest in child.items():
                    _source(sources, source, digest)
            _provenance(child, inputs, sources)
    elif isinstance(value, list):
        for child in value:
            _provenance(child, inputs, sources)


def audit_summary(path, *, max_seconds=120.):
    """Audit a complete saved band plus its context and all node checkpoints."""
    budget, pins, inputs = base.Budget(max_seconds), source_pins(), Inputs()
    check_source_pins(pins)
    record, summary_pin = inputs.load(path)
    previous._metadata(record)
    context = record['context']
    sources = dict(pins)
    _provenance(context, inputs, sources)
    for key, value in record.items():
        if key not in ('context', 'nodes'):
            _provenance({key: value}, inputs, sources)
    provenance = context.get('witness_source')
    require(isinstance(provenance, dict) and {'path', 'sha256'} <= provenance.keys(), 'Missing exact witness provenance')
    proposal_path = Path(provenance['path'])
    proposal, unused = inputs.load(proposal_path if proposal_path.is_absolute() else HERE/proposal_path, provenance['sha256'])
    for key in ('z', 'alpha', 'markers'):
        require(proposal['witness'][key] == context['witness'][key], 'Proposal/bank exact literal mismatch: '+key)
    directory = Path(summary_pin['path']).parent
    saved_context, context_pin = inputs.load(directory/'context.json')
    require(saved_context == dict(context=context, context_sha256=record['context_sha256'], directed=False,
                                  complete_requested_band=False, whole_code_certificate=False), 'Context artifact mismatch')
    checkpoint_count = record['checkpoints']
    require(type(checkpoint_count) is int and checkpoint_count == len(record['nodes'])+1, 'Missing node checkpoints')
    replay = {}
    for sequence in range(checkpoint_count):
        budget.check()
        checkpoint, unused = inputs.load(directory/f'node_{sequence:06d}.json')
        require(checkpoint['sequence'] == sequence and checkpoint['context_sha256'] == record['context_sha256'], 'Checkpoint linkage mismatch')
        for key in ('complete_requested_band', 'covers_all_supports_in_band', 'directed', 'batch_source_verified_at_exit', 'whole_code_certificate'):
            require(checkpoint[key] is False, 'Checkpoint overclaims completion')
        require(checkpoint['log2_upper_rational'] is None, 'Partial checkpoint has band bound')
        for row in checkpoint['changed_nodes']:
            replay[row['id']] = row
    require([replay[i] for i in sorted(replay)] == record['nodes'], 'Checkpoint replay differs from final tree')
    extra_path = directory/'driver.json'
    execution_pin = None
    if extra_path.exists():
        execution, execution_pin = inputs.load(extra_path)
        require(execution['requested_weight_range'] == record['requested_weight_range'],
                'Execution driver range mismatch')
        require(execution['witness_source'] == provenance, 'Execution witness source mismatch')
        require(execution['cell_failure_bits'] == record['cell_failure_bits']
                and execution['aggregate_failure_bits'] == record['requested_failure_bits'],
                'Execution budget mismatch')
        require(execution['numerical_provider'] == Path(producer.__file__).name,
                'Execution numerical provider mismatch')
        _provenance(execution, inputs, sources)
    numeric = previous.audit_numeric(record, 4096, budget)
    check_source_pins(sources)
    inputs.check()
    return dict(schema='gamma-midpoint-weight-band-typed-independent-audit-v1', audit_result='PASS',
        input_summary=summary_pin, context_file=context_pin,
        execution_driver_receipt=execution_pin, context_sha256=record['context_sha256'],
        construction=CONSTRUCTION, state_dimension=30, field_order=1 << 128, target_distance='2/5',
        requested_weight_range=record['requested_weight_range'], requested_weights=numeric['weights'],
        complete_requested_band=True, covers_all_matching_profiles=True, integer_partition_checked=True,
        every_original_weight_preparation_metadata_checked=True, normalized_histograms=numeric['normalized_histograms'],
        counts_recomputed=True, count_algorithm='independent five-prefix square/rectangle formula with exact array recurrence; one original weight at a time',
        cells_checked=numeric['cells_checked'], accepted_cells=numeric['accepted_cells'], checkpoints_checked=checkpoint_count,
        cell_log_enclosures_checked=True, positive_count_zero_cells=numeric['positive_count_zero_cells'],
        directed_union_precision=384, saved_band_log2_upper_rational=record['log2_upper_rational'],
        reaggregated_cell_log2_upper_rational=numeric['union'], reaggregated_cell_log2_upper_enclosure=numeric['union_enclosure'],
        log2_upper_rational=record['log2_upper_rational'], zero=numeric['zero'],
        source_sha256=pins, checked_source_sha256=dict(sorted(sources.items())), checked_file_pins=list(inputs.pins.values()),
        checked_non_json_file_pins=list(inputs.non_json_pins.values()),
        provenance_type_policy='JSON objects/lists recurse; UTF-8 Python/Markdown helper bytes are hashed without executing',
        typed_provenance_checked=True, frozen_numerical_auditor_sha256=PREVIOUS_SHA256,
        source_pin_conflicts=False, source_and_input_bytes_checked_at_exit=True,
        coefficient_replay_in_this_audit=False, normalized_regional_replay=False, cell_transfer_recomputed=False,
        trusted_numerical_producer='reviewed and byte-pinned complete band producer; stored cell caps are not authenticated by this audit',
        no_weight_interpolation=True, whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('summary', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-seconds', type=float, default=120.)
    args = parser.parse_args()
    result = audit_summary(args.summary, max_seconds=args.max_seconds)
    print(json.dumps(producer.profiles.write_new(args.output, result)))

