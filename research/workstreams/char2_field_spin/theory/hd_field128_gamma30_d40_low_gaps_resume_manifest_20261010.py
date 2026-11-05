"""Preserved low-gap proposals and non-admitted attempts after a resource stop."""

from hashlib import sha256
import json
from pathlib import Path
from shutil import disk_usage


def main():
    here, source = Path(__file__).resolve().parent, Path(__file__).resolve()
    source_hash = sha256(source.read_bytes()).hexdigest()
    inputs = {}
    def load(relative):
        path = here/relative
        raw = path.read_bytes()
        inputs[path] = raw
        return json.loads(raw), dict(path=str(path), sha256=sha256(raw).hexdigest())
    candidates = []
    for weight, name, expected_proposal, reviewed in (
            (960, 'r960_gap_refit', '80e0f727716f6c0044f81db404c7e80b121340d1f37c50d5ae3a4edfbab424b4', True),
            (1168, 'r1168_gap_calibration', 'fd0ac0024464d8d2f852dde4bb3bdc5566eac2bcbd1ab5f0c5577876da9bff81', True),
            (913, 'r913_gap_mark_search', 'fd99500cdaca668b4b59cc17f742f1c567ca1bd4376ad519fb04f9557239aa80', False)):
        folder = f'hd_field128_gamma30_d40_{name}_20261010'
        summary, summary_pin = load(folder+'/summary.json')
        proposal, proposal_pin = load(folder+'/proposal.json')
        replay, replay_pin = load(folder+'/replay.json')
        if (proposal_pin['sha256'] != expected_proposal or proposal['active_groups'] != weight
                or proposal['state_dimension'] != 30 or proposal['target_distance'] != '2/5'
                or summary['exact_replay_match'] is not True
                or replay['last_complete_step']['profile_pairs_evaluated'] != (weight//2+1)**2
                or replay['last_complete_step']['log2_upper_screen'] != summary['complete_score']):
            raise ValueError('Unexpected complete FLOAT proposal/replay linkage')
        audit_pin = None
        if reviewed:
            _, audit_pin = load(folder+'/independent_consistency_audit.json')
        candidates.append(dict(calibration_weight=weight, proposal=proposal_pin, summary=summary_pin,
            fresh_replay=replay_pin, independent_saved_lineage_audit=audit_pin,
            marker_count=len(proposal['witness']['markers']), complete_FLOAT_score=summary['complete_score'],
            directed_certificate=False, source_role='Literal proposal only; no numerical bound transfers'))
    partials = []
    for relative, expected in (
            ('hd_field128_gamma30_d40_r1168_gap_directed256_20261010/band_1145_1191/summary.json',
             'c105375d56ff3ff09f17a96acbbbd46c3e219561e5329d7dbf87bf015966dd61'),
            ('hd_field128_gamma30_d40_r1168_gap_directed256_20261010/retry2048/summary.json',
             '0175455552085790c8bac4121e00641ce5315848a8e158275b7d7580f59fa847')):
        row, pin = load(relative)
        if (pin['sha256'] != expected or row['complete_requested_band'] is not False
                or row['log2_upper_rational'] is not None or row['status'] != 'partial'
                or row['requested_weight_range'] != [1145, 1191]):
            raise ValueError('Unexpected preserved incomplete attempt')
        partials.append(dict(summary=pin, range=[1145, 1191], reason=row['stop_reason'],
            evaluated_nodes=row['evaluated_nodes'], tree_nodes=row['tree_nodes'], admitted=False))
    stop_root = here/'hd_field128_gamma30_d40_r1168_gap_subbands512_directed256_20261010'
    _, stop_bank = load(stop_root.name+'/bank.json')
    if list(stop_root.rglob('summary.json')):
        raise ValueError('Resource-stopped attempt unexpectedly has a summary')
    result = dict(schema='gamma-low-gaps-resource-stop-handoff-v1',
        assigned_ranges=[[913, 999], [1145, 1191]], pending_ranges=[[913, 999], [1145, 1191]],
        admitted_complete_events=[], admitted_weight_count=0, literal_candidates=candidates,
        excluded_partial_events=partials,
        resource_stopped_attempt=dict(path=str(stop_root), requested_first_range=[1176, 1191],
            bank=stop_bank, summary_present=False, admitted=False,
            observation='750MiB guard raised before checkpoint; failure also could not be written'),
        prepared_unrun_step=dict(path=str(here/'hd_field128_gamma30_d40_r913_gap_refit_20261010.py'),
            sha256=sha256((here/'hd_field128_gamma30_d40_r913_gap_refit_20261010.py').read_bytes()).hexdigest(),
            maximum_additional_fits=15, executed=False),
        runtime=dict(owned_bank_sessions_exited=[14874], owned_live_banks=[], all_numerical_processes_finished=True),
        source_sha256={source.name: source_hash}, whole_code_certificate=False,
        joint_probability_claim=False, no_weight_interpolation=True,
        scope='Pinned resumption record only; no admitted weight or aggregate probability claim')
    if sha256(source.read_bytes()).hexdigest() != source_hash or any(path.read_bytes() != raw for path, raw in inputs.items()):
        raise RuntimeError('Handoff input bytes changed')
    if disk_usage(here).free < 750 * 1024**2:
        raise OSError('Artifact guard: less than750MiB disk space remains')
    output = here/'hd_field128_gamma30_d40_low_gaps_resume_manifest_20261010.json'
    with output.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(dict(path=str(output), sha256=sha256(output.read_bytes()).hexdigest())), flush=True)


if __name__ == '__main__':
    main()
