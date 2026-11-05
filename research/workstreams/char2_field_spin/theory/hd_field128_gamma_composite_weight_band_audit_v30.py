"""Composite saved-data audit; no coefficient replay."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
import math
from pathlib import Path
import shutil

from flint import arb, ctx
import hd_field128_gamma_composite_weight_band_v30 as producer
import hd_field128_gamma_midpoint_weight_band_audit as numeric
import hd_field128_gamma_literal_midpoint_audit_v30 as literal

HERE, NAME = Path(__file__).resolve().parent, Path(__file__).name
PRODUCER_SHA = 'b37d1e3133a42675e911f8976077d2a2c061a6e1653d7cca038fa2ac9187aa87'
NUMERIC_SHA = 'b3099b35281ad95597de024912a7e52356b1df96b16c9af08c3c11a52d6cf94f'
LITERAL_SHA = '75f57e3ac55f247edf9a9e6e1fdee884cbb4a4a770bc46d77e21cd7c79df8e21'
SCHEMA = 'gamma-composite-weight-band-literal-only-remote-audit-v1'
MAX_OUTPUT_BYTES, RESERVE = 1 << 20, 20 << 30
MINIMUM_STAGE = 'complete frozen cell bounds, including each own outer constant and original-r scaling'
require, same, canonical = literal.require, literal.same, numeric.canonical


def source_pins(controller):
    for module, expected in ((producer, PRODUCER_SHA), (numeric, NUMERIC_SHA), (literal, LITERAL_SHA)):
        require(sha256(Path(module.__file__).read_bytes()).hexdigest() == expected, 'Reviewed helper changed')
    result = {}
    for pins in (producer.source_pins(controller), numeric.source_pins(), literal.source_pins()):
        for name, digest in pins.items():
            require(name not in result or result[name] == digest, 'Conflicting source pin')
            result[name] = digest
    result[NAME] = sha256((HERE/NAME).read_bytes()).hexdigest()
    return dict(sorted(result.items()))


class Inputs:
    def __init__(self):
        self.pins = {}

    def load(self, name, expected=None, *, limit=32 << 20):
        name = producer.selected.literals.proposal_name(name)
        path = HERE/name
        require(path.resolve().is_relative_to(HERE.resolve()) and path.stat().st_size <= limit,
                'Escaping or oversized input')
        raw = path.read_bytes()
        digest = sha256(raw).hexdigest()
        require(expected is None or expected == digest, 'Input byte pin mismatch')
        require(name not in self.pins or self.pins[name]['sha256'] == digest, 'Input changed')
        pin = dict(path=name, sha256=digest, bytes=len(raw))
        self.pins[name] = pin
        require(sum(p['bytes'] for p in self.pins.values()) <= 36 << 20, 'Input limit exceeded')
        value = json.loads(raw)
        require(isinstance(value, dict), 'Expected a JSON object')
        return value, pin

    def linked(self, pin, *, limit=32 << 20):
        require(isinstance(pin, dict) and set(pin) == {'path', 'sha256'}, 'Malformed input pin')
        digest(pin['sha256'])
        return self.load(pin['path'], pin['sha256'], limit=limit)

    def check(self):
        for pin in self.pins.values():
            require(sha256((HERE/pin['path']).read_bytes()).hexdigest() == pin['sha256'], 'Input changed at exit')


def digest(value):
    require(isinstance(value, str) and len(value) == 64 and all(c in '0123456789abcdef' for c in value),
            'Expected lowercase SHA256')
    return value


def seconds(value, maximum):
    require(type(value) in (int, float) and math.isfinite(value) and 0 <= value <= maximum,
            'Invalid reported time limit')


def log_bound(row):
    require(type(row.get('zero')) is bool, 'Missing explicit zero flag')
    if row['zero']:
        require(row['log2_upper_rational'] is None and row['log2_upper_enclosure'] is None,
                'Zero component has a logarithm')
        return None
    text = row['log2_upper_rational']
    require(isinstance(text, str) and str(F(text)) == text, 'Noncanonical ceiling')
    numeric._enclosure(row['log2_upper_enclosure'], text, 'Saved component ceiling')
    return F(text)


def component_minimum(receipt, count):
    parts = receipt['component_bounds']
    require(isinstance(parts, list) and len(parts) == count, 'Missing component bounds')
    values = []
    for i, part in enumerate(parts):
        same(part, dict(witness_index=i), 'ordered component')
        require(set(part) == {'witness_index', 'zero', 'log2_upper_rational', 'log2_upper_enclosure'}, 'Component fields')
        values.append(log_bound(part))
    winner = min(range(count), key=lambda i: (values[i] is not None, F(0) if values[i] is None else values[i]))
    same(receipt, dict(component_winner=winner, minimum_stage=MINIMUM_STAGE,
                      entrywise_cross_witness_minimum=False), 'complete-cell minimum')
    same(receipt, {key:parts[winner][key] for key in ('zero', 'log2_upper_rational', 'log2_upper_enclosure')},
         'chosen complete bound')
    return values[winner]


def audit_numeric(tree, count, budget, *, pairs=4096):
    first, last = tree['requested_weight_range']
    weights = list(numeric.producer.checked_weights(pairs, first, last))
    context = tree['context']
    digest_value = sha256(canonical(context)).hexdigest()
    require(tree['context_sha256'] == digest_value, 'Context digest mismatch')
    roots = numeric.expected_roots(pairs, weights)
    require(tree['requested_weights'] == context['requested_weights'] == weights,
            'Omitted weight')
    require(context['requested_weight_range'] == [first, last]
        and context['root_cells'] == tree['root_cells'] == roots, 'Incomplete original parity domain')
    for key in ('complete_requested_band', 'covers_all_supports_in_band', 'covers_all_matching_profiles',
                'every_integer_weight_prepared', 'no_weight_interpolation'):
        require(tree.get(key) is True, 'Incomplete generic tree: '+key)
    require(tree['status'] == 'complete' and tree['stop_reason'] is None, 'Partial tree')
    by_id, leaves = numeric._tree(tree, roots)
    counts = dict.fromkeys(by_id, 0)
    for r in weights:
        budget.check()
        table, u, eps = numeric.RectangleCounts(pairs, r, budget), r//2, r % 2
        limit = min(u, pairs-u-eps)
        for identity, node in by_id.items():
            (ul, uh), (dl, dh), (el, eh) = node['cell']['u_deficit_ranges']
            if node['cell']['parity'] == eps and ul <= u <= uh and dl <= limit and el <= limit:
                counts[identity] += table.rectangle(u-dh, u-dl, u-eh, u-el)
    cell_bits, target = tree['cell_failure_bits'], tree['requested_failure_bits']
    require(type(cell_bits) is int and type(target) is int and 40 <= target <= cell_bits <= 4096,
            'Wrong cell/event target')
    winners, leaf_ids = [0]*count, set(leaves)
    with ctx.workprec(384):
        for identity, node in by_id.items():
            budget.check()
            receipt = node['receipt']
            same(receipt, dict(cell=node['cell'], context_sha256=digest_value, complete_cell=True,
                status='complete', count_upper=counts[identity], count_method=numeric.COUNT_METHOD,
                failure_bits=cell_bits), 'complete shared-count cell')
            bound = component_minimum(receipt, count)
            require(bound is None or counts[identity] > 0, 'Nonzero bound for a zero-count cell')
            require(counts[identity] or all(part['zero'] for part in receipt['component_bounds']),
                    'Zero common count has a nonzero component')
            passing = bound is None or bound <= -cell_bits
            require(receipt['meets_cell_budget'] is passing and passing == (node['status'] == 'accepted'),
                    'Cell status mismatch')
            if identity in leaf_ids:winners[receipt['component_winner']] += 1
        total = sum((arb(2)**numeric.base.rational(F(by_id[i]['receipt']['log2_upper_rational']))
                     for i in leaves if not by_id[i]['receipt']['zero']), arb(0))
        zero = total.is_zero()
        bound = log_bound(tree)
        require(zero == (bound is None), 'Band/leaf zero mismatch')
        value = None if zero else numeric.tail.log2arb(total)
        union = None if zero else numeric.base.rational_ceiling(value)
        require(zero or bound >= union, 'Understated384 union')
        admitted = zero or bound <= -target
        require(tree['meets_requested_band_target'] is admitted, 'Wrong aggregate admission')
    return dict(weights=weights, roots=roots, cells_checked=len(by_id), accepted_cells=len(leaves),
        normalized_histograms=sum(r//2-max(0, r-pairs)+1 for r in weights), zero=zero,
        union=None if zero else str(union), union_enclosure=None if zero else str(value),
        positive_count_zero_cells=sum(node['receipt']['zero'] and counts[i] > 0 for i,node in by_id.items()),
        accepted_witness_counts=winners, meets_requested_band_target=admitted)


def component_context(bank, context, proposal, first, last, controller):
    weights = list(range(first, last+1))
    roots = numeric.expected_roots(4096, weights)
    bank_pins=controller.cert.source_pins()
    same(bank,dict(source_sha256=bank_pins,bank_source_sha256=bank_pins),'bank source maps')
    same(bank, dict(exact_input_witness=proposal['witness'], proposal_source=proposal['source'],
        calibration_weight=proposal['calibration_weight'], requested_weight_range=[first,last], fresh_bank_rebuilt=True,
        literal_source_role='calibration literals only; no calibration bound inherited',
        directed=True, whole_code_certificate=False, scope=literal.BANK_SCOPE), 'fresh component bank')
    same(bank['witness'], dict(proposal['witness'], precision=256, packet_mantissa_bits=80,
        state_dimension=30, inner_width=64, gamma_extension_degree=16, gamma_base_field_order=256,
        joint_rank_return_cap=True, all_pair_before_alpha=True), 'exact component witness')
    seconds(context['preparation_seconds'], producer.COMPUTE_SECONDS)
    excluded = {'exact_input_witness','proposal_source','calibration_weight','requested_weight_range',
                'fresh_bank_rebuilt','literal_source_role'}
    expected = {key:value for key,value in bank.items() if key not in excluded}
    projected = ('schema pairs repetitions paired_region_power normalized_every_original_weight scaling max_domain '
        'count_method count_pointwise_minimum_used count_only_complement_above_pairs no_profile_pair_table '
        'no_weight_interpolation complete_requested_band directed whole_code_certificate').split()
    expected.update({key:numeric.EXPECTED[key] for key in projected})
    expected.update(source_sha256=controller.bands.source_pins(), bank_source_sha256=bank_pins,
        witness_source=proposal['source'], requested_weight_range=[first,last], requested_weights=weights,
        parity_weights={str(eps):[r for r in weights if r%2==eps] for eps in (0,1) if any(r%2==eps for r in weights)},
        regional_histograms_prepared=sum(r//2-max(0,r-4096)+1 for r in weights), root_cells=roots,
        preparation_seconds=context['preparation_seconds'])
    require(context == expected, 'Bank/context mismatch')
    same(context,numeric.EXPECTED,'component construction')


def audit_summary(name, expected_sha256, *, max_seconds=120.):
    require(type(max_seconds) in (int,float) and 0 < max_seconds <= 120, 'Audit budget must be at most120 seconds')
    digest(expected_sha256)
    controller = producer.helpers()
    pins, inputs, budget = source_pins(controller), Inputs(), numeric.base.Budget(max_seconds)
    numeric.producer.check_source_pins(pins)
    record, summary_pin = inputs.load(name, expected_sha256)
    directory = Path(name).parent.as_posix()
    try:
        tree = record['partition_result']
        first,last = record['requested_weight_range']
        producer.validate_request(first,last,tree['max_nodes'],tree['cell_failure_bits'],record['requested_failure_bits'])
        weights = list(range(first,last+1))
        execution = producer.source_pins(controller)
        same(record, dict(schema=producer.SCHEMA,construction=producer.selected.literals.CONSTRUCTION,
            field_order=1<<128,state_dimension=30,target_distance='2/5',requested_weights=weights,
            composite_source_sha256=execution,status='complete',complete_requested_band=True,covers_all_matching_profiles=True,
            directed=True,numerical_backend='fresh-actual-all-pair-banks-complete-cell-minimum',
            generic_partition_flags_preserved=True,whole_code_certificate=False,historical_origin_verified=False,
            old_numerical_bounds_reused=False,no_weight_interpolation=True,
            runner=dict(path=producer.NAME,sha256=PRODUCER_SHA),
            scope='Separate fresh composite-band producer; nested generic engine is not the production credential'), 'production envelope')
        same(tree, dict(schema='gamma-midpoint-weight-band-cover-v1',directed=False,numerical_backend='synthetic-test-only',
            whole_code_certificate=False,batch_source_verified_at_entry=True,batch_source_verified_at_exit=True,
            requested_weight_range=[first,last],requested_weights=weights,construction_changed=False,
            controller=dict(path=controller.NAME,sha256=producer.CONTROLLER_SHA),controller_source_sha256=controller.source_pins(),
            partition_policy=controller.POLICY,checkpoint_format=controller.CHECKPOINT_FORMAT), 'unchanged generic engine')
        for key in ('context_sha256','zero','log2_upper_rational','log2_upper_enclosure','requested_failure_bits',
                    'meets_requested_band_target'):
            require(type(record[key]) is type(tree[key]) and record[key] == tree[key], 'Envelope differs: '+key)
        seconds(record['total_compute_seconds'],producer.COMPUTE_SECONDS)
        seconds(tree['max_seconds'],producer.COMPUTE_SECONDS)
        seconds(tree['elapsed_seconds'],producer.COMPUTE_SECONDS)
        context = tree['context'];proposals=context['proposals'];components=context['components']
        require(isinstance(proposals,list) and 2 <= len(proposals) <= 4 and len(components)==len(proposals), 'Wrong component count')
        proposal_names=set()
        for proposal in proposals:
            budget.check()
            original, pin = inputs.linked(proposal['source'],limit=1<<20)
            require(pin['path'] not in proposal_names,'Duplicate proposal');proposal_names.add(pin['path'])
            calibration=producer.selected.literals.validate_target(original['active_groups'])
            witness,key=producer.selected.literals.project_literals(original,calibration)
            require(proposal == dict(source=proposal['source'],witness=witness,source_marker_key=key,calibration_weight=calibration,
                source_role='Exact literals only; every requested original weight is normalized afresh'), 'Exact literal projection mismatch')
        expected_context=dict(schema='gamma-composite-weight-band-context-v30',construction=producer.selected.literals.CONSTRUCTION,
            field_order=1<<128,state_dimension=30,target_distance='2/5',requested_weight_range=[first,last],requested_weights=weights,
            proposals=proposals,components=components,source_sha256=execution,pairs=4096,precision=256,repetitions=256,
            paired_region_power=128,common_integer_count_domain=True,original_weights_normalized_individually=True,
            component_constants_and_slopes_kept_separate=True,entrywise_cross_witness_minimum=False,
            historical_origin_verified=False,old_numerical_bounds_reused=False,root_cells=numeric.expected_roots(4096,weights),
            complete_requested_band=False,directed=False,whole_code_certificate=False)
        require(context == expected_context,'Context metadata differs')
        require(sha256(canonical(context)).hexdigest() == record['context_sha256'],'Context hash mismatch')
        driver,driver_pin=inputs.load(directory+'/driver.json')
        expected_driver=dict(schema='gamma-composite-weight-band-driver-v30',source_sha256=execution,proposals=proposals,
            requested_weight_range=[first,last],maximum_nodes=tree['max_nodes'],cell_failure_bits=tree['cell_failure_bits'],
            requested_failure_bits=record['requested_failure_bits'],hard_wall_seconds=600,cooperative_seconds=producer.COMPUTE_SECONDS,
            memory_bytes=4<<30,table_bytes=64<<20,maximum_output_bytes=32<<20,minimum_free_bytes=20<<30,
            maximum_simultaneous_banks=1,fresh_banks_required=True,directed=False,whole_code_certificate=False)
        require(driver == expected_driver,'Exact composite driver mismatch')
        packed,context_pin=inputs.load(directory+'/context.json')
        require(packed == dict(context=context,context_sha256=record['context_sha256'],directed=False,
            complete_requested_band=False,whole_code_certificate=False),'External context mismatch')
        for i,(component,proposal) in enumerate(zip(components,proposals)):
            budget.check()
            same(component,dict(index=i),'component ordering')
            require(set(component)=={'index','context','bank_metadata'} and component['context']['path']==directory+f'/component_{i:02d}_context.json'
                and component['bank_metadata']['path']==directory+f'/bank_{i:02d}.json','Wrong component paths')
            saved,unused=inputs.linked(component['context']);bank,unused=inputs.linked(component['bank_metadata'])
            component_context(bank,saved,proposal,first,last,controller)
            literal.fresh_source_maps(saved,pins);literal.fresh_source_maps(bank,pins)
        count=tree['checkpoints']
        require(type(count) is int and count==len(tree['nodes'])+1 and len(tree['nodes'])<=tree['max_nodes'],
                'Checkpoint count/limit')
        replay={}
        for sequence in range(count):
            budget.check();checkpoint,unused=inputs.load(directory+f'/node_{sequence:06d}.json')
            same(checkpoint,dict(sequence=sequence,context_sha256=record['context_sha256'],complete_requested_band=False,
                covers_all_supports_in_band=False,directed=False,batch_source_verified_at_exit=False,
                log2_upper_rational=None,whole_code_certificate=False),'partial checkpoint scope')
            for node in checkpoint['changed_nodes']:replay[node['id']]=node
        require([replay[i] for i in sorted(replay)]==tree['nodes'],'Checkpoint tree replay differs')
        expected_nodes={HERE/directory/f'node_{i:06d}.json' for i in range(count)}
        require(set((HERE/directory).glob('node_*.json'))==expected_nodes,'Unexpected checkpoint files')
        literal.fresh_source_maps(context,pins);literal.fresh_source_maps(driver,pins)
        result=audit_numeric(tree,len(proposals),budget)
        output_bytes=sum(pin['bytes'] for path,pin in inputs.pins.items() if path not in proposal_names)
        require(output_bytes<=32<<20 and record['bytes_before_summary']==output_bytes-summary_pin['bytes'], 'Output accounting differs')
        budget.check()
        return dict(schema=SCHEMA,audit_result='PASS',input_summary=summary_pin,execution_driver_receipt=driver_pin,
            context_file=context_pin,context_sha256=record['context_sha256'],construction=record['construction'],field_order=1<<128,
            state_dimension=30,target_distance='2/5',requested_weight_range=[first,last],requested_weights=weights,
            complete_requested_band=True,covers_all_matching_profiles=True,integer_partition_checked=True,
            counts_recomputed=True,count_algorithm='frozen independent raw-midpoint rectangle prefix recurrence',
            component_count=len(proposals),component_literal_sources=[p['source'] for p in proposals],
            every_original_weight_preparation_metadata_checked=True,normalized_histograms=result['normalized_histograms'],
            exact_literal_bank_context_matches=True,generic_partition_flags_preserved=True,component_minima_recomputed=True,
            complete_component_scalar_enclosures_checked=True,component_outer_prefactors_already_included=True,
            accepted_witness_counts=result['accepted_witness_counts'],cells_checked=result['cells_checked'],accepted_cells=result['accepted_cells'],
            checkpoints_checked=count,cell_log_enclosures_checked=True,positive_count_zero_cells=result['positive_count_zero_cells'],
            directed_union_precision=384,requested_failure_bits=record['requested_failure_bits'],cell_failure_bits=tree['cell_failure_bits'],
            saved_band_log2_upper_rational=record['log2_upper_rational'],reaggregated_cell_log2_upper_rational=result['union'],
            log2_upper_rational=record['log2_upper_rational'],zero=result['zero'],meets_requested_band_target=result['meets_requested_band_target'],
            source_sha256=pins,producer_sha256=PRODUCER_SHA,numerical_auditor_sha256=NUMERIC_SHA,
            checked_file_pins=list(inputs.pins.values()),source_and_input_bytes_checked_at_exit=True,
            producer_json_bytes=output_bytes,producer_json_budget_bytes=32<<20,historical_origin_verified=False,
            coefficient_replay=False,normalized_regional_replay=False,cell_transfer_recomputed=False,
            no_weight_interpolation=True,whole_code_certificate=False,
            audit_limit='Saved bounds/minima, counts, tree and384 union only; no coefficient/region/transfer/history replay. Provenance is not authentication.')
    finally:
        numeric.producer.check_source_pins(pins)
        inputs.check()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary',required=True);parser.add_argument('--expected-sha256',required=True)
    parser.add_argument('--output',required=True);args=parser.parse_args()
    require(sha256(Path(literal.limits.__file__).read_bytes()).hexdigest()==literal.LIMIT_SHA,'Hard-limit helper changed')
    with literal.limits.remote_limits():
        output=literal.portable_path(args.output)
        require(not output.exists() and shutil.disk_usage(output.parent).free-MAX_OUTPUT_BYTES>=RESERVE,'Audit output/reserve preflight failed')
        result=audit_summary(args.summary,args.expected_sha256)
        raw=(json.dumps(result,indent=2,allow_nan=False)+'\n').encode()
        require(len(raw)<=MAX_OUTPUT_BYTES and shutil.disk_usage(output.parent).free-len(raw)>=RESERVE,'Audit output budget/reserve exceeded')
        pin=numeric.producer.profiles.write_new(output,result)
    print(json.dumps(dict(artifact=pin,audit_result=result['audit_result'],log2_upper_rational=result['log2_upper_rational'])),flush=True)


if __name__=='__main__':main()
