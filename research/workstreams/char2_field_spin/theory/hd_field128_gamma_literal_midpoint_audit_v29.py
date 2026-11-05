"""Literal-only saved-row audit using frozen counts/prefix arithmetic.
No coefficient, regional or transfer replay; no optimizer history reads.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath, PureWindowsPath
import shutil
from time import monotonic

import hd_field128_gamma_midpoint_receipt_audit as frozen
import hd_field128_gamma_portable_midpoint_directed_v29 as runner
import hd_field128_gamma_portable_neighbors_v27 as limits

HERE = Path(__file__).resolve().parent
NAME = Path(__file__).name
FROZEN_SHA = 'cd54e300893a3c657b7b95780523abd37b791032ea691a38b44e389c00129108'
RUNNER_SHA = 'd18e28ea1a504d592ed1c2670cc129158de528dd0fbcd7522b538d3d05baff96'
LIMIT_SHA = '169ee5c7719873ce4528ec9e9a5eea1e87e8c866ae32e713abde0d3678339a94'
SCHEMA = 'gamma-midpoint-selected-literal-only-remote-audit-v1'
MAX_OUTPUT_BYTES, RESERVE = 256 << 10, 20 << 30
require, same, fraction = frozen.require, frozen.same, frozen.fraction


def source_pins():
    for module, expected in ((frozen,FROZEN_SHA),(runner,RUNNER_SHA),(limits,LIMIT_SHA)):
        require(sha256(Path(module.__file__).read_bytes()).hexdigest()==expected,'Reviewed audit/runner changed')
    pins=frozen.source_pins()
    for name,digest in runner.source_pins().items():
        require(name not in pins or pins[name]==digest,'Conflicting source pin')
        pins[name]=digest
    pins[NAME]=sha256((HERE/NAME).read_bytes()).hexdigest()
    pins[limits.NAME]=LIMIT_SHA
    return dict(sorted(pins.items()))


def portable_path(name):
    require(isinstance(name,str) and '\\' not in name and ':' not in name
        and not PurePosixPath(name).is_absolute() and not PureWindowsPath(name).is_absolute()
        and all(part not in ('','.','..') for part in name.split('/')),'Expected theory-relative path')
    path=(HERE/name).resolve()
    path.relative_to(HERE.resolve())
    return path


def fresh_source_maps(value,pins):
    def mapping(value):
        require(isinstance(value,dict),'Source map must be an object')
        for name,digest in value.items():
            if isinstance(digest,dict):mapping(digest)
            else:require(name in pins and digest==pins[name],'Unknown or changed fresh source pin')
    if isinstance(value,dict):
        for key,child in value.items():
            if key=='source_sha256' or key.endswith('_source_sha256'):mapping(child)
            else:fresh_source_maps(child,pins)
    elif isinstance(value,list):
        for child in value:fresh_source_maps(child,pins)


def public_metadata(record):
    weight=record['active_groups']; first,last,rows,profiles=runner.geometry(weight)
    target=runner.target_bits(record['requested_failure_bits'])
    bound=fraction(record['log2_upper_rational']); meets=bound<=-target
    same(record,dict(frozen.previous.EXPECTED,status='complete',failure=None,
        complete_selected_weight=True,covers_all_supports_at_weight=True,covers_all_even_odd_splits=True,
        repetitions=256,paired_region_power=128,common_potential_used=False,
        source_pins_unchanged_at_exit=True,source_sha256=frozen.producer.source_pins(),
        coverage_driver_sha256=frozen.producer.source_pins(),
        bank_source_sha256=frozen.producer.provider.source_pins(),
        expected_even_rows=rows,expected_profile_pairs=profiles,completed_even_rows=rows,
        completed_even_full_range=[first,last],pending_even_full_range=None,profile_pairs_evaluated=profiles,
        required_normal_process_exit=True,meets_requested_selected_weight_target=meets,
        wrapper_admitted=meets,literal_seed_only=True,historical_origin_verified=False,
        calibration_history_replayed=False,old_numerical_bounds_reused=False,
        construction_unchanged=True,no_coefficient_truncation=True),'complete selected receipt')
    same(record['witness'],dict(state_dimension=30,inner_width=64,gamma_extension_degree=16,
        gamma_base_field_order=256,joint_rank_return_cap=True,all_pair_before_alpha=True,
        precision=256,grid=128,packet_mantissa_bits=80),'bank arithmetic')
    require(type(record['count_provenance']['maximum_refined_singletons']) is int
        and record['count_provenance']['maximum_refined_singletons']==64,'Wrong count cutoff')
    for key,cap in (('max_seconds',595),('total_job_seconds',600)):
        require(type(record[key]) in (int,float) and 0<record[key]<=cap,'Wrong reported time limit')
    return weight,target,meets


def expected_context(record,bank):
    result={key:value for key,value in bank.items()
        if key not in ('witness_source','runner_source_sha256','runtime','fresh_bank_rebuilt')}
    weight=record['active_groups']; unused,unused,rows,profiles=runner.geometry(weight)
    result.update(source_sha256=frozen.producer.source_pins(),
        coverage_driver_sha256=frozen.producer.source_pins(),witness_source=record['witness_source'],
        active_groups=weight,common_potential_used=False,
        transfer_method='frozen positive Arb row-zero sum of (R_even R_odd)^128',
        regional_normalization='frozen RegionalCache, entrywise minimum over all bank markers',
        repetitions=256,paired_region_power=128,requested_failure_bits=record['requested_failure_bits'],
        expected_even_rows=rows,expected_profile_pairs=profiles,max_seconds=record['max_seconds'],
        old_numerical_bounds_reused=False,directed=True,whole_code_certificate=False,
        count_provenance=record['count_provenance'])
    return result


def audit_summary(summary_name,*,expected_summary_sha256,expected_seed_sha256,max_seconds=120.):
    require(type(max_seconds) in (int,float) and 0<max_seconds<=120,'Expected at most120 seconds')
    previous=runner.baseline(); limits=previous.validation()
    limits.digest(expected_summary_sha256);limits.digest(expected_seed_sha256)
    portable=previous.pinned_module(runner.PORTABLE,runner.PORTABLE_SHA)
    started,pins=monotonic(),source_pins()
    frozen.producer.check_source_pins(pins)
    loader=frozen.index.Validator(F(2,5))
    def check_time():require(monotonic()-started<max_seconds,'Saved-data audit budget exhausted')
    record,summary_pin=loader.load(portable_path(summary_name),expected=expected_summary_sha256,role='summary')
    weight,target,meets=public_metadata(record)
    directory=Path(summary_pin['path']).parent
    link=record['witness_source']
    require(isinstance(link,dict) and set(link)=={'path','sha256'}
        and link['sha256']==expected_seed_sha256,'Wrong exact portable seed pin')
    name=portable.relative_name(link['path'],suffix='.json')
    seed,seed_pin=loader.load(portable_path(name),expected=expected_seed_sha256,role='literal-only-seed')
    portable.validate_seed(seed)
    for key in ('z','alpha','markers','grid','state_dimension'):
        require(record['witness'][key]==seed['witness'][key],'Seed/receipt exact literal mismatch')
    same(record,dict(source_calibration_weight=seed['active_groups']),'seed calibration weight')
    bank,bank_pin=loader.load(directory/'bank.json',role='fresh-bank')
    driver,driver_pin=loader.load(directory/'driver.json',role='execution-driver')
    saved_context,context_pin=loader.load(directory/'context.json',role='external-context')
    execution=runner.source_pins()
    require(driver==runner.driver_record(execution,link,seed,weight,runner.MEMORY_BYTES,target),
        'Exact v29 driver mismatch')
    same(record,dict(execution_driver=driver['driver'],execution_source_sha256=execution),'summary execution')
    same(bank,dict(witness_source=link,runner_source_sha256=execution,fresh_bank_rebuilt=True,
        directed=False,whole_code_certificate=False),'fresh bank scope')
    require(bank['witness']==record['witness'],'Bank/receipt witness mismatch')
    require(bank['runtime']['OPENBLAS_NUM_THREADS']=='1','Reported BLAS thread count differs')
    context=expected_context(record,bank)
    same(record,context,'inline metadata context')
    digest=sha256(json.dumps(context,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
    require(record['context_sha256']==digest,'Reconstructed context digest mismatch')
    require(saved_context==dict(context=context,context_sha256=digest,complete_selected_weight=False,
        covers_all_supports_at_weight=False,whole_code_certificate=False,scope='metadata only; no numerical bound'),
        'External context mismatch')
    first,last,rows,profiles=runner.geometry(weight)
    require(isinstance(record['rows'],list) and len(record['rows'])==rows,'Missing complete inline rows')
    paths=[directory/name for name in ('driver.json','bank.json','context.json','summary.json')]
    for position,row in zip(range(first,last+1),record['rows']):
        check_time();require(row['even_full']==position,'Wrong original-r row order')
        path=directory/f'row_{position:04d}.json';paths.append(path)
        checkpoint,unused=loader.load(path,role='external-row')
        require(checkpoint==dict(row,context_sha256=digest,active_groups=weight,directed=True,
            whole_code_certificate=False,
            scope='completed even-profile row only; requires midpoint-count metadata context'),'External row mismatch')
    require(set(directory.glob('row_*.json'))==set(paths[4:]),'Unexpected external row files')
    output_bytes=sum(path.stat().st_size for path in paths)
    require(output_bytes<=runner.MAX_OUTPUT_BYTES,'Producer exceeded2MiB JSON budget')
    for value in (record,bank,driver,context):fresh_source_maps(value,pins)
    check_time()
    report=frozen.audit_record(record,pairs=4096,precision=384)
    check_time()
    report.update(receipt=summary_pin,portable_seed=seed_pin,requested_failure_bits=target,
        meets_requested_target=meets)
    union=frozen.index.union_ceiling([fraction(report['saved_log2_upper_rational'])],384)
    require(source_pins()==pins,'Audit closure changed')
    frozen.producer.check_source_pins(pins)
    frozen.previous.band._check_links(tuple(loader.files.values()))
    return dict(schema=SCHEMA,audit_result='PASS',input_summary=summary_pin,portable_seed=seed_pin,
        execution_driver_receipt=driver_pin,bank_metadata=bank_pin,context_file=context_pin,
        construction=record['construction'],field_order=1<<128,state_dimension=30,target_distance='2/5',
        selected_weights=[weight],weight_audits=[report],
        requested_failure_bits=target,meets_requested_selected_weight_target=meets,
        complete_selected_weight=True,covers_all_matching_profiles=True,
        counts_recomputed=True,exact_seed_literals_checked=True,exact_driver_record_checked=True,
        fresh_bank_metadata_checked=True,external_checkpoint_context_reconstructed=True,
        external_rows_checked=rows,external_rows_equal_complete_inline_rows=True,context_sha256=digest,
        precision_bits=384,saved_log2_upper_rational=record['log2_upper_rational'],
        selected_events_union_log2_upper_rational=str(union),source_sha256=pins,
        numerical_auditor_sha256=FROZEN_SHA,execution_driver=driver['driver'],execution_source_sha256=execution,
        checked_file_pins=list(loader.files.values()),source_and_input_bytes_checked_at_exit=True,
        producer_json_bytes=output_bytes,producer_json_budget_bytes=runner.MAX_OUTPUT_BYTES,
        source_calibration_weight=seed['active_groups'],historical_origin_verified=False,
        calibration_history_replayed=False,historical_bounds_used=False,
        coefficient_replay=False,normalized_regional_replay=False,transfer_recomputed=False,
        coverage_index_modified=False,whole_code_certificate=False,
        normal_process_exit_evidence='External evidence required; JSON is not process authentication',
        audit_elapsed_seconds=monotonic()-started,
        audit_limit='Counts, saved rows and original-r prefix at384 bits only. Producer arithmetic and midpoint theorem are inherited. Fresh execution and runtime are unauthenticated provenance.')


def write_remote(name,report):
    path=portable_path(name)
    raw=(json.dumps(report,indent=2,allow_nan=False)+'\n').encode()
    require(len(raw)<=MAX_OUTPUT_BYTES,'Audit exceeds256KiB')
    require(shutil.disk_usage(path.parent).free-len(raw)>=RESERVE,'Audit breaches20GiB reserve')
    return frozen.producer.write_new(path,report)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('summary')
    parser.add_argument('--expected-summary-sha256',required=True)
    parser.add_argument('--expected-seed-sha256',required=True)
    parser.add_argument('--output',required=True)
    parser.add_argument('--max-seconds',type=float,default=120.)
    args=parser.parse_args()
    require(sha256(Path(limits.__file__).read_bytes()).hexdigest()==LIMIT_SHA,'Resource helper changed')
    with limits.remote_limits():
        output=portable_path(args.output)
        require(not output.exists(),'Audit output already exists')
        require(shutil.disk_usage(output.parent).free-MAX_OUTPUT_BYTES>=RESERVE,'Audit reserve preflight failed')
        result=audit_summary(args.summary,expected_summary_sha256=args.expected_summary_sha256,
            expected_seed_sha256=args.expected_seed_sha256,max_seconds=args.max_seconds)
        pin=write_remote(args.output,result)
    print(json.dumps(dict(artifact=pin,selected_weights=result['selected_weights'],
        log2_upper_rational=result['saved_log2_upper_rational'],historical_origin_verified=False)),flush=True)


if __name__=='__main__':main()
