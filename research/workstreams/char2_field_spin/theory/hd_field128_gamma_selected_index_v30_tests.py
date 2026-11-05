"""Memory-only index controls."""
from contextlib import ExitStack
from fractions import Fraction as F
from hashlib import sha256
import json
import unittest
from unittest.mock import patch

import hd_field128_gamma_selected_index_v30 as mode


def encode(value):
    return (json.dumps(value,sort_keys=True,allow_nan=False)+'\n').encode()


SOURCES = {'single':{'single.py':'a'*64}, 'multiple':{'multiple.py':'b'*64}}
EXECUTION = {'runner.py':'c'*64}


def audit(kind, r, ceiling='-49', target=48):
    directory = ('/home/prindal/job/theory/' if kind=='single' else '')+f'r{r}'
    first,last=max(0,r-4096),r//2; rows=last-first+1
    files = {}
    def pin(name):
        path=directory+'/'+name; digest=sha256(path.encode()).hexdigest()
        files[path]=digest
        return dict(path=path,sha256=digest)
    summary=pin('summary.json');driver=pin('driver.json');context=pin('context.json')
    base=dict(schema=mode.single.SCHEMA if kind=='single' else mode.multiple.SCHEMA,
        audit_result='PASS',construction=mode.CONSTRUCTION,field_order=1<<128,state_dimension=30,
        target_distance='2/5',complete_selected_weight=True,covers_all_matching_profiles=True,
        historical_origin_verified=False,normalized_regional_replay=False,whole_code_certificate=False,
        meets_requested_selected_weight_target=True,requested_failure_bits=target,
        source_sha256=SOURCES[kind],input_summary=summary,context_sha256='d'*64)
    if kind=='single':
        bank=pin('bank.json');seed=dict(path='/home/prindal/job/theory/seed.json',sha256='e'*64)
        files[seed['path']]=seed['sha256']
        for full in range(first,last+1):pin(f'row_{full:04d}.json')
        base.update(selected_weights=[r],counts_recomputed=True,exact_seed_literals_checked=True,
            exact_driver_record_checked=True,fresh_bank_metadata_checked=True,
            bank_context_scope_checked=True,selected_receipt_scope_checked=True,
            context_projection_overwritten_keys=['scope'],
            external_checkpoint_context_reconstructed=True,external_rows_checked=rows,
            external_rows_equal_complete_inline_rows=True,precision_bits=384,
            numerical_auditor_sha256=mode.single.FROZEN_SHA,execution_source_sha256=EXECUTION,
            execution_driver=dict(path=mode.single.runner.NAME,sha256=mode.single.RUNNER_SHA),
            source_and_input_bytes_checked_at_exit=True,producer_json_budget_bytes=2<<20,
            producer_json_bytes=1000,source_calibration_weight=913,calibration_history_replayed=False,
            historical_bounds_used=False,coefficient_replay=False,transfer_recomputed=False,
            coverage_index_modified=False,saved_log2_upper_rational=ceiling,
            selected_events_union_log2_upper_rational=str(mode.index.union_ceiling([F(ceiling)],384)),
            weight_audits=[dict(active_groups=r,even_rows=rows,profile_pairs=rows*rows,
                all_producer_midpoint_and_final_counts_match=True,saved_ceiling_dominates_row_enclosures=True,
                requested_failure_bits=target,meets_requested_target=True,exact_count_table_sha256='f'*64,
                saved_log2_upper_rational=ceiling,saved_row_reaggregation_log2_upper_rational=ceiling)],
            execution_driver_receipt=driver,bank_metadata=bank,context_file=context,portable_seed=seed)
    else:
        for i in range(2):pin(f'bank_{i:02d}.json');pin(f'proposal_{i}.json')
        for i in range(rows):pin(f'row_{i:06d}.json')
        base.update(active_groups=r,producer_sha256=mode.multiple.PRODUCER_SHA,even_rows=rows,
            profile_pairs=rows*rows,exact_literal_bank_context_matches=True,original_r_prefixes_recomputed=True,
            coefficient_replay_in_this_audit=False,transfer_or_winner_replay=False,
            row_partition_and_count_policy_checked=True,directed_union_precision=384,second_outer_prefix_applied=False,
            witness_winning_profile_counts=[rows*rows,0],positive_count_profiles=rows*rows,
            log2_upper_rational=ceiling,fresh_saved_row_union_384=dict(zero=False,log2_upper_rational=ceiling))
    base['checked_file_pins']=[dict(path=p,sha256=d) for p,d in files.items()]
    return base


class Fixture:
    def __init__(self):
        self.data={};self.specs=[]
        self.parent_name='hd_field128_gamma30_d40_partial_cover_v29_20261010.json'
        self.check_name='hd_field128_gamma30_d40_partial_cover_v29_20261010_independent_check.json'
        for name in (self.parent_name,self.check_name):self.data[name]=(mode.HERE/name).read_bytes()
        self.parent=json.loads(self.data[self.parent_name])
        self.links_checked=False;self.sources_checked=False

    def put(self,name,value):
        raw=encode(value);self.data[name]=raw
        return dict(path=name,sha256=sha256(raw).hexdigest())

    def add(self,kind,r,ceiling='-49',target=48):
        report=audit(kind,r,ceiling,target)
        auditpin=self.put(f'compact/audit_{r}.json',report)
        relative=f'r{r}/summary.json';auditrelative=f'r{r}/independent_audit.json'
        archive=dict(schema=mode.ARCHIVE_SCHEMA,host=mode.prior.HOST,zip_crc_checked=True,
            all_member_sha256_checked=True,exact_manifest_membership_checked=True,
            all_audit_input_files_included=True,audit_source_closure_included=True,
            storage_scope='remote durable archive; no local tree import',
            archive=dict(path=str(mode.ARCHIVE_ROOT/f'r{r}.zip'),sha256='1'*64,bytes=100000),
            manifest=dict(path=str(mode.ARCHIVE_ROOT/f'r{r}.json'),sha256='2'*64),
            relative_summary_path=relative,relative_audit_path=auditrelative,
            summary=dict(member='theory/'+relative,sha256=report['input_summary']['sha256']),
            audit=dict(member='theory/'+auditrelative,sha256=auditpin['sha256']),
            file_count=len(report['checked_file_pins'])+10,logical_bytes=1000000,restore='Restore verified archive.',
            audit_source_closure_sha256=sha256(mode.canonical(report['source_sha256'])).hexdigest())
        archivepin=self.put(f'compact/archive_{r}.json',archive)
        spec=dict(audit=auditpin,archive_receipt=archivepin);self.specs.append(spec)
        return spec

    def change(self,spec,kind,mutate):
        old=spec[kind];obj=json.loads(self.data[old['path']]);mutate(obj)
        spec[kind]=self.put(old['path'],obj)
        if kind=='audit':
            link=spec['archive_receipt'];archive=json.loads(self.data[link['path']])
            archive['audit']['sha256']=spec['audit']['sha256']
            spec['archive_receipt']=self.put(link['path'],archive)

    def run(self):
        fixture=self
        class Loader:
            def __init__(self,unused):self.files={}
        def load(loader,path,expected,role):
            name=str(path).replace('\\','/')
            if name.startswith(mode.HERE.as_posix()+'/'):name=name[len(mode.HERE.as_posix())+1:]
            raw=fixture.data[name];digest=sha256(raw).hexdigest()
            mode.require(digest==expected,'Memory input pin mismatch')
            pin=dict(path=name,sha256=digest,bytes=len(raw),roles=[role]);loader.files[name]=pin
            return json.loads(raw),pin
        def links(pins):
            for pin in pins:mode.require(sha256(fixture.data[pin['path']]).hexdigest()==pin['sha256'],'Exit changed')
            fixture.links_checked=True
        def sources(pins):fixture.sources_checked=True
        with ExitStack() as stack:
            for obj,name,value in ((mode,'source_pins',lambda:{}),(mode,'audit_sources',lambda kind:SOURCES[kind]),
                (mode.single.runner,'source_pins',lambda:EXECUTION),(mode.index,'Validator',Loader),
                (mode.prior,'_load',load),(mode.prior.prior.band,'_check_links',links),
                (mode.frozen.producer,'check_source_pins',sources)):
                stack.enter_context(patch.object(obj,name,value))
            return mode.extend_v29(self.parent_name,self.check_name,self.specs)


class SelectedIndexTests(unittest.TestCase):
    def test_bounds_and_parent(self):
        f=Fixture();f.add('single',914,'-83',54);f.add('multiple',548,'-464',54);f.add('single',1145,'-49',48)
        result=f.run()
        self.assertEqual(result['entries'][:260],f.parent['entries'])
        self.assertEqual((result['receipt_count'],result['covered_weight_count'],result['pending_weight_count']),(263,8129,63))
        self.assertFalse(result['whole_code_certificate']);self.assertTrue(f.links_checked and f.sources_checked)

    def test_pin_and_overlap(self):
        f=Fixture();f.add('single',914);f.data[f.parent_name]+=b' '
        with self.assertRaisesRegex(ValueError,'pin mismatch'):f.run()
        f=Fixture();spec=f.add('single',914);f.specs.append(spec)
        with self.assertRaisesRegex(ValueError,'overlap'):f.run()
        f=Fixture();f.add('single',913)
        with self.assertRaisesRegex(ValueError,'overlap'):f.run()

    def test_global_union(self):
        f=Fixture();f.add('single',914,'-40',40)
        with self.assertRaisesRegex(ValueError,'union misses'):f.run()

    def test_invalid_reports(self):
        changes=[lambda d:d.update(whole_code_certificate=True),lambda d:d.update(audit_result='FAIL'),
            lambda d:d.update(schema='gamma-midpoint-selected-literal-only-remote-audit-v1'),
            lambda d:d.update(bank_context_scope_checked=False),
            lambda d:d.update(context_projection_overwritten_keys=['scope','witness']),
            lambda d:d.update(requested_failure_bits=39),
            lambda d:d.update(meets_requested_selected_weight_target=False),
            lambda d:d.update(external_rows_checked=1),lambda d:d['weight_audits'][0].update(profile_pairs=1),
            lambda d:d.update(source_sha256={}),lambda d:d['checked_file_pins'].pop()]
        for change in changes:
            f=Fixture();s=f.add('single',914);f.change(s,'audit',change)
            with self.assertRaises((ValueError,KeyError)):f.run()
        f=Fixture();s=f.add('multiple',548)
        f.change(s,'audit',lambda d:d.update(second_outer_prefix_applied=True))
        with self.assertRaises(ValueError):f.run()
        f=Fixture();s=f.add('multiple',548)
        f.change(s,'audit',lambda d:d.update(log2_upper_rational='-50'))
        with self.assertRaisesRegex(ValueError,'understated'):f.run()

    def test_archive_bindings(self):
        changes=[lambda d:d['archive'].update(path='/tmp/r.zip'),
            lambda d:d.update(relative_summary_path='../summary.json'),
            lambda d:d.update(all_audit_input_files_included=False),
            lambda d:d.update(audit_source_closure_sha256='0'*64),
            lambda d:d['summary'].update(sha256='0'*64),lambda d:d['audit'].update(member='theory/other/audit.json')]
        for change in changes:
            f=Fixture();s=f.add('multiple',548);f.change(s,'archive_receipt',change)
            with self.assertRaises(ValueError):f.run()

if __name__=='__main__':unittest.main()
