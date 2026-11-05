from contextlib import ExitStack
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import hd_field128_gamma_composite_weight_band_audit_v30 as mode
from hd_field128_gamma_composite_weight_band_v30_tests import table, controller


def partition(pairs=4, first=3, last=5):
    parts=[table(pairs,first,last,i,constant=-200,slope=0) for i in range(2)]
    tables=mode.producer.CompositeTables(parts,controller.bands)
    context=dict(requested_weight_range=[first,last],requested_weights=list(range(first,last+1)),
                 root_cells=[c.descriptor() for c in tables.roots])
    evaluator=SimpleNamespace(tables=tables,_active=False,check_pins=lambda:None,context=lambda:deepcopy(context))
    checkpoints=[]
    record=controller.cover(evaluator,max_nodes=64,cell_failure_bits=54,aggregate_bits=48,checkpoint=checkpoints.append)
    return record,checkpoints


class Fixture:
    def __init__(self):
        self.files={};self.execution={'producer.py':'a'*64};self.pins=dict(self.execution,**{'audit.py':'b'*64})
        self.component_pins={'component.py':'c'*64};self.bank_pins={'bank.py':'d'*64}
        self.proposals=[];prepared=[]
        for i in range(2):
            witness=dict(z='1/4',alpha='1/8',markers=[['1',str(i+1)]],grid=128,state_dimension=30)
            original=dict(active_groups=2,construction=mode.producer.selected.literals.CONSTRUCTION,witness=witness,
                state_dimension=30,target_distance='2/5')
            source=self.put(f'p{i}.json',original)
            proposal=dict(source=source,witness=witness,source_marker_key='markers',calibration_weight=2,
                source_role='Exact literals only; every requested original weight is normalized afresh')
            self.proposals.append(proposal)
            bank=object.__new__(controller.cert.AliasBank)
            for key,value in dict(construction=mode.producer.selected.literals.CONSTRUCTION,state=30,q=1<<128,
                field_bits=128,forbidden=838860,distance=F(2,5),outer_measure='invertible-toeplitz',precision=256,
                bits=80,z=F('1/4'),alpha=F('1/8'),markers=((F(1),F(i+1)),),grid=128,diagnostics=[],pins=self.bank_pins).items():
                setattr(bank,key,value)
            metadata=bank.metadata();bank.constant=bank.slope=0
            tiny=table(4096,2,3,i,constant=-200,slope=0)
            with patch.object(controller.bands,'source_pins',return_value=self.component_pins),\
                 patch.object(controller.bands.Evaluator,'check_pins'),\
                 patch.object(controller.bands.restored,'RegionalCache',return_value=SimpleNamespace(regional=None)),\
                 patch.object(controller.bands,'BandTables',return_value=tiny):
                saved=controller.bands.Evaluator(bank,2,3,witness_source=source).context()
            metadata.update(exact_input_witness=witness,proposal_source=source,calibration_weight=2,
                requested_weight_range=[2,3],fresh_bank_rebuilt=True,
                literal_source_role='calibration literals only; no calibration bound inherited')
            prepared.append((tiny,saved,metadata))
        def publish(i,context,metadata):
            return dict(index=i,context=self.put(f'run/component_{i:02d}_context.json',context),
                bank_metadata=self.put(f'run/bank_{i:02d}.json',metadata))
        nodes=[]
        with patch.object(mode.producer,'prepare_component',side_effect=prepared),\
             patch.object(mode.producer,'source_pins',return_value=self.execution),\
             patch.object(mode.producer.selected,'check_pins'):
            evaluator=mode.producer.CompositeEvaluator(self.proposals,2,3,controller,controller.base.Budget(),publish)
            part=controller.cover(evaluator,max_nodes=64,max_seconds=540,cell_failure_bits=54,aggregate_bits=48,checkpoint=nodes.append)
            self.record=mode.producer.production_envelope(evaluator,part,controller)
        self.part=part;context=part['context'];digest=part['context_sha256']
        for checkpoint in nodes:
            self.put(f"run/node_{checkpoint['sequence']:06d}.json",checkpoint)
        self.put('run/context.json',dict(context=context,context_sha256=digest,directed=False,
            complete_requested_band=False,whole_code_certificate=False))
        self.put('run/driver.json',dict(schema='gamma-composite-weight-band-driver-v30',source_sha256=self.execution,
            proposals=self.proposals,requested_weight_range=[2,3],maximum_nodes=64,cell_failure_bits=54,requested_failure_bits=48,
            hard_wall_seconds=600,cooperative_seconds=540.,memory_bytes=4<<30,table_bytes=64<<20,maximum_output_bytes=32<<20,
            minimum_free_bytes=20<<30,maximum_simultaneous_banks=1,fresh_banks_required=True,directed=False,whole_code_certificate=False))
        self.record.update(total_compute_seconds=1.,
            runner=dict(path=mode.producer.NAME,sha256=mode.PRODUCER_SHA))
        self.record['bytes_before_summary']=sum(len(raw) for name,raw in self.files.items() if name.startswith('run/'))
        self.summary=self.put('run/summary.json',self.record)

    def put(self,name,value):
        raw=(json.dumps(value,indent=2,allow_nan=False)+'\n').encode();self.files[name]=raw
        return dict(path=name,sha256=sha256(raw).hexdigest())

    def run(self):
        fixture=self
        class Inputs:
            def __init__(self):self.pins={}
            def load(self,name,expected=None,**unused):
                raw=fixture.files[name];actual=sha256(raw).hexdigest()
                mode.require(expected is None or expected==actual,'Input mismatch')
                pin=dict(path=name,sha256=actual,bytes=len(raw));self.pins[name]=pin
                return json.loads(raw),pin
            def linked(self,pin,**options):return self.load(pin['path'],pin['sha256'],**options)
            def check(self):
                for name,pin in self.pins.items():mode.require(sha256(fixture.files[name]).hexdigest()==pin['sha256'],'Exit mismatch')
        pins=dict(self.pins,**self.component_pins,**self.bank_pins)
        with ExitStack() as stack:
            for obj,name,value in ((mode,'Inputs',Inputs),(mode,'source_pins',lambda c:pins),
                (mode.producer,'helpers',lambda:controller),(mode.producer,'source_pins',lambda c:self.execution),
                (controller,'source_pins',lambda:self.part['controller_source_sha256']),
                (controller.bands,'source_pins',lambda:self.component_pins),(controller.cert,'source_pins',lambda:self.bank_pins),
                (mode.numeric.producer,'check_source_pins',Mock()),
                (Path,'glob',lambda p,pattern:(mode.HERE/name for name in self.files if name.startswith('run/node_')))):
                stack.enter_context(patch.object(obj,name,value))
            return mode.audit_summary('run/summary.json',self.summary['sha256'])


class AuditTests(unittest.TestCase):
    def test_tiny_domains_counts_and_union(self):
        for first,last in ((1,3),(3,5),(5,7)):
            part,unused=partition(4,first,last)
            result=mode.audit_numeric(part,2,mode.numeric.base.Budget(),pairs=4)
            self.assertEqual(result['weights'],list(range(first,last+1)))
            self.assertTrue(result['meets_requested_band_target'])
            self.assertEqual(sum(result['accepted_witness_counts']),result['accepted_cells'])

    def test_memory_envelope(self):
        fixture=Fixture();result=fixture.run()
        self.assertEqual(result['schema'],mode.SCHEMA);self.assertEqual(result['requested_weights'],[2,3])
        self.assertTrue(result['component_minima_recomputed']);self.assertFalse(result['coefficient_replay'])
        self.assertFalse(result['historical_origin_verified']);self.assertFalse(result['whole_code_certificate'])
        self.assertFalse(fixture.record['partition_result']['directed'])

    def test_original_bank_maps_cannot_be_hidden_by_context_projection(self):
        fixture=Fixture()
        bank=json.loads(fixture.files['run/bank_00.json'])
        context=json.loads(fixture.files['run/component_00_context.json'])
        with patch.object(controller.cert,'source_pins',return_value=fixture.bank_pins),\
             patch.object(controller.bands,'source_pins',return_value=fixture.component_pins):
            def check(value):mode.component_context(value,context,fixture.proposals[0],2,3,controller)
            check(bank)
            for key in ('source_sha256','bank_source_sha256'):
                for value in ({},None):
                    changed=deepcopy(bank)
                    if value is None:del changed[key]
                    else:changed[key]=value
                    with self.subTest(key=key,value=value),self.assertRaises(ValueError):check(changed)
            for key in ('construction','field_order'):
                changed=deepcopy(bank);del changed[key]
                with self.subTest(key=key),self.assertRaises(ValueError):check(changed)

    def test_partition_count_minimum_mutations(self):
        changes=(lambda p:p['nodes'][0]['receipt'].__setitem__('count_upper',0),
            lambda p:p['nodes'][0]['receipt'].__setitem__('component_winner',True),
            lambda p:p['nodes'][0]['receipt']['component_bounds'].pop(),
            lambda p:p['nodes'][0]['receipt'].__setitem__('entrywise_cross_witness_minimum',True),
            lambda p:p['nodes'][0]['receipt']['component_bounds'][0].__setitem__('log2_upper_enclosure','100'),
            lambda p:p.__setitem__('pending_frontier',[0]),lambda p:p['nodes'].pop(),
            lambda p:p.__setitem__('log2_upper_rational','-1000'))
        for change in changes:
            part,unused=partition();change(part)
            with self.assertRaises((ValueError,KeyError)):mode.audit_numeric(part,2,mode.numeric.base.Budget(),pairs=4)

    def test_envelope_context_checkpoint_mutations(self):
        for filename,change in (('run/summary.json',lambda v:v['partition_result'].__setitem__('directed',True)),
            ('run/component_00_context.json',lambda v:v.__setitem__('requested_weights',[2])),
            ('run/bank_00.json',lambda v:v['witness'].__setitem__('alpha','1/2')),
            ('run/node_000001.json',lambda v:v['changed_nodes'][0]['receipt'].__setitem__('count_upper',0))):
            fixture=Fixture();value=json.loads(fixture.files[filename]);change(value);pin=fixture.put(filename,value)
            if filename=='run/summary.json':fixture.summary=pin
            with self.subTest(filename=filename),self.assertRaises(ValueError):fixture.run()

    def test_zero_and_ties(self):
        row=dict(component_bounds=[dict(witness_index=i,zero=True,log2_upper_rational=None,log2_upper_enclosure=None) for i in range(2)],
            component_winner=0,minimum_stage=mode.MINIMUM_STAGE,entrywise_cross_witness_minimum=False,
            zero=True,log2_upper_rational=None,log2_upper_enclosure=None)
        self.assertIsNone(mode.component_minimum(row,2))
        row['component_winner']=1
        with self.assertRaises(ValueError):mode.component_minimum(row,2)

    def test_main_hard_guard(self):
        with patch.object(mode.argparse.ArgumentParser,'parse_args',return_value=None),\
             patch.object(mode.literal.limits,'remote_limits',side_effect=RuntimeError('hard guard')):
            with self.assertRaisesRegex(RuntimeError,'hard guard'):mode.main()

if __name__=='__main__':unittest.main()
