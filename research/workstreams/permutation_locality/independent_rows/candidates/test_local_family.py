"""Memoized proof bounds must retain exact values and all dependencies."""
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from flint import arb,ctx

import local_family as local
from mass_density import blend
from occupancy_memory import Z,C
from test_mass_density import operators


class LocalMemo(unittest.TestCase):
    def setUp(self):
        self.precision=ctx.prec; ctx.prec=192
        self.directory=tempfile.TemporaryDirectory()
        self.key=local.parameters('.072','.9',192,16,10,0,True,10,'test-sources')
        self.family=(operators(32),{j:(arb(j)*arb(2)**-95,arb(3)*arb(2)**-230,arb(0)) for j in range(1,17)})

    def tearDown(self):
        self.directory.cleanup(); ctx.prec=self.precision

    def test_exact_roundtrip_and_same_coupled_replay(self):
        local.save(self.directory.name,self.key,self.family)
        loaded=local.load(self.directory.name,self.key)
        self.assertEqual(loaded,self.family)
        for base,coefficients in (self.family,loaded):
            changed=blend(base,coefficients,{6:Q(1,3),7:Q(1)},target=Z)
            changed=blend(changed,coefficients,{8:Q(2,5)},target=C)
            if base is self.family[0]: original=changed
            else: self.assertEqual(changed,original)
        self.assertFalse(list(Path(self.directory.name).glob('*.tmp')))

    def test_each_bound_option_and_source_changes_key(self):
        options=dict(self.key)
        variants=dict(precision=256,maximum=17,exact_feedback=9,spectral_feedback=8,
                      joint_four=False,density_through=9,sources='new-source',tilt='1/20',
                      penalty='1',rounds=1,input_weight='7/8',full_feedback=5,
                      window_histogram=7,joint_cancellation=False,column_density=False,
                      feedback_density=5,mature_tail=56,windows=31,dimension=12,
                      stage='other',schema=2)
        local.save(self.directory.name,self.key,self.family)
        for name,value in variants.items():
            changed=dict(options,**{name:value})
            self.assertNotEqual(local.cache_path('.',changed),local.cache_path('.',self.key))
            self.assertIsNone(local.load(self.directory.name,changed))
        canonical=local.parameters('9/125','9/10',192,16,10,0,True,10,'test-sources')
        self.assertEqual(canonical,self.key)

    def test_reject_corruption_even_with_matching_path(self):
        local.save(self.directory.name,self.key,self.family)
        path=local.cache_path(self.directory.name,self.key)
        original=json.loads(path.read_bytes())
        original['body']['base'][1][0][0][0]+=1
        path.write_bytes(local.operator_cache.encoded(original))
        with self.assertRaisesRegex(ValueError,'checksum'): local.load(self.directory.name,self.key)
        original['checksum']=hashlib.sha256(local.operator_cache.encoded(original['body'])).hexdigest()
        original['body']['key']['rounds']=1
        path.write_bytes(local.operator_cache.encoded(original))
        with self.assertRaisesRegex(ValueError,'key'): local.load(self.directory.name,self.key)

    def test_reject_precision_dimensions_and_noninteger_endpoints(self):
        local.save(self.directory.name,self.key,self.family)
        ctx.prec=128
        with self.assertRaisesRegex(ValueError,'precision'): local.load(self.directory.name,self.key)
        ctx.prec=192
        with self.assertRaisesRegex(ValueError,'dimensions'):
            local.save(self.directory.name,self.key,(self.family[0][:-1],self.family[1]))
        for pair in ([True,0],[1.0,0],[-1,0],[1,0,0],['1',0]):
            with self.assertRaises(ValueError): local.unpack(pair)
        ctx.prec=16
        with self.assertRaisesRegex(ValueError,'endpoint'): local.unpack([2**100+1,-100])

    def test_hit_skips_builder_and_support_search_does_not_bind_cache(self):
        key=local.parameters('.072','.9',192,16,0,0,False,0,'fixed')
        local.save(self.directory.name,key,self.family)
        with patch('local_family.source_digest',return_value='fixed'),patch('local_family.epoch_grid') as build:
            result=local.build(['.072'],'.9',directory=self.directory.name)
        build.assert_not_called()
        self.assertEqual(result['.072'],self.family)

    def test_partial_hit_builds_only_missing_tilts(self):
        key=local.parameters('.072','.9',192,16,0,0,False,0,'fixed')
        local.save(self.directory.name,key,self.family)
        with patch('local_family.source_digest',return_value='fixed'),\
                patch('local_family.epoch_grid',return_value=({'.076':self.family[0]},{})) as build,\
                patch('local_family.conditioned_coefficients',return_value=self.family[1]):
            result=local.build(['.072','.076'],'.9',directory=self.directory.name)
        self.assertEqual(set(result),{'.072','.076'})
        self.assertEqual(build.call_args.args,(['.076'],'.9',192))
        self.assertEqual(len(list(Path(self.directory.name).glob('*.json'))),2)

    def test_source_change_during_calculation_cannot_create_stale_memo(self):
        with patch('local_family.source_digest',side_effect=['before','after']),\
                patch('local_family.epoch_grid',return_value=({'.072':self.family[0]},{})),\
                patch('local_family.conditioned_coefficients',return_value=self.family[1]),\
                patch('builtins.print'):
            result=local.build(['.072'],'.9',directory=self.directory.name)
        self.assertEqual(result['.072'],self.family)
        self.assertFalse(list(Path(self.directory.name).glob('*.json')))

    def test_no_cache_is_fresh_and_no_digest_needed(self):
        with patch('local_family.source_digest') as digest,\
                patch('local_family.epoch_grid',return_value=({'.072':self.family[0]},{})) as build,\
                patch('local_family.conditioned_coefficients',return_value=self.family[1]):
            local.build(['.072'],'.9')
        build.assert_called_once(); digest.assert_not_called()


if __name__=='__main__':
    unittest.main()
