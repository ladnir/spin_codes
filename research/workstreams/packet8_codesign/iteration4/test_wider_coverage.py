import copy
import unittest
import wider_coverage as wc


class CoverageTests(unittest.TestCase):
    def test_wrong_shape_or_map_rejected(self):
        shape=wc.wf.wider.geometry(symbol_bits=32)
        shape['state_bits']=24
        saved=dict(schema='packet8-iteration4-wider-actual24-fractional-proposal-1',
            geometry=shape,map_record={'state_bits':24},group_steps=4,
            source_pins_verified_at_finish=True,source_sha256={})
        wc.validate_receipt(saved,shape,{'state_bits':24})
        for key,value in (('geometry',{}),('map_record',{'state_bits':16}),
                          ('group_steps',8),('source_pins_verified_at_finish',False)):
            changed=copy.deepcopy(saved)
            changed[key]=value
            with self.assertRaises(ValueError):
                wc.validate_receipt(changed,shape,{'state_bits':24})

    def test_all_q1_shell_mass_is_nonzero_message_mass(self):
        shape=wc.wf.wider.geometry(symbol_bits=32)
        _,counts=wc.wf.wider.outer(shape)
        self.assertEqual(len(counts),65)
        self.assertEqual(counts[0],0)
        self.assertEqual(sum(counts),2**256-1)


if __name__=='__main__': unittest.main()
