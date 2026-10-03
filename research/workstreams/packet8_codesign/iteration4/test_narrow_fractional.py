import unittest
import narrow_fractional as nf


class NarrowTests(unittest.TestCase):
    def test_full_geometry(self):
        shape=nf.geometry()
        self.assertEqual(shape['regions']*shape['physical_steps_per_region'],4096)
        self.assertEqual(shape['physical_steps_total']*shape['physical_t'],131072)
        self.assertEqual(shape['outer_groups']*shape['group_dimension'],65536)
        self.assertEqual(shape['physical_packet_slots'],4)

    def test_literal_map_matches_geometry(self):
        data,_=nf.wf.birth.maps.prepare('byte_native_t32')
        self.assertEqual((data['bits'],data['windows'],data['width']),(16,4,32))


if __name__=='__main__': unittest.main()
