"""Check the no-constant nested maps and isolated model selection."""
import unittest

import parameter_no_constant as maps


class NoConstantTests(unittest.TestCase):
    def test_nested_maps_and_spectra(self):
        for t in (64,128,256):
            records = [maps.inner(t,s) for s in (t.bit_length(),t.bit_length()+1,20)]
            for record in records:
                self.assertNotIn(t,record['spectrum'])
                self.assertEqual(sum(record['spectrum'].values()),(1 << record['s'])-1)
                self.assertNotIn(0,record['expansion_columns'])
                self.assertEqual(record['feedback_columns'],maps.grid.inner(t,record['s'])['feedback_columns'])
            a,b = records[:2]
            mask = (1 << a['s'])-1
            self.assertEqual(a['expansion_columns'],[v & mask for v in b['expansion_columns']])

    def test_model_selection_restores_even_on_error(self):
        original = maps.grid.inner
        with self.assertRaisesRegex(RuntimeError,'test exit'):
            with maps.use():
                self.assertIs(maps.grid.inner,maps.inner)
                self.assertNotIn(64,maps.full.prepare(64,7)[0]['spectrum'])
                raise RuntimeError('test exit')
        self.assertIs(maps.grid.inner,original)
        self.assertIn(64,maps.full.prepare(64,7)[0]['spectrum'])


if __name__ == '__main__':
    unittest.main()
