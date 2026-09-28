"""Occupancy batching must neither interpolate nor count failed covers."""
import unittest
from types import SimpleNamespace
from unittest.mock import patch
from flint import arb

from mass_verify import requested_occupancies,run_covers


class MassBatch(unittest.TestCase):
    def test_requested_domains(self):
        self.assertEqual(requested_occupancies(64,None,False),(64,))
        self.assertEqual(requested_occupancies(64,[59,60,63],True),(59,60,63))
        for groups,qs,full in ((64,[59,59],True),(64,[59],False),(64,[],True),
                               (0,None,True),(True,None,True),(64,[2049],True)):
            with self.assertRaises(ValueError): requested_occupancies(groups,qs,full)

    def test_each_cover_is_called_and_failure_stays_missing(self):
        seen=[]
        def fake_cover(args,operators,counts,terminal,shells,**options):
            seen.append(args.groups)
            self.assertTrue(options['prefix_rank'])
            return None if args.groups==60 else arb(2)**-80
        with patch('mass_verify.cover',side_effect=fake_cover),patch('builtins.print'):
            results=run_covers(SimpleNamespace(groups=64),(59,60,63),{}, {}, {})
        self.assertEqual(seen,[59,60,63])
        self.assertEqual(set(results),{59,60,63})
        self.assertIsNone(results[60])
        self.assertTrue(results[59]>0 and results[63]>0)


if __name__=='__main__':
    unittest.main()
