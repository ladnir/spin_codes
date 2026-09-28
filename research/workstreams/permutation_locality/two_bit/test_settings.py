import unittest

from model import resolve_updates
from full_cover import saved_updates,ENSEMBLE_PREFIX


class SettingsTests(unittest.TestCase):
    def test_update_count_is_part_of_the_ensemble(self):
        self.assertEqual(resolve_updates(None),2)
        self.assertEqual(resolve_updates(None,{}),2)
        self.assertEqual(resolve_updates(None,{'updates':3}),3)
        self.assertEqual(resolve_updates(3),3)
        self.assertEqual(resolve_updates(3,{'updates':3}),3)
        self.assertEqual(resolve_updates(3,{},retarget=True),3)
        self.assertEqual(resolve_updates(2,{'updates':3},retarget=True),2)
        for requested,record in ((3,{}),(2,{'updates':3})):
            with self.assertRaises(ValueError):resolve_updates(requested,record)
        with self.assertRaises(ValueError):resolve_updates(3,retarget=True)
        for value in (0,33,True,2.0,'2'):
            with self.assertRaises(ValueError):resolve_updates(value)
            with self.assertRaises(ValueError):resolve_updates(None,{'updates':value})

    def test_sparse_replay_checks_the_whole_ensemble(self):
        record={'schema':'two-bit-support-cover-1','ensemble':[*ENSEMBLE_PREFIX,3]}
        self.assertEqual(saved_updates(record),3)
        for changed in ([*ENSEMBLE_PREFIX,0],ENSEMBLE_PREFIX,[*ENSEMBLE_PREFIX,3,0],
                        [*ENSEMBLE_PREFIX[:-1],20,3],[*ENSEMBLE_PREFIX,True]):
            with self.assertRaises(ValueError):saved_updates({**record,'ensemble':changed})
        with self.assertRaises(ValueError):saved_updates({**record,'schema':'unknown'})


if __name__=='__main__':unittest.main()
