import copy
import unittest
import shared_relaxed_warm_start as w


class WarmStartTests(unittest.TestCase):
    def records(self):
        source = dict(schema=w.proof.DENSE_SCHEMA, updates=2, K=1 << 20, N=1 << 21,
            minimum_groups=33, zero_bits=64, variance_bins=16, base_tilt='3/16', cost_tilt='1/4',
            mixture=[dict(mass='8', activity='1/2')],
            results=[dict(distance='3/50', threshold=125829, root=['0', '1'], cover=dict(
                leaves={'':dict(witness=dict(parameters=['1/100','0','0']), upper=[1,-999])}, unresolved={}))])
        comparison = copy.deepcopy(source)
        comparison['mixture'][0]['mass'] = '4'
        comparison['results'][0].update(distance='7/100', threshold=146800)
        comparison['results'][0].pop('cover')
        return source, comparison

    def test_all_cells_become_unresolved_and_no_numeric_bounds_transfer(self):
        source, comparison = self.records()
        originals = copy.deepcopy((source, comparison))
        result = w.seed(source, comparison)
        self.assertEqual((source, comparison), originals)
        cover = result['results'][0]['cover']
        self.assertEqual(cover['leaves'], {})
        self.assertEqual(set(cover['unresolved']), {''})
        self.assertNotIn('upper', cover['unresolved'][''])
        self.assertEqual(cover['unresolved']['']['witness'], source['results'][0]['cover']['leaves']['']['witness'])
        self.assertEqual(result['mixture'], comparison['mixture'])
        with self.assertRaises(ValueError): w.proof.validate_record(result, complete=True)

    def test_geometry_changes_and_mass_increases_rejected(self):
        for key, value in (('minimum_groups', 65), ('base_tilt', '1/8'), ('variance_bins', 32)):
            source, comparison = self.records(); comparison[key] = value
            with self.assertRaises(ValueError): w.seed(source, comparison)
        for mixture in ([dict(mass='9', activity='1/2')], [dict(mass='4', activity='1/3')],
                        [dict(mass='4', activity='1/2')]*2):
            source, comparison = self.records(); comparison['mixture'] = mixture
            with self.assertRaises(ValueError): w.seed(source, comparison)
        source, comparison = self.records(); comparison['results'][0]['root'] = ['1/10','1']
        with self.assertRaises(ValueError): w.seed(source, comparison)

    def test_nonuniform_partition_preserved_but_every_leaf_requires_replay(self):
        source, comparison = self.records()
        source['results'][0]['root'] = comparison['results'][0]['root'] = ['1/4', '3/4']
        template = source['results'][0]['cover']['leaves']['']
        source['results'][0]['cover']['leaves'] = {
            path: copy.deepcopy(template) for path in ('0', '10', '110', '111')}
        result = w.seed(source, comparison)
        expected = {'0': ['1/4', '1/2'], '10': ['1/2', '5/8'],
                    '110': ['5/8', '11/16'], '111': ['11/16', '3/4']}
        cells = result['results'][0]['cover']['unresolved']
        self.assertEqual({path: row['cell'] for path, row in cells.items()}, expected)
        self.assertTrue(all(set(row) == {'cell', 'witness'} for row in cells.values()))
        cells['0']['witness']['parameters'][0] = '1/200'
        self.assertEqual(source['results'][0]['cover']['leaves']['0']['witness']['parameters'][0], '1/100')

    def test_missing_or_overlapping_source_cells_rejected(self):
        for paths in (('0',), ('0', '00', '1'), ('0', '10')):
            source, comparison = self.records()
            leaf = source['results'][0]['cover']['leaves']['']
            source['results'][0]['cover']['leaves'] = {p: copy.deepcopy(leaf) for p in paths}
            with self.assertRaises(ValueError): w.seed(source, comparison)

    def test_append_points_preserves_partial_partition_and_ignores_bounds(self):
        source, _ = self.records()
        source['results'][0]['cover']['unresolved'] = {'1': {}}
        source['results'][0]['cover']['leaves']['0'] = source['results'][0]['cover']['leaves'].pop('')
        comparison = copy.deepcopy(source)
        comparison['results'][0].pop('cover')
        comparison['results'][0]['probes'] = [dict(mean='1/3', witness=dict(parameters=['1/100','0','0']),
                                                   upper=[1,-999999], proposal=-999999.)]
        result = w.append_points(source, comparison)
        self.assertEqual(result['results'][0]['cover'], source['results'][0]['cover'])
        self.assertEqual(set(result['results'][0]['probes'][0]), {'mean', 'witness'})
        self.assertNotIn('probes', source['results'][0])
        changed = copy.deepcopy(comparison); changed['results'][0]['threshold'] += 1
        with self.assertRaises(ValueError): w.append_points(source, changed)
        changed = copy.deepcopy(comparison); changed['mixture'][0]['mass'] = '9'
        with self.assertRaises(ValueError): w.append_points(source, changed)
        with self.assertRaises(ValueError): w.append_points(result, comparison)


if __name__ == '__main__': unittest.main()
