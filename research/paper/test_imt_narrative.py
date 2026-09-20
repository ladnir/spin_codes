"""Keep the paper's current inner distinct from frozen evidence filenames."""
import unittest

from check_finite_integration import check_inner_terminology


class IMTNarrativeTests(unittest.TestCase):
    def test_old_names_rejected(self):
        for name in ('RM2Sub', 'MR2Sub', 'rm2-sub', 'RM2 Sub'):
            with self.subTest(name=name), self.assertRaises(ValueError):
                check_inner_terminology('The selected inner is ' + name)

    def test_frozen_url_preserved(self):
        check_inner_terminology(
            r'\artifactlink{https://example.org/single_sampled_ba_rm2sub/README.md}'
            r'{outer certificate snapshot}')

    def test_visible_link_label_checked(self):
        with self.assertRaises(ValueError):
            check_inner_terminology(r'\artifactlink{https://example.org/imt}{RM2Sub inner}')

    def test_reed_muller_math_preserved(self):
        check_inner_terminology(r'The IMT expansion uses $\operatorname{RM}(2,m)$.')


if __name__ == '__main__':
    unittest.main()
