import unittest

import extrapolate_parameters


class ParameterExtrapolationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.payload = extrapolate_parameters.build(extrapolate_parameters.DEFAULT_DB)

    def test_fit_is_bound_to_expected_diagnostic_slice(self):
        fit = self.payload["fit"]
        self.assertEqual(fit["observation_count"], 15)
        self.assertLess(fit["rmse_bits"], 0.35)
        self.assertLess(fit["maximum_absolute_residual_bits"], 0.87)

    def test_certificate_anchored_staircase_starts_at_known_point(self):
        first = self.payload["rm2sub_staircases"]["certificate_anchored_offset_4"][0]
        self.assertEqual(
            (first["message_exponent"], first["step_bits"], first["state_bits"]),
            (16, 64, 14),
        )

    def test_cross_validation_covers_every_source_row(self):
        cross_validation = self.payload["cross_validation"]
        self.assertEqual(
            len(cross_validation["leave_one_message_exponent_out"]["predictions"]), 15
        )
        self.assertEqual(len(cross_validation["leave_one_block_size_out"]["predictions"]), 15)

    def test_staircase_respects_rm2_capacity(self):
        for schedule in self.payload["rm2sub_staircases"].values():
            for row in schedule:
                self.assertLessEqual(row["state_bits"], row["rm2_capacity"])


if __name__ == "__main__":
    unittest.main()
