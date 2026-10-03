from fractions import Fraction as F
from unittest import TestCase,main
import summarize_screen as summary


class SummaryTests(TestCase):
    def receipt(self,value,margin):
        return dict(whole_code_certificate=False,all_occupancies_checked=False,
            source_pins_verified_at_finish=True,results={'18':dict(geometry={'groups':1024},
                selected_endpoints={'2':summary.screen.cert.dyadic_record(value)},
                selected_raw_margins_display={'2':margin})})

    def test_minimum_never_claims_full_coverage(self):
        value=summary.merge_results([self.receipt(F(1,8),3),self.receipt(F(1,16),4)])['18']
        self.assertEqual(summary.screen.cert.read_dyadic(value['endpoints']['2']),F(1,16))
        self.assertEqual(value['untested_q_count'],1023)
        self.assertFalse(value['whole_code_certificate'])

    def test_reject_geometry_mismatch(self):
        a,b=self.receipt(F(1,8),3),self.receipt(F(1,16),4)
        b['results']['18']['geometry']['groups']=4096
        with self.assertRaises(ValueError):summary.merge_results([a,b])


if __name__=='__main__':main()
