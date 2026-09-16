import json
import unittest
import search_basis as s


def multiply(a,b):
    v=0
    for _ in range(19):
        if b&1:v^=a
        b>>=1;carry=a>>18;a=(a<<1)&((1<<19)-1)
        if carry:a^=0x27
    return v


class BasisTests(unittest.TestCase):
    def candidates(self):
        for path in s.HERE.glob('BASIS_SEARCH*.json'):
            data=json.loads(path.read_text())
            if 'monomials' in data: yield from data['candidates']

    def test_all_candidates(self):
        _,a,b,_=s.load()
        for r in self.candidates():
            aa,bb=s.audit(r['V'],a,b)
            self.assertEqual(aa,[int(x,16) for x in r['emission_rows_hex']])
            self.assertEqual(bb,[int(x,16) for x in r['feedback_monomials_hex']])

    def test_transition_bilinear_basis(self):
        # Linearity in both coefficient and state proves all coefficient/state pairs.
        for r in self.candidates():
            v=r['V'];inv=r['V_inverse']
            for k in range(19):
                transition=[multiply(1<<j,1<<k) for j in range(19)]
                transformed=s.compose(v,s.compose(transition,inv))
                for j in range(19):
                    self.assertEqual(s.apply(transformed,s.apply(v,1<<j)),s.apply(v,s.apply(transition,1<<j)))

    def test_feedback_and_emission_pairing(self):
        record,a,b,_=s.load()
        for r in self.candidates():
            aa,_=s.audit(r['V'],a,b)
            for p,col in enumerate(record['columns']):
                new_col=sum(((row>>p)&1)<<j for j,row in enumerate(aa))
                for j in range(19):
                    self.assertEqual((new_col&s.apply(r['V'],1<<j)).bit_count()%2,(col>>j)&1)


if __name__=='__main__':unittest.main()
