"""Sensitivity probe for a new literature input, NOT yet an application certificate.

The source table was returned by web search, but the primary PDF page has not
yet been inspected. Keep the resulting objective provisional until that audit.
"""
import json
from pathlib import Path
from prepare_bch_hull_dual18_probe import build as previous_build
from audit_bch_closure_envelope import kraw_table
from prepare_bch_hull_probe import export
from bch_hull_cut_iteration import warm_basis,solve
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/schaubplus_provisional'
PREVIOUS=ROOT/'generated/oa21_hull_dual18_probe'


def build():
    model,scales=previous_build()
    kt=kraw_table()
    def coeffs(degree,multipliers):
        return {f'{a}_{w}':str(mult*(kt[degree][w]+(kt[degree][256-w] if w!=128 else 0)))
                for a,mult in multipliers for w in range(0,129,2)}
    for prefix in ('q','h'):
        model['constraints'].append(dict(name='PROVISIONAL_SchaubPlus_Qdual24_'+prefix,
            coeffs=coeffs(22,((prefix,1),)),sense='eq',rhs='0'))
    model['constraints'].append(dict(name='PROVISIONAL_SchaubPlus_Pdual26',
        coeffs=coeffs(24,(('q',1),('h',255))),sense='eq',rhs='0'))
    model['metadata']['classification']='PROVISIONAL literature-input sensitivity probe; NOT a certified BCH bound'
    model['metadata']['provisional_input']={
        'authors':'Federico Ponchio and Massimiliano Sala',
        'title':'A lower bound on the distance of cyclic codes',
        'date_on_indexed_manuscript':'2003-02-25',
        'source_url':'https://citeseerx.ist.psu.edu/document?doi=f4cc68b3184f5a89579bf7375f463e0ae7d4e4ea&repid=rep1&type=pdf',
        'table':1,'manuscript_page':16,'column':'d-prime SchaubPlus',
        'Q0_dimension':123,'Q0_designed_distance':39,'Q0_dual_distance_claimed_lower':24,
        'P0_dimension':131,'P0_designed_distance':37,'P0_dual_distance_claimed_lower':26,
        'primary_PDF_page_visually_verified':False,'rank_algorithm_independently_reproduced':False,
        'extension_argument':'Translate a zero coordinate of a putative low-weight extended-dual word to the parity position, then puncture, as in BCH_OA21_REFINEMENT.md',
        'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),ROOT/'code/prepare_bch_hull_dual18_probe.py')}}
    return model,scales


if __name__=='__main__':
    model,scales=build()
    lp,count=export(model,scales)
    FOLDER.mkdir(exist_ok=False)
    write_new(FOLDER/'model.json',model)
    write_new(FOLDER/'scales.json',scales)
    old=json.loads((PREVIOUS/'model.json').read_text())
    basis,mapping=warm_basis(old,scales,model,scales,PREVIOUS/'h_38.bas')
    for name,text in (('h_38.lp',lp),('warm.bas',basis)):
        with (FOLDER/name).open('x') as f:
            f.write(text)
    write_new(FOLDER/'warm_mapping.json',mapping)
    print(json.dumps(dict(classification='PROVISIONAL ONLY',variables=len(scales),rows=count)),flush=True)
    print(solve(FOLDER,90,True),flush=True)
