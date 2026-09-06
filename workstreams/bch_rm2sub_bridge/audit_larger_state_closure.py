"""Final deterministic audit: full coverage, exact outer caps, and arithmetic tests.

No optimization or encoder benchmark is run. All subprocess checks are serial.
The producer/replay receipts for the inner bounds must already exist.
"""
import math
import subprocess
import sys
from pathlib import Path
from fractions import Fraction as F
import bridge as base
import larger_state_maps as maps
import exact_joint_shell_caps as exact
import christoffel_caps as christoffel
from bch_joint_shell_objective import audit as joint_audit


def run():
    output=base.HERE/'generated/larger_t64_s20_full_closure_audit.json'
    assert not output.exists()
    ledger_path=base.HERE/'generated/larger_coverage_full.json'
    ledger=base.read(ledger_path)
    assert ledger['configuration']=='t64_s20'
    assert ledger['covered_intervals']==[[1,8192]] and ledger['covered_occupancies']==8192
    assert ledger['all_occupations_certified'] and ledger['full_target_proved']
    hashed={str(ledger_path.relative_to(base.HERE)):base.sha(ledger_path)}
    for name,digest in ledger['local_sha256'].items():
        assert base.sha(base.HERE/name)==digest;hashed[name]=digest
    best={}
    for name in ledger['local_sha256']:
        if not name.endswith('_outward.json'):continue
        receipt=base.read(base.HERE/name)
        assert receipt['configuration']=='t64_s20'
        assert receipt['parameters']==dict(message_bits=1<<20,output_bits=1<<21,step_bits=64,state_bits=20,distance_cutoff=209716)
        for source,digest in receipt['local_sha256'].items():
            assert base.sha(base.HERE/source)==digest;hashed[source]=digest
        if 'Q1_upper' in receipt:
            best[1]=base.decode(receipt['Q1_upper'])
            for source,digest in receipt['outer_sha256'].items():assert base.sha(base.BCH/source)==digest
            coefficients={int(w):base.decode(v) for w,v in receipt['coefficient_upper'].items()}
            assert base.bch_bound(coefficients)==tuple(base.decode(receipt[k]) for k in ('Q1_upper','dual_domination_factor','remaining_shell_upper'))
        else:
            lo,hi=receipt['occupancy_range']
            assert [r['occupation'] for r in receipt['rows']]==list(range(lo,hi+1))
            for row in receipt['rows']:
                q=row['occupation'];bound=base.decode(row['upper'])
                assert bound>0
                if bound<F(1,1<<60):best[q]=min(best.get(q,bound),bound)
    assert sorted(best)==list(range(1,8193))
    total=sum(best.values(),F(0))
    assert total==base.decode(ledger['covered_union_upper']) and total<F(1,1<<50)
    higher=total-best[1]
    print('All8192 occupancies reaggregated exactly; union <2^-50',flush=True)
    t,s,spectrum,kernel=maps.load('t64_s20')
    assert (t,s)==(64,20) and min(spectrum)==16
    assert min(w for w,c in kernel.items() if w and c)==8
    reference=base.BCH/'generated/shift_rank_oa29_joint'
    rebuilt=joint_audit('prepare_bch_shift_rank_oa29_probe',reference)
    assert rebuilt==base.read(reference/'audit.json')
    assert rebuilt['all_primal_dual_checks_passed'] and rebuilt['rows_checked']==1163
    print('Exact outer joint LP reconstructed:1163 rows',flush=True)
    model,scales=exact.model_build()
    assert model==base.read(reference/'model.json') and scales==base.read(reference/'scales.json')
    cap_files=sorted((base.HERE/'generated').glob('joint_shell_*/cap.json'))
    for index,path in enumerate(cap_files):
        saved=base.read(path)
        for source,digest in saved['local_sha256'].items():assert base.sha(base.HERE/source)==digest
        for source,digest in saved['outer_sha256'].items():assert base.sha(base.BCH/source)==digest
        rebuilt_cap=exact.audit(path.parent,saved['weight'],model,scales)
        assert all(saved[k]==v for k,v in rebuilt_cap.items())
        hashed[str(path.relative_to(base.HERE))]=base.sha(path)
        if index%10==0:print('Exact shell cap replay',index+1,'of',len(cap_files),flush=True)
    christoffel.build(verify=True)
    tests=[]
    for name in ('test_larger_state_constant_band.py','test_polynomial_regions.py','test_scaled_adaptive.py',
                 'test_shared_range.py','test_gap_dyadics.py','test_constant_box_envelope.py'):
        path=base.HERE/name
        result=subprocess.run([sys.executable,'-B',str(path)],capture_output=True,text=True,check=True)
        tests.append(dict(file=name,output=result.stdout+result.stderr))
        hashed[name]=base.sha(path);print('Passed',name,flush=True)
    for path in (Path(__file__),Path(exact.__file__),Path(christoffel.__file__)):
        hashed[str(path.relative_to(base.HERE))]=base.sha(path)
    base.write_new(output,dict(status='FULL_T64_S20_DISTANCE_CLOSURE_AUDIT_PASSED',
        parameters=dict(message_bits=1<<20,output_bits=1<<21,outer_rows=8192,outer_length=256,
                        inner_step_bits=64,inner_state_bits=20,distance_cutoff=209716),
        all_8192_occupancies_covered=True,full_first_moment_upper=base.encode(total),
        Q1_upper=base.encode(best[1]),Q2_through_Q8192_upper=base.encode(higher),
        full_failure_below_2_to_minus_50=True,requested_2_to_minus_40_target_closed=True,
        margin_bits_diagnostic=math.log2(total.denominator)-math.log2(total.numerator),
        exact_outer_joint_rows=1163,exact_shell_caps_replayed=len(cap_files),tests=tests,
        assumption_scope='Fixed snapshotted RM2Sub map; uniform independent row and region permutations; fresh independent nonzero GF(2^20) multipliers each epoch; zero initial state; output before update; no final flush',
        prg_replacement_established=False,full_spin_protocol_security_established=False,performance_claimed=False,
        local_sha256=hashed,outer_sha256={str(p.relative_to(base.BCH)):base.sha(p) for p in
            (reference/'audit.json',reference/'model.json',reference/'scales.json',reference/'objective.json')}))
    print('FULL CLOSURE AUDIT PASSED; margin',math.log2(total.denominator)-math.log2(total.numerator),flush=True)


if __name__=='__main__':run()
