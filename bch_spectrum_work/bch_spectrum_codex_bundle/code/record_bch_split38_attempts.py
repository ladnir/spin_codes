"""Record terminal, inconclusive solver attempts observed by the root agent."""
from pathlib import Path
from run_higher_endpoint_preflight import sha,write_new
ROOT=Path(__file__).resolve().parents[1]
paths=[ROOT/name for name in (
    'code/prepare_bch_split38_probe.py','code/wrap_bch_split38_lp.py',
    'code/prepare_bch_split38_local_probe.py','code/verify_bch_split38_solution.py',
    'generated/split38_probe/model.json','generated/split38_probe/h38_wrapped.lp',
    'generated/split38_local_probe/model.json','generated/split38_local_probe/h38.lp',
    'generated/bch256_wambach_shortening.json','generated/bch_split38_identity_toy.json')]
solutions=[ROOT/'generated/split38_probe/h38_wrapped.sol',ROOT/'generated/split38_probe/h38_240s.sol',
           ROOT/'generated/split38_local_probe/h38_512.sol']
assert not any(path.exists() for path in solutions)
write_new(ROOT/'generated/bch_split38_attempts.json',dict(
    classification='Observed solver failures/time limits, not mathematical infeasibility certificates',
    observations=[dict(model='split38',precision_bits=256,time_limit_seconds=50,tool_session=14972,exit_code=1,result='TIME_LIMIT_REACHED'),
                  dict(model='split38',precision_bits=256,time_limit_seconds=240,tool_session=63869,exit_code=1,result='TIME_LIMIT_REACHED'),
                  dict(model='split38_local',precision_bits=512,time_limit_seconds=180,tool_session=17294,exit_code=1,result='TIME_LIMIT_REACHED')],
    no_solution_files_present=True,new_A38_cap_proved=False,
    preliminary_parse_issue='Long lines in the first unwrapped LP exceeded the reader limit; wrapping preserved the token stream',
    exact_projection_and_shortening_evidence='generated/bch256_wambach_shortening.json',
    all_solver_processes_terminal=True,source_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths+[Path(__file__)]}))
print('Recorded three terminal, inconclusive attempts; no new cap')
