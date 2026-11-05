"""Fresh serial256/384 replay of the frozen final-gap box proposal."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from hd_shifted_swap_ledger import replay_swap_plan


def main():
    path = Path(__file__).with_name('hd_shifted_swap_4068_4095_frozen_plan3.json')
    raw = path.read_bytes()
    proposal_hash = sha256(raw).hexdigest()
    plan = json.loads(raw)
    signatures = []
    for precision in (256, 384):
        output = path.with_name(f'hd_shifted_swap_4068_4095_directed{precision}.json')
        if output.exists():
            raise FileExistsError(output)
        ledger, result = replay_swap_plan(plan, precision=precision, threshold=-60,
            progress=lambda row: print(json.dumps(dict(row, precision=precision)), flush=True))
        if result['failed'] or not result['summary']['complete_local_cover']:
            raise ArithmeticError('The complete final-gap proposal did not replay')
        result.update(kind='swap-plan', precision=precision, construction=ledger.construction,
                      accepted=ledger.accepted, proposal_sha256=proposal_hash,
                      fresh_arithmetic=True, whole_code_certificate=False)
        signatures.append(sorted((tuple(map(tuple, row['ranges'])), row['witness_index'],
                                  F(row['log2_upper_rational'])) for row in ledger.accepted))
        with output.open('x', encoding='utf-8') as target:
            json.dump(result, target, separators=(',', ':'), allow_nan=False)
            target.write('\n')
        print(json.dumps(dict(path=str(output), sha256=sha256(output.read_bytes()).hexdigest(),
                              precision=precision, summary=result['summary'], failed=0)), flush=True)
    if signatures[0] != signatures[1]:
        raise ArithmeticError('Dual replay leaf ceilings disagree')
    print(json.dumps(dict(dual_leaf_ceiling_agreement=True, leaves=len(signatures[0]),
                          proposal_sha256=proposal_hash)), flush=True)


if __name__ == '__main__':
    main()
