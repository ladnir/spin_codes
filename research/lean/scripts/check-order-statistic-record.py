"""Record the unconditional canonical order-statistic integral proof."""
import hashlib
import json
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'scripts/map_data'
STEMS=['ConcreteFixedOrderStatisticVolume','ConcreteFixedOrderStatisticSplit',
       'ConcreteFixedOrderStatistic']

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    records={}
    audits=[]
    for stem in STEMS:
        src=ROOT/f'SpinCodes/Structured/{stem}.lean'
        obj=ROOT/f'.lake/build/lib/lean/SpinCodes/Structured/{stem}.olean'
        log=DATA/f'peach_{stem}.log'
        assert obj.exists() and obj.stat().st_mtime>=src.stat().st_mtime,stem
        output=log.read_text(encoding='utf-8')
        assert ': error:' not in output,stem
        assert not re.search(r'\b(sorry|admit|axiom|native_decide|trustCompiler)\b',src.read_text(encoding='utf-8')),stem
        rows=re.findall(r"'([^']+)' depends on axioms: \[([^]]*)\]",output)
        assert all(set(a.strip() for a in ax.split(','))<={'propext','Classical.choice','Quot.sound'} for _,ax in rows),stem
        audits.extend(t for t,_ in rows)
        records[src.relative_to(ROOT).as_posix()]={'source_sha256':sha(src),'object_sha256':sha(obj),'log':log.name}
    assert len(audits)==2,len(audits)
    (DATA/'fixed_order_statistic_verification.json').write_text(json.dumps({
        'status':'PASS','scope':'Three modules individually compiled on Peach against checked cached dependencies. '
        'Exact scaled ordered-simplex volume, measure-preserving coordinate split and Fubini proof of the '
        'canonical order-statistic integral for every continuous real F. No probabilistic spacing premise.',
        'modules':records,'axiom_audits':len(audits),'axiom_theorems':audits,
    },indent=2)+'\n',encoding='utf-8')
    print('PASS: canonical order-statistic integral; three modules, two standard-only audits.')

if __name__=='__main__':
    main()
