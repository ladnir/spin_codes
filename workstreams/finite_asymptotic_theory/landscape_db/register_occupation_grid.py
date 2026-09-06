"""Append verified occupation batches to the explicit source catalog."""
import argparse
import json
from pathlib import Path

import run_occupation_grid as producer


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    args=parser.parse_args()
    directory=args.directory.resolve()
    directory.relative_to(producer.grid.HERE)
    path=producer.grid.HERE/'catalog.json'
    catalog=json.loads(path.read_text())
    known={s['path'] for s in catalog['activation_sources']}
    count=0
    for manifest in sorted(directory.glob('t*/manifest.json')):
        receipt=json.loads(manifest.read_text())
        producer.verify(manifest.parent,receipt['arguments'])
        csv_path=(manifest.parent/'occupations.csv').relative_to(producer.grid.HERE).as_posix()
        if csv_path in known: continue
        catalog['activation_sources'].append(dict(path=csv_path,base='landscape_db',
            manifest=manifest.relative_to(producer.grid.HERE).as_posix(),
            study=directory.name+'_'+manifest.parent.name,map_tag='nested-pilot',
            transfer_review_status='activation_aware'))
        count+=1
    if not list(directory.glob('t*/manifest.json')):
        raise ValueError('no completed occupation batches')
    path.write_text(json.dumps(catalog,indent=2)+'\n')
    print(f'appended {count} verified occupation sources')


if __name__=='__main__':
    main()
