"""Pack or restore reproducible compressed landscape snapshots."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3

HERE=Path(__file__).resolve().parent
FILES=('spin_landscape.sqlite3','landscape_export.csv')


def sha(path):
    digest=hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda:handle.read(1<<20),b''):digest.update(chunk)
    return digest.hexdigest()


def compress(source, destination):
    pending=destination.with_name(destination.name+'.pending')
    digest=hashlib.sha256();size=0
    try:
        with source.open('rb') as incoming,pending.open('wb') as raw:
            with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0,compresslevel=6) as outgoing:
                for chunk in iter(lambda:incoming.read(1<<20),b''):
                    digest.update(chunk);size+=len(chunk);outgoing.write(chunk)
        pending.replace(destination)
    finally:
        pending.unlink(missing_ok=True)
    return dict(file=source.name,plain_sha256=digest.hexdigest(),plain_bytes=size,
                compressed_file=destination.name,compressed_sha256=sha(destination),compressed_bytes=destination.stat().st_size)


def restore_one(directory, record):
    if record['file'] not in FILES or record['compressed_file']!=record['file']+'.gz':
        raise ValueError('unknown snapshot path')
    source=directory/record['compressed_file'];target=directory/record['file']
    if sha(source)!=record['compressed_sha256']:raise ValueError('compressed snapshot hash mismatch')
    if target.exists():
        if sha(target)!=record['plain_sha256']:raise ValueError(f'{target.name} already exists with different content; preserve it before restoring')
        return
    pending=target.with_name(target.name+'.pending')
    digest=hashlib.sha256();size=0
    try:
        with gzip.open(source,'rb') as incoming,pending.open('wb') as outgoing:
            for chunk in iter(lambda:incoming.read(1<<20),b''):
                digest.update(chunk);size+=len(chunk);outgoing.write(chunk)
        if size!=record['plain_bytes'] or digest.hexdigest()!=record['plain_sha256']:
            raise ValueError('expanded snapshot hash mismatch')
        pending.replace(target)
    finally:
        pending.unlink(missing_ok=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=('pack','restore'))
    args=parser.parse_args();manifest=HERE/'landscape_snapshot.json'
    if args.action=='pack':
        db=sqlite3.connect(f'{(HERE/FILES[0]).as_uri()}?mode=ro',uri=True)
        try:
            if db.execute('PRAGMA quick_check').fetchone()[0]!='ok':raise ValueError('database quick check failed')
            if db.execute("SELECT value FROM schema_metadata WHERE key='catalog_sha256'").fetchone()[0]!=sha(HERE/'catalog.json'):
                raise ValueError('database was built from a different catalog')
            rows=db.execute('SELECT COUNT(*) FROM landscape').fetchone()[0]
        finally:db.close()
        with (HERE/FILES[1]).open(newline='',encoding='utf-8') as handle:
            if sum(1 for _ in csv.DictReader(handle))!=rows:raise ValueError('database/export row counts differ')
        records=[compress(HERE/name,HERE/(name+'.gz')) for name in FILES]
        payload=dict(schema='finite-landscape-compressed-snapshot-v1',row_count=rows,
                     catalog_sha256=sha(HERE/'catalog.json'),files=records)
        manifest.write_text(json.dumps(payload,indent=2)+'\n')
    else:
        payload=json.loads(manifest.read_text())
        if payload['schema']!='finite-landscape-compressed-snapshot-v1':raise ValueError('unknown snapshot schema')
        if payload['catalog_sha256']!=sha(HERE/'catalog.json'):raise ValueError('snapshot catalog differs from current catalog; rebuild instead')
        if {r['file'] for r in payload['files']}!=set(FILES):raise ValueError('incomplete snapshot')
        for record in payload['files']:restore_one(HERE,record)
    print(json.dumps(payload,indent=2))


if __name__=='__main__':main()
