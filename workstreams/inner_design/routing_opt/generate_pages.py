"""Optional Linux huge-page advice on owned workspace; no system setting changes."""
import hashlib
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
PARENT=HERE.parent
SOURCE=PARENT/'generated/asymmetric_greedy3_2_sparse/CandidateSpin.cpp'
variants=[]
for mode in ('tile','both','bothwrite'):
    name=f'routeopt_pages_{mode}'
    source=SOURCE if mode!='bothwrite' else PARENT/'generated/routeopt_scalarwrite32/CandidateSpin.cpp'
    code=source.read_text().replace('#include <algorithm>','''#include <algorithm>
#ifdef __linux__
#include <sys/mman.h>
#include <unistd.h>
#include <cstdio>
#include <cerrno>
#endif''')
    helper='''namespace {
void adviseWorkspace(block* p,std::size_t count) {
#ifdef __linux__
    const auto page=static_cast<std::uintptr_t>(sysconf(_SC_PAGESIZE));
    if(!page || page==static_cast<std::uintptr_t>(-1)) return;
    // Only whole pages contained in this owned allocation, never neighboring objects.
    const auto begin=(reinterpret_cast<std::uintptr_t>(p)+page-1)&~(page-1);
    const auto end=(reinterpret_cast<std::uintptr_t>(p)+count*sizeof(block))&~(page-1);
    if(end<=begin) return;
    const int hint=madvise(reinterpret_cast<void*>(begin),end-begin,MADV_HUGEPAGE);
    const int hintError=hint?errno:0;
    int collapse=-1,collapseError=0;
#ifdef MADV_COLLAPSE
    collapse=madvise(reinterpret_cast<void*>(begin),end-begin,MADV_COLLAPSE);
    collapseError=collapse?errno:0;
#endif
    std::fprintf(stderr,"WORKSPACE_PAGES bytes=%zu hint=%d errno=%d collapse=%d errno=%d\\n",
        static_cast<std::size_t>(end-begin),hint,hintError,collapse,collapseError);
#else
    (void)p;(void)count;
#endif
}
}
'''
    old='Spin::Workspace::Workspace(const Spin& code):buckets(code.codeBlocks()),tile(code.tileBlocks()) {}'
    body='adviseWorkspace(tile.data(),tile.size());'
    if mode.startswith('both'):body+='adviseWorkspace(buckets.data(),buckets.size());'
    assert old in code
    code=code.replace(old,helper+'Spin::Workspace::Workspace(const Spin& code):buckets(code.codeBlocks()),tile(code.tileBlocks()) {'+body+'}')
    code=code.replace('"t128_s19_asymmetric_greedy3_2_r1"',f'"{name}"')
    code=code.replace('"routeopt_scalarwrite32"',f'"{name}"')
    directory=PARENT/'generated'/name
    directory.mkdir(exist_ok=True)
    (directory/'CandidateSpin.cpp').write_text(code,encoding='utf-8',newline='\n')
    (directory/'AsymmetricMap.h').write_text((SOURCE.parent/'AsymmetricMap.h').read_text(),encoding='utf-8',newline='\n')
    variants.append(dict(name=name,sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in directory.iterdir() if p.is_file()}))
(HERE/'PAGES.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),candidates=variants),indent=2)+'\n')
