"""Figure 1: one-round IMT with fixed versus length-dependent step size."""
import argparse
import math
from pathlib import Path

import imt_results as evidence

GRID_PIN = ('finite_migration/ADAPTIVE_LENGTH_v1.json',
            'e5fbfaace653ef5597817c4241ff31789fa82d4b57a5b4426d824424aa747de0')
REPLAY_PIN = ('finite_migration/ADAPTIVE_LENGTH_VERIFIED_v1.json',
              'd5d5a536e8fc1ab86b6987c713281252199cc3744d73ec638c934553b16d1d17')
OUT = Path(__file__).resolve().parent / 'figures/imt_parameter_k_b.tex'
require = evidence.require
TOLERANCE = .1


def baseline_t(b):
    return 128 if b == 256 else 64


def exponents(b):
    return tuple(range(14,25,2) if b == 256 else range(12,27,2))


def validate(grid,replay,selected_inner=None):
    require(grid['status'] == 'BINARY64_IMT_ADAPTIVE_LENGTH_STUDY','Wrong length study')
    require(replay['status'] == 'VERIFIED_BINARY64_IMT_ADAPTIVE_LENGTH_STUDY','Wrong length replay')
    require(replay['producer_sha256'] == GRID_PIN[1],'Replay binds another grid')
    require(not grid['full_distance_proved'] and not replay['full_distance_proved'],'Wrong Q1 scope')
    expected = {(b,t,m,r) for b in (64,128,256) for m in exponents(b)
                for t in ((16,32,64,128) if b==256 else (8,16,32,64))
                for r in ([False,True] if t==baseline_t(b) else [False])}
    require(replay['cells_checked'] == replay['log_domain_cells_checked'] == len(expected),
            'Incomplete log-domain replay')
    require(replay['maps_checked'] == 12 and 0 <= replay['largest_difference'] < 1e-7,
            'Wrong map coverage or numerical disagreement')
    rows = {}
    for row in grid['cells']:
        b,t,m,r = key = row['b'],row['t'],row['exponent'],row['refresh']
        require(key not in rows and key in expected,'Duplicate or unexpected cell')
        require(row['sharp'] and not row['full_distance_proved'],'Wrong transfer or evidence scope')
        require(math.isfinite(row['q1_margin_bits']) and not row['dominant_grid_edge'],'Invalid margin/witness')
        length = 2**(m+1)//b
        require((row['outer_rows'],row['epochs_per_region']) == (length,length//t),'Wrong geometry')
        s = min(19 if b==256 else 20,(t.bit_length()-1)*t.bit_length()//2)
        record = grid['maps'][f'b{b}_t{t}']
        require(row['s'] == record['s'] == s and record['t'] == t,'Wrong state size')
        require(record['transvection_rounds'] == 1,'Extra mixing rounds in IMT curve')
        require(sum(record['spectrum'].values()) == 2**s-1 and str(t) not in record['spectrum'],
                'Wrong expansion spectrum')
        for field in ('expansion_columns','feedback_columns'):
            require(len(record[field]) == len(set(record[field])) == t and
                    all(0 < v < 2**s for v in record[field]),'Invalid fixed map')
        rows[key] = row
    require(set(rows) == expected and len(grid['maps'])==12,'Missing grid geometry')
    if selected_inner:
        record = grid['maps']['b256_t128']
        require(all(record[f] == selected_inner[f] for f in ('expansion_columns','feedback_columns')),
                'BCH-256 baseline differs from implemented maps')
    return rows


def load(selected=None):
    seen = {}
    inner = selected['half'][16]['instance']['inner'] if selected else None
    return validate(evidence.authenticate(GRID_PIN,seen),evidence.authenticate(REPLAY_PIN,seen),inner)


def adaptive(rows,b,m):
    choices = [r for (bb,t,mm,refresh),r in rows.items() if (bb,mm,refresh)==(b,m,False)]
    best = max(r['q1_margin_bits'] for r in choices)
    return max((r for r in choices if r['q1_margin_bits'] >= best-TOLERANCE),key=lambda r:r['t'])


def panel(rows,b,wide=False):
    points = exponents(b)
    fixed = [rows[b,baseline_t(b),m,False] for m in points]
    tuned = [adaptive(rows,b,m) for m in points]
    ideal = [rows[b,baseline_t(b),m,True] for m in points]
    options = (r'width=.95\linewidth,height=5.4cm,' if wide else r'width=\linewidth,height=5.0cm,')
    options += (r'title={BCH-'+str(b)+r'},xlabel={$\log_2 K$},ylabel={Q1 margin (bits)},'
                r'tick label style={font=\footnotesize},label style={font=\small},title style={font=\small},'
                r'grid=major,grid style={black!10},'
                +f'xmin={points[0]-.4},xmax={points[-1]+.4},'
                +('ymin=28,ymax=58,xtick={14,16,18,20,22,24},' if wide else
                  'ymin=0,ymax='+('20' if b==64 else '42')+',xtick={12,16,20,24,26},'))
    if wide:
        options += r'legend style={font=\footnotesize,draw=none,at={(.5,1.17)},anchor=south,legend columns=3},'
    text = '\\begin{tikzpicture}\n\\begin{axis}['+options+']\n'
    for records,style,label in ((fixed,'blue!70!black,dashed,thick,mark=o','Fixed step'),
                                (tuned,'orange!85!black,thick,mark=*','Adaptive step'),
                                (ideal,'black,densely dotted,thick,mark=none','Full refresh (fixed maps)')):
        coords = ' '.join(f"({r['exponent']},{r['q1_margin_bits']:.10f})" for r in records)
        text += f'\\addplot[{style}] coordinates {{{coords}}};\n'
        if wide: text += '\\addlegendentry{'+label+'}\n'
    previous = None
    for row in tuned:
        if wide and row['t'] != previous:
            anchor = 'south west' if row['exponent']==points[0] else 'south'
            text += (f"\\node[font=\\scriptsize,anchor={anchor},fill=white,inner sep=1pt,yshift=6pt] at (axis cs:{row['exponent']},"
                     f"{row['q1_margin_bits']:.10f}) {{$t={row['t']}$}};\n")
        previous = row['t']
    return text+'\\end{axis}\n\\end{tikzpicture}\n'


def render(rows):
    text = '% Generated by paper/build_imt_length_figure.py; do not hand edit.\n'
    text += '% Producer SHA-256: '+GRID_PIN[1]+'\n\\begin{figure}[!tbp]\n\\centering\n'
    text += panel(rows,256,True)+'\\par\\medskip\n'
    for b in (64,128):
        text += '\\begin{minipage}[t]{.485\\linewidth}\n\\centering\n'+panel(rows,b)+'\\end{minipage}\n'
        if b==64: text += '\\hfill\n'
    return text+r'''\caption{Length dependence with one transvection per IMT update.
The fixed curves use $(t,s)=(128,19)$ for BCH-256 and $(64,20)$ for the
smaller outers. The adaptive curves choose the largest tested step within
$0.1$ bits of the best tested Q1 margin. Top labels mark step changes;
the smaller outers use $t=16,32,32,64,\ldots$ at successive markers.
The full-refresh references retain the fixed curves' maps.
All points are replayed binary64 Q1 evaluations, not full certificates;
segments only connect tested lengths.}
\label{fig:finite-k-b}
\end{figure}
'''


def check(selected=None):
    rows = load(selected)
    require(OUT.read_text() == render(rows),'Stale adaptive length figure')
    return dict(cells=110,maps=12,transvections_per_step=1,full_distance_proved=False)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--check',action='store_true')
    a = p.parse_args()
    if a.check: print(check())
    else: OUT.write_text(render(load()),encoding='utf-8',newline='\n')
