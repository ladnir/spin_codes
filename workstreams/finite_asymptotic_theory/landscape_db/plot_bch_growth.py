"""Plot continuous BCH margins and their explanatory row-count adjustment."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
COLORS={8:'#787878',32:'#bf7720',64:'#226f96',128:'#643e91'}


def main():
    rows=list(csv.DictReader((HERE/'bch_growth.csv').open(newline='')))
    payload=json.loads((HERE/'bch_growth.json').read_text())
    if hashlib.sha256((HERE/'bch_growth.csv').read_bytes()).hexdigest()!=payload['csv_sha256']:
        raise ValueError('growth CSV does not match study receipt')
    fig,axes=plt.subplots(1,2,figsize=(12,4.8))
    for b,color in COLORS.items():
        chosen=sorted((r for r in rows if int(r['block_bits'])==b and r['step_bits']=='64' and r['state_bits']=='20'),
                      key=lambda r:int(r['message_exponent']))
        x=[int(r['message_exponent']) for r in chosen]
        axes[0].plot(x,[float(r['margin_bits']) for r in chosen],color=color,marker='o',markersize=3,label=f'BCH {b}')
        axes[1].plot(x,[float(r['intercept_bits']) for r in chosen],color=color,marker='o',markersize=3,label=f'BCH {b}')
        estimate=payload['models'][str(b)]['chernoff_intercept_bits']
        if estimate is not None:
            axes[1].axhline(estimate,color=color,linestyle='--',alpha=.65)
    axes[0].axhline(0,color='#bbb',linewidth=.8)
    axes[0].set_title('At fixed constituent size, growing k consumes margin')
    axes[0].set_ylabel('Q1 diagnostic margin (bits)')
    axes[1].set_title('Removing the row count exposes the BCH contribution')
    axes[1].set_ylabel('Q1 margin + log2(number of outer rows)')
    for ax in axes:
        ax.set_xlabel('Message size log2(k)');ax.grid(alpha=.2);ax.set_xticks(range(12,27,2))
    axes[0].legend(frameon=False,ncol=2)
    axes[1].text(.04,.64,'Dashed: persistent-state Chernoff model\n(BCH 64 and 128 only)',transform=axes[1].transAxes,va='top',fontsize=9)
    fig.suptitle('BCH growth at fixed t=64, s=20',fontsize=14)
    fig.text(.5,.015,'Q1 only; four-state uniform-refresh bounds. Models omit cancellations and finite-epoch fluctuations. No full-distance claim.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.05,1,.94));fig.savefig(HERE/'bch_growth_curves.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(10,7),sharex=True,sharey=True)
    colors={64:'#226f96',128:'#bf7720',256:'#643e91'}
    for ax,b in zip(axes.flat,COLORS):
        for t,color in colors.items():
            chosen=sorted((r for r in rows if int(r['block_bits'])==b and int(r['step_bits'])==t and r['message_exponent']=='20'),
                          key=lambda r:int(r['state_bits']))
            ax.plot([int(r['state_bits']) for r in chosen],[float(r['intercept_bits']) for r in chosen],
                    color=color,marker='o',markersize=3,label=f't={t}')
        ax.set_title(f'BCH {b}');ax.grid(alpha=.2)
        estimate=payload['models'][str(b)]['chernoff_intercept_bits']
        if estimate is not None:ax.axhline(estimate,color='#777',linestyle='--',linewidth=1)
    for ax in axes[1]:ax.set_xlabel('State dimension s')
    for ax in axes[:,0]:ax.set_ylabel('Q1 margin + log2(outer rows)')
    axes[0,0].legend(frameon=False)
    fig.suptitle('State size improves the bound until the outer-code contribution dominates (k=2^20)',fontsize=12)
    fig.text(.5,.015,'The three t curves nearly coincide. Dashed: persistent-state model. Q1 diagnostics only; no full-distance claim.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.05,1,.95));fig.savefig(HERE/'bch_state_tradeoff.png',dpi=160);plt.close(fig)


if __name__=='__main__':main()
