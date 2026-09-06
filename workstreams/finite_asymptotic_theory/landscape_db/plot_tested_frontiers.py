"""Plot tested state frontiers while retaining their occupation scope."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

HERE=Path(__file__).resolve().parent


def main():
    with (HERE/'tested_parameter_frontiers.csv').open(newline='') as handle:rows=list(csv.DictReader(handle))
    selected=[('extended BCH [128,64,22] exact','BCH [128,64,22]'),
              ('RM(4,9) [512,256,32] exact','RM(4,9) [512,256,32]'),
              ('random full-rank [512,256] reused','Random [512,256]')]
    colors={64:'#176b87',128:'#b65325',256:'#7451a8'}
    fig,axes=plt.subplots(2,3,figsize=(12,7),sharex=True,sharey=True)
    exponents=[16,18,20,22,24]
    for col,(series,title) in enumerate(selected):
        for ri,target in enumerate((20,40)):
            ax=axes[ri,col];present=False
            for t,color in colors.items():
                for coverage,style in (('q1_q4','--'),('q1_q64','-')):
                    matching={int(r['message_exponent']):r for r in rows if r['series']==series
                              and int(r['step_bits'])==t and int(r['target_margin_bits'])==target
                              and r['occupation_coverage']==coverage}
                    values=[float(matching[e]['smallest_tested_passing_state'])
                            if matching[e]['smallest_tested_passing_state'] else np.nan for e in exponents]
                    if np.isfinite(values).any():
                        present=True;ax.plot(exponents,values,color=color,linestyle=style,marker='o',markersize=4)
                full=[r for r in rows if r['series']==series and int(r['step_bits'])==t
                      and int(r['target_margin_bits'])==target and r['occupation_coverage']=='full'
                      and r['smallest_tested_passing_state']]
                if full:
                    present=True;ax.scatter([int(r['message_exponent']) for r in full],
                        [int(r['smallest_tested_passing_state']) for r in full],color=color,marker='*',s=160,zorder=5)
            ax.set_title(f'{title}\n{target}-bit diagnostic threshold',fontsize=10)
            if not present:ax.text(.5,.5,'No passing result\nin this snapshot',ha='center',va='center',transform=ax.transAxes,color='#666')
            ax.grid(alpha=.2);ax.set_xticks(exponents);ax.set_ylim(6.5,20.5)
            if ri==1:ax.set_xlabel('Message exponent log2(k)')
            if col==0:ax.set_ylabel('Smallest observed passing state s')
    handles=[Line2D([],[],color=c,label=f't={t}') for t,c in colors.items()]
    handles += [Line2D([],[],color='#333',linestyle='--',label='Q1..4 partial'),
                Line2D([],[],color='#333',label='Q1..64 partial'),
                Line2D([],[],color='#333',marker='*',linestyle='None',markersize=11,label='Full range diagnostic')]
    fig.legend(handles=handles,loc='lower center',ncol=6,bbox_to_anchor=(.5,.045),frameon=False,fontsize=9)
    fig.suptitle('Tested RM2Sub state choices by message length and epoch size',fontsize=14)
    fig.text(.5,.025,'Missing points may reflect absent coverage or a weak bound. All values are binary64 diagnostics; no runtime or outward-certificate claim.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.11,1,.95));fig.savefig(HERE/'tested_parameter_frontiers.png',dpi=160);plt.close(fig)


if __name__=='__main__':main()
