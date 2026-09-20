"""Static scientific figure: certified frontier versus reference, and backtests."""
import argparse
import io
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import bridge as base


def run(directory):
    directory.mkdir(exist_ok=True)
    for suffix in ('svg','png'):assert not (directory/f'curve_assessment.{suffix}').exists()
    data=base.read(base.HERE/'generated/curve_spectrum_calibration_v1.json')
    curve=base.read(base.HERE/'generated/frontier_curve_v3.json')
    full=base.read(base.HERE/'generated/curve_k28_full_retained_v1.json')
    fig,(ax,bx)=plt.subplots(1,2,figsize=(12,4.4),layout='constrained')
    rows=[r for r in curve['rows'] if r['full_certificate_margin_bits'] is not None]
    x=[r['message_exponent'] for r in rows]
    certified=[full['margin_bits'] if r['message_exponent']==28 else r['full_certificate_margin_bits'] for r in rows]
    ax.plot(x,[r['modeled_q1_bits'] for r in rows],'o--',color='#8052a4',label='Reference-spectrum Q1 (model)')
    ax.plot(x,certified,'s-',color='#176b87',label='Full distance certificate')
    ax.scatter([28],[data['conditional_endpoint']['cases'][0]['conditional_margin_bits']],marker='D',
        facecolors='white',edgecolors='#8052a4',zorder=4,label='Full bound if weighted score ≤ reference')
    ax.axhline(40,color='#777777',linewidth=.8,linestyle=':')
    ax.set(xlabel='log₂ K (message bits)',ylabel='Margin (bits)',title='BCH-256 × selected RM2Sub t64/s20',xticks=x,ylim=(38,74))
    ax.legend(fontsize=8,loc='upper right');ax.grid(axis='y',alpha=.2)
    known=[r for r in data['rows'] if r['message_exponent']==20 and r['exact_spectrum_q1_bits'] is not None]
    bx.bar([str(r['block_bits']) for r in known],[r['weighted_inflation_bits'] for r in known],color='#176b87',width=.55)
    bx.axhline(0,color='#555555',linewidth=.8)
    bx.set(xlabel='Known constituent length B',ylabel='Model margin − exact-spectrum margin (bits)',
        title='Unfitted spectrum backtest at K = 2²⁰',ylim=(-.65,1.15))
    for i,r in enumerate(known):
        v=r['weighted_inflation_bits'];bx.text(i,v+(.04 if v>=0 else -.08),f'{v:+.3f}',ha='center',fontsize=10)
    bx.text(.03,.03,'Positive = optimistic model\nB64 includes odd weights',transform=bx.transAxes,fontsize=9)
    bx.grid(axis='y',alpha=.2)
    fig.suptitle('Model agreement is evidence, not a confidence interval',fontsize=12)
    svg=io.StringIO();fig.savefig(svg,format='svg')
    with (directory/'curve_assessment.svg').open('x',encoding='utf-8',newline='\n') as handle:
        handle.write('\n'.join(line.rstrip() for line in svg.getvalue().splitlines())+'\n')
    fig.savefig(directory/'curve_assessment.png',dpi=160);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--directory',type=Path,required=True)
    a=p.parse_args();run(a.directory)
