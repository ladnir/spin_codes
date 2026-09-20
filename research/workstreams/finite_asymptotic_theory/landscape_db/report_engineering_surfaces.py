"""Continuous surface summaries and scientific plots; all outputs stay local."""
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
FAMILIES=('bch','rm','random_mean')
LABELS={'bch':'BCH (exact spectrum)','rm':'RM (exact spectrum)',
        'random_mean':'Random (ensemble average)','random_caps60':'Random (conditional caps)'}
COLORS={'bch':'#643e91','rm':'#226f96','random_mean':'#bd721d'}


def read():
    path=HERE/'engineering_surfaces.csv'
    payload=json.loads((HERE/'engineering_surfaces.json').read_text())
    if hashlib.sha256(path.read_bytes()).hexdigest()!=payload['csv_sha256']:
        raise ValueError('surface CSV does not match receipt')
    with path.open(newline='') as handle:
        rows=list(csv.DictReader(handle))
    for row in rows:
        for key in ('block_bits','dimension','step_bits','state_bits','message_exponent','outer_rows',
                    'epochs_per_region','dominant_weight','witness_at_edge'):
            row[key]=int(row[key])
        for key in ('margin_bits','intercept_bits','dominant_log_tilt','dominant_share'):
            row[key]=float(row[key])
    return rows,payload


def select(rows,**where):
    return [r for r in rows if all(r[k]==v for k,v in where.items())]


def main():
    rows,payload=read();summary=[]
    for family in FAMILIES:
        blocks=sorted({r['block_bits'] for r in rows if r['family']==family})
        for b in blocks:
            item=dict(family=family,block_bits=b)
            anchor=select(rows,family=family,block_bits=b,step_bits=64,state_bits=20,message_exponent=20)[0]
            item.update({k:anchor[k] for k in ('margin_bits','intercept_bits','dominant_weight','dominant_share')})
            slopes={};spans={}
            for s in (10,12,16,20):
                selected=[r for r in select(rows,family=family,block_bits=b,step_bits=64,state_bits=s)
                          if 16<=r['message_exponent']<=24]
                selected.sort(key=lambda r:r['message_exponent'])
                slopes[s]=float(np.polyfit([r['message_exponent'] for r in selected],
                                           [r['margin_bits'] for r in selected],1)[0])
                spans[s]=max(r['intercept_bits'] for r in selected)-min(r['intercept_bits'] for r in selected)
            item['k_slopes_by_s']=slopes;item['adjusted_span_by_s']=spans
            state_rows=select(rows,family=family,block_bits=b,step_bits=64,message_exponent=20)
            item['smallest_s_within_bits_of_s20']={str(tolerance):min(r['state_bits'] for r in state_rows
                if r['margin_bits']>=anchor['margin_bits']-tolerance) for tolerance in (.25,1.)}
            late_slopes={};late_spans={}
            for s in (10,12,16,20):
                selected=[r for r in select(rows,family=family,block_bits=b,step_bits=64,state_bits=s)
                          if 20<=r['message_exponent']<=24]
                late_slopes[s]=float(np.polyfit([r['message_exponent'] for r in selected],
                                                [r['margin_bits'] for r in selected],1)[0])
                late_spans[s]=max(r['intercept_bits'] for r in selected)-min(r['intercept_bits'] for r in selected)
            item['late_k_slopes_by_s']=late_slopes;item['late_adjusted_span_by_s']=late_spans
            item['maximum_t_spread_bits']=max(
                max(r['margin_bits'] for r in select(rows,family=family,block_bits=b,state_bits=s,message_exponent=20))-
                min(r['margin_bits'] for r in select(rows,family=family,block_bits=b,state_bits=s,message_exponent=20))
                for s in range(9,21))
            model=payload['models'][f'{family}_{b}']
            item['onset_model_intercept']=model['onset_intercept_bits']
            item['chernoff_model_intercept']=model['chernoff_intercept_bits']
            finite=[sh for sh in model['shells'] if sh['log_upper'] is not None]
            if finite:
                contributions=[math.log(sh['count'])+sh['log_upper'] for sh in finite]
                minimum=finite[0]
                item['minimum_shell_model_intercept']=-(math.log(minimum['count'])+minimum['log_upper'])/math.log(2)
                item['minimum_shell_model_share']=math.exp(contributions[0]-float(np.logaddexp.reduce(contributions)))
            summary.append(item)
    shared=0
    for b in (8,32):
        for row in select(rows,family='bch',block_bits=b):
            other=select(rows,family='rm',block_bits=b,step_bits=row['step_bits'],state_bits=row['state_bits'],message_exponent=row['message_exponent'])[0]
            if row['margin_bits']!=other['margin_bits']: raise ArithmeticError('shared BCH/RM anchors differ')
            shared+=1
    old_path=HERE/'bch_growth.csv'
    old_receipt=json.loads((HERE/'bch_growth.json').read_text())
    if hashlib.sha256(old_path.read_bytes()).hexdigest()!=old_receipt['csv_sha256']:
        raise ValueError('BCH reference CSV changed')
    differences=[]
    with old_path.open(newline='') as handle:
        for old in csv.DictReader(handle):
            matched=select(rows,family='bch',**{k:int(old[k]) for k in
                ('block_bits','step_bits','state_bits','message_exponent')})[0]
            differences.append(abs(matched['margin_bits']-float(old['margin_bits'])))
    if max(differences)>.01: raise ArithmeticError('BCH reference changed by more than the witness-grid allowance')
    result=dict(anchors=summary,identical_shared_anchor_pairs=shared,
                bch_reference_rows_checked=len(differences),bch_reference_maximum_difference_bits=max(differences),
                edge_witnesses=sum(r['witness_at_edge'] for r in rows),
                report_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                study_sha256=hashlib.sha256((HERE/'engineering_surfaces.json').read_bytes()).hexdigest())
    (HERE/'engineering_surface_summary.json').write_text(json.dumps(result,indent=2)+'\n')

    # Show the familiar RM engineering slices in the same form as BCH.
    fig,axes=plt.subplots(1,2,figsize=(12,4.8))
    palette={8:'#888888',32:'#bf7720',128:'#226f96',512:'#643e91'}
    for b,color in palette.items():
        selected=sorted(select(rows,family='rm',block_bits=b,step_bits=64,state_bits=20),key=lambda r:r['message_exponent'])
        for ax,field in zip(axes,('margin_bits','intercept_bits')):
            ax.plot([r['message_exponent'] for r in selected],[r[field] for r in selected],
                    color=color,marker='o',markersize=3,label=f'RM {b}')
        model=payload['models'][f'rm_{b}']['chernoff_intercept_bits']
        if model is not None: axes[1].axhline(model,color=color,linestyle='--',alpha=.65)
    axes[0].set_title('Growing K consumes margin after the finite-epoch transition')
    axes[1].set_title('Removing the row count reveals the constituent contribution')
    axes[0].set_ylabel('Q1 margin (bits)');axes[1].set_ylabel('Q1 margin + log2(outer rows)')
    for ax in axes: ax.set_xlabel('Message size log2(K)');ax.grid(alpha=.2)
    axes[0].legend(frameon=False,ncol=2)
    fig.suptitle('RM growth at t=64, s=20')
    fig.text(.5,.01,'Q1 diagnostics only. Dashed: persistent-state Chernoff model; cancellations and finite-epoch fluctuations omitted.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.05,1,.94));fig.savefig(HERE/'rm_growth_curves.png',dpi=160);plt.close(fig)

    fig,axes=plt.subplots(2,2,figsize=(10,7),sharex=True,sharey=True)
    for ax,b in zip(axes.flat,palette):
        for t,color in ((64,'#226f96'),(128,'#bf7720'),(256,'#643e91')):
            selected=sorted(select(rows,family='rm',block_bits=b,step_bits=t,message_exponent=20),key=lambda r:r['state_bits'])
            ax.plot([r['state_bits'] for r in selected],[r['intercept_bits'] for r in selected],color=color,marker='o',markersize=3,label=f't={t}')
        model=payload['models'][f'rm_{b}']['chernoff_intercept_bits']
        if model is not None: ax.axhline(model,color='#888',linestyle='--')
        ax.set_title(f'RM {b}');ax.grid(alpha=.2)
    axes[0,0].legend(frameon=False)
    for ax in axes[1]: ax.set_xlabel('State dimension s')
    for ax in axes[:,0]: ax.set_ylabel('Q1 margin + log2(outer rows)')
    fig.suptitle('RM state-size tradeoffs at K=2^20')
    fig.text(.5,.01,'Exact outer spectra; corrected four-state transfer. Dashed: onset-model Chernoff level. Q1 only.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.05,1,.94));fig.savefig(HERE/'rm_state_tradeoff.png',dpi=160);plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(12,4.8))
    random_colors={8:'#888888',32:'#bf7720',64:'#a69b22',128:'#226f96',256:'#2a9d8f',512:'#643e91'}
    for b,color in random_colors.items():
        selected=sorted(select(rows,family='random_mean',block_bits=b,step_bits=64,state_bits=20),key=lambda r:r['message_exponent'])
        for ax,field in zip(axes,('margin_bits','intercept_bits')):
            ax.plot([r['message_exponent'] for r in selected],[r[field] for r in selected],
                    color=color,marker='o',markersize=3,label=f'B={b}')
        model=payload['models'][f'random_mean_{b}']['chernoff_intercept_bits']
        if model is not None: axes[1].axhline(model,color=color,linestyle='--',alpha=.65)
    axes[0].set_title('Message growth consumes the ensemble margin')
    axes[1].set_title('Larger blocks raise the constituent contribution')
    axes[0].set_ylabel('Q1 ensemble margin (bits)');axes[1].set_ylabel('Q1 ensemble margin + log2(outer rows)')
    for ax in axes: ax.set_xlabel('Message size log2(K)');ax.grid(alpha=.2)
    axes[0].legend(frameon=False,ncol=3,loc='upper left',bbox_to_anchor=(.02,.80))
    fig.suptitle('Random ensemble growth at t=64, s=20')
    fig.text(.5,.01,'Average over one reused random outer and encoder setup. Dashed: onset-model Chernoff level. Q1 only.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.05,1,.94));fig.savefig(HERE/'random_growth_curves.png',dpi=160);plt.close(fig)

    fig,axes=plt.subplots(3,2,figsize=(10,10),sharex=True)
    for ax,b in zip(axes.flat,random_colors):
        for t,color in ((64,'#226f96'),(128,'#bf7720'),(256,'#643e91')):
            selected=sorted(select(rows,family='random_mean',block_bits=b,step_bits=t,message_exponent=20),key=lambda r:r['state_bits'])
            ax.plot([r['state_bits'] for r in selected],[r['intercept_bits'] for r in selected],color=color,marker='o',markersize=3,label=f't={t}')
        model=payload['models'][f'random_mean_{b}']['chernoff_intercept_bits']
        if model is not None: ax.axhline(model,color='#888',linestyle='--')
        ax.set_title(f'Random B={b}');ax.grid(alpha=.2)
    axes[0,0].legend(frameon=False)
    for ax in axes[2]: ax.set_xlabel('State dimension s')
    for ax in axes[:,0]: ax.set_ylabel('Q1 margin + log2(outer rows)')
    fig.suptitle('Random ensemble state-size tradeoffs at K=2^20')
    fig.text(.5,.01,'Panels use different vertical ranges. Dashed: onset-model Chernoff level. Ensemble Q1 diagnostics, not fixed-code certificates.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.04,1,.95));fig.savefig(HERE/'random_state_tradeoff.png',dpi=160);plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(12,4.8))
    for family in FAMILIES:
        chosen=sorted([r for r in summary if r['family']==family],key=lambda r:r['block_bits'])
        axes[0].plot([r['block_bits'] for r in chosen],[r['intercept_bits'] for r in chosen],
                     color=COLORS[family],marker='o',label=LABELS[family])
        b=128 if family=='bch' else 512
        selected=sorted(select(rows,family=family,block_bits=b,step_bits=64,message_exponent=20),key=lambda r:r['state_bits'])
        axes[1].plot([r['state_bits'] for r in selected],[r['margin_bits'] for r in selected],
                     color=COLORS[family],marker='o',markersize=3,label=f'{LABELS[family]}, B={b}')
    axes[0].set_xscale('log',base=2);axes[0].set_xticks([8,32,64,128,256,512],labels=['8','32','64','128','256','512'])
    axes[0].set_xlabel('Outer block size B');axes[0].set_ylabel('Q1 margin + log2(outer rows)')
    axes[0].set_title('How much contribution does a larger outer buy? (s=20)')
    axes[1].set_xlabel('State dimension s');axes[1].set_ylabel('Q1 margin (bits)')
    axes[1].set_title('How much margin does a larger state buy?')
    for ax in axes: ax.grid(alpha=.2);ax.legend(frameon=False,fontsize=8)
    fig.suptitle('Three constituent families at K=2^20, t=64')
    fig.text(.5,.01,'Q1 diagnostics. Random averages over one reused random outer; it does not certify a particular sampled code.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.05,1,.94));fig.savefig(HERE/'constituent_engineering_comparison.png',dpi=160);plt.close(fig)

    fig,axes=plt.subplots(2,2,figsize=(11,7))
    exponents=(16,18,20,22,24);states=(10,12,16,20)
    for ax,(family,b) in zip(axes.flat,(('bch',128),('rm',128),('rm',512),('random_mean',512))):
        values=np.array([[select(rows,family=family,block_bits=b,step_bits=64,state_bits=s,message_exponent=e)[0]['margin_bits']
                          for e in exponents] for s in states])
        # Preserve actual numeric spacing: the s samples are not equally spaced.
        im=ax.pcolormesh([15,17,19,21,23,25],[9,11,14,18,22],values,cmap='viridis',shading='flat')
        for i in range(len(states)):
            for j in range(len(exponents)):
                fraction=(values[i,j]-values.min())/(values.max()-values.min())
                ax.text(exponents[j],states[i],f'{values[i,j]:.1f}',ha='center',va='center',color='black' if fraction>.55 else 'white')
        ax.set_xticks(exponents);ax.set_yticks(states)
        ax.set_xlabel('Message size log2(K)');ax.set_ylabel('State dimension s')
        ax.set_title(f'{LABELS[family]}, B={b}')
        fig.colorbar(im,ax=ax,label='Q1 margin (bits)',fraction=.046,pad=.04)
    fig.suptitle('Measured K/s surfaces at t=64')
    fig.text(.5,.01,'Numbers are measured diagnostics; colors scale separately in each panel. Distance cutoff is fixed at 10%. Q1 only.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.05,1,.94));fig.savefig(HERE/'engineering_surface_slices.png',dpi=160);plt.close(fig)


if __name__=='__main__': main()
