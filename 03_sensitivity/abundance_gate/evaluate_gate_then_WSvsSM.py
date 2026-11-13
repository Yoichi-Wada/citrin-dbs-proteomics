#!/usr/bin/env python3 -I
"""Exploratory three-step design (2026-11-07): SLC25A13 absence -> LoF/LoF; SLC25A13/ACTB ratio gate (family-LOFO cut) -> Ref/Ref;
remaining specimens scored by a nested family-LOFO ensemble trained on carriers vs LoF/Mis only.
Usage: python3 -I evaluate_gate_then_WSvsSM.py GATE.csv STAGE2_outer_predictions.csv PRODUCTION_outer_predictions.csv OUT_DIR"""
import sys, pandas as pd, numpy as np
from scipy.stats import beta
gate=pd.read_csv(sys.argv[1]); s2=pd.read_csv(sys.argv[2]); prod=pd.read_csv(sys.argv[3]); out=sys.argv[4]
def cp(k,n):
    lo=beta.ppf(0.025,k,n-k+1) if k>0 else 0.0; hi=beta.ppf(0.975,k+1,n-k) if k<n else 1.0; return f"{k}/{n} ({lo:.3f}-{hi:.3f})"
def auc(pos,neg): pos=np.asarray(pos); neg=np.asarray(neg); return float(np.mean([(a>b)+0.5*(a==b) for a in pos for b in neg]))
rows=[]
for strat in ['consensus30','diverse30','adaptive']:
    for pol in sorted(s2.policy.unique()):
        q=s2[(s2.strategy==strat)&(s2.policy==pol)][['SampleID','truth','probability','threshold','predicted_class']]
        if q.empty: continue
        d=gate.merge(q,on='SampleID',how='left')
        def final(r):
            if not r.quantified: return 'LoF/LoF'
            if r.gate_call_midpoint.startswith(('Ref/Ref','Negative branch')): return 'Ref/Ref'   # negative branch of the abundance gate
            return 'LoF/Mis' if r.predicted_class==1 else 'carrier'
        d['final_call']=d.apply(final,axis=1)
        truth_map={'SS':'LoF/LoF','SM':'LoF/Mis','WS':'carrier','DD':'Ref/Ref','WW':'Ref/Ref'}
        d['truth_class']=d.group.map(truth_map)
        d['correct4']=d.final_call==d.truth_class
        d['correct3']=(d.final_call.replace({'carrier':'other','Ref/Ref':'other'})==d.truth_class.replace({'carrier':'other','Ref/Ref':'other'}))
        sm_pos=(d.group=='SM')&(d.final_call=='LoF/Mis'); fp=(d.group!='SM')&(d.group!='SS')&(d.final_call=='LoF/Mis')
        st=d[d.SampleID.str.match(r'^(WS|SM)_')&d.probability.notna()]
        r3=st[st.batch=='Run3']
        rows.append(dict(strategy=strat,policy=pol,
            LoFMis_sensitivity=cp(int(sm_pos.sum()),7), LoFMis_specificity_vs37=cp(int(37-fp.sum()),37),
            FP_specimens=';'.join(d.SampleID[fp]), accuracy_3class=cp(int(d.correct3.sum()),60), accuracy_4class=cp(int(d.correct4.sum()),60),
            carriers_called_RefRef=';'.join(d.SampleID[(d.group=='WS')&(d.final_call=='Ref/Ref')]),
            stage2_AUC_WS_vs_SM=round(auc(st.probability[st.truth==1],st.probability[st.truth==0]),3),
            stage2_AUC_Run3carriers_vs_SM=round(auc(r3.probability[r3.truth==1],r3.probability[r3.truth==0]),3),
            n_stage2=len(st)))
        if strat=='consensus30' and pol=='full_inner_recall':
            d[[c for c in ['SampleID','group','batch','family_id','quantified','ratio_x1e4','gate_call_midpoint','probability','threshold','final_call','truth_class','correct4'] if c in d.columns]].to_csv(f"{out}/per_specimen_consensus30_full_inner_recall.csv",index=False)
tab=pd.DataFrame(rows); tab.to_csv(f"{out}/gate_then_WSvsSM_summary.csv",index=False)
pd.set_option('display.width',250); print(tab.to_string(index=False))
# reference: current production 44-classifier on the same 17
q=prod[(prod.strategy=='consensus30')&(prod.policy=='full_inner_recall')&prod.SampleID.str.match(r'^(WS|SM)_')]
print('\nreference (current 44-classifier OOF): AUC WS vs SM',round(auc(q.probability[q.truth==1],q.probability[q.truth==0]),3),
      '; Run3 carriers vs SM',round(auc(q[q.batch=='Run3'].probability[q.truth==1],q[(q.batch=='Run3')&(q.truth==0)].probability),3))
