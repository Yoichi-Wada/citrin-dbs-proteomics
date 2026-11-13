"""For each sample, locate each target peptide near an expected RT (aligned via LQVAGEITTGPR + YEGFFGLYR),
report apex intensity of the diagnostic fragments (max within +/-0.3 min of the best co-elution point)."""
import pandas as pd, numpy as np, sys, warnings, json; warnings.filterwarnings('ignore')
REF={'LQVAGEITTGPR':63.99,'YEGFFGLYR':52.83}   # SM_1 (Run3) anchors
# peptide: (diagnostic fragments for scoring, RT in SM_1)
TARGETS={
 'INLPAPNPDHVGGYK':(['y8','y10','y8++','y10++','y12++'],70.73),
 'EEFVLAAQK':(['y7','y8','b5'],47.70),
 'TVELLSGVVDQTK':(['y7','y8','y9'],57.53),
 'LAVATFAGIENK':(['b3','y11','y10'],72.91),
 'LATATFAGIENK':(['b3','y3','y6'],73.03),
 'LQVAGEITTGPR':(['y11','y6','b4'],63.99),
 'YEGFFGLYR':(['y5','y7'],52.83),
 'KIYSTLAGTR':(['y3','b3','b6'],51.40),
 'GLLPQLLGVAPEK':(['y7','y5','y6','y8'],68.07),
}
def smoothed(d,pep,frags):
    g=d[(d.peptide==pep)&(d.fragment.isin(frags))]
    w=g.pivot_table(index='rt_min',columns='fragment',values='intensity',aggfunc='sum').fillna(0).sort_index()
    return w
def find_anchor(d,pep,frags,rt_guess,win=12):
    w=smoothed(d,pep,frags); sm=w.rolling(3,center=True,min_periods=1).mean()
    sub=sm[(sm.index>rt_guess-win)&(sm.index<rt_guess+win)]
    sc=np.sort(sub.values,axis=1)[:,::-1][:,1] if sub.shape[1]>1 else sub.values[:,0]  # 2nd-largest
    return sub.index[int(np.argmax(sc))]
def main(samples):
    rows=[]
    for s in samples:
        d=pd.read_csv(f'/tmp/slc25a13_mzml/xic/xic_{s}.csv')
        a1=find_anchor(d,'LQVAGEITTGPR',['y11','y6'],REF['LQVAGEITTGPR']); a2=find_anchor(d,'YEGFFGLYR',['y5','y7'],REF['YEGFFGLYR'])
        slope=(a1-a2)/(REF['LQVAGEITTGPR']-REF['YEGFFGLYR']); 
        def align(rt): return a2+(rt-REF['YEGFFGLYR'])*slope
        for pep,(frags,rt_sm) in TARGETS.items():
            exp=align(rt_sm); w=smoothed(d,pep,frags)
            sub=w[(w.index>exp-1.0)&(w.index<exp+1.0)]
            if sub.empty: continue
            sm=sub.rolling(3,center=True,min_periods=1).mean()
            sc=np.sort(sm.values,axis=1)[:,::-1][:,min(1,sm.shape[1]-1)]
            apex=sm.index[int(np.argmax(sc))]
            pk=sub[(sub.index>=apex-0.3)&(sub.index<=apex+0.3)].max()
            rows.append(dict(sample=s,peptide=pep,expected_rt=round(exp,2),apex_rt=round(apex,2),**{f:int(pk.get(f,0)) for f in frags},
                             anchor_LQVAG=round(a1,2),anchor_YEGFF=round(a2,2)))
    out=pd.DataFrame(rows); pd.set_option('display.width',250); pd.set_option('display.max_columns',30)
    return out
if __name__=='__main__':
    out=main(sys.argv[1:]); print(out.fillna('').to_string(index=False))
    out.to_csv('/tmp/slc25a13_mzml/apex_table.csv',index=False)
