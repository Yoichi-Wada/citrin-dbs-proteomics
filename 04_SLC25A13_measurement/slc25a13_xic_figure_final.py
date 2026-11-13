#!/usr/bin/env python3 -I
"""Final supplementary figure: MS2 fragment-ion XICs for 3 SLC25A13-specific peptides, 1 shared peptide and
1 SLC25A12-specific peptide in one representative specimen per group (all Run3), columns ordered as in Fig. 2A.
Reference fragment pattern (top-6) = DIA-NN XIC apex in the disease-control run. Expected RT for runs without a
DIA-NN identification = Run3 median DIA-NN RT + offset estimated from shared/SLC25A12 peptides in that run.
Usage: python3 -I slc25a13_xic_figure_final.py --xicdir DIR --report diann_report_SLC25A12_13_mapped.csv --xicparquet DBS-Tohoku_xic --out FIG.pdf --table T.csv
"""
import argparse, re, pandas as pd, numpy as np, pyarrow.parquet as pq, pyarrow.compute as pc, pyarrow as pa
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
PAL=["#0072B2","#E69F00","#009E73","#D55E00","#CC79A7","#56B4E9"]
# Column order and labels follow Fig. 2A (healthy adult, carrier, LoF/LoF, LoF/Mis, disease control; 2026-11-07).
# All five specimens were acquired in the same batch (Run 3, July 2025; WW_6 = 2025-07-20).
SAMPLES=[("WW_6","C25006-47_NANDA_e480_DBS_ET_C006","Healthy adult\n(Ref/Ref)"),
         ("WS_10","C25006-14_NANDA_e480_DBS_ET_P020M","Carrier\n(Ref/LoF)"),
         ("SS_13","C25006-09_NANDA_e480_DBS_ET_P019","Citrin deficiency\nLoF/LoF"),
         ("SM_5","C25006-08_NANDA_e480_DBS_ET_P018","Citrin deficiency\nLoF/Mis"),
         ("DD_17","C25006-37_NANDA_e480_DBS_ET_D017","Disease control\n(Ref/Ref)")]
REF_RUN="C25006-37_NANDA_e480_DBS_ET_D017"   # reference fragment pattern = DIA-NN XIC apex in the disease-control run (DD_17)
PEPS=[("EEFVLAAQK",2,"SLC25A13-specific"),("INLPAPNPDHVGGYK",3,"SLC25A13-specific"),("TVELLSGVVDQTK",2,"SLC25A13-specific"),
      ("YEGFFGLYR",2,"shared SLC25A12 / SLC25A13"),("DLGIFGLYK",2,"SLC25A12-specific")]
ALIGN=["YEGFFGLYR","LQVAGEITTGPR","IVQLLAGVADQTK","LATATFAGIENK","DLGIFGLYK","EEGPSAFWK"]
WIN=0.8
ap=argparse.ArgumentParser(); ap.add_argument("--xicdir",required=True); ap.add_argument("--report",required=True); ap.add_argument("--xicparquet",required=True)
ap.add_argument("--out",required=True); ap.add_argument("--table",required=True); ap.add_argument("--peptides",default="",help="comma-separated subset of peptides (default all)")
ap.add_argument("--compact",action="store_true",help="main-figure style: no sample IDs, short annotations, 7.2-inch width")
ap.add_argument("--width",type=float,default=None); a=ap.parse_args()
if a.peptides:
    keep=[x.strip() for x in a.peptides.split(",")]; PEPS=[t for t in PEPS if t[0] in keep]
rep=pd.read_csv(a.report); run3=rep[rep['Run'].str.startswith('C25006')]
med=run3.groupby('Stripped.Sequence')['RT'].median()
# reference patterns from the disease-control DIA-NN XIC
t=pq.read_table(f"{a.xicparquet}/{REF_RUN}.xic.parquet"); t=t.filter(pc.is_in(t['pr'],pa.array([p+str(z) for p,z,_ in PEPS]))).to_pandas()
ref={}
for pr,g in t.groupby('pr'):
    g=g[g.feature.str.match(r'[by]\d+\^')]; apx=g.groupby('feature')['value'].max()
    apx.index=[re.sub(r'\^1$','',f).replace('^2','++') for f in apx.index]; apx=apx[apx>0].sort_values(ascending=False).head(6); ref[re.sub(r'\d$','',pr)]=apx/apx.sum()
def own(s,pep):
    r=rep[(rep['sample']==s)&(rep['Stripped.Sequence']==pep)]; return (float(r['RT'].iloc[0]), float(r['Q.Value'].iloc[0]), float(r['Precursor.Quantity'].iloc[0])) if len(r) else (None,None,None)
rows=[]; data={s:pd.read_csv(f"{a.xicdir}/xic_{s}.csv") for s,_,_ in SAMPLES}
W=a.width or (7.2 if a.compact else 3.3*len(SAMPLES)); H=(1.75 if a.compact else 2.3)*len(PEPS)+(0.3 if a.compact else 0)
fig,axes=plt.subplots(len(PEPS),len(SAMPLES),figsize=(W,H),squeeze=False)
for i,(pep,z,cls) in enumerate(PEPS):
    pat=ref[pep]; frs=list(pat.index); p=pat.values/np.linalg.norm(pat.values)
    for j,(s,run,slab) in enumerate(SAMPLES):
        rt,q,quant=own(s,pep)
        if rt is None:
            offs=[own(s,x)[0]-med[x] for x in ALIGN if own(s,x)[0] is not None and x in med]
            rt=float(med[pep])+float(np.median(offs)); src="aligned"
        else: src="DIA-NN"
        d=data[s]; g=d[(d.peptide==pep)&(d.fragment.isin(frs))&(d.rt_min>=rt-WIN)&(d.rt_min<=rt+WIN)]
        w=g.pivot_table(index="rt_min",columns="fragment",values="intensity",aggfunc="sum").fillna(0).sort_index()
        for f in frs:
            if f not in w: w[f]=0.0
        w=w[frs]; sub=w[(w.index>=rt-0.3)&(w.index<=rt+0.3)]; X=sub.values
        cos=np.array([(x@p)/np.linalg.norm(x) if np.linalg.norm(x)>0 else 0 for x in X]); wint=X@p; b=int(np.argmax(wint*cos))
        # co-eluting = fragment > 0 at the best scan and showing a local maximum within +/-1 cycle of it
        def local_max(col,b):
            for c in (b-1,b,b+1):
                if 0<=c<len(col) and col[c]>0 and col[c]>=(col[c-1] if c>0 else 0) and col[c]>=(col[c+1] if c<len(col)-1 else 0): return True
            return False
        co=sum(1 for k,f in enumerate(frs) if X[b,k]>0 and local_max(sub[f].values,b))
        rows.append(dict(sample=s,group=slab,peptide=pep,charge=z,cls=cls,rt_used=round(rt,2),rt_source=src,diann_q=q,diann_quantity=quant,
                         top6_fragments=";".join(frs),n_coeluting=co,cosine=round(float(cos[b]),2),pattern_weighted_intensity=int(wint[b])))
        ax=axes[i][j]
        for k,f in enumerate(frs): ax.plot(w.index,w[f]/1e3,color=PAL[k],lw=1.2,label=f)
        ax.axvline(rt,color="0.55",ls="--",lw=0.8,zorder=0); ax.set_xlim(rt-WIN,rt+WIN); ax.set_ylim(bottom=0)
        ax.spines[["top","right"]].set_visible(False); ax.tick_params(labelsize=7)
        if a.compact:
            tag=("identified (q = %.0e)"%q if q is not None else "not identified")
            ax.text(0.03,0.96,f"{tag}\n{co}/6 fragments co-eluting",transform=ax.transAxes,fontsize=5.8,va="top",color=("0.2" if q is not None else "#B22222"))
            if i==0: ax.set_title(slab,fontsize=7.5,pad=4)
            if j==0: ax.set_ylabel(f"{pep}\n{cls}\nintensity (x10$^3$)",fontsize=6.3)
            if i==len(PEPS)-1: ax.set_xlabel("Retention time (min)",fontsize=7)
            if j==len(SAMPLES)-1: ax.legend(fontsize=5,frameon=False,loc="center right",bbox_to_anchor=(1.02,0.42),handlelength=1.2,labelspacing=0.2)
            ax.tick_params(labelsize=6); ax.locator_params(axis="x",nbins=4); ax.locator_params(axis="y",nbins=4)
            ax.set_ylim(top=ax.get_ylim()[1]*1.45)
        else:
            tag=(f"DIA-NN q = {q:.1e}" if q is not None else "not identified by DIA-NN")
            ax.text(0.02,0.97,f"{tag}\n{co}/6 library fragments with a local maximum here",transform=ax.transAxes,fontsize=6.5,va="top",color=("0.15" if q is not None else "#B22222"))
            if i==0: ax.set_title(f"{slab.replace(chr(10),' ')}\n{s}",fontsize=8.5)
            if j==0: ax.set_ylabel(f"{pep} ({z}+)\n{cls}\nintensity (x10$^3$)",fontsize=7)
            if i==len(PEPS)-1: ax.set_xlabel("retention time (min)",fontsize=8)
            if j==len(SAMPLES)-1: ax.legend(fontsize=6,frameon=False,loc="upper right")
if a.compact:   # common y-scale within each peptide row so that absence is visible as a flat trace
    for i in range(len(PEPS)):
        top=max(axes[i][j].get_ylim()[1] for j in range(len(SAMPLES)))
        for j in range(len(SAMPLES)):
            axes[i][j].set_ylim(0,top)
            if j>0: axes[i][j].tick_params(labelleft=False)
fig.tight_layout(pad=0.5 if a.compact else 1.08); fig.savefig(a.out); pd.DataFrame(rows).to_csv(a.table,index=False); print("wrote",a.out,a.table)
