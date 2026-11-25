#!/usr/bin/env python3 -I
"""MS2 fragment-ion XICs for 3 SLC25A13-specific peptides, 1 shared peptide and 1 SLC25A12-specific peptide in one
representative specimen per group (all Run3), columns ordered as in Fig. 2A. Supplementary figure (Fig. S4) by default;
--compact gives the main-figure panel (Fig. 2C).
Reference fragment pattern (top-6) = DIA-NN XIC apex in the disease-control run. Expected RT for runs without a
DIA-NN identification = Run3 median DIA-NN RT + offset estimated from shared/SLC25A12 peptides in that run.
Usage: python3 -I slc25a13_xic_figure_final.py --xicdir DIR --report diann_report_SLC25A12_13_mapped.csv --xicparquet DBS-Tohoku_xic --out FIG.pdf --table T.csv
       (Fig. 2C: add --compact --peptides EEFVLAAQK,DLGIFGLYK)
"""
import argparse, re, pandas as pd, numpy as np, pyarrow.parquet as pq, pyarrow.compute as pc, pyarrow as pa
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
PAL=["#0072B2","#E69F00","#009E73","#D55E00","#CC79A7","#56B4E9"]
# Column order follows Fig. 2A (healthy adult, carrier, LoF/LoF, LoF/Mis, disease control; 2026-11-07).
# All five specimens were acquired in the same batch (Run 3, July 2025; WW_6 = 2025-07-20).
SAMPLES=[("WW_6","C25006-47_NANDA_e480_DBS_ET_C006","Healthy adult"),
         ("WS_10","C25006-14_NANDA_e480_DBS_ET_P020M","Carrier"),
         ("SS_13","C25006-09_NANDA_e480_DBS_ET_P019","Citrin deficiency\nLoF/LoF"),
         ("SM_5","C25006-08_NANDA_e480_DBS_ET_P018","Citrin deficiency\nLoF/Mis"),
         ("DD_17","C25006-37_NANDA_e480_DBS_ET_D017","Disease control")]
# Main-figure (--compact) column headings: the group names of Fig. 2A and 2B, with LoF/LoF and LoF/Mis under one heading.
HEAD=["Healthy controls","Carriers","LoF/LoF","LoF/Mis","Disease controls"]
SPAN=("Citrin deficiency",2,3)   # heading, first column, last column
REF_RUN="C25006-37_NANDA_e480_DBS_ET_D017"   # reference fragment pattern = DIA-NN XIC apex in the disease-control run (DD_17)
PEPS=[("EEFVLAAQK",2,"SLC25A13-specific"),("INLPAPNPDHVGGYK",3,"SLC25A13-specific"),("TVELLSGVVDQTK",2,"SLC25A13-specific"),
      ("YEGFFGLYR",2,"shared SLC25A12 / SLC25A13"),("DLGIFGLYK",2,"SLC25A12-specific")]
ALIGN=["YEGFFGLYR","LQVAGEITTGPR","IVQLLAGVADQTK","LATATFAGIENK","DLGIFGLYK","EEGPSAFWK"]
WIN=0.8

def qtext(q):
    """q-value with one significant digit, as 2x10^-3."""
    m,e=("%.0e"%q).split("e"); return f"{m}x10$^{{{int(e)}}}$"

def draw(panels,peps,compact,width=None):
    """Draw the grid of chromatograms. panels[(i,j)] holds, for peptide i and specimen j: w (data frame indexed by
    retention time, one column per fragment in reference-pattern order), rt, q (None if not identified) and co."""
    nP,nS=len(peps),len(SAMPLES)
    if compact:   # main figure: Arial, text kept as text in PDF/SVG
        plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
                             "mathtext.default":"regular","pdf.fonttype":42,"svg.fonttype":"none"})
    W=width or (7.2 if compact else 3.3*nS); H=(1.75*nP+0.75) if compact else 2.3*nP
    fig,axes=plt.subplots(nP,nS,figsize=(W,H),squeeze=False)
    for i,(pep,z,cls) in enumerate(peps):
        for j,(s,run,slab) in enumerate(SAMPLES):
            d=panels[(i,j)]; w,rt,q,co=d["w"],d["rt"],d["q"],d["co"]; frs=list(w.columns)
            ax=axes[i][j]
            for k,f in enumerate(frs): ax.plot(w.index,w[f]/1e3,color=PAL[k],lw=1.2,label=f)
            ax.axvline(rt,color="0.55",ls="--",lw=0.8,zorder=0); ax.set_xlim(rt-WIN,rt+WIN); ax.set_ylim(bottom=0)
            ax.spines[["top","right"]].set_visible(False); ax.tick_params(labelsize=7)
            if compact:
                tag=(f"identified (q = {qtext(q)})" if q is not None else "not identified")
                for txt,dy in ((tag,11.5),(f"{co}/6 fragments co-eluting",4)):   # two lines above the panel, on fixed baselines
                    ax.annotate(txt,xy=(0.5,1),xycoords="axes fraction",xytext=(0,dy),textcoords="offset points",ha="center",va="baseline",fontsize=6.5)
                if i==0: ax.annotate(HEAD[j],xy=(0.5,1),xycoords="axes fraction",xytext=(0,25),textcoords="offset points",ha="center",va="baseline",fontsize=7,fontweight="bold")
                if j==0:
                    ax.set_ylabel("intensity (x10$^3$)",fontsize=6.5,labelpad=2)
                    ax.annotate(f"{pep}\n{cls}",xy=(0,0.5),xycoords="axes fraction",xytext=(-38,0),textcoords="offset points",rotation=90,
                                ha="center",va="center",fontsize=7,fontweight="bold",linespacing=1.1)
                if i==nP-1: ax.set_xlabel("Retention time (min)",fontsize=6.5)
                if j==nS-1: ax.legend(fontsize=6.5,frameon=False,loc="upper left",bbox_to_anchor=(1.06,0.93),borderaxespad=0,handlelength=1.0,labelspacing=0.15)
                ax.tick_params(labelsize=6.5); ax.locator_params(axis="x",nbins=4); ax.locator_params(axis="y",nbins=4)
                ax.set_ylim(top=ax.get_ylim()[1]*1.45)
            else:
                tag=(f"DIA-NN q = {q:.1e}" if q is not None else "not identified by DIA-NN")
                ax.text(0.02,0.97,f"{tag}\n{co}/6 library fragments with a local maximum here",transform=ax.transAxes,fontsize=6.5,va="top",color=("0.15" if q is not None else "#B22222"))
                if i==0: ax.set_title(f"{slab.replace(chr(10),' ')}\n{s}",fontsize=8.5)
                if j==0: ax.set_ylabel(f"{pep} ({z}+)\n{cls}\nintensity (x10$^3$)",fontsize=7)
                if i==nP-1: ax.set_xlabel("retention time (min)",fontsize=8)
                if j==nS-1: ax.legend(fontsize=6,frameon=False,loc="upper right")
    if compact:   # common y-scale within each peptide row so that absence is visible as a flat trace
        for i in range(nP):
            top=max(axes[i][j].get_ylim()[1] for j in range(nS))
            for j in range(nS):
                axes[i][j].set_ylim(0,top)
                if j>0: axes[i][j].tick_params(labelleft=False)
        fig.tight_layout(pad=0.5,w_pad=0.9,h_pad=0.8,rect=(0,0,1,1-15/72/H))   # 15 pt kept free at the top for the spanning heading
        lab,c0,c1=SPAN; b0=axes[0][c0].get_position(); b1=axes[0][c1].get_position()
        x0=b0.x0+0.33*b0.width; x1=b1.x0+0.67*b1.width; y=b0.y1+36/72/H
        fig.add_artist(Line2D([x0,x1],[y,y],color="black",lw=0.6)); fig.text((x0+x1)/2,y+3/72/H,lab,ha="center",va="bottom",fontsize=7,fontweight="bold")
    else: fig.tight_layout(pad=1.08)
    return fig

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--xicdir",required=True); ap.add_argument("--report",required=True); ap.add_argument("--xicparquet",required=True)
    ap.add_argument("--out",required=True); ap.add_argument("--table",required=True); ap.add_argument("--peptides",default="",help="comma-separated subset of peptides (default all)")
    ap.add_argument("--compact",action="store_true",help="main-figure style (Fig. 2C): group headings, no sample IDs, annotations above the panels, Arial, 7.2-inch width")
    ap.add_argument("--width",type=float,default=None); a=ap.parse_args()
    peps=PEPS
    if a.peptides:
        keep=[x.strip() for x in a.peptides.split(",")]; peps=[t for t in PEPS if t[0] in keep]
    rep=pd.read_csv(a.report); run3=rep[rep['Run'].str.startswith('C25006')]
    med=run3.groupby('Stripped.Sequence')['RT'].median()
    # reference patterns from the disease-control DIA-NN XIC
    t=pq.read_table(f"{a.xicparquet}/{REF_RUN}.xic.parquet"); t=t.filter(pc.is_in(t['pr'],pa.array([p+str(z) for p,z,_ in peps]))).to_pandas()
    ref={}
    for pr,g in t.groupby('pr'):
        g=g[g.feature.str.match(r'[by]\d+\^')]; apx=g.groupby('feature')['value'].max()
        apx.index=[re.sub(r'\^1$','',f).replace('^2','++') for f in apx.index]; apx=apx[apx>0].sort_values(ascending=False).head(6); ref[re.sub(r'\d$','',pr)]=apx/apx.sum()
    def own(s,pep):
        r=rep[(rep['sample']==s)&(rep['Stripped.Sequence']==pep)]; return (float(r['RT'].iloc[0]), float(r['Q.Value'].iloc[0]), float(r['Precursor.Quantity'].iloc[0])) if len(r) else (None,None,None)
    rows=[]; panels={}; data={s:pd.read_csv(f"{a.xicdir}/xic_{s}.csv") for s,_,_ in SAMPLES}
    for i,(pep,z,cls) in enumerate(peps):
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
            panels[(i,j)]=dict(w=w,rt=rt,q=q,co=co)
    fig=draw(panels,peps,a.compact,a.width)
    fig.savefig(a.out); pd.DataFrame(rows).to_csv(a.table,index=False); print("wrote",a.out,a.table)

if __name__=="__main__": main()
