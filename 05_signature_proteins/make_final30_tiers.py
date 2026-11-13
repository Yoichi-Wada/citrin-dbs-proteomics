#!/usr/bin/env python3
"""Selection-stability table for the 30 proteins of the final consensus-30 panel (WW-fixed production ensemble).
Tier rule for the 30-protein panel: A = in the consensus panel of all 35 outer folds; B = 30-34 folds; C = fewer than 30 folds (matches the Results sentence '6 in every partition, 18 in at least 30').
Adds final-model descriptors: members (of 10) with a non-zero coefficient and the mean coefficient in the frozen model."""
import csv, os, statistics
from collections import defaultdict
REPO=os.path.abspath(os.environ.get("CITRIN_REPO","."))   # repository root (public release)
P=os.path.join(REPO,"05_signature_proteins","results")+os.sep
E=os.path.join(REPO,"01_nested_ensemble","results")+os.sep
panel=list(csv.DictReader(open(P+"final30_panel.csv")))
feat=list(csv.DictReader(open(E+"outer_consensus30_features.csv")))
coef=list(csv.DictReader(open(E+"outer_consensus_member_coefficients.csv")))
fcoef=list(csv.DictReader(open(E+"final_consensus_member_coefficients.csv")))
fams=sorted({r['held_family'] for r in feat}, key=int); n_outer=len(fams); assert n_outer==35
in_panel=defaultdict(set)
for r in feat: in_panel[r['FeatureID']].add(r['held_family'])
nz_folds=defaultdict(set); nz_members=defaultdict(lambda: defaultdict(int))
for r in coef:
    if r['nonzero'] in ('TRUE','1','True'): nz_folds[r['FeatureID']].add(r['held_family']); nz_members[r['FeatureID']][r['held_family']]+=1
fin_nz=defaultdict(int); fin_coef=defaultdict(list)
for r in fcoef:
    if r['FeatureID']=='(Intercept)': continue
    fin_coef[r['FeatureID']].append(float(r['coefficient']))
    if r['nonzero'] in ('TRUE','1','True'): fin_nz[r['FeatureID']]+=1
pos_fams={'13','14','17','21','22'}
out=[]
for p in panel:
    a=p['FeatureID']; n_in=len(in_panel[a]); nzf=len(nz_folds[a])
    med=statistics.median([nz_members[a].get(f,0) for f in in_panel[a]]) if n_in else 0
    miss=[f for f in fams if f not in in_panel[a]]; miss_pos=[f for f in miss if f in pos_fams]
    tier='A' if n_in==n_outer else ('B' if n_in>=30 else 'C')
    absent=';'.join(miss_pos) if tier in ('B','C') else ''
    out.append({'rank_final_panel':int(p['rank']),'gene':p['Gene'],'accession':a,'protein':p['Protein'],'unique_peptides':p['UniquePeptides'],
                'final_members_top30_of_10':p['members_top30'],'final_members_nonzero_of_10':fin_nz[a],
                'final_mean_coefficient':round(statistics.mean(fin_coef[a]),4) if fin_coef[a] else '',
                'in_candidate_panel_of_35':n_in,'folds_nonzero_any_member_of_35':nzf,
                'median_members_nonzero_of_10':int(med) if float(med).is_integer() else med,
                'absent_when_family_held_out':absent,'n_outer_folds_missing':len(miss),'tier':tier,
                'tier_label':f"{tier} · {n_in}/{n_outer}"+(f" (fam {absent})" if absent else ''),'in_original15':p['in_original15']})
order={'A':0,'B':1,'C':2}
out.sort(key=lambda r:(order[r['tier']],-r['in_candidate_panel_of_35'],-r['final_members_nonzero_of_10'],r['gene']))
with open(P+"final30_ensemble_tiers.csv","w",newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
for r in out: print(f"{r['tier_label']:<22} {r['gene']:<9} fin_nz {r['final_members_nonzero_of_10']:>2} coef {r['final_mean_coefficient']:>8} nzfolds {r['folds_nonzero_any_member_of_35']:>2} med {r['median_members_nonzero_of_10']:>3} orig15 {r['in_original15']}")
from collections import Counter; print(Counter(r['tier'] for r in out))
