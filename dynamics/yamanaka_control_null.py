"""Stage 2.11C: preregistered control and null analysis for Yamanaka signals."""
from __future__ import annotations
import hashlib
import itertools
import json
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
from .yamanaka_poc import condition_table, load_expression, _signature
from .yamanaka_trajectory_poc import RDS_FILES, _convert_rds, _log_cpm
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "Data"
OUT = ROOT / "results" / "Dynamics" / "stage2_11c_control_null"
SEED = 20260911
N_MONOTONIC = 10_000
N_BOOTSTRAP = 1_000
EXPECTED_DAYS = (0.0, 3.0, 7.0, 10.0)
CONFOUNDER_PROGENY = {"JAK-STAT"}
CONFOUNDER_TF = {"IRF1", "IRF2", "IRF9", "STAT1", "STAT2"}
PERTURBATION_FILES = {"GSE297233": DATA / "GSE297233_raw_counts_matrix.csv.gz", "GSE304042": DATA / "GSE304042_5_ARPE_single_triple_OSK_30-1011800743.csv.gz"}

def _protocol_hash():
    p = ROOT / "dynamics" / "STAGE_2_11C_CONTROL_AND_NULL_PROTOCOL.md"
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else "missing"

def _activity_networks():
    try:
        import decoupler as dc
    except ImportError as exc:
        raise RuntimeError("Stage 2.11C requires decoupler. Install requirements-ai.txt.") from exc
    return dc.op.progeny(organism="human", top=100), dc.op.dorothea(organism="human", levels=["A", "B", "C"])

def _score_time_series(X, meta, net):
    import decoupler as dc
    meta = meta.copy(); meta["day"] = pd.to_numeric(meta["day"], errors="coerce"); meta["group"] = meta["group"].astype(str)
    meta = meta.dropna(subset=["day"]); meta = meta[meta.day.isin(EXPECTED_DAYS)]; meta = meta[meta.group.isin(X.columns)].sort_values(["day", "group"])
    if meta.empty: return pd.DataFrame()
    samples = _log_cpm(X.loc[:, meta.group.tolist()]).T; samples.index = meta.group.tolist()
    shared = set(samples.columns) & set(net["target"].astype(str))
    if not shared: return pd.DataFrame()
    net_use = net[net["target"].astype(str).isin(shared)].copy()
    if net_use.empty or net_use["source"].nunique() == 0: return pd.DataFrame()
    acts, _ = dc.mt.ulm(data=samples.loc[:, sorted(shared)], net=net_use, tmin=1)
    acts = acts.apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)
    acts["day"] = meta.set_index("group").loc[acts.index, "day"].astype(float)
    return acts

def _day_centroids(frame): return frame.groupby("day")[[c for c in frame.columns if c != "day"]].mean().sort_index()

def _trajectory_rows(scored):
    names = sorted(scored)
    if len(names) < 2: return pd.DataFrame()
    ca, cb = _day_centroids(scored[names[0]]), _day_centroids(scored[names[1]]); common = sorted(set(ca.index) & set(cb.index))
    if len(common) != len(EXPECTED_DAYS): return pd.DataFrame()
    rows=[]
    for feature in sorted(set(ca.columns)&set(cb.columns)):
        va, vb = ca.loc[common,feature].to_numpy(float), cb.loc[common,feature].to_numpy(float)
        if np.std(va)==0 or np.std(vb)==0: continue
        ra = pd.Series(common).corr(pd.Series(va),method="spearman"); rb = pd.Series(common).corr(pd.Series(vb),method="spearman")
        rows.append({"feature":feature,"observed_corr":float(np.corrcoef(va,vb)[0,1]),"rho_day_a":float(ra),"rho_day_b":float(rb),"conserved_direction":bool(np.sign(va[-1]-va[0])==np.sign(vb[-1]-vb[0]) and va[-1]!=va[0] and vb[-1]!=vb[0]),"abs_temporal_signal":float((abs(ra)+abs(rb))/2),"a_values":va.tolist(),"b_values":vb.tolist()})
    return pd.DataFrame(rows)

def _temporal_null(rows):
    perms=list(itertools.permutations(range(len(EXPECTED_DAYS)))); identity=tuple(range(len(EXPECTED_DAYS))); null_perms=[p for p in perms if p != identity]; out=[]
    for r in rows.itertuples(index=False):
        va,vb=np.asarray(r.a_values,float),np.asarray(r.b_values,float); null=np.array([np.corrcoef(va,vb[list(p)])[0,1] for p in null_perms],float); obs=abs(float(r.observed_corr)); ge=int(np.sum(np.abs(null)>=obs-1e-12)); p=(ge+1)/(len(null)+1)
        out.append({"feature":r.feature,"n_exact_null_permutations":len(null),"null_q95":float(np.quantile(null,.95)),"null_q99":float(np.quantile(np.abs(null),.99)),"empirical_two_sided_p":float(p),"passes_null1":bool(obs>np.quantile(np.abs(null),.99) and p<=.05)})
    return pd.DataFrame(out)

def _monotonic_series(rng,start,end):
    if start < end: interior=np.sort(rng.uniform(start,end,2))
    else: interior=np.sort(rng.uniform(end,start,2))[::-1]
    return np.array([start,*interior,end],float)

def _monotonic_null(rows,seed):
    rng=np.random.default_rng(seed); out=[]
    for r in rows.itertuples(index=False):
        va,vb=np.asarray(r.a_values,float),np.asarray(r.b_values,float); null=np.empty(N_MONOTONIC)
        for i in range(N_MONOTONIC): null[i]=np.corrcoef(_monotonic_series(rng,va[0],va[-1]),_monotonic_series(rng,vb[0],vb[-1]))[0,1]
        out.append({"feature":r.feature,"n_monotonic_controls":N_MONOTONIC,"control_q95":float(np.quantile(null,.95)),"control_q99":float(np.quantile(null,.99)),"passes_null2":bool(r.observed_corr>np.quantile(null,.95))})
    return pd.DataFrame(out)

def _perturbation_activity(network,condition):
    import decoupler as dc
    result={}
    audit=[]
    for ds,path in PERTURBATION_FILES.items():
        row={"dataset":ds,"condition":condition,"path":str(path.relative_to(ROOT)),"status":"skipped","skip_reason":""}
        try:
            X,_=load_expression(path)
        except Exception as exc:
            row["skip_reason"]="load_expression_failed:" + type(exc).__name__
            row["error"] = str(exc)
            audit.append(row); continue
        try:
            labels=condition_table(X,ds).set_index("sample")["condition"]
            row["n_samples"] = int(X.shape[1])
            row["condition_counts"] = json.dumps(labels.value_counts().to_dict(), sort_keys=True)
            sig=_signature(X,labels,condition)
        except Exception as exc:
            row["skip_reason"]="signature_failed:" + type(exc).__name__
            row["error"] = str(exc)
            audit.append(row); continue
        if sig is None:
            row["skip_reason"]="signature_unavailable"
            audit.append(row); continue
        available=set(sig.index.astype(str)) & set(network["target"].astype(str))
        row["n_signature_genes"] = int(len(sig))
        row["n_shared_network_genes"] = int(len(available))
        if not available:
            row["skip_reason"]="no_shared_network_genes"
            audit.append(row); continue
        net_use=network[network["target"].astype(str).isin(available)].copy()
        row["n_network_edges"] = int(len(net_use))
        row["n_network_sources"] = int(net_use["source"].nunique())
        if net_use.empty or net_use["source"].nunique()==0:
            row["skip_reason"]="empty_activity_network"
            audit.append(row); continue
        cols=sorted(available); sample=pd.DataFrame([sig.reindex(cols).to_numpy(float)],columns=cols,index=[f"{ds}__{condition}"])
        try:
            acts,_=dc.mt.ulm(data=sample,net=net_use,tmin=1)
        except Exception as exc:
            row["skip_reason"]="decoupler_failed:" + type(exc).__name__
            row["error"] = str(exc)
            audit.append(row); continue
        if acts.empty:
            row["skip_reason"]="decoupler_empty_result"
            audit.append(row); continue
        result[ds]=acts.iloc[0].apply(float)
        row["status"]="scored"
        row["n_activity_features"] = int(len(result[ds]))
        audit.append(row)
    return result,audit

def _context_control(progeny,dorothea,candidates):
    progeny_activity, progeny_audit = _perturbation_activity(progeny,"OSK")
    dorothea_activity, dorothea_audit = _perturbation_activity(dorothea,"OSK")
    activity={"PROGENy":progeny_activity,"DoRothEA":dorothea_activity}; audits=progeny_audit+dorothea_audit; rows=[]; datasets={}
    for rep,acts in activity.items():
        datasets[rep]=sorted(acts); n_features=max((len(s) for s in acts.values()),default=0)
        for feature in candidates.get(rep,[]):
            if len(acts)<2 or any(feature not in a.index for a in acts.values()): rows.append({"representation":rep,"feature":feature,"status":"inconclusive","survives":False,"known_delivery_confounded":False}); continue
            names=sorted(acts); a,b=acts[names[0]][feature],acts[names[1]][feature]; direction=np.sign(a)==np.sign(b) and a!=0 and b!=0; ranks=[s.abs().sort_values(ascending=False).index.tolist().index(feature)+1 for s in acts.values()]; frac=max(ranks)/max(n_features,1); confound=(rep=="PROGENy" and feature in CONFOUNDER_PROGENY) or (rep=="DoRothEA" and feature in CONFOUNDER_TF)
            rows.append({"representation":rep,"feature":feature,"activity_a":float(a),"activity_b":float(b),"concordant_direction":bool(direction),"worst_rank_fraction":float(frac),"known_delivery_confounded":bool(confound),"survives":bool(direction and frac<=.25 and not confound),"status":"ok"})
    columns=["representation","feature","activity_a","activity_b","concordant_direction","worst_rank_fraction","known_delivery_confounded","survives","status"]
    return pd.DataFrame(rows,columns=columns),datasets,audits

def _bootstrap_stability(scored,rep,features,seed):
    rng=np.random.default_rng(seed); rows=[]
    for feature in features:
        same=total=0
        for _ in range(N_BOOTSTRAP):
            vals=[]
            for ds in sorted(scored):
                f=scored[ds]; parts=[]
                for day in EXPECTED_DAYS:
                    g=f[f.day==day][feature].to_numpy(float) if feature in f.columns else np.array([])
                    if len(g): parts.append(float(rng.choice(g)))
                if len(parts)==4: vals.append(parts)
            if len(vals)==2:
                a,b=map(np.asarray,vals)
                if np.std(a) and np.std(b): same+=int(np.sign(a[-1]-a[0])==np.sign(b[-1]-b[0])); total+=1
        rows.append({"representation":rep,"feature":feature,"bootstrap_replicates":N_BOOTSTRAP,"stable_direction_fraction":float(same/total) if total else np.nan,"stable":bool(total and same/total>=.8)})
    return pd.DataFrame(rows)

def _assign_tiers(detail,context,stability):
    context_cols=["representation","feature","survives","known_delivery_confounded"]
    context_safe=context.reindex(columns=context_cols,fill_value=False)
    stability_safe=stability.reindex(columns=["representation","feature","stable"],fill_value=False)
    out=detail.merge(context_safe,on=["representation","feature"],how="left").merge(stability_safe,on=["representation","feature"],how="left")
    out["tier"]="NONE"
    a=out.passes_null1.fillna(False)&out.passes_null2.fillna(False)&out.conserved_direction.fillna(False)&out.stable.fillna(False)
    out.loc[a,"tier"]="A"
    b=a&out.survives.fillna(False)&~out.known_delivery_confounded.fillna(False)
    out.loc[b,"tier"]="B"
    return out

def _decision(tiered):
    a=tiered[tiered.tier=="A"]
    if a.empty: return {"decision":"CONFOUNDED_UNSUPPORTED","tier_a_n":0,"tier_b_n":0,"survival_rate":0.0,"reason":"No candidate passed temporal-order, monotonic-shape and stability gates."}
    b=tiered[tiered.tier=="B"]; rate=len(b)/len(a); decision="PROCEED" if rate>=.60 else "MIXED" if rate>=.30 else "CONFOUNDED_UNSUPPORTED"; return {"decision":decision,"tier_a_n":int(len(a)),"tier_b_n":int(len(b)),"survival_rate":float(rate),"reason":"Decision follows the fixed Stage 2.11C thresholds."}

def run():
    OUT.mkdir(parents=True,exist_ok=True); progeny,dorothea=_activity_networks(); scored={"PROGENy":{},"DoRothEA":{}}; audits=[]
    with tempfile.TemporaryDirectory(prefix="stage2_11c_",dir=str(ROOT)) as td:
        tmp=Path(td)
        for ds,rds in RDS_FILES.items():
            X,meta=_convert_rds(ds,rds,tmp); audits.append({"dataset":ds,"n_genes":int(X.shape[0]),"n_groups":int(X.shape[1]),"days":sorted(pd.to_numeric(meta.day,errors="coerce").dropna().unique().tolist())}); scored["PROGENy"][ds]=_score_time_series(X,meta,progeny); scored["DoRothEA"][ds]=_score_time_series(X,meta,dorothea)
    candidates={}; frames=[]
    for rep,values in scored.items():
        r=_trajectory_rows(values)
        if not r.empty: r["representation"]=rep; frames.append(r); candidates[rep]=r.feature.tolist()
    trajectory=pd.concat(frames,ignore_index=True) if frames else pd.DataFrame(); null1=_temporal_null(trajectory) if not trajectory.empty else pd.DataFrame(); null2=_monotonic_null(trajectory,SEED+1) if not trajectory.empty else pd.DataFrame(); detail=trajectory.merge(null1,on="feature",how="left").merge(null2,on="feature",how="left"); context,context_datasets,perturbation_audit=_context_control(progeny,dorothea,candidates); stability=pd.concat([_bootstrap_stability(scored["PROGENy"],"PROGENy",candidates.get("PROGENy",[]),SEED+2),_bootstrap_stability(scored["DoRothEA"],"DoRothEA",candidates.get("DoRothEA",[]),SEED+3)],ignore_index=True); tiered=_assign_tiers(detail,context,stability) if not detail.empty else pd.DataFrame(); decision=_decision(tiered) if not tiered.empty else {"decision":"INCONCLUSIVE","reason":"No complete two-dataset temporal activity matrix was available."}
    trajectory.to_csv(OUT/"01_candidate_trajectories.csv",index=False); null1.to_csv(OUT/"02_null1_temporal_order.csv",index=False); null2.to_csv(OUT/"03_null2_monotonic_shape.csv",index=False); context.to_csv(OUT/"04_context_control_OSK.csv",index=False); stability.to_csv(OUT/"05_candidate_stability.csv",index=False); tiered.to_csv(OUT/"06_candidate_tiers.csv",index=False); pd.DataFrame(perturbation_audit).to_csv(OUT/"07_perturbation_activity_audit.csv",index=False)
    summary={"status":"ok","protocol_hash":_protocol_hash(),"seed":SEED,"n_monotonic_controls":N_MONOTONIC,"n_bootstrap":N_BOOTSTRAP,"dataset_audit":audits,"context_control_datasets":context_datasets,"perturbation_activity_audit":perturbation_audit,"decision":decision,"interpretation_guardrail":"Stage 2.11C is a falsification gate. Passing supports candidate temporal regulators; it does not establish causality, lineage, or a molecular mechanism.","important_limitations":["The temporal-order null excludes the observed ordering from the exact permutation reference; with four days there are 23 null permutations and the minimum +1-corrected empirical p-value is 1/24 ≈ 0.0417.","GSE304042 is a heterologous biological/context control, not a clean Sendai-only control because cell type and experimental design differ.","Pseudobulk cannot establish single-cell lineage or cell-level transition dynamics."]}
    (OUT/"stage2_11c_summary.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False,default=str),encoding="utf-8"); return summary

def main(): print(json.dumps(run(),indent=2,ensure_ascii=False,default=str))
if __name__=="__main__": main()
