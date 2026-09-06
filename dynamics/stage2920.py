"""Stage 2.9.20: controlled repair of the common biological state space."""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_9_20"
OUT.mkdir(parents=True, exist_ok=True)
TARGET = ["GSE67462", "GSE28688", "GSE297234"]

def log(x): print(f"Stage 2.9.20: {x}", flush=True)

def load_space():
    from dynamics.validation import _load_common_space
    matrix, meta = _load_common_space()
    meta = meta[meta.dataset.isin(TARGET)].copy()
    cols = [c for c in meta.matrix_column if c in matrix.columns]
    meta = meta[meta.matrix_column.isin(cols)].copy()
    meta = meta.set_index("matrix_column").loc[cols].reset_index()
    return matrix.loc[:, cols], meta

def variance_table(matrix, meta):
    dummies = pd.get_dummies(meta.dataset, dtype=float).to_numpy()
    t = meta.time_hours.to_numpy(float); ok_t = np.isfinite(t)
    tn = np.zeros_like(t)
    if ok_t.any() and np.ptp(t[ok_t]) > 0: tn[ok_t] = (t[ok_t]-np.mean(t[ok_t]))/np.std(t[ok_t])
    Xd=np.column_stack([np.ones(len(meta)),dummies[:,1:]]); Xt=np.column_stack([np.ones(len(meta)),tn]); Xdt=np.column_stack([Xd,tn])
    rows=[]
    for gene,s in matrix.loc[:,meta.matrix_column].iterrows():
        y=s.to_numpy(float); good=np.isfinite(y)&ok_t
        if good.sum()<6 or np.var(y[good])<=1e-12: continue
        yy=y[good]
        def r2(X):
            xx=X[good]; beta=np.linalg.lstsq(xx,yy,rcond=None)[0]; pred=xx@beta; den=np.sum((yy-yy.mean())**2)
            return max(0.,1.-np.sum((yy-pred)**2)/den) if den>0 else 0.
        rd,rt,rdt=r2(Xd),r2(Xt),r2(Xdt)
        rows.append({"gene":gene,"r2_dataset":rd,"r2_time":rt,"r2_dataset_time":rdt,"dataset_excess_r2":max(0.,rdt-rt),"time_minus_dataset_r2":rt-rd})
    out=pd.DataFrame(rows);out.to_csv(OUT/"01_gene_variance_diagnostics.csv",index=False);return out

def feature_sets(matrix,var):
    genes=set(matrix.index.astype(str)); sets={"baseline_all":sorted(genes)}
    if len(var):
        v=var.set_index("gene")
        sets["low_dataset_excess_q80"]=sorted(set(v.index[v.dataset_excess_r2<=v.dataset_excess_r2.quantile(.80)]).intersection(genes))
        sets["low_dataset_excess_q70"]=sorted(set(v.index[v.dataset_excess_r2<=v.dataset_excess_r2.quantile(.70)]).intersection(genes))
        sets["time_dominant"]=sorted(set(v.index[(v.r2_time>=v.r2_dataset)&(v.r2_time>=v.r2_time.quantile(.50))]).intersection(genes))
    for k,g in list(sets.items()):
        if len(g)<100: sets[k]=sorted(genes)
    pd.DataFrame([{"variant":k,"n_genes":len(g)} for k,g in sets.items()]).to_csv(OUT/"02_feature_set_summary.csv",index=False);return sets

def aggregate_vectors(matrix,meta,genes):
    m=matrix.loc[genes,meta.matrix_column]; rows=[]
    for ds,g in meta.groupby("dataset"):
        for t,gt in g.groupby("time_hours"):
            cols=gt.matrix_column.tolist()
            if cols: rows.append({"dataset":ds,"time_hours":float(t),"vector":m[cols].mean(axis=1).to_numpy(float)})
    return rows

def evaluate_variant(matrix,meta,genes,variant):
    m=matrix.loc[genes,meta.matrix_column]; X=np.array(m.T.to_numpy(float),copy=True)
    med=np.nanmedian(X,axis=0); inds=np.where(~np.isfinite(X))
    if len(inds[0]): X[inds]=np.take(med,inds[1])
    X=StandardScaler().fit_transform(X)
    ncomp=min(10,X.shape[0]-1,X.shape[1]); coords=PCA(n_components=ncomp).fit_transform(X); ds=meta.dataset.to_numpy(); k=min(5,len(ds)-1)
    idx=NearestNeighbors(n_neighbors=k+1).fit(coords).kneighbors(return_distance=False)[:,1:]
    same=np.mean(np.array([[ds[j]==ds[i] for j in row] for i,row in enumerate(idx)])); props=pd.Series(ds).value_counts(normalize=True).to_numpy(); expected=float(np.sum(props**2))
    agg=aggregate_vectors(matrix,meta,genes); by={(r["dataset"],r["time_hours"]):r["vector"] for r in agg}; pairs=[]
    for i,a in enumerate(TARGET):
        for b in TARGET[i+1:]:
            common=sorted(set(t for d,t in by if d==a)&set(t for d,t in by if d==b));cors=[]
            for t in common:
                x,y=by[(a,t)],by[(b,t)]; good=np.isfinite(x)&np.isfinite(y)
                if good.sum()>2 and np.std(x[good])>0 and np.std(y[good])>0: cors.append(float(np.corrcoef(x[good],y[good])[0,1]))
            if cors:pairs.append({"dataset_a":a,"dataset_b":b,"n_common_times":len(cors),"mean_pearson":np.mean(cors)})
    return {"variant":variant,"n_genes":len(genes),"same_dataset_knn":same,"expected_same_dataset_knn":expected,"knn_excess":same-expected,"mean_matched_time_pearson":np.mean([p["mean_pearson"] for p in pairs]) if pairs else np.nan,"n_pairs":len(pairs)},pairs

def program_variant(matrix,meta):
    """Evaluate fixed Stage 2.9.14 programs using the same within-sample rank rule."""
    from dynamics.stage2914 import PROGRAMS
    ranks=matrix.rank(axis=0,method="average",pct=True); lookup={str(g).upper():g for g in matrix.index}; rows=[]
    for _,r in meta.iterrows():
        activities={}
        for pid,spec in PROGRAMS.items():
            pos=[lookup[g.upper()] for g in spec["positive"] if g.upper() in lookup]; neg=[lookup[g.upper()] for g in spec["negative"] if g.upper() in lookup]
            if len(pos)<3: continue
            score=float(ranks.loc[pos,r.matrix_column].mean())
            if neg: score-=float(ranks.loc[neg,r.matrix_column].mean())
            activities[pid]=score
        rows.append({"dataset":r.dataset,"sample":r["sample"],"time_hours":r.time_hours,**activities})
    a=pd.DataFrame(rows);a.to_csv(OUT/"05_fixed_program_activity.csv",index=False);return a

def run():
    log("loading common gene space");matrix,meta=load_space();log(f"datasets={','.join(TARGET)}; genes={matrix.shape[0]}; samples={matrix.shape[1]}")
    var=variance_table(matrix,meta);sets=feature_sets(matrix,var);results=[];pair_rows=[]
    for variant,genes in sets.items():
        log(f"evaluating {variant}: {len(genes)} genes");r,pairs=evaluate_variant(matrix,meta,genes,variant);results.append(r);pair_rows.extend([{**p,"variant":variant} for p in pairs])
    summary=pd.DataFrame(results);summary.to_csv(OUT/"03_harmonization_variant_comparison.csv",index=False);pd.DataFrame(pair_rows).to_csv(OUT/"04_matched_time_by_variant.csv",index=False);program_variant(matrix,meta)
    base=summary.loc[summary.variant=="baseline_all"].iloc[0];best=summary.sort_values(["mean_matched_time_pearson","knn_excess"],ascending=[False,True]).iloc[0]
    overall=pd.DataFrame([{"n_variants":len(summary),"baseline_matched_time_pearson":float(base.mean_matched_time_pearson),"best_variant":str(best.variant),"best_matched_time_pearson":float(best.mean_matched_time_pearson),"baseline_knn_excess":float(base.knn_excess),"best_knn_excess":float(best.knn_excess),"repair_improves_matched_time":bool(best.mean_matched_time_pearson>base.mean_matched_time_pearson+.05),"repair_reduces_dataset_knn_enrichment":bool(best.knn_excess<base.knn_excess-.02),"stage3_readiness":False,"interpretation":"compare variants; no automatic acceptance of a repaired state"}]);overall.to_csv(OUT/"07_STAGE2_9_20_SUMMARY.csv",index=False)
    log("complete");print("\nStage 2.9.20 variant comparison:");print(summary.to_string(index=False));print("\nStage 2.9.20 summary:");print(overall.to_string(index=False));return overall

if __name__=="__main__":run()
