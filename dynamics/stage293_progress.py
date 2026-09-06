"""Fast Stage 2.9.3 diagnostics for latent-progress stability."""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_9_3"
OUT.mkdir(parents=True, exist_ok=True)
SRC = ROOT / "results" / "Dynamics" / "stage2_9_2" / "02_latent_progress_trajectories.csv"

def corr(a,b,method="spearman"):
    a=np.asarray(a,float); b=np.asarray(b,float); ok=np.isfinite(a)&np.isfinite(b)
    if ok.sum()<3:return np.nan
    return float(pd.Series(a[ok]).corr(pd.Series(b[ok]),method=method))

def run(bootstrap_replicates=200,permutations=1000,seed=42):
    if not SRC.exists(): raise RuntimeError("Run Stage 2.9.2 first.")
    df=pd.read_csv(SRC)
    rng=np.random.default_rng(seed); boots=[]
    groups=list(df.groupby("held_out_dataset",sort=True))
    print(f"Stage 2.9.3: found {len(groups)} held-out datasets",flush=True)
    for di,(ds,g) in enumerate(groups,1):
        g=g.dropna(subset=["normalized_time","latent_progress"]); a=g.normalized_time.to_numpy(float); z=g.latent_progress.to_numpy(float)
        print(f"Stage 2.9.3: bootstrap dataset {di}/{len(groups)}: {ds} ({bootstrap_replicates})",flush=True)
        for b in range(bootstrap_replicates):
            idx=rng.integers(0,len(g),len(g)); boots.append({"held_out_dataset":ds,"bootstrap":b,"spearman":corr(a[idx],z[idx]),"pearson":corr(a[idx],z[idx],"pearson")})
    boot=pd.DataFrame(boots); boot.to_csv(OUT/"02_bootstrap_stability.csv",index=False)
    sm=boot.groupby("held_out_dataset").agg(n_bootstrap=("bootstrap","size"),bootstrap_mean_spearman=("spearman","mean"),bootstrap_p05_spearman=("spearman",lambda x:x.quantile(.05)),bootstrap_p95_spearman=("spearman",lambda x:x.quantile(.95)),bootstrap_mean_pearson=("pearson","mean"),bootstrap_p05_pearson=("pearson",lambda x:x.quantile(.05)),bootstrap_p95_pearson=("pearson",lambda x:x.quantile(.95))).reset_index()
    sm.to_csv(OUT/"03_bootstrap_stability_summary.csv",index=False); print(sm.to_string(index=False),flush=True)
    print(f"Stage 2.9.3: permutation null {permutations} replicates",flush=True); null=[]
    for p in range(permutations):
        ss=[];pp=[]
        for _,g in groups:
            a=g.normalized_time.to_numpy(float);z=g.latent_progress.to_numpy(float); perm=rng.permutation(len(z));ss.append(corr(a,z[perm]));pp.append(corr(a,z[perm],"pearson"))
        null.append({"permutation":p,"mean_spearman":np.nanmean(ss),"mean_pearson":np.nanmean(pp)})
        if (p+1)%25==0 or p==0: print(f"Stage 2.9.3: permutation {p+1}/{permutations}",flush=True)
    null=pd.DataFrame(null);null.to_csv(OUT/"04_permutation_null.csv",index=False)
    obs_s=float(np.nanmean([corr(g.normalized_time,g.latent_progress) for _,g in groups]));obs_p=float(np.nanmean([corr(g.normalized_time,g.latent_progress,"pearson") for _,g in groups]))
    ps=(1+(null.mean_spearman>=obs_s).sum())/(len(null)+1);pp=(1+(null.mean_pearson>=obs_p).sum())/(len(null)+1)
    out=pd.DataFrame([{"n_datasets":len(groups),"n_cases":len(df),"observed_mean_spearman":obs_s,"observed_mean_pearson":obs_p,"permutation_p_spearman":ps,"permutation_p_pearson":pp,"bootstrap_replicates":bootstrap_replicates,"permutations":len(null)}]);out.to_csv(OUT/"01_overall_summary.csv",index=False)
    print(f"Observed mean Spearman={obs_s:.3f}, Pearson={obs_p:.3f}",flush=True);print(f"Permutation p: Spearman={ps:.4f}, Pearson={pp:.4f}",flush=True);print("Stage 2.9.3 complete.",flush=True)
    return out

if __name__=="__main__":run()
