from pathlib import Path
import gzip
import tarfile

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
ARCHIVE = ROOT / "Data" / "GSE52052_RAW.tar"
OUT = ROOT / "results" / "GSE52052"
OUT.mkdir(parents=True, exist_ok=True)

extract = ROOT / "GSE52052_extracted"
extract.mkdir(exist_ok=True)
if not any(extract.rglob("*.txt.gz")):
    with tarfile.open(ARCHIVE, "r") as tar: tar.extractall(extract)

signals = []
for path in sorted(extract.rglob("*.txt.gz")):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()
    header_i = next(i for i, line in enumerate(lines) if line.startswith("FEATURES\t"))
    cols = lines[header_i].split("\t")[1:]
    rows = []
    for line in lines[header_i + 1:]:
        if line.startswith("*"): break
        if line.startswith("DATA\t"):
            fields = line.split("\t")[1:]
            if len(fields) == len(cols): rows.append(fields)
    df = pd.DataFrame(rows, columns=cols)
    ct = pd.to_numeric(df["ControlType"], errors="coerce")
    df = df[ct.fillna(0) == 0].copy()
    df["gProcessedSignal"] = pd.to_numeric(df["gProcessedSignal"], errors="coerce")
    df = df.dropna(subset=["gProcessedSignal"])
    df["ProbeName"] = df["ProbeName"].astype(str).str.strip()
    df = df[(df["ProbeName"] != "") & (df["ProbeName"] != "nan")]
    df.loc[df["gProcessedSignal"] < 0, "gProcessedSignal"] = 0
    sample = path.name[:-7]
    signals.append(df.groupby("ProbeName")["gProcessedSignal"].mean().rename(sample))

expr = pd.concat(signals, axis=1, join="inner")
expr.to_csv(OUT / "expression_raw_processed.csv")
log = np.log2(expr + 1)
log.to_csv(OUT / "expression_log2.csv")

fig, ax = plt.subplots(figsize=(12, 6)); ax.boxplot([log[c].dropna() for c in log], tick_labels=log.columns, showfliers=False); ax.set_title("GSE52052 - log2 processed signal"); ax.set_ylabel("log2(signal + 1)"); ax.tick_params(axis="x",rotation=45); fig.tight_layout(); fig.savefig(OUT/"01_boxplot.png",dpi=250); plt.close(fig)
qc = pd.DataFrame({"mean":expr.mean(),"median":expr.median(),"sd":expr.std(),"missing":expr.isna().sum()}); qc.to_csv(OUT/"02_sample_QC.csv")

x = log.to_numpy(float); order=np.argsort(x,axis=0); means=np.sort(x,axis=0).mean(axis=1); norm=np.empty_like(x)
for j in range(x.shape[1]): norm[order[:,j],j]=means
norm=pd.DataFrame(norm,index=log.index,columns=log.columns); norm.to_csv(OUT/"expression_quantile_normalized.csv")

corr=norm.corr(); corr.to_csv(OUT/"03_sample_correlation.csv")
fig,ax=plt.subplots(figsize=(8,7)); im=ax.imshow(corr,vmin=-1,vmax=1); ax.set_xticks(range(len(corr))); ax.set_xticklabels(corr.columns,rotation=45,ha="right"); ax.set_yticks(range(len(corr))); ax.set_yticklabels(corr.index); ax.set_title("GSE52052 - sample correlation"); fig.colorbar(im,ax=ax,label="Pearson r"); fig.tight_layout(); fig.savefig(OUT/"04_sample_correlation.png",dpi=250); plt.close(fig)

var=norm.var(axis=1).sort_values(ascending=False); var.to_csv(OUT/"05_probe_variance.csv",header=["variance"])
top=var.head(min(2000,len(var))).index; X=norm.loc[top].T; C=X.to_numpy()-X.to_numpy().mean(axis=0); U,S,_=np.linalg.svd(C,full_matrices=False); PC=U*S; EV=S**2/np.sum(S**2)
pd.DataFrame(PC[:,:min(5,PC.shape[1])],index=X.index,columns=[f"PC{i+1}" for i in range(min(5,PC.shape[1]))]).to_csv(OUT/"06_PCA_coordinates.csv")
if len(EV)>=2:
    fig,ax=plt.subplots(figsize=(8,7)); ax.scatter(PC[:,0],PC[:,1],s=90)
    for i,s in enumerate(X.index): ax.annotate(s,(PC[i,0],PC[i,1]),xytext=(6,6),textcoords="offset points",fontsize=8)
    ax.set_xlabel(f"PC1 ({EV[0]*100:.1f}%)"); ax.set_ylabel(f"PC2 ({EV[1]*100:.1f}%)"); ax.set_title("GSE52052 - PCA"); fig.tight_layout(); fig.savefig(OUT/"07_PCA.png",dpi=250); plt.close(fig)

features=var.head(min(50,len(var))).index; z=norm.loc[features]; z=z.sub(z.mean(axis=1),axis=0).div(z.std(axis=1).replace(0,np.nan),axis=0).fillna(0)
fig,ax=plt.subplots(figsize=(10,max(6,len(features)*.18))); im=ax.imshow(z,aspect="auto"); ax.set_xticks(range(len(z.columns))); ax.set_xticklabels(z.columns,rotation=45,ha="right"); ax.set_yticks(range(len(z.index))); ax.set_yticklabels(z.index,fontsize=6); ax.set_title("GSE52052 - top variable probes"); fig.colorbar(im,ax=ax,label="row z-score"); fig.tight_layout(); fig.savefig(OUT/"08_top_variable_probes.png",dpi=250); plt.close(fig)

report=["Dataset: GSE52052",f"Probes: {expr.shape[0]:,}",f"Samples: {expr.shape[1]:,}","Signal: Agilent gProcessedSignal","Control probes: removed","Duplicate ProbeName: averaged","Exploratory quantile normalization: applied"]
if len(EV)>=2: report += [f"PC1: {EV[0]*100:.2f}%",f"PC2: {EV[1]*100:.2f}%",f"PC1+PC2: {(EV[0]+EV[1])*100:.2f}%"]
(OUT/"REPORT.txt").write_text("\n".join(report),encoding="utf-8")
print("GSE52052 exploration complete")
