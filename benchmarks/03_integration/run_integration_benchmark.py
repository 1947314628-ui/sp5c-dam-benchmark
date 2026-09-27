#!/usr/bin/env python
# BIB 任务③：GSE285766 六样本整合方式比较（G26 云端运行）
# 方法：无校正基线(PCA) / Harmony / BBKNN / scVI（哪个可用跑哪个）
# 输入：gse285766_annotated.h5ad（36,290 细胞）
# 输出（下载回本地）：
#   integration_metrics.csv            —— 各方法 runtime/batch_ASW/ARI/NMI
#   integration_umap_{method}.csv.gz   —— UMAP 坐标+注释（FigB3 用）
import time, sys, os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

sc.settings.verbosity = 1
OUT = "integration_out"
os.makedirs(OUT, exist_ok=True)

adata = sc.read_h5ad("gse285766_annotated.h5ad")
print("obs columns:", list(adata.obs.columns), file=sys.stderr)
print("layers:", list(adata.layers.keys()), "| raw:", adata.raw is not None, file=sys.stderr)

def pick_col(cands):
    for c in cands:
        for col in adata.obs.columns:
            if col.lower() == c or c in col.lower():
                return col
    raise KeyError(f"no obs column for {cands}")

cond_col = pick_col(["condition", "treatment", "group"])
smpl_col = pick_col(["sample", "orig.ident", "batch"])
type_col = pick_col(["celltype", "annotation", "cell_type", "label"])
print(f"mapped: condition={cond_col} sample={smpl_col} celltype={type_col}", file=sys.stderr)

has_counts = "counts" in adata.layers or (adata.raw is not None)
if "counts" in adata.layers:
    counts = adata.layers["counts"]
elif adata.raw is not None:
    counts = adata.raw.to_adata().X
else:
    counts = None  # 只有 lognorm X 时 scVI 跳过

# 通用预处理：log1p + HVG + scale + PCA（无校正基线即 PCA 后直接聚类/UMAP）
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3", batch_key=smpl_col)
sc.tl.pca(adata, svd_solver="arpack", n_comps=50)
print(f"PCA done, HVG: {adata.var['highly_variable'].sum()}", file=sys.stderr)

def metric_block(label, X_emb, runtime, cluster_key):
    """batch ASW（0=混得好, 1=分得开）+ ARI/NMI vs celltype，存 CSV 行"""
    from sklearn.metrics import silhouette_score, adjusted_rand_score, normalized_mutual_info_score
    idx = np.random.RandomState(0).choice(len(adata), min(5000, len(adata)), replace=False)
    batch_ASW = silhouette_score(X_emb[idx], adata.obs[smpl_col].astype("category").cat.codes[idx])
    ari = adjusted_rand_score(adata.obs[type_col], adata.obs[cluster_key])
    nmi = normalized_mutual_info_score(adata.obs[type_col], adata.obs[cluster_key])
    print(f"[{label}] runtime={runtime:.1f}s batch_ASW={batch_ASW:.3f} ARI={ari:.3f} NMI={nmi:.3f}", file=sys.stderr)
    return dict(method=label, runtime_s=round(runtime, 1),
                batch_ASW=round(batch_ASW, 4), ARI=round(ari, 4), NMI=round(nmi, 4))

def finish(method, X_emb, runtime):
    key = f"leiden_{method}"
    sc.pp.neighbors(adata, use_rep=X_emb)
    sc.tl.leiden(adata, resolution=0.5, key_added=key)
    sc.tl.umap(adata, min_dist=0.3, spread=1.0)
    out = pd.DataFrame({
        "umap1": adata.obsm["X_umap"][:, 0], "umap2": adata.obsm["X_umap"][:, 1],
        "condition": adata.obs[cond_col].values, "sample": adata.obs[smpl_col].values,
        "celltype": adata.obs[type_col].values, "cluster": adata.obs[key].values})
    out.to_csv(f"{OUT}/integration_umap_{method}.csv.gz", index=False)
    return metric_block(method, adata.obsm[X_emb], runtime, key)

metrics = []
# 0) 无校正基线
t0 = time.time()
finish("PCA", "X_pca", time.time() - t0)

# 1) Harmony
try:
    import harmonypy
    t0 = time.time()
    hm = harmonypy.run_harmony(adata.obsm["X_pca"][:, :30], adata.obs, smpl_col)
    adata.obsm["X_harmony"] = np.array(hm.Z_corr.T)
    metrics.append(finish("Harmony", "X_harmony", time.time() - t0))
except Exception as e:
    print("Harmony FAILED:", e, file=sys.stderr)

# 2) BBKNN（图整合，neighbors 层面）
try:
    import bbknn
    t0 = time.time()
    sc.pp.neighbors(adata, use_rep="X_pca", n_pcs=30)
    bbknn.bbknn(adata, batch_key=smpl_col, n_pcs=30)
    sc.tl.leiden(adata, resolution=0.5, key_added="leiden_bbknn")
    sc.tl.umap(adata, min_dist=0.3, spread=1.0)
    out = pd.DataFrame({
        "umap1": adata.obsm["X_umap"][:, 0], "umap2": adata.obsm["X_umap"][:, 1],
        "condition": adata.obs[cond_col].values, "sample": adata.obs[smpl_col].values,
        "celltype": adata.obs[type_col].values, "cluster": adata.obs["leiden_bbknn"].values})
    out.to_csv(f"{OUT}/integration_umap_bbknn.csv.gz", index=False)
    # BBKNN 输出邻居图而非低维嵌入，ASW 用 PCA 坐标（与聚类同一邻居空间）近似
    metrics.append(metric_block("BBKNN", "X_pca", time.time() - t0, "leiden_bbknn"))
except Exception as e:
    print("BBKNN FAILED:", e, file=sys.stderr)

# 3) scVI（需要原始 counts：优先 layers["counts"]，否则 adata.raw）
try:
    if "counts" in adata.layers:
        target, setup = adata, dict(layer="counts")
    elif adata.raw is not None:
        target = ad.AnnData(X=adata.raw.X, obs=adata.obs, var=adata.raw.var)
        setup = dict(layer=None)
    else:
        raise ValueError("no counts layer nor raw")
    import scvi
    t0 = time.time()
    scvi.model.SCVI.setup_anndata(target, batch_key=smpl_col, **setup)
    model = scvi.model.SCVI(target, n_latent=30, n_layers=2)
    model.train(max_epochs=100, early_stopping=True)
    adata.obsm["X_scvi"] = model.get_latent_representation()
    metrics.append(finish("scVI", "X_scvi", time.time() - t0))
except Exception as e:
    print("scVI FAILED:", e, file=sys.stderr)

pd.DataFrame(metrics).to_csv(f"{OUT}/integration_metrics.csv", index=False)
print("INTEGRATION DONE")
print("metrics:"); print(pd.DataFrame(metrics).to_string())
