#!/usr/bin/env python
# liana-py 三种口径在 CFA 子集六型上的通讯打分：
#   cellphonedb / cellchat（liana 重实现）/ rank_aggregate（共识）
# 输入：export/（build_seurat_sp5c.R 导出）
# 输出：liana_py_{cellphonedb|cellchat|consensus}_cfa.csv
import sys
import pandas as pd
import scanpy as sc
import liana as li
from liana.method import cellphonedb, cellchat, rank_aggregate

sc.settings.verbosity = 1
adata = sc.read_mtx("export/counts.mtx").T
adata.var_names = pd.read_csv("export/genes.csv", header=None)[0].astype(str).values
adata.obs_names = pd.read_csv("export/barcodes.txt", header=None)[0].values
meta = pd.read_csv("export/meta.csv", index_col=0)
adata.obs = meta.loc[adata.obs_names]
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
cfa = adata[adata.obs["condition"] == "CFA"].copy()
cfa.obs["label"] = cfa.obs["celltype"].astype(str)
cfa.raw = cfa  # liana 方法需要 .raw 取平均表达
print(f"CFA cells: {cfa.n_obs}, types: {sorted(cfa.obs['label'].unique())}", file=sys.stderr)

runs = [("cellphonedb", cellphonedb), ("cellchat", cellchat), ("consensus", rank_aggregate)]
for name, fn in runs:
    print(f"running {name} ...", file=sys.stderr, flush=True)
    try:
        res = fn(cfa, groupby="label", resource_name="mouseconsensus", expr_prop=0.1,
                 inplace=False, n_perms=100)
        res.to_csv(f"liana_py_{name}_cfa.csv", index=False)
        print(f"{name}: {res.shape[0]} rows -> liana_py_{name}_cfa.csv", file=sys.stderr)
    except Exception as e:
        print(f"{name} FAILED: {e}", file=sys.stderr)
print("LIANA-PY DONE")
