# -*- coding: utf-8 -*-
# 评审修复计算：表1 SD 列、DEG 效应量分布、LIANA Apoe 轴聚合排名、表S2 26 样本全表
# 数据只读：D:/01_研究项目/公共单细胞_SP5C/投稿包/拟投稿_SP5C/data_r/ 与 tables/
# 输出：本目录 评审补缺数值_2026-08-31.txt + 表S2_26样本_markdown.txt
import pandas as pd, os

DATA = r"D:/01_研究项目/公共单细胞_SP5C/投稿包/拟投稿_SP5C/data_r"
TAB  = r"D:/01_研究项目/公共单细胞_SP5C/投稿包/拟投稿_SP5C/tables"
out  = []

def p(*a): out.append(" ".join(str(x) for x in a))

# ---- 1. 表1：每样本 DAM 评分 SD（细胞间） ----
mg = pd.read_csv(os.path.join(DATA, "gse285766_mg_cells.csv"))
p("== 1. 表1 补列：每样本 dam_score SD ==")
sds = mg.groupby("sample")["dam_score"].agg(["mean", "std", "count"]).round(4)
p(sds.to_string())
p("(对照 pseudobulk CSV 核对：)\n")
pb = pd.read_csv(os.path.join(TAB, "gse285766_mg_pseudobulk.csv"))
p(pb.to_string())

# ---- 2. DEG 效应量分布 ----
deg = pd.read_csv(os.path.join(DATA, "deg_damlike_vs_homeostatic.csv"))
up = deg[(deg.padj < 0.05) & (deg.log2fc > 0)]
dn = deg[(deg.padj < 0.05) & (deg.log2fc < 0)]
sig = deg[deg.padj < 0.05]
p("\n== 2. DEG 效应量分布（278 显著基因）==")
p("n=", len(sig), "上调", len(up), "下调", len(dn))
for nm, d in [("全部显著", sig), ("上调", up), ("下调", dn)]:
    a = d.log2fc.abs()
    p(f"{nm}: |log2FC| 中位 {a.median():.2f}  IQR {a.quantile(.25):.2f}-{a.quantile(.75):.2f}  范围 {a.min():.2f}-{a.max():.2f}")

# ---- 3. LIANA：Apoe 轴在 rank_aggregate 中的位置 ----
li = pd.read_csv(os.path.join(DATA, "liana_mouseconsensus_6types.csv"))
p("\n== 3. LIANA mouseconsensus CSV 概览 ==")
p("行数:", len(li), "列:", list(li.columns))
p("source 取值:", sorted(li.source.unique()))
p("target 取值:", sorted(li.target.unique()))
p("ligand 取值数:", li.ligand_complex.nunique())
apo = li[li.ligand_complex.str.contains("Apoe", na=False)]
p("Apoe 配体对行数:", len(apo))
rec = ["Lrp1", "Sorl1", "Ldlr", "Vldlr"]
axis = apo[apo.receptor_complex.isin(rec)]
p("\nApoe->Lrp1/Sorl1/Ldlr/Vldlr 各行（合并口径）:")
cols = ["source", "target", "ligand_complex", "receptor_complex", "lr_means", "expr_prod", "magnitude_rank"]
p(axis[cols].sort_values(["receptor_complex", "source"]).to_string())
p("\nmagnitude_rank 分布（全部对）: 最小", li.magnitude_rank.min(), "中位", li.magnitude_rank.median())
p("Apoe 轴 magnitude_rank 分位数（越小越好）:")
p(axis.magnitude_rank.rank(pct=True).round(3).to_string())
p("全部对中 magnitude_rank 小于 Apoe 轴最小值的对数:",
  int((li.magnitude_rank < axis.magnitude_rank.min()).sum()))
p("Apoe 轴 magnitude_rank 最小值:", axis.magnitude_rank.min(),
  "最大值:", axis.magnitude_rank.max())

# ---- 4. 表S2：26 样本全表 markdown ----
per = pd.read_csv(os.path.join(DATA, "gse162807_dam_per_sample.csv"))
p("\n== 4. gse162807 每样本 ==")
p("行数:", len(per), "列:", list(per.columns))
md = ["| 样本 | 性别 | 手术 | 时点 | n | DAM 评分均值 | DAM-like% | 稳态% | Apoe 均值 |",
      "|---|---|---|---|---|---|---|---|---|"]
for _, r in per.iterrows():
    md.append(f"| {r['sample']} | {r.sex} | {r.surgery} | {r.tp} | {int(r.n)} | {r.dam_mean:.3f} | {r['damlike%']:.1f} | {r['homeo%']:.1f} | {r.apoe_mean:.3f} |")
with open(os.path.join(TAB, "表S2_26样本_markdown.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(md))
p("表S2 26 样本 markdown 已写 tables/表S2_26样本_markdown.txt")

with open(os.path.join(TAB, "评审补缺数值_2026-08-31.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(out))
print("done")
