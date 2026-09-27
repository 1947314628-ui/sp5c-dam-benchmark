# 任务③ 四法整合对比 —— 输出说明

生成时间: 2026-09-02 11:24:14（G26，192 核 CPU，无 GPU）
输入: `/mnt/DATA/home/s20221110329/SP5C/analysis/gse285766_clustered.h5ad`（36,290 细胞 × 25,644 基因，X=1e4+ln(1+x)）

## 统一下游协议（四法一致）
kNN k=15 → UMAP min_dist=0.3 → Leiden res=0.8（随机种子 42）；与 08-31 手稿重跑协议一致。
差异仅在上游表示: harmony=X_pca_harmony / baseline=X_pca(50PC) / bbknn=X_pca+批量kNN / scvi=scVI 隐空间(n_latent=30)。

## 各方法要点
- **harmony**: 复用 08-31 权威结果（38 簇），未重跑。
- **baseline**: 直接对 50PC 建图，无批次校正（展示批次效应基线）。
- **bbknn**: batch=sample, neighbors_within_batch=3 → 总 k=18（方法固有，非 15）。
- **scvi**: 3,000 HVG raw counts 训练（batch=sample），n_latent=20/n_layers=1（16 核 CPU 配额下取常规轻量配置），max_epochs=200/早停 patience=25/train_size=0.9/batch_size=2048；最终训练指标见 `scvi_train_history.csv`（scvi-tools 1.4 仅保留最终值）。

## 指标解读
- `ari_vs_sample`: Leiden 簇与 6 个样本标签的调整兰德指数 —— **越低批次混合越好**。
- `ari_vs_condition`: 簇与 CON/CFA/TB 的调整兰德指数 —— 越高生物信号保留越好。
- `entropy_sample`: 每簇内 6 样本归一化香农熵（1=完全混合，0=纯单样本簇）；`wmean_entropy_sample` 为按簇大小加权均值。
- `top_sample_share`: 每簇最大样本占比。

## 环境
- scanpy: 1.11.5
- anndata: 0.12.19
- scvi-tools: 1.4.2
- torch: 2.13.0+cu130
- bbknn: 1.6.0
- python: 3.11.2

## 文件
- `cells_umap_cluster.csv.gz`: 4 方法 × 36,290 细胞长表（UMAP 坐标/簇/样本/条件）
- `clusters_summary.csv`: 每方法每簇（大小/条件构成/批次熵/最大样本占比）
- `methods_summary.csv`: 方法级汇总（簇数/aRI/熵/耗时）
- `scvi_train_history.csv`: scVI 训练损失曲线
