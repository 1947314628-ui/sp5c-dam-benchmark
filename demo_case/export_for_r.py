#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 任务 #19：从幸存 h5ad 导出 R 绘图 tidy CSV + 重算 XGBoost/SHAP/Scn9a
import numpy as np, pandas as pd, anndata as ad, os, shutil
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.metrics import roc_curve, roc_auc_score
import xgboost as xgb
OUT = '/mnt/DATA/home/s20221110329/拟投稿_SP5C/data_r'
SRC = '/mnt/DATA/home/s20221110329/jiaoyan/公共单细胞/analysis'
os.makedirs(OUT, exist_ok=True)

a = ad.read_h5ad(f'{SRC}/gse285766_annotated.h5ad', backed='r')
X = a.X
um = a.obsm['X_umap']
obs = a.obs

# ---- 1. 全细胞表（Fig1/S1/S2）----
ac = pd.DataFrame({
    'cell': obs.index, 'sample': obs['sample'].values, 'condition': obs['condition'].values,
    'celltype': obs['celltype'].values, 'leiden': obs['leiden'].values,
    'n_genes': obs['n_genes'].values, 'n_counts': obs['n_counts'].values,
    'pct_mt': obs['pct_counts_mt'].values, 'umap1': um[:,0], 'umap2': um[:,1]})
ac.to_csv(f'{OUT}/gse285766_all_cells.csv', index=False)
print('all_cells', ac.shape, '| celltype:', sorted(ac.celltype.unique()))

# ---- 2. 小胶质子集 + DAM 状态（Fig2）----
mg = ac[ac.celltype == 'Microglia'].copy()
genes = ['Apoe','Lpl','Itgax','Clec7a','Axl','Trem2','Cst7','Ccl6','Cd9','Spp1','Csf1','Tyrobp','B2m','Ctsb',
         'Tmem119','P2ry12','Cx3cr1','Hexb','C1qa','Olfml3','Siglech','Selplg','Scn9a','Csf1r','Lrp1','Sorl1','Ldlr','Vldlr']
avail = [g for g in genes if g in a.var_names]
pos = {c: i for i, c in enumerate(obs.index)}
mg_idx = [pos[c] for c in mg.cell.values]
G = X[mg_idx][:, [a.var_names.get_loc(g) for g in avail]].toarray()
for i, g in enumerate(avail):
    mg[g] = G[:, i]
mg['dam_score'] = obs.dam_score.values[mg_idx]
q30, q70 = np.quantile(mg.dam_score.dropna(), [0.30, 0.70])
mg['dam_state'] = pd.cut(mg.dam_score, [-np.inf, q30, q70, np.inf], labels=['Homeostatic','Intermediate','DAM-like'])
mg.to_csv(f'{OUT}/gse285766_mg_cells.csv', index=False)
print('MG', mg.shape, '| 状态比例:', (mg.dam_state.value_counts(normalize=True)*100).round(1).to_dict())
print('dam_score q30/q70:', round(q30,4), round(q70,4))
print('缺失基因:', [g for g in genes if g not in a.var_names])

# ---- 3. 神经元 Scn9a（Fig5）----
neu = ac[ac.celltype == 'Neurons']
neu_idx = [pos[c] for c in neu.cell.values]
scn = X[neu_idx][:, a.var_names.get_loc('Scn9a')].toarray().ravel()
neu['Scn9a'] = scn
scn_sum = neu.groupby('condition', observed=True).agg(scn9a_mean=('Scn9a','mean'), scn9a_pct=('Scn9a', lambda s:(s>0).mean()*100), n=('cell','size'))
scn_sum.round(3).to_csv(f'{OUT}/gse285766_scn9a_by_condition.csv')
print('Scn9a per condition:\n', scn_sum.round(3).to_string())

# ---- 4. Apoe 受体表达（神经元）与配体（胶质）Fig3 ----
rec = ['Lrp1','Sorl1','Ldlr','Vldlr']
rec_df = pd.DataFrame({r: X[neu_idx][:, a.var_names.get_loc(r)].toarray().ravel() for r in rec})
rec_df['condition'] = neu.condition.values
rec_sum = rec_df.groupby('condition', observed=True).agg(**{f'{r}_mean':(r,'mean') for r in rec}, **{f'{r}_pct':(r, lambda s:(s>0).mean()*100) for r in rec})
rec_sum.round(3).to_csv(f'{OUT}/gse285766_neuron_receptors.csv')
print('neuron receptors:\n', rec_sum.round(2).to_string())
lig_df = ac[ac.celltype.isin(['Astrocytes','OPCs','Microglia','Oligodendrocytes'])].copy()
lig_idx = [pos[c] for c in lig_df.cell.values]
lig_df['Apoe'] = X[lig_idx][:, a.var_names.get_loc('Apoe')].toarray().ravel()
lig_sum = lig_df.groupby(['celltype','condition'], observed=True).agg(apoe_mean=('Apoe','mean'), apoe_pct=('Apoe',lambda s:(s>0).mean()*100), n=('cell','size'))
lig_sum.round(3).to_csv(f'{OUT}/gse285766_apoe_by_celltype.csv')
print('Apoe by celltype:\n', lig_sum.round(2).to_string())

# ---- 5. XGBoost + SHAP（Fig4）----
up13 = ['Apoe','Lpl','Itgax','Axl','Trem2','Cst7','Ccl6','Cd9','Spp1','Csf1','Tyrobp','B2m','Ctsb']
ho8  = ['Tmem119','P2ry12','Cx3cr1','Hexb','C1qa','Olfml3','Siglech','Selplg']
feat = [g for g in up13+ho8 if g in a.var_names]
lab = mg[mg.dam_state != 'Intermediate'].copy()
y = (lab.dam_state == 'DAM-like').astype(int).values
F = lab[feat].values
print('XGBoost 样本:', len(lab), '阳性:', y.sum(), '特征:', len(feat))
skf = StratifiedKFold(5, shuffle=True, random_state=42)
aucs = []
for tr, te in skf.split(F, y):
    m = xgb.XGBClassifier(n_estimators=200, max_depth=6, learning_rate=0.1, subsample=0.8,
                          colsample_bytree=0.8, eval_metric='auc', n_jobs=64, random_state=42)
    m.fit(F[tr], y[tr]); aucs.append(roc_auc_score(y[te], m.predict_proba(F[te])[:,1]))
print('5折 AUROC:', [round(v,4) for v in aucs], 'mean', round(np.mean(aucs),4))
m = xgb.XGBClassifier(n_estimators=200, max_depth=6, learning_rate=0.1, subsample=0.8,
                      colsample_bytree=0.8, eval_metric='auc', n_jobs=64, random_state=42)
Ftr, Fte, ytr, yte = train_test_split(F, y, test_size=0.3, stratify=y, random_state=42)
m.fit(Ftr, ytr)
fpr, tpr, _ = roc_curve(yte, m.predict_proba(Fte)[:,1])
pd.DataFrame({'fpr':fpr,'tpr':tpr}).to_csv(f'{OUT}/xgboost_roc.csv', index=False)
import shap
ex = shap.TreeExplainer(m)
sv = np.asarray(ex(Fte).values)
if sv.ndim == 3: sv = sv[:, :, 1]
imp = pd.DataFrame({'gene':feat,'shap_importance':np.abs(sv).mean(0)}).sort_values('shap_importance', ascending=False)
imp.to_csv(f'{OUT}/xgboost_shap.csv', index=False)
print('SHAP top10:\n', imp.head(10).round(4).to_string(index=False))

# ---- 6. 复制既有 CSV ----
for f in ['gse172167_dam_scores.csv','gse285766_mg_pseudobulk.csv','gse162807_dam_per_sample.csv',
          'gse162807_dam_summary.csv','liana_mouseconsensus_6types.csv','deg_damlike_vs_homeostatic.csv',
          'deg_pure_mg_damlike_vs_homeostatic.csv','gse162807_dam_scores.csv']:
    p = f'{SRC}/{f}'
    if os.path.exists(p): shutil.copy(p, f'{OUT}/{f}')
# gse172167 条件汇总（Apoe/模块，前次已算）
g17 = pd.DataFrame({
 'condition':['A_Uninj','B_1dpi','C_1wpi','D_3wpi','E_6wpi'],
 'apoe_mean':[0.914,1.072,4.734,4.525,4.376],
 'up13_mean':[0.546,0.689,1.184,0.988,1.116], 'ho8_mean':[1.759,1.230,1.349,1.059,1.287]})
g17.to_csv(f'{OUT}/gse172167_condition_summary.csv', index=False)
print('\n导出完成 →', OUT, '| 文件数:', len(os.listdir(OUT)))
