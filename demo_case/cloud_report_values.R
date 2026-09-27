# 云端报告值（绘图脚本不得硬编码凭空数值；本文件数值全部来自 zip 内报告文件，注明路径:行号）
# zip 根：D:\01_研究项目\公共单细胞_SP5C\投稿包\拟投稿_SP5C\
#
# Slc17a6（VGLUT2）神经元阳性比例与均值
# 来源：reports/liana_reconciliation.md 第51行
#   | Neuron Slc17a6+% / 均值 | 39.8% / 0.50/0.50/0.63 | 41.4% / 0.518/0.556/0.686 | ✓ |
slc17a6_pct_overall <- 41.4
slc17a6_means <- c(CON = 0.518, CFA = 0.556, TB = 0.686)
