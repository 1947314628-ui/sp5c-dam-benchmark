#!/usr/bin/env Rscript
# 标准 CellChat（jinworks/CellChat master）在 CFA 子集六型上的通讯分析
# 输入：sp5c_seurat.rds（build_seurat_sp5c.R 产出）
# 输出：cellchat_cfa_pairs.csv + cellchat_cfa.rds
suppressMessages({library(Seurat); library(CellChat); library(dplyr); library(readr)})
seu <- readRDS("sp5c_cfa_seurat.rds")
Idents(seu) <- "celltype"
cat("CFA 细胞:", ncol(seu), "| 类型:", paste(unique(seu$celltype), collapse = ","), "\n")

cc <- createCellChat(seu, group.by = "celltype")
cc@DB <- CellChatDB.mouse
cc <- subsetData(cc)
cc <- identifyOverExpressedGenes(cc)
cc <- identifyOverExpressedInteractions(cc)
cat("计算通讯概率（triMean）...\n")
cc <- computeCommunProb(cc, type = "triMean")
cc <- filterCommunication(cc, min.cells = 10)
df <- subsetCommunication(cc)
write_csv(df, "cellchat_cfa_pairs.csv")
saveRDS(cc, "cellchat_cfa.rds")
cat("CellChat 完成:", nrow(df), "对 | 正交互:", sum(df$pval < 0.05), "\n")
cat("显著对按 pval 前 5:\n")
print(head(df %>% arrange(pval) %>% select(source, target, ligand, receptor, pval), 5))
