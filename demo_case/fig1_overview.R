#!/usr/bin/env Rscript
# 图1 SP5C 单核细胞图谱（A–G；G 20 基因点图）
# 数据来源（只读）：
#   ../data_r/gse285766_all_cells.csv   （36,290 细胞，UMAP 坐标/注释，由云端 gse285766_annotated.h5ad 经 export_for_r.py §1 导出）
#   ../data_r/gse285766_mg_cells.csv    （870 小胶质，含 dam_state，同 §2）
#   ../data_r/gse285766_markers_summary.csv（G 面板：云端 export_markers.py 2026-08-31 补导出，6 细胞型 × 20 基因 mean/pct）
suppressMessages({library(ggplot2); library(dplyr); library(readr); library(patchwork); library(scales); library(tidyr)})
source("theme_jtm.R"); set_figdir()
DATA <- normalizePath("../data_r")

ac <- read_csv(file.path(DATA, "gse285766_all_cells.csv"), show_col_types = FALSE)
ac$condition <- factor(ac$condition, levels = names(col_condition))
ac$celltype  <- factor(ac$celltype, levels = names(col_celltype))
mg <- read_csv(file.path(DATA, "gse285766_mg_cells.csv"), show_col_types = FALSE)
mg$condition <- factor(mg$condition, levels = names(col_condition))
mg$dam_state <- factor(mg$dam_state, levels = names(col_dam))

# 图内数值核对（对账用）
cat("cells:", nrow(ac), "| leiden clusters:", length(unique(ac$leiden)), "\n")
print(table(ac$condition)); print(table(ac$celltype))
cat("MG:", nrow(mg), "\n"); print(table(mg$condition, mg$dam_state))

# A UMAP by celltype / B by condition / C by sample
pA <- ggplot(ac, aes(umap1, umap2, colour = celltype)) +
  geom_point(size = 0.08) + scale_colour_manual(values = col_celltype, name = NULL) +
  labs(x = "UMAP1", y = "UMAP2", title = "A  Cell types (n = 36,290)") +
  theme_jtm() + theme(legend.position = "bottom", legend.key.size = unit(2.5, "mm")) +
  guides(colour = guide_legend(nrow = 3))
pB <- ggplot(ac, aes(umap1, umap2, colour = condition)) +
  geom_point(size = 0.08) + scale_colour_manual(values = col_condition, name = NULL) +
  labs(x = "UMAP1", y = "UMAP2", title = "B  Conditions") +
  theme_jtm() + theme(legend.position = "bottom", legend.key.size = unit(2.5, "mm"))
pC <- ggplot(ac, aes(umap1, umap2, colour = sample)) +
  geom_point(size = 0.08) + scale_colour_manual(values = col_sample, name = NULL) +
  labs(x = "UMAP1", y = "UMAP2", title = "C  Samples") +
  theme_jtm() + theme(legend.position = "bottom", legend.key.size = unit(2.5, "mm")) +
  guides(colour = guide_legend(nrow = 2))

# D 各条件细胞类型构成（比例堆积）
comp <- ac %>% count(condition, celltype) %>% group_by(condition) %>% mutate(pct = n / sum(n) * 100)
pD <- ggplot(comp, aes(condition, pct, fill = celltype)) +
  geom_col(width = 0.6) + scale_fill_manual(values = col_celltype, name = NULL) +
  labs(x = NULL, y = "Proportion (%)", title = "D  Cell-type composition") +
  theme_jtm() + theme(legend.position = "bottom", legend.key.size = unit(2.5, "mm")) +
  guides(fill = guide_legend(nrow = 3))

# E 小胶质 UMAP by DAM 状态 / F 星形胶质亚群 UMAP by leiden
pE <- ggplot(mg, aes(umap1, umap2, colour = dam_state)) +
  geom_point(size = 0.5) + scale_colour_manual(values = col_dam, name = NULL) +
  labs(x = "UMAP1", y = "UMAP2", title = "E  MG DAM states (n = 870)") +
  theme_jtm() + theme(legend.position = "bottom", legend.key.size = unit(2, "mm")) +
  guides(colour = guide_legend(nrow = 2))
astro <- ac %>% filter(celltype == "Astrocytes") %>% mutate(leiden = factor(leiden))
pF <- ggplot(astro, aes(umap1, umap2, colour = leiden)) +
  geom_point(size = 0.6) +
  scale_colour_manual(values = hue_pal()(length(unique(astro$leiden))), name = "Leiden") +
  guides(colour = guide_legend(nrow = 4)) +
  labs(x = "UMAP1", y = "UMAP2", title = sprintf("F  Astro subclusters (n = %d)", nrow(astro))) +
  theme_jtm() + theme(legend.position = "bottom", legend.key.size = unit(2, "mm"))

# G 20 个标记基因点图（每细胞型：均值=颜色，阳性%=点大小）
mk <- read_csv(file.path(DATA, "gse285766_markers_summary.csv"), show_col_types = FALSE)
genes <- c("Mbp","Mog","Plp1","Aqp4","Gfap","Slc1a3","Aldh1l1",
           "Pdgfra","Cspg4","Vcan","Cx3cr1","Tmem119","P2ry12","Aif1","Csf1r",
           "Snap25","Syt1","Rbfox3","Pdgfrb","Vtn")
mkL <- mk %>% pivot_longer(-celltype, names_to = "k", values_to = "v") %>%
  mutate(metric = ifelse(grepl("_pct$", k), "pct", "mean"),
         gene = sub("_(pct|mean)$", "", k)) %>%
  select(-k) %>%
  pivot_wider(names_from = metric, values_from = v)
mkL$gene <- factor(mkL$gene, levels = genes)
mkL$celltype <- factor(mkL$celltype, levels = names(col_celltype))
cat("G 面板核验：\n")
print(mkL %>% filter(gene %in% c("Cx3cr1","Tmem119","Snap25","Aqp4","Mbp")) %>%
        select(celltype, gene, mean, pct), n = 30)
pG <- ggplot(mkL, aes(celltype, gene, colour = mean, size = pct)) +
  geom_point() +
  scale_colour_gradient(low = "#f2f2f2", high = "#0072B2", name = "Mean expr.") +
  scale_size(range = c(0.4, 5), limits = c(0, 100), name = "% positive") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL, title = "G  Marker genes (20)") +
  theme_jtm() + theme(axis.text.x = element_text(angle = 45, hjust = 0),
                      legend.key.size = unit(3, "mm"))

Fig1 <- (pA + pB + pC) / (pD + pE + pF) / pG + plot_layout(heights = c(1, 1.02, 1.5))
save_fig(Fig1, "Fig1_overview", w = 185, h = 250)
cat("图1 A–G 完成\n")
