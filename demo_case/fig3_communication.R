#!/usr/bin/env Rscript
# 图3 ApoE 为中心通讯网络（A CFA 条件弦图 / B CFA-CON 折变 / C 配体-受体表达跨条件）
# 数据来源（只读）：
#   ../data_r/gse285766_apoe_by_celltype.csv  （Apoe 四胶质型 × 条件均值，云端 export_for_r.py §4）
#   ../data_r/gse285766_neuron_receptors.csv  （Lrp1/Sorl1/Ldlr/Vldlr × 条件均值，同 §4）
#   ../data_r/gse285766_mg_cells.csv          （MG Trem2/Csf1/Tyrobp/Csf1r × 条件）
#   ../tables/Table_communication_scores.csv  （折变 2.89/5.16/3.22/3.21/3.34）
suppressMessages({library(ggplot2); library(dplyr); library(readr); library(tidyr); library(circlize)})
source("theme_jtm.R"); set_figdir()
ZIP <- normalizePath("..")

apo <- read_csv(file.path(ZIP, "data_r/gse285766_apoe_by_celltype.csv"), show_col_types = FALSE)
rec <- read_csv(file.path(ZIP, "data_r/gse285766_neuron_receptors.csv"), show_col_types = FALSE)
mg  <- read_csv(file.path(ZIP, "data_r/gse285766_mg_cells.csv"), show_col_types = FALSE)
mg$condition <- factor(mg$condition, levels = names(col_condition))
receptors <- c("Lrp1", "Sorl1", "Ldlr", "Vldlr")
glia <- c("Astrocytes", "OPCs", "Microglia", "Oligodendrocytes")
abb <- c(Oligodendrocytes = "Oligo", Astrocytes = "Astro", OPCs = "OPC",
         Microglia = "MG", Neurons = "Neuron")   # 扇区/行名缩写（图注说明全称）

# ---- A：CFA 条件弦图（分数 = mean(L)×mean(R)，每条件各自均值）----
apo_cfa <- apo %>% filter(condition == "CFA") %>% select(celltype, apoe_mean)
recv <- rec %>% filter(condition == "CFA") %>%
  select(all_of(paste0(receptors, "_mean"))) %>% unlist() %>% as.numeric()
names(recv) <- receptors
links <- do.call(rbind, lapply(glia, function(g) {
  data.frame(from = unname(abb[g]), to = "Neuron", value = apo_cfa$apoe_mean[apo_cfa$celltype == g] * recv)
}))
mgc <- mg %>% filter(condition == "CFA") %>% summarise(Trem2 = mean(Trem2), Tyrobp = mean(Tyrobp),
                                                       Csf1 = mean(Csf1), Csf1r = mean(Csf1r))
links <- rbind(links,
  data.frame(from = "MG", to = "MG",
             value = c(Trem2 = mgc$Trem2 * mgc$Tyrobp, Csf1 = mgc$Csf1 * mgc$Csf1r)))
cat("CFA 条件通讯分数（与 liana_reconciliation §2 本机 CFA 列对账）:\n")
print(round(links$value, 3))

grid_col <- c(col_celltype[glia], Neurons = "#3498DB"); names(grid_col) <- abb[names(grid_col)]
sector_order <- c(abb[glia], "Neuron")
chord_title <- "A  CFA communication (mean(L) × mean(R))"
pdf(file.path(figdir, "Fig3A_chord_cfa.pdf"), width = 90/25.4, height = 90/25.4)
par(mar = c(2, 2, 2, 2))   # 边距 2 行：标题有上缘空间，圆不顶满画布（扇区名不被裁）
circos.clear()
chordDiagramFromDataFrame(links, order = sector_order, self.link = 1,
  grid.col = grid_col, annotationTrack = c("grid", "name"))
par(cex = 1)
mtext(chord_title, side = 3, line = 0.5, cex = 0.75, font = 2)  # 12pt×0.75=9pt 粗体，与 B/C 一致
circos.clear(); dev.off()
png(file.path(figdir, "Fig3A_chord_cfa.png"), width = 90, height = 90, units = "mm", res = 300)
par(mar = c(2, 2, 2, 2))
circos.clear()
chordDiagramFromDataFrame(links, order = sector_order, self.link = 1,
  grid.col = grid_col, annotationTrack = c("grid", "name"))
par(cex = 1)
mtext(chord_title, side = 3, line = 0.5, cex = 0.75, font = 2)  # 12pt×0.75=9pt 粗体，与 B/C 一致
circos.clear(); dev.off()
cat("saved Fig3A_chord_cfa\n")

# ---- B：CFA/CON 折变（来自 Table_communication_scores.csv 的 fold 行）----
tbl <- read_csv(file.path(ZIP, "tables/Table_communication_scores.csv"), show_col_types = FALSE)
fold <- tbl %>% filter(source == "fold") %>%
  mutate(axis = gsub("\n", " ", receptor), value = score) %>%
  arrange(desc(value)) %>% mutate(axis = factor(axis, levels = axis))
pB <- ggplot(fold, aes(value, axis)) +
  geom_col(fill = "#D55E00", width = 0.6) +
  geom_text(aes(label = sprintf("%.2f×", value)), hjust = -0.2, size = 3.2) +
  coord_cartesian(xlim = c(0, 6.2)) +
  labs(x = "CFA / CON fold change", y = NULL, title = "B  CFA/CON fold changes") +
  theme_jtm()

# ---- C：配体与受体表达跨条件（行 z 标准化）----
lig <- apo %>% mutate(row = paste0(abb[celltype], " Apoe")) %>% select(row, condition, value = apoe_mean)
mgl <- mg %>% group_by(condition) %>%
  summarise(Trem2 = mean(Trem2), Csf1 = mean(Csf1), Tyrobp = mean(Tyrobp), Csf1r = mean(Csf1r),
            .groups = "drop") %>%
  pivot_longer(-condition, names_to = "gene", values_to = "value") %>%
  mutate(row = paste0("MG ", gene))
rec_l <- rec %>% pivot_longer(-condition, names_to = "k", values_to = "value") %>%
  filter(k %in% paste0(receptors, "_mean")) %>%
  mutate(row = paste0("Neuron ", sub("_mean$", "", k))) %>% select(row, condition, value)
expr <- bind_rows(lig, mgl, rec_l)
expr <- expr %>% group_by(row) %>% mutate(z = as.numeric(scale(value))) %>% ungroup()
ord <- c(paste0(abb[glia], " Apoe"), "MG Trem2", "MG Csf1", "MG Tyrobp", "MG Csf1r",
         paste0("Neuron ", receptors))
expr$row <- factor(expr$row, levels = rev(ord))
expr$grp <- ifelse(grepl("^Neuron", expr$row), "Receptors", "Ligands")
pC <- ggplot(expr, aes(condition, row, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0,
                       breaks = pretty(range(expr$z), n = 3), name = "z-score") +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL, title = "C  Ligand & receptor expression") +
  theme_jtm() + theme(legend.position = "bottom") +
  guides(fill = guide_colorbar(barwidth = unit(40, "mm"), barheight = unit(2, "mm"),
                               title.position = "top"))

# 弦图为独立图形文件（Fig3A_chord_cfa），B/C 上下堆叠同图输出（与 A 并排拼接；
# B 55mm 5 行横条足够，C 12 行热图需 75mm 防行名重叠；全图文字统一 9pt）
save_fig(pB / pC + plot_layout(heights = c(55, 75)), "Fig3_BC_communication", w = 95, h = 130)
cat("图3B/C 完成（图3A 见 Fig3A_chord_cfa）\n")
