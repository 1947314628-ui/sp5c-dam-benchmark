#!/usr/bin/env Rscript
# 图S3 DEG 火山图 + SHAP 重要性 + DAM 基因热图（原图2 热图面板并入此图）
# 数据来源（只读）：
#   ../data_r/deg_damlike_vs_homeostatic.csv   （DEG 278，云端复制）
#   ../data_r/deg_pure_mg_damlike_vs_homeostatic.csv （纯小胶质 DEG，用于 Apoe 3.03 标注核验）
#   ../data_r/xgboost_shap.csv                 （云端 export_for_r.py §5 重算）
#   ../data_r/gse285766_mg_cells.csv           （21 DAM 基因 × 870 细胞 → 6 样本 pseudobulk）
# C 面板热图用 ggplot2 geom_tile 绘制（ComplexHeatmap 在本机与 patchwork 组合时行名丢失/布局错乱，弃用；
# 左侧类别条与右侧热图为一行两格、共享 y 因子级，顶部条件色带 = bar 行 + annotate rect，全部确定性渲染）
suppressMessages({library(ggplot2); library(dplyr); library(readr); library(tidyr)})
source("theme_jtm.R"); set_figdir()
DATA <- normalizePath("../data_r")

deg <- read_csv(file.path(DATA, "deg_damlike_vs_homeostatic.csv"), show_col_types = FALSE)
degp <- read_csv(file.path(DATA, "deg_pure_mg_damlike_vs_homeostatic.csv"), show_col_types = FALSE)
shap <- read_csv(file.path(DATA, "xgboost_shap.csv"), show_col_types = FALSE)
mg <- read_csv(file.path(DATA, "gse285766_mg_cells.csv"), show_col_types = FALSE)

cat("DEG 行数:", nrow(deg), "| 上调:", sum(deg$log2fc > 0 & deg$padj < 0.05),
    "| 下调:", sum(deg$log2fc < 0 & deg$padj < 0.05), "\n")
cat("纯小胶质 DEG 行数:", nrow(degp), "\n")
if ("Apoe" %in% degp$gene) print(degp %>% filter(gene == "Apoe"))
cat("SHAP top8:\n"); print(head(shap, 8))

# ---- 火山图 ----
deg <- deg %>% mutate(
  neglogp = -log10(padj),
  cl = ifelse(padj < 0.05 & log2fc > 0, "up",
              ifelse(padj < 0.05 & log2fc < 0, "down", "ns")))
apoe_row <- deg %>% filter(gene == "Apoe")
pVol <- ggplot(deg, aes(log2fc, neglogp, colour = cl)) +
  geom_point(size = 0.7) +
  scale_colour_manual(values = c(up = "#D73027", down = "#4575B4", ns = "grey80"), guide = "none") +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3) +
  {if (nrow(apoe_row) > 0) geom_text(data = apoe_row, aes(label = sprintf("Apoe\n(%.2f)", log2fc)),
                                     colour = "black", size = 3.2, vjust = 1.6) else NULL} +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(x = "log2 fold change (DAM-like vs homeostatic)",
       y = "-log10(adjusted P)",
       title = "A") +
  theme_jtm()

# ---- SHAP top15 ----
sh <- shap %>% head(15) %>% mutate(gene = factor(gene, levels = rev(gene)),
                                   is_apoe = gene == "Apoe")
pSh <- ggplot(sh, aes(shap_importance, gene, fill = is_apoe)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "#D55E00"), guide = "none") +
  labs(x = "Mean |SHAP value|", y = NULL,
       title = "B") +
  theme_jtm()

# ---- C：DAM 基因热图：21 基因 × 6 样本 pseudobulk z 分 ----
up13 <- c("Apoe","Lpl","Itgax","Axl","Trem2","Cst7","Ccl6","Cd9","Spp1","Csf1","Tyrobp","B2m","Ctsb")
ho8  <- c("Tmem119","P2ry12","Cx3cr1","Hexb","C1qa","Olfml3","Siglech","Selplg")
genes21 <- c(up13, ho8)
samples6 <- c("CON1","CON2","CFA1","CFA2","TB1","TB2")
pbmat <- sapply(genes21, function(g) {
  sapply(split(mg[[g]], factor(mg$sample, levels = samples6)), mean)
})
zmat <- t(apply(pbmat, 2, scale))
rownames(zmat) <- genes21
colnames(zmat) <- samples6
cat("z 分数范围:", round(range(zmat), 3), "| max|z|:", round(max(abs(zmat)), 3), "\n")

ylevels <- c(rev(genes21), "bar")   # "bar" 行在顶部（条件色带）
zdf <- as.data.frame(zmat) %>%
  mutate(gene = factor(genes21, levels = rev(genes21)),
         Class = factor(c(rep("Upregulated", length(up13)), rep("Homeostatic", length(ho8))),
                        levels = c("Upregulated", "Homeostatic")))
zlong <- zdf %>%
  pivot_longer(-c(gene, Class), names_to = "sample", values_to = "z") %>%
  mutate(sample = factor(sample, levels = samples6), y = as.character(gene))
barrow <- data.frame(sample = factor(samples6, levels = samples6),
                     y = "bar", z = NA_real_, Class = NA_character_)
zlong <- bind_rows(zlong, barrow) %>%
  mutate(y = factor(y, levels = ylevels),
         Class = factor(Class, levels = c("Upregulated", "Homeostatic")))

# 左缘类别条（红=DAM 上调 / 蓝=稳态），标题 C 在条左上方（与 A/B 同款粗体左对齐）
pC_cls <- ggplot(zlong, aes(factor(""), y, fill = Class)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_manual(values = c(Upregulated = "#D73027", Homeostatic = "#4575B4"),
                    name = NULL, na.value = "white") +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "C") +
  theme_jtm() + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                      axis.line = element_blank(),
                      legend.key.size = unit(2.5, "mm"),
                      plot.margin = margin(4, 0, 0, 0, "mm"))

# 热图主体：行名 = y 轴文字（8pt），列名 = x 轴文字（9pt），顶部条件色带 = bar 行 annotate
pC_body <- ggplot(zlong, aes(sample, y, fill = z)) +
  geom_tile(width = 1, height = 1) +
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = 21.5, ymax = 22.5,
           fill = col_condition[["CON"]]) +
  annotate("rect", xmin = 2.5, xmax = 4.5, ymin = 21.5, ymax = 22.5,
           fill = col_condition[["CFA"]]) +
  annotate("rect", xmin = 4.5, xmax = 6.5, ymin = 21.5, ymax = 22.5,
           fill = col_condition[["TB"]]) +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027",
                       midpoint = 0, limits = c(-1.5, 1.5), oob = scales::squish,
                       na.value = "white",
                       name = "z", breaks = c(-1.5, 0, 1.5),
                       guide = guide_colorbar(barwidth = unit(30, "mm"), barheight = unit(3, "mm"),
                                              direction = "horizontal", title.position = "top")) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(breaks = ylevels,
                   labels = c(bar = "", setNames(rev(genes21), rev(genes21)))) +
  labs(x = NULL, y = NULL) +
  theme_jtm() + theme(axis.text.y = element_text(size = 8),
                      axis.text.x = element_text(size = 9),
                      axis.ticks = element_blank(),
                      plot.margin = margin(0, 4, 4, 0, "mm"))

pC <- (pC_cls | pC_body) + plot_layout(widths = c(1, 60), guides = "collect") &
  theme(legend.position = "bottom", legend.box.margin = margin(0, 0, 0, 0, "mm"))
FigS3 <- (pVol + pSh) / pC + plot_layout(heights = c(72, 148))
save_fig(FigS3, "FigS3_deg_shap_heatmap", w = 185, h = 220)
cat("图S3 A/B/C 合并组图完成\n")
