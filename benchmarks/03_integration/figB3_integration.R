#!/usr/bin/env Rscript
# FigB3：四种整合方式（harmony/baseline/bbknn/scvi）下游稳定性比较
#   输入：benchmark_cloud/回传_2026-09-02/analysis/integration_compare/ 长表+指标表
#         + 本地权威注释 投稿包/拟投稿_SP5C/data_r/（celltype / MG dam_state）
#   输出：FigB3_integration.png/pdf + integration_long_annotated.csv.gz（图数据溯源）
suppressMessages({library(ggplot2); library(dplyr); library(readr); library(patchwork)})
DIR <- "."
META <- "../../data_r/gse285766_all_cells.csv"
MG   <- "../../data_r/gse285766_mg_cells.csv"

lt <- read_csv(file.path(DIR, "cells_umap_cluster.csv.gz"), show_col_types = FALSE)
ms <- read_csv(file.path(DIR, "methods_summary.csv"), show_col_types = FALSE)
meta <- read_csv(META, show_col_types = FALSE) %>% select(cell, celltype)
mg   <- read_csv(MG, show_col_types = FALSE) %>% select(cell, dam_state)

lt <- lt %>%
  left_join(meta, by = c(barcode = "cell")) %>%
  left_join(mg, by = c(barcode = "cell")) %>%
  mutate(method = factor(method, levels = c("harmony", "baseline", "bbknn", "scvi")))
stopifnot(sum(!is.na(lt$celltype)) == 145160, sum(!is.na(lt$dam_state)) == 3480)
write_csv(lt, "integration_long_annotated.csv.gz")   # 图数据溯源

METHOD_COLS <- c(harmony = "#377EB8", baseline = "#999999", bbknn = "#4DAF4A", scvi = "#E41A1C")
STATE_COLS  <- c(`DAM-like` = "#E41A1C", Intermediate = "#FF7F00", Homeostatic = "#377EB8")

thm <- theme_bw(base_size = 9) +
  theme(legend.position = "bottom", legend.key.size = unit(2.5, "mm"),
        plot.title = element_text(size = 9, face = "bold"),
        axis.text = element_text(size = 7), axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 7), legend.title = element_text(size = 8))

# ---- A：四方法 UMAP，灰色全部细胞 + MG 按 DAM 状态着色 ----
pA <- ggplot(lt, aes(UMAP1, UMAP2)) +
  geom_point(data = filter(lt, is.na(dam_state)), color = "#E8E8E8", size = 0.05) +
  geom_point(data = filter(lt, !is.na(dam_state)), aes(color = dam_state), size = 0.3) +
  facet_wrap(~method, nrow = 2) +
  scale_color_manual(values = STATE_COLS) +
  labs(title = "A  MG DAM-like states in each integration space (n = 870 per method)",
       color = "MG DAM state") +
  thm + theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 2.5)))

# ---- B：三指标（aRI vs 样本 / 加权簇熵 / 运行时间）----
ms <- mutate(ms, method = factor(method, levels = c("harmony", "baseline", "bbknn", "scvi")))
mkbar <- function(df, yvar, ylab, lab_fmt, title, note = NULL) {
  p <- ggplot(df, aes(method, .data[[yvar]], fill = method)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf(lab_fmt, .data[[yvar]])), vjust = -0.4, size = 3.2) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
    scale_fill_manual(values = METHOD_COLS, guide = "none") +
    labs(x = NULL, y = ylab, title = title) + thm +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  if (!is.null(note)) p <- p + annotate("text", x = 2.5, y = Inf,
        label = note, size = 2.5, vjust = 2, color = "grey35")
  p
}
pB1 <- mkbar(ms, "ari_vs_sample", "aRI (cluster vs sample)", "%.4f",
             "B  Batch effect (aRI vs sample; lower = better mixed)",
             "all methods near 0: batch effect negligible even uncorrected")
pB2 <- mkbar(ms, "wmean_entropy_sample", "Weighted cluster entropy", "%.3f",
             "C  Sample mixing entropy (higher = better)", "y-axis zoomed") +
  coord_cartesian(ylim = c(0.9, 0.93))
pB3 <- mkbar(ms, "runtime_s", "Runtime (s)", "%.0f",
             "D  Runtime", "harmony reused 08-31 embedding (excluded)")

fig <- pA / (pB1 | pB2 | pB3)
ggsave("FigB3_integration.png", fig, width = 185, height = 165, units = "mm", dpi = 300)
ggsave("FigB3_integration.pdf", fig, width = 185, height = 165, units = "mm")
cat("FigB3 done | clusters:", paste(paste0(ms$method, "=", ms$n_clusters), collapse = " "), "\n")
