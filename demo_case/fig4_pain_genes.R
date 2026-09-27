#!/usr/bin/env Rscript
# 图4 神经元疼痛相关基因表达（A 阳性比例 / B 跨条件表达）
# 数据来源（只读）：
#   ../data_r/gse285766_scn9a_by_condition.csv  （云端 export_for_r.py §3）
#   ../data_r/gse285766_neurons_scn9a.csv       （同 §3，细胞级）
#   cloud_report_values.R（同目录）→ Slc17a6 41.4% / 0.518 / 0.556 / 0.686
#       来源：../../投稿包/拟投稿_SP5C/reports/liana_reconciliation.md 第51行
suppressMessages({library(ggplot2); library(dplyr); library(readr)})
source("theme_jtm.R"); set_figdir(); source("cloud_report_values.R")
DATA <- normalizePath("../data_r")

scn_sum <- read_csv(file.path(DATA, "gse285766_scn9a_by_condition.csv"), show_col_types = FALSE)
scn_sum$condition <- factor(scn_sum$condition, levels = names(col_condition))
scn <- read_csv(file.path(DATA, "gse285766_neurons_scn9a.csv"), show_col_types = FALSE)
scn$condition <- factor(scn$condition, levels = names(col_condition))

# 总体 Scn9a+ 比例 = 条件比例按细胞数加权（计算，非硬编码）
scn_overall <- sum(scn_sum$scn9a_pct * scn_sum$n) / sum(scn_sum$n)
cat("Scn9a 条件比例/均值:\n"); print(scn_sum)
cat(sprintf("Scn9a 总体阳性（加权）: %.1f%%；手稿锚点 43.9%%\n", scn_overall))
cat("Slc17a6（报告值）: 41.4%；", paste(sprintf("%s=%.3f", names(slc17a6_means), slc17a6_means), collapse=" "), "\n")

# A 阳性比例
pA_scn <- scn_sum %>% mutate(row = "Scn9a (Nav1.7)")
pA1 <- ggplot(scn_sum, aes(condition, scn9a_pct, fill = condition)) +
  geom_col(width = 0.55) + scale_fill_manual(values = col_condition, guide = "none") +
  geom_text(aes(label = sprintf("%.1f%%", scn9a_pct)), vjust = -0.4, size = 3.2) +
  geom_hline(yintercept = scn_overall, linetype = 2, linewidth = 0.5) +
  annotate("text", x = 1, y = scn_overall + 2, size = 3.2, hjust = 0,
           label = sprintf("overall %.1f%%", scn_overall)) +
  coord_cartesian(ylim = c(0, 60)) +
  labs(x = NULL, y = "Neuron positive (%)",
       title = "A  Scn9a positive fraction") + theme_jtm()
pA2 <- ggplot(data.frame(g = "Slc17a6\n(VGLUT2)", v = slc17a6_pct_overall), aes(g, v)) +
  geom_col(fill = "grey50", width = 0.4) +
  geom_text(aes(label = sprintf("%.1f%%", v)), vjust = -0.4, size = 3.2) +
  coord_cartesian(ylim = c(0, 60)) +
  labs(x = NULL, y = NULL, title = "B  Slc17a6 positive fraction") + theme_jtm()

# B 跨条件均值：Scn9a 小提琴 + Slc17a6 报告值柱
pB1 <- ggplot(scn, aes(condition, Scn9a, fill = condition)) +
  geom_violin(scale = "width", linewidth = 0.2) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2.5, colour = "black") +
  scale_fill_manual(values = col_condition, guide = "none") +
  labs(x = NULL, y = "Scn9a expression (log-norm)",
       title = "C  Scn9a expression across conditions") + theme_jtm()
slc <- data.frame(condition = factor(names(slc17a6_means), levels = names(col_condition)),
                  v = as.numeric(slc17a6_means))
pB2 <- ggplot(slc, aes(condition, v, fill = condition)) +
  geom_col(width = 0.55) + scale_fill_manual(values = col_condition, guide = "none") +
  geom_text(aes(label = sprintf("%.3f", v)), vjust = -0.4, size = 3.2) +
  coord_cartesian(ylim = c(0, 0.8)) +
  labs(x = NULL, y = "Slc17a6 expression (log-norm)",
       title = "D  Slc17a6 expression across conditions") + theme_jtm()

Fig4 <- (pA1 + pA2) / (pB1 + pB2) + plot_layout(heights = c(1, 1.15))
save_fig(Fig4, "Fig4_neuron_pain_genes", w = 185, h = 130)
cat("图4 完成\n")
