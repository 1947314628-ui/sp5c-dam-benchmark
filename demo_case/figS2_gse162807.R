#!/usr/bin/env Rscript
# 图S2 GSE162807 每样本 DAM-like 比例（26 样本，a priori 协议；替代丢失的旧 Figure_S7，统一口径）
# 数据来源（只读）：
#   ../data_r/gse162807_dam_per_sample.csv
suppressMessages({library(ggplot2); library(dplyr); library(readr)})
source("theme_jtm.R"); set_figdir()
DATA <- normalizePath("../data_r")

pb <- read_csv(file.path(DATA, "gse162807_dam_per_sample.csv"), show_col_types = FALSE)
pb$tp <- factor(pb$tp, levels = c("naive", "3d", "14d", "5mo"))
pb$surgery <- factor(pb$surgery, levels = names(col_surgery))
s8 <- read_csv(file.path(DATA, "gse162807_dam_scores.csv"), show_col_types = FALSE)  # 细胞级
cat("样本数:", nrow(pb), "\n"); print(table(pb$tp, pb$surgery))
g <- s8 %>% group_by(timepoint, surgery) %>%
  summarise(m = mean(dam_state == "DAM-like") * 100, .groups = "drop") %>%
  rename(tp = timepoint)
cat("5mo SNI 细胞级比例（对账 37.3 口径）:", round(g$m[g$tp == "5mo" & g$surgery == "SNI"], 2), "\n")
print(round(g$m, 2))

p <- ggplot(pb, aes(tp, `damlike%`, colour = surgery)) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.75) +
  geom_line(data = g, aes(tp, m, group = surgery, colour = surgery), linewidth = 0.7) +
  geom_point(data = g, aes(tp, m, colour = surgery), shape = 17, size = 2.6) +
  scale_colour_manual(values = col_surgery, name = NULL) +
  annotate("text", x = Inf, y = 52, hjust = 1.05, size = 3.2,
           label = sprintf("SNI 5mo cell-level %.1f%%", g$m[g$tp == "5mo" & g$surgery == "SNI"])) +
  coord_cartesian(ylim = c(0, 55), clip = "off") +
  labs(x = NULL, y = "DAM-like (%)", title = NULL) +
  theme_jtm() + theme(legend.position = c(0.15, 0.12))
save_fig(p, "FigS2_gse162807_per_sample", w = 150, h = 90)
cat("图S2 完成\n")
