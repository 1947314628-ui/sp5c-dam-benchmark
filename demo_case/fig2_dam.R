#!/usr/bin/env Rscript
# 图2 CFA 诱导 DAM-like 状态转变（A/B/C）
# 数据来源（只读）：
#   ../data_r/gse285766_mg_pseudobulk.csv  （每样本 DAM 均值/比例，云端复制）
#   ../data_r/gse285766_mg_cells.csv       （870 细胞级 dam_score/dam_state）
suppressMessages({library(ggplot2); library(dplyr); library(readr); library(tidyr)})
source("theme_jtm.R"); set_figdir()
DATA <- normalizePath("../data_r")

pb <- read_csv(file.path(DATA, "gse285766_mg_pseudobulk.csv"), show_col_types = FALSE)
pb$cond <- factor(pb$cond, levels = names(col_condition))
pb$sample <- factor(pb$sample, levels = c("CON1","CON2","CFA1","CFA2","TB1","TB2"))
mg <- read_csv(file.path(DATA, "gse285766_mg_cells.csv"), show_col_types = FALSE)
mg$condition <- factor(mg$condition, levels = names(col_condition))
mg$dam_state <- factor(mg$dam_state, levels = names(col_dam))

# 条件均值：细胞级（水平线）与样本级（折变计算）
cell_means <- mg %>% group_by(condition) %>% summarise(m = mean(dam_score), .groups = "drop")
samp_means <- pb %>% group_by(cond) %>% summarise(ms = mean(dam_mean), pp = mean(`damlike%`), .groups = "drop")
d_s <- round(samp_means$ms[2] - samp_means$ms[1], 3)   # CFA-CON
d_stb <- round(samp_means$ms[3] - samp_means$ms[2], 3) # TB-CFA
d_pp <- round(samp_means$pp[2] - samp_means$pp[1], 1)  # CFA-CON pp
cat("细胞级条件均值:"); print(round(cell_means$m, 3))
cat("样本级均值:"); print(samp_means)
cat("Δscore CFA-CON:", d_s, "| TB-CFA:", d_stb, "| Δpp CFA-CON:", d_pp, "\n")

# A 每样本 DAM 评分柱 + 细胞级条件均值水平线
pA <- ggplot(pb, aes(sample, dam_mean, fill = cond)) +
  geom_col(width = 0.55) +
  geom_hline(data = cell_means, aes(yintercept = m, colour = condition), linewidth = 0.6, linetype = 2) +
  scale_fill_manual(values = col_condition, guide = "none") +
  scale_colour_manual(values = col_condition, name = "Condition mean\n(cell-level)") +
  labs(x = NULL, y = "DAM score (mean)",
       title = "A  Per-sample DAM score\n(dashed = cell-level means)") +
  theme_jtm() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# B 状态构成
prop <- mg %>% count(condition, dam_state) %>% group_by(condition) %>% mutate(pct = n / sum(n) * 100)
pB <- ggplot(prop, aes(condition, pct, fill = dam_state)) +
  geom_col(width = 0.55) + scale_fill_manual(values = col_dam, name = NULL) +
  geom_text(data = prop %>% filter(dam_state == "DAM-like"),
            aes(label = sprintf("%.1f%%", pct)), vjust = -0.4, size = 3.2) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = NULL, y = "Proportion (%)", title = "B  DAM state composition") +
  theme_jtm() + theme(legend.position = "bottom")

# C 方向性折变（样本层，与 2026-08-31 权威值 +0.222/+11.1pp/−0.236 对账）
fc <- data.frame(term = c("CFA−CON ΔDAM score", "TB−CFA ΔDAM score", "CFA−CON ΔDAM-like (pp)"),
                 value = c(d_s, d_stb, d_pp),
                 sign = factor(c("up", "down", "up"), levels = c("up", "down")),
                 lab = c(sprintf("%+.3f", d_s), sprintf("%+.3f", d_stb), sprintf("%+.1f pp", d_pp)))
pC <- ggplot(fc, aes(term, value, fill = sign)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = lab, vjust = ifelse(value >= 0, -0.4, 1.4)), size = 3.2) +
  scale_fill_manual(values = c(up = "#D73027", down = "#4575B4"), guide = "none") +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_y_continuous(expand = expansion(mult = c(0.3, 0.30))) +
  labs(x = NULL, y = "Fold change (sample-level)", title = "C  Directional changes") +
  theme_jtm() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

Fig2 <- (pA + pB + pC) + plot_layout(widths = c(1.15, 1, 1.25))
save_fig(Fig2, "Fig2_DAM_state", w = 185, h = 80)
cat("图2 完成\n")
