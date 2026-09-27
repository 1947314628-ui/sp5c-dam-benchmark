#!/usr/bin/env Rscript
# 图5 三数据集 DAM 轨迹与 Apoe 表达（A/B，三面板，统一 a priori 13+8 基因协议）
# 数据来源（只读）：
#   ../data_r/gse285766_mg_cells.csv       （CFA 柱：条件 DAM-like % 与细胞级 Apoe 均值）
#   ../data_r/gse172167_dam_scores.csv     （SCI 折线：每条件状态计数）
#   ../data_r/gse172167_condition_summary.csv （SCI Apoe 均值）
#   ../data_r/gse162807_dam_per_sample.csv （SNI 折线：每样本 damlike%/Apoe，26 样本）
suppressMessages({library(ggplot2); library(dplyr); library(readr); library(tidyr); library(patchwork)})
source("theme_jtm.R"); set_figdir()
DATA <- normalizePath("../data_r")

mg <- read_csv(file.path(DATA, "gse285766_mg_cells.csv"), show_col_types = FALSE)
mg$condition <- factor(mg$condition, levels = names(col_condition))
s17 <- read_csv(file.path(DATA, "gse172167_dam_scores.csv"), show_col_types = FALSE)
s17sum <- read_csv(file.path(DATA, "gse172167_condition_summary.csv"), show_col_types = FALSE)
pb8 <- read_csv(file.path(DATA, "gse162807_dam_per_sample.csv"), show_col_types = FALSE)
pb8$tp <- factor(pb8$tp, levels = c("naive", "3d", "14d", "5mo"))
s8 <- read_csv(file.path(DATA, "gse162807_dam_scores.csv"), show_col_types = FALSE)  # 细胞级（与手稿正文口径一致）
tp17 <- c(A_Uninj = "Uninj", B_1dpi = "1dpi", C_1wpi = "1wpi", D_3wpi = "3wpi", E_6wpi = "6wpi")

# ---- A：DAM-like 比例 ----
g766 <- mg %>% count(condition, dam_state) %>% group_by(condition) %>%
  mutate(pct = n / sum(n) * 100) %>% filter(dam_state == "DAM-like")
g172 <- s17 %>% count(condition, state) %>% group_by(condition) %>%
  mutate(pct = n / sum(n) * 100) %>% filter(state == "DAM-like") %>%
  mutate(group = factor(tp17[condition], levels = unname(tp17)))
g807 <- s8 %>% filter(surgery != "Naive") %>%
  group_by(timepoint, surgery) %>% summarise(pct = mean(dam_state == "DAM-like") * 100, .groups = "drop") %>%
  rename(tp = timepoint)
cat("GSE285766 DAM-like%:"); print(round(g766$pct, 2))
cat("GSE172167 DAM-like%:"); print(round(g172$pct, 2))
cat("GSE162807 DAM-like%（细胞级，与手稿锚点 30.3/34.4/37.3 及 24.4/24.8/25.4 对账）:\n"); print(round(g807$pct, 2))

pA1 <- ggplot(g766, aes(condition, pct, fill = condition)) +
  geom_col(width = 0.55) + scale_fill_manual(values = col_condition, guide = "none") +
  geom_text(aes(label = sprintf("%.1f", pct)), vjust = -0.4, size = 3.2) +
  coord_cartesian(ylim = c(0, 45)) +
  labs(x = NULL, y = "DAM-like (%)", title = "A  GSE285766 (CFA)") + theme_jtm()
pA2 <- ggplot(g172, aes(group, pct, group = 1)) +
  geom_line(colour = "#0072B2", linewidth = 0.7) + geom_point(colour = "#0072B2", size = 2) +
  coord_cartesian(ylim = c(0, 45)) +
  labs(x = NULL, y = "DAM-like (%)", title = "B  GSE172167 (SCI)") + theme_jtm()
pA3 <- ggplot(g807, aes(tp, pct, group = surgery, colour = surgery)) +
  geom_line(linewidth = 0.7) + geom_point(aes(shape = surgery), size = 2) +
  scale_colour_manual(values = col_surgery, name = NULL) +
  scale_shape_manual(values = c(Sham = 16, SNI = 15), name = NULL) +
  coord_cartesian(ylim = c(0, 45)) +
  labs(x = NULL, y = "DAM-like (%)", title = "C  GSE162807 (SNI)") +
  theme_jtm() + theme(legend.position = c(0.12, 0.85))

# ---- B：Apoe 均值 ----
apo766 <- mg %>% group_by(condition) %>% summarise(apoe = mean(Apoe), .groups = "drop")
apo172 <- s17sum %>% mutate(group = factor(tp17[condition], levels = unname(tp17)),
                            apoe = apoe_mean) %>% select(group, apoe)
apo807 <- s8 %>% filter(surgery != "Naive") %>%
  group_by(timepoint, surgery) %>% summarise(apoe = mean(apoe), .groups = "drop") %>%
  rename(tp = timepoint)
cat("GSE285766 Apoe 细胞级均值（手稿锚点 0.40/1.36/0.91）:"); print(round(apo766$apoe, 3))
cat("GSE172167 Apoe:"); print(round(apo172$apoe, 3))
cat("GSE162807 Apoe 细胞级（手稿锚点 0.991/2.004/1.766 与 1.224/1.214/0.829）:"); print(round(apo807$apoe, 3))

pB1 <- ggplot(apo766, aes(condition, apoe, fill = condition)) +
  geom_col(width = 0.55) + scale_fill_manual(values = col_condition, guide = "none") +
  geom_text(aes(label = sprintf("%.2f", apoe)), vjust = -0.4, size = 3.2) +
  coord_cartesian(ylim = c(0, 1.6)) +
  labs(x = NULL, y = "Apoe mean (log-norm)", title = "D  GSE285766 (CFA)") + theme_jtm()
pB2 <- ggplot(apo172, aes(group, apoe, group = 1)) +
  geom_line(colour = "#0072B2", linewidth = 0.7) + geom_point(colour = "#0072B2", size = 2) +
  coord_cartesian(ylim = c(0, 5.2)) + labs(x = NULL, y = NULL, title = "E  GSE172167 (SCI)") + theme_jtm()
pB3 <- ggplot(apo807, aes(tp, apoe, group = surgery, colour = surgery)) +
  geom_line(linewidth = 0.7) + geom_point(aes(shape = surgery), size = 2) +
  scale_colour_manual(values = col_surgery, guide = "none") +
  scale_shape_manual(values = c(Sham = 16, SNI = 15), guide = "none") +
  coord_cartesian(ylim = c(0, 2.5)) + labs(x = NULL, y = NULL, title = "F  GSE162807 (SNI)") + theme_jtm()

Fig5 <- (pA1 + pA2 + pA3) / (pB1 + pB2 + pB3) + plot_layout(heights = c(1, 0.8))
save_fig(Fig5, "Fig5_three_datasets", w = 185, h = 130)
cat("图5 完成\n")
