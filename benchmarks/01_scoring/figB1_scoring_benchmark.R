#!/usr/bin/env Rscript
# 图B1（BIB Method Review 改造新增）DAM-like 比例定义方法基准：
#   M-A 云端 state 标签（现行口径）
#   M-B 每数据集 dam_score top-25% 分位
#   M-C 每数据集 dam_score 2-均值 k-means 聚类（高中心簇 = DAM-like，stats::kmeans 免依赖）
# 数据来源（只读）：
#   ../../data_r/gse285766_mg_cells.csv   （870 细胞 × 21 基因 + dam_score/dam_state）
#   ../../data_r/gse162807_dam_scores.csv （26 样本细胞级 dam_score/dam_state）
#   ../../data_r/gse172167_dam_scores.csv （细胞级 dam_score/state）
# B 面板（仅 GSE285766，有基因矩阵）：打分函数比较
#   21 基因净评分（up13−ho8 均值）vs Keren-Shaul 2017 DAM 上调子集（矩阵可得的 11 基因：Apoe,Axl,Cst7,Itgax,Lpl,Tyrobp,Trem2,Cd9,Spp1,Ctsb,B2m）
#   vs 双阳性标记法（Itgax>0 且 Apoe>0）；前两者取 top-25% 分位定义 DAM-like，后者为阳性率
suppressMessages({library(ggplot2); library(dplyr); library(readr); library(tidyr)})
source("theme_jtm.R"); set_figdir()
DATA <- normalizePath("../../data_r")

mg  <- read_csv(file.path(DATA, "gse285766_mg_cells.csv"), show_col_types = FALSE)
s8  <- read_csv(file.path(DATA, "gse162807_dam_scores.csv"), show_col_types = FALSE)
s17 <- read_csv(file.path(DATA, "gse172167_dam_scores.csv"), show_col_types = FALSE)
cat("细胞数: 285766 =", nrow(mg), "| 162807 =", nrow(s8), "| 172167 =", nrow(s17), "\n")
cat("172167 state 取值:\n"); print(table(s17$state, useNA = "ifany"))

# ---- 三方法比例（每样本）----
three_props <- function(d, score, label, group) {
  q75 <- quantile(d[[score]], 0.75)
  km  <- kmeans(d[[score]], centers = 2, nstart = 10)
  high_cl <- which.max(km$centers)
  data.frame(score = d[[score]], label = d[[label]], group = d[[group]]) %>%
    mutate(mA = label == "DAM-like", mB = score >= q75, mC = km$cluster == high_cl) %>%
    group_by(group) %>%
    summarise(A = mean(mA) * 100, B = mean(mB) * 100, C = mean(mC) * 100,
              n = n(), .groups = "drop")
}

p_gs  <- three_props(mg,  "dam_score", "dam_state", "sample")
p_s8  <- three_props(s8,  "dam_score", "dam_state", "sample")
p_s17 <- three_props(s17, "dam_score", "state",     "condition")

# ---- 锚点核对（M-A 必须复现权威值）----
condA <- mg %>% mutate(condition = factor(condition, levels = names(col_condition))) %>%
  group_by(condition) %>%
  summarise(A = mean(dam_state == "DAM-like") * 100, .groups = "drop")
cat("285766 M-A 条件比例（对账 28.51/37.86/26.11）:\n"); print(condA)
g8 <- s8 %>% group_by(timepoint, surgery) %>%
  summarise(A = mean(dam_state == "DAM-like") * 100, .groups = "drop")
cat("162807 M-A 5mo SNI（对账 37.3）:",
    round(g8$A[g8$timepoint == "5mo" & g8$surgery == "SNI"], 2), "\n")
cat("172167 M-A 条件比例:\n"); print(p_s17 %>% select(group, A, n) %>% arrange(group))
cat("方法间均值（样本级）:\n")
cat("285766: A/B/C =", round(colMeans(p_gs[, c("A", "B", "C")]), 2), "\n")
cat("162807: A/B/C =", round(colMeans(p_s8[, c("A", "B", "C")]), 2), "\n")
cat("172167: A/B/C =", round(colMeans(p_s17[, c("A", "B", "C")]), 2), "\n")

# ---- A：三数据集 × 三方法 箱线 + 样本点 ----
longA <- bind_rows(
  p_gs  %>% select(-n) %>% mutate(dataset = "GSE285766 (CFA, trigeminal)"),
  p_s8  %>% select(-n) %>% mutate(dataset = "GSE162807 (SNI, spinal)"),
  p_s17 %>% select(-n) %>% mutate(dataset = "GSE172167")
) %>% pivot_longer(c(A, B, C), names_to = "method", values_to = "pct") %>%
  mutate(method = factor(method, levels = c("A", "B", "C"),
                         labels = c("Label", "Top 25%", "k-means 2C")),
         dataset = factor(dataset, levels = unique(dataset)))
pA <- ggplot(longA, aes(method, pct)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 0.8, alpha = 0.6) +
  facet_wrap(~dataset, nrow = 1) +
  labs(x = NULL, y = "DAM-like (%) per sample",
       title = "A  Proportion-definition methods across datasets") +
  theme_jtm()

# ---- B：GSE285766 打分函数比较 ----
ks11 <- c("Apoe", "Axl", "Cst7", "Itgax", "Lpl", "Tyrobp", "Trem2", "Cd9", "Spp1", "Ctsb", "B2m")
up13 <- c("Apoe","Lpl","Itgax","Axl","Trem2","Cst7","Ccl6","Cd9","Spp1","Csf1","Tyrobp","B2m","Ctsb")
ho8  <- c("Tmem119","P2ry12","Cx3cr1","Hexb","C1qa","Olfml3","Siglech","Selplg")
mg2 <- mg %>% mutate(
  s_net = rowMeans(.[up13]) - rowMeans(.[ho8]),
  s_ks  = rowMeans(.[ks11]),
  s_dp  = as.numeric(Itgax > 0 & Apoe > 0),
  dam_state = factor(dam_state, levels = names(col_dam)))
up <- rowMeans(mg[up13]); ho <- rowMeans(mg[ho8])
cat("sd(up13)=", round(sd(up), 3), "| sd(ho8)=", round(sd(ho), 3),
    "| cor(up,ho)=", round(cor(up, ho), 3),
    "| cor(dam_score,s_net)=", round(cor(mg$dam_score, mg2$s_net), 3), "\n")
cat("s_net vs s_ks Spearman:", round(cor(mg2$s_net, mg2$s_ks, method = "spearman"), 3),
    "| Itgax+Apoe 双阳性率(%):", round(mean(mg2$s_dp) * 100, 2), "\n")
st <- mg2 %>% group_by(dam_state) %>%
  summarise(n = n(), net = mean(s_net), ks = mean(s_ks), .groups = "drop")
cat("285766 各状态评分均值:\n"); print(st)

rho <- round(cor(mg2$s_net, mg2$s_ks, method = "spearman"), 2)
pB1 <- ggplot(mg2, aes(s_net, s_ks, colour = dam_state)) +
  geom_point(size = 0.8, alpha = 0.55) +
  scale_colour_manual(values = col_dam, name = NULL) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.6, size = 3.2,
           label = sprintf("Spearman rho = %.2f", rho)) +
  labs(x = "21-gene net score", y = "Keren-Shaul subset score",
       title = "B  Scoring functions rank cells differently") +
  theme_jtm() + theme(legend.position = "bottom")

stlong <- st %>% pivot_longer(c(net, ks), names_to = "score", values_to = "m") %>%
  mutate(score = factor(score, levels = c("net", "ks"),
                        labels = c("21-gene net", "Keren-Shaul subset")))
pB2 <- ggplot(stlong, aes(dam_state, m, fill = dam_state)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = sprintf("%.2f", m), vjust = ifelse(m >= 0, -0.4, 1.4)), size = 3.2) +
  scale_fill_manual(values = col_dam, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.3, 0.30))) +
  facet_wrap(~score, nrow = 1) +
  labs(x = NULL, y = "Mean score per state",
       title = "C  State separation by score") +
  theme_jtm()

FigB1 <- pA / (pB1 | pB2) + plot_layout(heights = c(1.15, 1))
save_fig(FigB1, "FigB1_scoring_benchmark", w = 185, h = 165)

# ---- 输出长表 CSV（正文+代码仓用）----
out <- longA %>% rename(sample = group)
write_csv(out, file.path(DATA, "scoring_method_concordance.csv"))
cat("图B1 + scoring_method_concordance.csv 完成\n")
