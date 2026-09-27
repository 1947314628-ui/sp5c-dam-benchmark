#!/usr/bin/env Rscript
# FigB2：四种通讯打分口径在 CFA 子集六型上的一致性比较
#   输入：benchmark_R/cellchat_cfa_pairs.csv（CellChat 规范版）
#         benchmark_R/liana_py_{cellphonedb|cellchat|consensus}_cfa.csv
#   输出：FigB2_comm_tools.png/pdf + comm_tool_shared_pairs.csv + top30_overlap.csv + mg_top_pairs.csv
suppressMessages({library(ggplot2); library(dplyr); library(readr); library(tidyr); library(patchwork)})
setwd(".")   # CSV 在 benchmark_R 下；图输出回 figures 目录
OUT <- "."

cc  <- read_csv("cellchat_cfa_pairs.csv", show_col_types = FALSE) %>%
  transmute(key = paste(source, target, ligand, receptor, sep = "|"),
            cellchat_prob = prob, cellchat_pval = pval,
            source, target, ligand, receptor)
cp  <- read_csv("liana_py_cellphonedb_cfa.csv", show_col_types = FALSE) %>%
  transmute(key = paste(source, target, ligand_complex, receptor_complex, sep = "|"),
            cpdb_means = lr_means, cpdb_pval = cellphone_pvals,
            source, target, ligand = ligand_complex, receptor = receptor_complex)
lc  <- read_csv("liana_py_cellchat_cfa.csv", show_col_types = FALSE) %>%
  transmute(key = paste(source, target, ligand_complex, receptor_complex, sep = "|"),
            lianaCC_probs = lr_probs, lianaCC_pval = cellchat_pvals,
            source, target, ligand = ligand_complex, receptor = receptor_complex)
cs  <- read_csv("liana_py_consensus_cfa.csv", show_col_types = FALSE) %>%
  transmute(key = paste(source, target, ligand_complex, receptor_complex, sep = "|"),
            cs_mag_rank = magnitude_rank, cs_spec_rank = specificity_rank,
            source, target, ligand = ligand_complex, receptor = receptor_complex)

cat("对空间规模 | CellChat:", nrow(cc), "| liana-CPDB:", nrow(cp),
    "| liana-CellChat:", nrow(lc), "| consensus:", nrow(cs), "\n")
shared_cc_lc <- inner_join(cc, lc, by = "key", suffix = c(".cc", ".lc"))
cat("CellChat ∩ liana-CellChat 共同键:", nrow(shared_cc_lc), "\n")
rho1 <- cor.test(shared_cc_lc$cellchat_prob, shared_cc_lc$lianaCC_probs, method = "spearman")
cat("ρ(CellChat prob, lianaCC probs) =", round(rho1$estimate, 3), "p =", format.pval(rho1$p.value), "\n")
rho2 <- cor.test(cp$cpdb_means, lc$lianaCC_probs[match(cp$key, lc$key)], method = "spearman")
cat("ρ(cpdb means, lianaCC probs) =", round(rho2$estimate, 3), "p =", format.pval(rho2$p.value), "\n")

# ---- top-30 Jaccard ----
top30 <- list(CellChat       = head(arrange(cc, desc(cellchat_prob)), 30)$key,
              `LIANA-CellChat`    = head(arrange(lc, desc(lianaCC_probs)), 30)$key,
              `LIANA-CellPhoneDB` = head(arrange(cp, desc(cpdb_means)), 30)$key,
              `LIANA-consensus`   = head(arrange(cs, cs_mag_rank), 30)$key)
jacc <- function(a, b) length(intersect(a, b)) / length(union(a, b))
J <- outer(seq_along(top30), seq_along(top30),
           Vectorize(function(i, j) jacc(top30[[i]], top30[[j]])))
dimnames(J) <- list(names(top30), names(top30))
cat("top-30 Jaccard:\n"); print(round(J, 2))

# ---- 小胶质相关对 ----
top30_all <- bind_rows(
  head(arrange(cc, desc(cellchat_prob)), 30) %>% mutate(tool = "CellChat"),
  head(arrange(lc, desc(lianaCC_probs)), 30) %>% mutate(tool = "LIANA-CellChat"),
  head(arrange(cp, desc(cpdb_means)), 30)    %>% mutate(tool = "LIANA-CellPhoneDB"),
  head(arrange(cs, cs_mag_rank), 30)         %>% mutate(tool = "LIANA-consensus"))
mg_top <- top30_all %>%
  mutate(mg_role = case_when(source == "Microglia" & target == "Microglia" ~ "MG→MG",
                             source == "Microglia" ~ "MG as source",
                             target == "Microglia" ~ "MG as target",
                             TRUE ~ NA_character_)) %>%
  filter(!is.na(mg_role)) %>%
  mutate(pair = paste(ligand, receptor, sep = "→"))
cat("top-30 中 MG 相关对:\n"); print(mg_top %>% select(tool, pair, source, target))
write_csv(mg_top %>% select(tool, source, target, pair, mg_role),
          file.path(OUT, "mg_top_pairs.csv"))

# ---- 图 ----
thm <- theme_bw(base_size = 9) + theme(legend.position = "bottom",
  legend.key.size = unit(2.5, "mm"), plot.title = element_text(size = 9, face = "bold"))

# A：共享键上的散点（CellChat×lianaCC、CellChat×cpdb、lianaCC×cpdb）
shared_cc_cp <- inner_join(cc, cp, by = "key")
sc <- bind_rows(
  inner_join(cc, lc, by = "key") %>% transmute(x = cellchat_prob, y = lianaCC_probs,
        cmp = paste0("CellChat vs LIANA-CellChat   (rho=",
                     round(rho1$estimate, 2), ")")),
  shared_cc_cp %>% transmute(x = cellchat_prob, y = cpdb_means,
        cmp = paste0("CellChat vs LIANA-CellPhoneDB   (rho=",
                     round(cor.test(shared_cc_cp$cellchat_prob, shared_cc_cp$cpdb_means, method = "spearman")$estimate, 2), ")")),
  transmute(cp, x = lc$lianaCC_probs[match(cp$key, lc$key)], y = cpdb_means,
        cmp = paste0("LIANA-CellChat vs LIANA-CellPhoneDB   (rho=",
                     round(rho2$estimate, 2), ")"))) %>%
  mutate(cmp = factor(cmp, levels = unique(cmp)))
pA <- ggplot(sc, aes(x, y)) + geom_point(size = 0.4, alpha = 0.35) +
  facet_wrap(~cmp, nrow = 1) + labs(x = "score (x)", y = "score (y)", title = "A  Shared-pair score agreement") +
  thm + theme(strip.text = element_text(size = 9))

# B：top-30 Jaccard 热图
Jdf <- as.data.frame(as.table(J)) %>% setNames(c("t1", "t2", "jacc")) %>%
  mutate(across(c(t1, t2), ~factor(., levels = names(top30))))
pB <- ggplot(Jdf, aes(t1, t2, fill = jacc)) + geom_tile(color = "white") +
  geom_text(aes(label = round(jacc, 2)), size = 3.2) +
  scale_fill_gradient(low = "white", high = "#2171b5", limits = c(0, 1)) +
  labs(x = NULL, y = NULL, title = "B  Top-30 pair overlap (Jaccard)") +
  thm + theme(axis.text.x = element_text(angle = 30, hjust = 1))

# C：各工具 top-30 中 MG 相关对方向计数
cnt <- mg_top %>% count(tool, mg_role) %>%
  mutate(tool = factor(tool, levels = names(top30)))
pC <- ggplot(cnt, aes(tool, n, fill = mg_role)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = n), position = position_dodge(0.8), vjust = -0.4, size = 3.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = NULL, y = "MG-involving pairs in top-30", title = "C  Microglia-related top pairs") +
  thm + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  guides(fill = guide_legend(nrow = 1))

fig <- pA / (pB | pC)
ggsave(file.path(OUT, "FigB2_comm_tools.png"), fig, width = 185, height = 165, units = "mm", dpi = 300)
ggsave(file.path(OUT, "FigB2_comm_tools.pdf"), fig, width = 185, height = 165, units = "mm")
# 共享对明细（复查用）
write_csv(full_join(shared_cc_lc %>% select(key, cellchat_prob, lianaCC_probs),
                    cp %>% select(key, cpdb_means), by = "key"),
          file.path(OUT, "comm_tool_shared_pairs.csv"))
cat("FigB2 完成\n")
