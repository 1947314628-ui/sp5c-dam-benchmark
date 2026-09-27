# 统一绘图主题与调色板 —— 自云端 zip（投稿包/拟投稿_SP5C/R/theme_jtm.R）适配本机路径
# 所有图脚本 source 本文件；figdir 由各脚本设为脚本所在目录
suppressMessages({library(ggplot2); library(patchwork)})

# ---- 调色板 ----
col_condition <- c(CON = "#0072B2", CFA = "#D55E00", TB = "#009E73")       # 条件（GSE285766）
col_celltype  <- c(Oligodendrocytes = "#95A5A6", Astrocytes = "#E74C3C",
                   OPCs = "#F39C12", Microglia = "#8E44AD",
                   Neurons = "#3498DB", Pericytes = "#16A085")            # 六型
col_dam       <- c(Homeostatic = "#4575B4", Intermediate = "#FEE090",
                   `DAM-like` = "#D73027")                                # DAM 状态
col_surgery   <- c(Naive = "#91BFDB", Sham = "#4575B4", SNI = "#D73027")  # GSE162807
col_sample    <- c(CON1 = "#99C9E8", CON2 = "#0072B2", CFA1 = "#F2A77E", CFA2 = "#D55E00",
                   TB1 = "#7FCDBB", TB2 = "#009E73")                      # 六样本
col_dataset   <- c(`GSE285766 (CFA)` = "#D55E00", `GSE172167 (SCI)` = "#0072B2",
                   `GSE162807 (SNI)` = "#009E73")                         # 三数据集

# ---- 统一主题 ----
theme_jtm <- function(base = 9) {
  theme_classic(base_size = base) +
    theme(
      text = element_text(family = "sans", colour = "black"),
      axis.text = element_text(colour = "black", size = base),
      axis.title = element_text(size = base),
      axis.line = element_line(colour = "black", linewidth = 0.4),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = unit(1.2, "mm"),
      legend.text = element_text(size = base), legend.title = element_text(size = base),
      legend.key.size = unit(3.5, "mm"),
      strip.background = element_blank(), strip.text = element_text(size = base),
      panel.spacing = unit(4, "mm"),
      plot.margin = margin(4, 4, 4, 4, "mm"),
      plot.title = element_text(size = base, face = "bold", hjust = 0)
    )
}

# ---- 统一导出：PDF（矢量）+ PNG（300dpi）；figdir 由调用脚本定义 ----
save_fig <- function(p, name, w = 90, h = 80, units = "mm") {
  if (capabilities("cairo")) {
    ggsave(file.path(figdir, paste0(name, ".pdf")), p, width = w, height = h, units = units, device = cairo_pdf)
  } else {
    ggsave(file.path(figdir, paste0(name, ".pdf")), p, width = w, height = h, units = units)
  }
  ggsave(file.path(figdir, paste0(name, ".png")), p, width = w, height = h, units = units, dpi = 300)
  cat("saved", name, "\n")
}

# 各脚本共用：定位脚本所在目录为 figdir（脚本文件与输出同目录）
set_figdir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  figdir <<- normalizePath(dirname(f))
}
