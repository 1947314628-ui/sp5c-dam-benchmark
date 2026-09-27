#!/usr/bin/env Rscript
# 组装 GSE285766 CFA 两样本 GEX 矩阵 → Seurat 对象（仅 CFA 子集，通讯比较用）
# 关键修正：features 是 6 列 TSV（含 Peaks 行），只保留 Gene Expression 行
# 输入：../../downloads/gex/（CFA1/CFA2 三件套）
#       ../../投稿包/拟投稿_SP5C/data_r/gse285766_all_cells.csv（细胞 ID 权威源）
# 输出：sp5c_cfa_seurat.rds + export/（counts.mtx 等，liana-py 用）
suppressMessages({library(Matrix); library(readr); library(dplyr); library(Seurat)})
GEX  <- normalizePath(Sys.getenv("SP5C_GEX", "../../downloads/gex"))
META <- read_csv("../../投稿包/拟投稿_SP5C/data_r/gse285766_all_cells.csv", show_col_types = FALSE)
samples <- c("CFA1", "CFA2")
dir.create("mats", showWarnings = FALSE)

for (s in samples) {
  feat <- read_tsv(file.path(GEX, paste0(s, "-features_gex.csv.gz")), col_names = FALSE, show_col_types = FALSE)
  keep_gex <- feat[[3]] == "Gene Expression"
  genes   <- make.unique(feat[[2]][keep_gex])
  mtx <- readMM(gzfile(file.path(GEX, paste0(s, "-matrix_gex.mtx.gz"))))
  mtx <- mtx[keep_gex, , drop = FALSE]                     # 先去峰行（行序与 features 一致）
  bc  <- read_csv(file.path(GEX, paste0(s, "-barcodes_gex.csv.gz")), col_names = FALSE, show_col_types = FALSE)[[1]]
  bc  <- bc[!grepl("^barcode", bc, ignore.case = TRUE)]
  suffix <- if (any(grepl("-1$", head(bc, 5)))) "" else "-1"
  cells <- paste0(s, "_", bc, suffix)
  keep <- cells %in% META$cell
  mtx <- mtx[, keep]
  rownames(mtx) <- genes; colnames(mtx) <- cells[keep]
  mtx <- mtx[Matrix::rowSums(mtx) > 0, , drop = FALSE]
  cat(s, ": 特征", nrow(feat), "→基因行", sum(keep_gex), "→保留", nrow(mtx),
      "| 细胞", length(cells), "→命中", sum(keep), "\n")
  saveRDS(mtx, file.path("mats", paste0(s, ".rds")))
  rm(mtx, feat, bc); gc()
}

m <- readRDS(file.path("mats", paste0(samples[1], ".rds")))
for (s in samples[-1]) {
  nxt <- readRDS(file.path("mats", paste0(s, ".rds")))
  add <- setdiff(rownames(nxt), rownames(m))
  if (length(add)) m <- rbind(m, Matrix(0, length(add), ncol(m), dimnames = list(add, colnames(m)), sparse = TRUE))
  add2 <- setdiff(rownames(m), rownames(nxt))
  if (length(add2)) nxt <- rbind(nxt, Matrix(0, length(add2), ncol(nxt), dimnames = list(add2, colnames(nxt)), sparse = TRUE))
  m <- cbind(m, nxt[rownames(m), ])
  rm(nxt); gc()
}
meta <- META[match(colnames(m), META$cell), ] %>% as.data.frame()
rownames(meta) <- meta$cell
seu <- CreateSeuratObject(counts = m, meta.data = meta)
seu <- NormalizeData(seu)
saveRDS(seu, "sp5c_cfa_seurat.rds")
cat("Seurat CFA:", ncol(seu), "细胞 ×", nrow(seu), "基因\n")
print(table(seu$celltype))

dir.create("export", showWarnings = FALSE)
Matrix::writeMM(seu[["RNA"]]$counts, "export/counts.mtx")
writeLines(colnames(seu), "export/barcodes.txt")
writeLines(rownames(seu), "export/genes.csv")
write_csv(data.frame(cell = colnames(seu),
                     condition = seu$condition, celltype = seu$celltype,
                     leiden = seu$leiden, sample = seu$sample),
          "export/meta.csv")
cat("export/ 导出完成\n")
