# ============================================
# 04_prepare_IPA.R
# 整理 DESeq2 输出数据供 IPA 使用
# ============================================

library(data.table)

# 输入输出目录
deg_dir <- "results/deseq2/"
ipa_dir <- "results/ipa_input/"
dir.create(ipa_dir, recursive = TRUE, showWarnings = FALSE)

# 加载所有 DEG 文件
deg_files <- list.files(deg_dir, pattern = "^DEG_.*\\.csv$", full.names = TRUE)

# IPA 需要的列：GeneName | log2FC | p-value | FDR
for (file in deg_files) {
  df <- fread(file)

  # 确保列名一致
  required_cols <- c("Gene", "log2FoldChange", "pvalue", "padj")
  if (!all(required_cols %in% colnames(df))) {
    stop("错误：找不到必要列，请确认输入文件格式")
  }

  # IPA 输入要求简洁清晰
  ipa_df <- df[, .(Gene, log2FoldChange, pvalue, padj)]
  setnames(ipa_df,
           old = c("Gene", "log2FoldChange", "pvalue", "padj"),
           new = c("Gene", "log2FC", "p-value", "FDR"))

  # 导出为 .txt 文件（IPA 建议）
  outname <- paste0("IPA_", gsub("DEG_", "", basename(file)))
  outname <- sub("\\.csv$", ".txt", outname)
  fwrite(ipa_df, file = file.path(ipa_dir, outname), sep = "\t", quote = FALSE)
}

message("✅ IPA 输入文件已生成：", ipa_dir)