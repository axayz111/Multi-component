# ============================
# 02_enrichment_analysis.R
# GO + KEGG 富集分析 + 可视化
# ============================

# 📦 依赖加载
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(DOSE)
  library(dplyr)
  library(stringr)
  library(readr)
})

# 🛠 设置路径
input_dir <- "results/deseq2/"
output_dir <- "results/enrichment/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 🔍 查找所有差异基因文件
files <- list.files(input_dir, pattern = "^sig_genes_.*\\.csv$", full.names = TRUE)

# 🔁 多组循环处理
for (f in files) {
  group_name <- str_extract(basename(f), "sig_genes_(.*)\\.csv") %>% str_replace("sig_genes_", "") %>% str_replace(".csv", "")
  message("📊 Analyzing enrichment for: ", group_name)
  
  degs <- read.csv(f)
  
  # 🧬 转换为 Entrez ID
  gene_symbols <- degs$Gene
  entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL",
                     toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>% na.omit()

  gene_list <- entrez_ids$ENTREZID

  # ======= GO 分析 =======
  ego <- enrichGO(gene = gene_list,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)

  go_out <- paste0(output_dir, "GO_", group_name, ".csv")
  write.csv(as.data.frame(ego), go_out, row.names = FALSE)

  # 🎨 GO 可视化
  pdf(paste0(output_dir, "GO_barplot_", group_name, ".pdf"), width = 8, height = 6)
  barplot(ego, showCategory = 20, title = paste("GO Enrichment -", group_name))
  dev.off()

  pdf(paste0(output_dir, "GO_dotplot_", group_name, ".pdf"), width = 8, height = 6)
  dotplot(ego, showCategory = 20, title = paste("GO Dotplot -", group_name))
  dev.off()

  # ======= KEGG 分析 =======
  ekegg <- enrichKEGG(gene = gene_list,
                      organism = 'hsa',
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05)

  kegg_out <- paste0(output_dir, "KEGG_", group_name, ".csv")
  write.csv(as.data.frame(ekegg), kegg_out, row.names = FALSE)

  # 🎨 KEGG 可视化
  pdf(paste0(output_dir, "KEGG_barplot_", group_name, ".pdf"), width = 8, height = 6)
  barplot(ekegg, showCategory = 20, title = paste("KEGG Enrichment -", group_name))
  dev.off()

  pdf(paste0(output_dir, "KEGG_dotplot_", group_name, ".pdf"), width = 8, height = 6)
  dotplot(ekegg, showCategory = 20, title = paste("KEGG Dotplot -", group_name))
  dev.off()
}

message("✅ 富集分析完成，结果保存在：", output_dir)