# ============================
# 03_visualization.R
# 差异表达可视化汇总
# ============================

# 📦 加载包
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(VennDiagram)
  library(EnhancedVolcano)
  library(DESeq2)
  library(pheatmap)
  library(umap)
  library(ggpubr)
  library(RColorBrewer)
  library(data.table)
  library(stringr)
})

# 🛠 路径设置
input_dir <- "results/deseq2/"
plot_dir <- "results/plots/"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# 📂 加载差异表达结果
files <- list.files(input_dir, pattern = "^DEG_.*\\.csv$", full.names = TRUE)

# 🎨 配色设置
cols <- brewer.pal(3, "Set1")

for (file in files) {
  # 📊 提取比较组名称
  group_name <- str_extract(basename(file), "DEG_(.*)\\.csv") %>% str_replace(".csv", "")
  df <- fread(file)

  # ✅ 火山图
  pdf(paste0(plot_dir, "volcano_", group_name, ".pdf"), width = 7, height = 6)
  EnhancedVolcano(df,
                  lab = df$Gene,
                  x = 'log2FoldChange',
                  y = 'padj',
                  title = paste("Volcano Plot:", group_name),
                  pCutoff = 0.05,
                  FCcutoff = 1.5,
                  pointSize = 2.5,
                  labSize = 3.0)
  dev.off()

  # ✅ MA 图
  pdf(paste0(plot_dir, "MAplot_", group_name, ".pdf"), width = 6, height = 5)
  with(df, plot(baseMean, log2FoldChange, log = "x", pch = 20, col = ifelse(padj < 0.05, "red", "gray"),
                main = paste("MA Plot:", group_name),
                xlab = "Mean Expression", ylab = "log2 Fold Change"))
  abline(h = c(-1, 1), col = "blue", lty = 2)
  dev.off()

  # ✅ 表达密度图（如有表达矩阵）
  expr_file <- "data/expr.csv"
  if (file.exists(expr_file)) {
    expr <- fread(expr_file)
    expr_melt <- reshape2::melt(expr)
    pdf(paste0(plot_dir, "density_", group_name, ".pdf"), width = 7, height = 5)
    ggplot(expr_melt, aes(x = value, color = variable)) +
      geom_density() +
      theme_minimal() +
      labs(title = "Expression Density", x = "Expression", y = "Density")
    dev.off()
  }

  # ✅ Q-Q 图
  pdf(paste0(plot_dir, "QQ_", group_name, ".pdf"), width = 5, height = 5)
  observed <- -log10(sort(df$pvalue))
  expected <- -log10(ppoints(length(df$pvalue)))
  plot(expected, observed, main = paste("Q-Q Plot:", group_name), xlab = "Expected", ylab = "Observed", pch = 20)
  abline(0, 1, col = "red")
  dev.off()

  # ✅ 箱线图（表达值）
  if (exists("expr")) {
    pdf(paste0(plot_dir, "boxplot_", group_name, ".pdf"), width = 8, height = 6)
    boxplot(as.data.frame(expr), las = 2, main = "Expression Value Boxplot")
    dev.off()
  }
}

# ✅ Venn 图（基因交集）
if (length(files) >= 2) {
  deg_sets <- list()
  for (f in files) {
    name <- str_extract(basename(f), "DEG_(.*)\\.csv") %>% str_replace(".csv", "")
    deg <- fread(f)
    sig <- deg[padj < 0.05, Gene]
    deg_sets[[name]] <- sig
  }

  venn_file <- paste0(plot_dir, "venn_genes.pdf")
  pdf(venn_file, width = 6, height = 6)
  venn.diagram(deg_sets, filename = NULL, fill = cols[1:length(deg_sets)],
               main = "Venn of DEGs", main.cex = 1.5)
  dev.off()
}

# ✅ UMAP（降维）
expr_file <- "data/expr.csv"
group_file <- "data/group.csv"
if (file.exists(expr_file) & file.exists(group_file)) {
  expr <- fread(expr_file)
  group <- fread(group_file)
  rownames(expr) <- expr$Gene
  expr <- expr[, -1, with = FALSE]

  # 计算样本数量
  n_samples <- ncol(expr)
  n_neighbors <- min(15, n_samples - 1)  # 避免邻居数超过样本数

  # 设置UMAP参数
  umap_config <- umap.defaults
  umap_config$n_neighbors <- n_neighbors

  # 执行UMAP
  umap_result <- umap(t(expr), config = umap_config)
  umap_df <- as.data.frame(umap_result$layout)
  umap_df$group <- group$Group

  # 保存UMAP图
  pdf(paste0(plot_dir, "UMAP.pdf"), width = 6, height = 6)
  ggplot(umap_df, aes(x = V1, y = V2, color = group)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = "UMAP of Samples", x = "UMAP1", y = "UMAP2")
  dev.off()
}

# ✅ 均值-方差趋势图
if (exists("expr")) {
  means <- rowMeans(as.matrix(expr))
  sds <- apply(as.matrix(expr), 1, sd)
  pdf(paste0(plot_dir, "mean_variance_trend.pdf"), width = 6, height = 5)
  plot(log10(means), log10(sds), pch = 20, col = "steelblue",
       xlab = "log10 Mean Expression", ylab = "log10 SD",
       main = "Mean-Variance Trend")
  dev.off()
}

message("✅ 所有可视化已完成，查看: ", plot_dir)