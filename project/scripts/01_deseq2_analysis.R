# ==============================================================================
# 01_deseq2_analysis_enhanced.R
# 差异表达分析模块（支持多组比较，增加绘图自定义选项）
# ==============================================================================

# 📦 加载必要的包
suppressPackageStartupMessages({
  library(DESeq2)
  library(tibble)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(matrixStats) # For rowVars used in heatmap
})

# ============================
# 🎯 参数设置区域 (可在此处修改)
# ============================

# --- 输入文件 ---
expr_file <- "data/expr.csv"        # 表达矩阵文件路径 (行:基因, 列:样本)
group_file <- "data/group.csv"      # 分组信息文件路径 (列: Sample, Group)

# --- 输出目录 ---
output_dir <- "results/deseq2_enhanced/" # 输出结果的总目录

# --- DESeq2 分析参数 ---
min_count_sum <- 10                 # 过滤低表达基因的阈值 (基因在所有样本中的总count)

# --- 火山图 (Volcano Plot) 参数 ---
volcano_padj_threshold <- 0.05      # 显著性 P 值阈值
volcano_lfc_threshold <- 1          # log2 Fold Change 阈值
volcano_top_n_genes <- 20           # 自动标注的最显著差异基因数量 (按 padj 排序)
volcano_sig_color <- "red"          # 显著差异基因点的颜色
volcano_nonsig_color <- "grey"      # 非显著差异基因点的颜色
volcano_point_alpha <- 0.6          # 点的透明度
volcano_label_size <- 3.0           # 标注基因名称的文字大小
volcano_title_prefix <- "Volcano Plot:" # 图标题的前缀
volcano_width <- 7                  # 输出 PNG/PDF 图片宽度 (英寸)
volcano_height <- 6                 # 输出 PNG/PDF 图片高度 (英寸)

# --- PCA 图参数 ---
pca_point_size <- 3                 # PCA 图中点的大小
pca_title <- "PCA Plot of Samples"  # PCA 图的标题
pca_width <- 6                      # 输出 PNG/PDF 图片宽度 (英寸)
pca_height <- 5                     # 输出 PNG/PDF 图片高度 (英寸)

# --- 热图 (Heatmap) 参数 ---
heatmap_top_n_genes <- 30           # 用于绘制热图的 Top N 变化最大基因数量 (基于 vst 转换后的方差)
heatmap_color_palette <- "RdBu"     # 热图颜色方案 (来自 RColorBrewer, 如 "RdBu", "YlGnBu", "Blues" 等)
heatmap_color_breaks <- 255         # 热图颜色分段数量
heatmap_show_rownames <- TRUE       # 是否在热图上显示基因名 (行名)
heatmap_cluster_rows <- TRUE        # 是否对行进行聚类
heatmap_cluster_cols <- TRUE        # 是否对列进行聚类
heatmap_scale <- "row"              # 热图数值缩放方式 ('row', 'column', 'none')
heatmap_title <- paste0("Heatmap of Top ", heatmap_top_n_genes, " Variable Genes") # 热图标题
heatmap_width <- 8                  # 输出 PNG 图片宽度 (英寸) - pheatmap 直接输出文件
heatmap_height <- 10                # 输出 PNG 图片高度 (英寸) - pheatmap 直接输出文件

# ============================
# 脚本主体部分 (一般无需修改)
# ============================

# --- 准备工作 ---
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message("Output directory: ", normalizePath(output_dir))

# --- 📥 读取数据 ---
message("Reading expression data from: ", expr_file)
count_data <- read.csv(expr_file, row.names = 1, check.names = FALSE)
message("Reading group information from: ", group_file)
group_info <- read.csv(group_file)

# --- 数据预处理 ---
# 清理样本名中的空格，确保与表达矩阵列名匹配
colnames(count_data) <- gsub("\\s+", "_", colnames(count_data))
group_info$Sample <- gsub("\\s+", "_", group_info$Sample)

# 检查分组文件中的样本是否存在于表达矩阵中
missing_samples <- setdiff(group_info$Sample, colnames(count_data))
if (length(missing_samples) > 0) {
  warning("Samples found in group file but not in expression data: ", paste(missing_samples, collapse=", "))
  group_info <- group_info[group_info$Sample %in% colnames(count_data), ]
  if (nrow(group_info) == 0) stop("No matching samples found between group file and expression data.")
}

# 检查表达矩阵中的样本是否存在于分组文件中
missing_groups <- setdiff(colnames(count_data), group_info$Sample)
if (length(missing_groups) > 0) {
  warning("Samples found in expression data but not in group file: ", paste(missing_groups, collapse=", "), ". These samples will be removed.")
  count_data <- count_data[, colnames(count_data) %in% group_info$Sample]
  if (ncol(count_data) == 0) stop("No matching samples found between expression data and group file.")
}

# 保证表达矩阵列序与分组文件行序一致
count_data <- count_data[, group_info$Sample]

message("Dimensions of count data (Genes x Samples): ", nrow(count_data), " x ", ncol(count_data))
message("Group information summary:")
print(table(group_info$Group))

# --- 🧪 构建 DESeq2 对象 ---
dds <- DESeqDataSetFromMatrix(
  countData = round(count_data), # DESeq2 需要整数 counts
  colData = group_info,
  design = ~ Group
)

# --- 🔍 过滤低表达基因 ---
message("Filtering genes with total counts <= ", min_count_sum)
keep <- rowSums(counts(dds)) > min_count_sum
dds <- dds[keep, ]
message("Number of genes after filtering: ", nrow(dds))

# --- 🚀 差异分析主流程 ---
message("Running DESeq2 analysis...")
dds <- DESeq(dds)
message("DESeq2 analysis complete.")

# --- 🌟 获取比较组 ---
group_levels <- levels(factor(group_info$Group))
if (length(group_levels) < 2) {
  stop("Need at least two groups for differential expression analysis. Found groups: ", paste(group_levels, collapse=", "))
}
combinations <- combn(group_levels, 2, simplify = FALSE)
message("Found ", length(combinations), " comparison pairs.")

# --- 🔁 多组比较 ---
for (pair in combinations) {
  group1 <- pair[1]
  group2 <- pair[2]
  comparison_name <- paste0(group2, "_vs_", group1)
  message("Comparing: ", group2, " vs ", group1)

  # 获取结果并进行 LFC shrinkage
  res <- results(dds, contrast = c("Group", group2, group1), alpha = volcano_padj_threshold)
  # 使用 ashr 进行 LFC shrinkage 可以提高准确性
  message("Performing LFC shrinkage using 'ashr'...")
  res_shrink <- lfcShrink(dds, contrast = c("Group", group2, group1), res = res, type = "ashr")

  # 转换为数据框并添加基因名
  res_df <- as.data.frame(res_shrink) %>%
    rownames_to_column("Gene") %>%
    arrange(padj) # 按调整后的 p 值排序

  # 输出差异表达结果表格
  out_file <- file.path(output_dir, paste0("DEG_", comparison_name, ".csv"))
  message("Writing results to: ", out_file)
  write.csv(res_df, out_file, row.names = FALSE)

  # --- 🔥 绘图 ---
  plot_dir <- file.path(output_dir, paste0("plots_", comparison_name))
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  message("Saving plots to: ", plot_dir)

  # -- 火山图 (Volcano Plot) --
  res_plot_df <- res_df %>%
    filter(!is.na(padj) & !is.na(log2FoldChange)) # 过滤掉NA值

  # 添加显著性标签
  res_plot_df$Significance <- "Not Sig"
  res_plot_df$Significance[res_plot_df$padj < volcano_padj_threshold & abs(res_plot_df$log2FoldChange) > volcano_lfc_threshold] <- "Significant"
  res_plot_df$Significance <- factor(res_plot_df$Significance, levels = c("Significant", "Not Sig")) # 控制图例顺序

  # 自动标注基因：选择最显著的 N 个基因
  top_genes <- res_plot_df %>%
    filter(Significance == "Significant") %>%
    arrange(padj) %>%
    head(volcano_top_n_genes) # 使用参数控制数量

  message("Labeling top ", min(nrow(top_genes), volcano_top_n_genes), " significant genes on volcano plot.")

  p_volcano <- ggplot(res_plot_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = Significance), alpha = volcano_point_alpha, size = 1.5) +
    scale_color_manual(values = c("Significant" = volcano_sig_color, "Not Sig" = volcano_nonsig_color)) +
    geom_text_repel(
        data = top_genes,
        aes(label = Gene),
        size = volcano_label_size,
        max.overlaps = Inf, # 尽量显示所有标签
        box.padding = 0.5,
        point.padding = 0.2,
        segment.color = 'grey50',
        segment.size = 0.5
    ) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste(volcano_title_prefix, comparison_name),
      x = "log2 Fold Change",
      y = "-log10 Adjusted p-value"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    ) +
    geom_hline(yintercept = -log10(volcano_padj_threshold), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = c(-volcano_lfc_threshold, volcano_lfc_threshold), linetype = "dashed", color = "grey50")

  # 保存火山图
  ggsave(filename = file.path(plot_dir, "volcano_plot.png"), plot = p_volcano, width = volcano_width, height = volcano_height, dpi = 300)
  ggsave(filename = file.path(plot_dir, "volcano_plot.pdf"), plot = p_volcano, width = volcano_width, height = volcano_height)

  # -- 保存显著差异基因列表 --
  sig_genes <- res_df %>%
    filter(padj < volcano_padj_threshold & abs(log2FoldChange) > volcano_lfc_threshold) %>%
    arrange(padj)
  sig_out_file <- file.path(output_dir, paste0("sig_genes_", comparison_name, ".csv"))
  message("Writing significant gene list to: ", sig_out_file)
  write.csv(sig_genes, sig_out_file, row.names = FALSE)

} # End of comparison loop

# --- 📊 PCA 可视化 (所有样本) ---
message("Generating PCA plot...")
# 使用 VST (Variance Stabilizing Transformation) 进行数据转换，适用于中大样本量
# 对于小样本量 (<30)，rlog 可能更好，但计算较慢
vsd <- vst(dds, blind = TRUE) # blind=TRUE 推荐用于可视化，忽略实验设计

pcaData <- plotPCA(vsd, intgroup = "Group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Group, shape = Group)) + # 可选：用形状区分组
  geom_point(size = pca_point_size) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(title = pca_title) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom") +
  coord_fixed() # 使 x 和 y 轴比例相同

# 保存 PCA 图
pca_plot_file_base <- file.path(output_dir, "PCA_plot_all_samples")
ggsave(paste0(pca_plot_file_base, ".png"), plot = p_pca, width = pca_width, height = pca_height, dpi = 300)
ggsave(paste0(pca_plot_file_base, ".pdf"), plot = p_pca, width = pca_width, height = pca_height)

# --- 🎨 热图 (Heatmap) ---
message("Generating heatmap of top ", heatmap_top_n_genes, " variable genes...")
# 计算 vst 转换后数据的行方差，找出变化最大的基因
topVarGenesIndices <- head(order(rowVars(assay(vsd)), decreasing = TRUE), heatmap_top_n_genes)
topVarGeneNames <- rownames(assay(vsd))[topVarGenesIndices] # 获取基因名

mat <- assay(vsd)[topVarGeneNames, ] # 提取 Top N 基因的表达数据

# 根据参数进行缩放
if (heatmap_scale == "row") {
  mat <- t(scale(t(mat)))
  message("Heatmap data scaled by row (gene).")
} else if (heatmap_scale == "column") {
  mat <- scale(mat)
  message("Heatmap data scaled by column (sample).")
} else {
   message("Heatmap data not scaled.")
}
# 处理缩放后可能产生的 NA/NaN/Inf 值 (如果某行方差为0)
mat[is.na(mat)] <- 0
mat[!is.finite(mat)] <- 0

# 准备列注释信息 (样本分组)
ann_col <- as.data.frame(colData(vsd)[, "Group", drop = FALSE])
rownames(ann_col) <- colnames(mat) # 确保行名匹配

# 定义颜色
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = heatmap_color_palette)))(heatmap_color_breaks)

# 绘制热图并直接保存到文件
heatmap_file_base <- file.path(output_dir, paste0("heatmap_top", heatmap_top_n_genes, "_variable_genes"))

# 使用 pheatmap 绘制热图
pheatmap(
  mat,
  annotation_col = ann_col,
  color = heatmap_colors,
  cluster_rows = heatmap_cluster_rows,
  cluster_cols = heatmap_cluster_cols,
  show_rownames = heatmap_show_rownames,
  show_colnames = FALSE, # 通常不显示样本名，除非样本少
  scale = "none", # 数据已手动缩放，此处设为 none
  main = heatmap_title,
  fontsize_row = ifelse(heatmap_show_rownames, 8, 0), # 如果显示行名，设置字体大小
  filename = paste0(heatmap_file_base, ".png"),
  width = heatmap_width,
  height = heatmap_height
)
# 如果需要 PDF 输出 (pheatmap 对 PDF 的直接支持有时效果不如 PNG)
pdf(paste0(heatmap_file_base, ".pdf"), width = heatmap_width, height = heatmap_height)
pheatmap(
  mat,
  annotation_col = ann_col,
  color = heatmap_colors,
  cluster_rows = heatmap_cluster_rows,
  cluster_cols = heatmap_cluster_cols,
  show_rownames = heatmap_show_rownames,
  show_colnames = FALSE,
  scale = "none",
  main = heatmap_title,
  fontsize_row = ifelse(heatmap_show_rownames, 8, 0)
)
dev.off()


message("✅ DESeq2 differential expression analysis and visualization complete.")
message("Results saved in: ", normalizePath(output_dir))
