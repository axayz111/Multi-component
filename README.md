# Multi-component
# 🧬 Multi-component hematopoietic stem cell differentiation analysis platform: 多组造血干细胞分化分析平台

基于 RNA-seq 数据，自动完成差异表达分析、GO/KEGG 富集、IPA 输入准备和可视化展示。支持 Shiny 图形化界面。

---

## 📦 项目结构
```
my_analysis/
├── my_analysis.Rproj           # RStudio 项目文件
├── app.R                       # Shiny 可视化主程序
├── scripts/                    # 分析脚本模块
│   ├── 01_deseq2_analysis.R
│   ├── 02_enrichment_analysis.R
│   ├── 03_visualization.R
│   └── 04_prepare_IPA.R
├── data/                       # 输入数据目录（表达矩阵 + 分组信息）
├── results/                    # 输出目录（表格 + 图片）
├── Makefile                    # 一键运行所有流程
└── README.md                   # 项目说明文件
```
---


🧪 分析功能包括
	•	✅ 差异表达分析（DESeq2）
	•	🎯 GO / KEGG 富集分析
	•	🧠 IPA 输入准备
	•	📊 可视化图表：
	•	火山图、箱线图、热图、UMAP、均值-方差图
	•	Q-Q图、密度分布图、Venn图等

---

🛠️ 依赖包（可通过 install.packages 安装）
```
library(shiny)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(ggpubr)
library(RColorBrewer)
library(umap)
```


---

📧 联系

作者：zzx
邮箱：lidachuilige123@gmail
项目仅供科研用途

---

