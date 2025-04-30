# ==============================
# 📦 Makefile: 一键执行分析流程
# ==============================

# Rscript 路径
RSCRIPT = Rscript

# 脚本路径
SCRIPT_DIR = scripts

# 默认目标：执行全部分析步骤
all: deseq2 enrichment visualize ipa

# 差异表达分析
deseq2:
	@echo "🔍 Step 1: Running DESeq2 analysis..."
	$(RSCRIPT) $(SCRIPT_DIR)/01_deseq2_analysis.R

# 富集分析
enrichment: deseq2
	@echo "📈 Step 2: Running enrichment analysis (GO + KEGG)..."
	$(RSCRIPT) $(SCRIPT_DIR)/02_enrichment_analysis.R

# 结果可视化
visualize: enrichment
	@echo "🖼️  Step 3: Generating visualizations..."
	$(RSCRIPT) $(SCRIPT_DIR)/03_visualization.R

# 准备 IPA 输入文件
ipa: deseq2
	@echo "🧠 Step 4: Preparing IPA input files..."
	$(RSCRIPT) $(SCRIPT_DIR)/04_prepare_IPA.R

# 启动 Shiny 界面（可选）
shiny:
	@echo "🚀 Launching Shiny interface..."
	$(RSCRIPT) app.R

# 清理结果（谨慎）
clean:
	@echo "🧹 Cleaning up result folders..."
	rm -rf results/*

.PHONY: all deseq2 enrichment visualize ipa shiny clean