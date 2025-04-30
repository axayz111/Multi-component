library(shiny)
library(shinythemes)

# 加载后端脚本（假设你已经把它们封装成函数）
source("scripts/01_deseq2_analysis.R")
source("scripts/02_enrichment_analysis.R")
source("scripts/03_visualization.R")
source("scripts/04_prepare_IPA.R")

ui <- fluidPage(
  theme = shinytheme("cosmo"),
  titlePanel("🧬 造血干细胞多组分析平台"),
  sidebarLayout(
    sidebarPanel(
      fileInput("exprFile", "上传表达矩阵（csv）", accept = ".csv"),
      fileInput("groupFile", "上传分组信息（csv）", accept = ".csv"),

      selectInput("filetype", "输出图片格式", choices = c("png", "pdf", "svg"), selected = "png"),
      sliderInput("fig_width", "图片宽度 (inch)", min = 4, max = 15, value = 8),
      sliderInput("fig_height", "图片高度 (inch)", min = 4, max = 15, value = 6),
      selectInput("color_scheme", "颜色主题", choices = c("viridis", "Reds", "Blues", "Greens", "custom")),

      actionButton("run", "🚀 开始分析", class = "btn-success")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("火山图", plotOutput("volcanoPlot")),
        tabPanel("箱线图", plotOutput("boxPlot")),
        tabPanel("热图", plotOutput("heatmapPlot")),
        tabPanel("UMAP", plotOutput("umapPlot")),
        tabPanel("GO/KEGG", plotOutput("barPlot")),
        tabPanel("日志信息", verbatimTextOutput("log"))
      )
    )
  )
)

server <- function(input, output, session) {
  observeEvent(input$run, {
    req(input$exprFile, input$groupFile)
    expr_path <- input$exprFile$datapath
    group_path <- input$groupFile$datapath

    output$log <- renderText("🚧 正在执行差异表达分析 ...")

    # Step 1: DESeq2
    res <- run_deseq2(expr_path, group_path)

    # Step 2: Enrichment
    enrich <- run_enrichment(res)

    # Step 3: Visualization
    plots <- generate_all_plots(
      res,
      enrich,
      filetype = input$filetype,
      width = input$fig_width,
      height = input$fig_height,
      color = input$color_scheme
    )

    # Step 4: IPA
    prepare_ipa_input(res, outfile = "results/tables/ipa_input.txt")

    # 显示图表
    output$volcanoPlot <- renderPlot(plots$volcano)
    output$boxPlot     <- renderPlot(plots$box)
    output$heatmapPlot <- renderPlot(plots$heatmap)
    output$umapPlot    <- renderPlot(plots$umap)
    output$barPlot     <- renderPlot(plots$bar)

    output$log <- renderText("✅ 所有分析完成，结果已保存至 results/ 文件夹！")
  })
}

shinyApp(ui, server)