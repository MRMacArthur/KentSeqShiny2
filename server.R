library(ggplot2)
library(plotly)
library(shiny)
library(DT)
library(edgeR)
library(statmod)
library(GO.db)
library(org.Mm.eg.db)

x <- readRDS("./Data/naive_LPS_WT_KO.RData")
group <- factor(c(rep(1:4, each = 3)),
                labels = c("WT_n", "WT_6h", "KO_n", "KO_6h"))
designMat <- model.matrix(~0+group)
colnames(designMat) <- levels(group)
fit <- glmQLFit(x, designMat, robust = T)
rpkmMat <- readRDS("Data/naiveRPKMmat.RData")

xT <- readRDS("Data/T_TL_DG_NAM.RData")
groupT <- factor(c(rep(1:6, each = 3)),
                 labels = c("WT_T", "WT_TL", "KO_T", "KO_TL",
                            "WT_TL_DG", "WT_TL_NAM"))
designMatT <- model.matrix(~0+groupT)
colnames(designMatT) <- levels(groupT)
fitT <- glmQLFit(xT, designMatT, robust = T)
rpkmMatT <- readRDS("Data/T_DG_NAM_RPKMmat.RData")


function(input, output, session) {
  
  makeDEtab <- reactive({
    
    comparison <- makeContrasts(contrasts = paste(as.character(input$compChoose)), levels = designMat)
    tr <- glmTreat(fit, contrast = comparison, 
                   lfc = log2(1.2), null = "worst.case")
    trTags <- topTags(tr, n = nrow(tr$table))
    trTab <- trTags$table
    trTab$Significant <- F
    trTab$Significant[trTab$FDR<0.05 & abs(trTab$logFC) > 1.1 ] <- T
    return(trTab)
    
  })
  
  makeDEtabT <- reactive({
    
    comparisonT <- makeContrasts(contrasts = paste(as.character(input$compChooseT)), levels = designMatT)
    trT <- glmTreat(fitT, contrast = comparisonT, 
                   lfc = log2(1.2), null = "worst.case")
    trTagsT <- topTags(trT, n = nrow(trT$table))
    trTabT <- trTagsT$table
    trTabT$Significant <- F
    trTabT$Significant[trTabT$FDR<0.05 & abs(trTabT$logFC) > 1.1 ] <- T
    return(trTabT)
    
  })
  
  makePathway <- reactive({
    
    comparison <- makeContrasts(contrasts = paste(as.character(input$compChoose)), levels = designMat)
    tr <- glmTreat(fit, contrast = comparison, 
                   lfc = log2(1.2), null = "worst.case")
    goann <- goana(tr, species = "Mm")
    keggan <- kegga(tr, species = "Mm")
    
    if (input$pathCheck == T){
    pathOut <- topKEGG(keggan, n = 30)
    } else {
    pathOut <- topGO(goann, n = 30)
    }
  
    return(pathOut)
      
  })
  
  makePathwayT <- reactive({
    
    comparisonT <- makeContrasts(contrasts = paste(as.character(input$compChooseT)), levels = designMatT)
    trT <- glmTreat(fitT, contrast = comparisonT, 
                   lfc = log2(1.2), null = "worst.case")
    goannT <- goana(trT, species = "Mm")
    kegganT <- kegga(trT, species = "Mm")
    
    if (input$pathCheckT == T){
      pathOutT <- topKEGG(kegganT, n = 30)
    } else {
      pathOutT <- topGO(goannT, n = 30)
    }
    
    return(pathOutT)
    
  })
  
  output$volcanoOut <- renderPlotly({
    
    p <- ggplot(data = makeDEtab(), 
                aes(x = logFC, y = -log10(FDR), 
                    color = Significant,
                    text = paste("Gene:", make.names(Symbol)),
                    key = make.names(Symbol))) +
      geom_point(alpha = (1/4)) +
      labs(x = "Log2 Fold Change", y = "-log10 q Value") +
      guides(color = guide_legend(title = "Significant"))
    
    p <- ggplotly(p, source = "nVolcano")
    
    print(p)
    
  })
  
  output$diffTable <- DT::renderDataTable(DT::datatable({
    
    makeDEtab()
    
  }))
  
  output$barPlotNaive <- renderPlotly({
    
    event.data <- event_data("plotly_click", source = "nVolcano")
    if (is.null(event.data)) {
      
      p <- ggplot(rpkmMat, aes(x = group, y = Cxcl1)) + 
        geom_boxplot()
      p <- ggplotly(p)
      
    } else {
      
      p <- ggplot(rpkmMat, aes_string(x = "group", y = event.data$key)) + 
        geom_boxplot()
      p <- ggplotly(p)
      
    }
    
  })
  
  output$barPlotNaive2 <- renderPlotly({
    
    p <- ggplot(rpkmMat, aes_string(x = "group", y = input$geneChoose)) + 
      geom_boxplot()
    p <- ggplotly(p)
    
  })
  
  output$pathPrint <- renderPrint({
    
    makePathway()
    
  })
  
  output$volcanoOutT <- renderPlotly({
    
    p <- ggplot(data = makeDEtabT(), 
                aes(x = logFC, y = -log10(FDR), 
                    color = Significant,
                    text = paste("Gene:", make.names(Symbol)),
                    key = make.names(Symbol))) +
      geom_point(alpha = (1/4)) +
      labs(x = "Log2 Fold Change", y = "-log10 q Value") +
      guides(color = guide_legend(title = "Significant"))
    
    p <- ggplotly(p, source = "nVolcanoT")
    
    print(p)
    
  })
  
  output$diffTableT <- DT::renderDataTable(DT::datatable({
    
    makeDEtabT()
    
  }))
  
  output$barPlotNaiveT <- renderPlotly({
    
    if (input$dgCheck == F) {
      rpkmMatTp <- rpkmMatT[1:12,]
    } else {
      rpkmMatTp <- rpkmMatT
    }
    
    event.data <- event_data("plotly_click", source = "nVolcanoT")
    if (is.null(event.data)) {
      
      p <- ggplot(rpkmMatTp, aes(x = group, y = Cxcl1)) + 
        geom_boxplot()
      p <- ggplotly(p)
      
    } else {
      
      p <- ggplot(rpkmMatTp, aes_string(x = "group", y = event.data$key)) + 
        geom_boxplot()
      p <- ggplotly(p)
      
    }
    
  })
  
  output$barPlotNaive2T <- renderPlotly({
    
    if (input$dgCheck == F) {
      rpkmMatTp <- rpkmMatT[1:12,]
    } else {
      rpkmMatTp <- rpkmMatT
    }
    
    
    p <- ggplot(rpkmMatTp, aes_string(x = "group", y = input$geneChooseT)) + 
      geom_boxplot()
    p <- ggplotly(p)
    
  })
  
  output$pathPrintT <- renderPrint({
    
    makePathwayT()
    
  })
  
}