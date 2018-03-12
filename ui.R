library(ggplot2)
library(plotly)
library(shiny)
library(DT)
library(edgeR)
library(statmod)
library(GO.db)
library(org.Mm.eg.db)

navbarPage("PKL Seq Data Jr",
           tabPanel("Naive v LPS",
                    sidebarLayout(
                          sidebarPanel(width = 2,
                                     selectInput("compChoose", label = "Choose Sample Comparison",
                                                 choices = c( "KO_n - WT_n", "KO_6h - WT_6h",
                                                              "KO_6h - WT_n", "WT_6h - KO_n"),
                                                 selected = "KO_n - WT_n"),
                                     selectInput("geneChoose", label = "Gene for DiffTab Barplot",
                                                 choices = colnames(readRDS("Data/naiveRPKMmat.RData"))),
                                     checkboxInput("pathCheck", "KEGG Pathway (uncheck for GO)", T)
                              ),
                            mainPanel(
                              tabsetPanel(type = "tabs",
                                          tabPanel(
                                              "Volano Naive", plotlyOutput('volcanoOut', height = 600),
                                              plotlyOutput("barPlotNaive"),
                                              verbatimTextOutput("pathPrint")),
                                          tabPanel("Diff Tables", DT::dataTableOutput("diffTable"),
                                                   plotlyOutput("barPlotNaive2"))
                                    )))),
           tabPanel("Tolerant v LPS",
                    sidebarLayout(
                      sidebarPanel(width = 2,
                                    selectInput("compChooseT", label = "Choose Sample Comparison",
                                                choices = c( "KO_T - WT_T", "KO_TL - WT_TL",
                                                             "WT_TL_DG - WT_TL", "WT_TL_NAM - WT_TL",
                                                             "WT_TL_DG - KO_TL", "WT_TL_NAM - KO_TL", 
                                                             "WT_TL_DG - WT_TL_NAM")
                                                             ),
                                    selectInput("geneChooseT", label = "Gene for DiffTab Barplot",
                                                choices = colnames(readRDS("Data/tolerantRPKMmat.RData"))),
                                    checkboxInput("pathCheckT", "KEGG Pathway (uncheck for GO)", T),
                                   checkboxInput("dgCheck", "Boxplot Include 2DG/NAM", F)
                                    ),
                      mainPanel(
                        tabsetPanel(type = "tabs",
                                    tabPanel(
                                      "Volano Tolerant", plotlyOutput('volcanoOutT', height = 600),
                                      plotlyOutput("barPlotNaiveT"),
                                      verbatimTextOutput("pathPrintT")),
                                    tabPanel("Diff Tables", DT::dataTableOutput("diffTableT"),
                                             plotlyOutput("barPlotNaive2T")
                                             )
                                    )
                        )
                      )
                    )
)
