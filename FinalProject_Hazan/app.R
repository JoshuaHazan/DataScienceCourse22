#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(stats)
library(BiocGenerics)
library(DESeq2)
library(pheatmap)
library(BASiCS)
library(SingleCellExperiment)
library(Rtsne)
library(corrplot)
library(gprofiler2)
library(gage)

set.seed(592)

# Raw data
normalized_counts <- read.table("C:/Users/besterlab/Dropbox (Technion Dropbox)/Josh/Genome Data Science Course/Final Project/normalized_data.txt", sep = "\t", header = TRUE)
metadata <- read.table("C:/Users/besterlab/Dropbox (Technion Dropbox)/Josh/Genome Data Science Course/Final Project/metadata.txt", sep = "\t", header = TRUE)

# PCA data
MACS <- metadata %>%
  dplyr::filter(Celltype %in% "MACS-purified Naive") %>%
  dplyr::filter(Strain %in% "Mus musculus domesticus")
PCA_Data <- normalized_counts %>%
  dplyr::select(MACS$CellName)
# Filtering for age
Young_B6 <- MACS %>%
  dplyr::filter(Age %in% "Young") %>%
  dplyr::pull(CellName)
Old_B6 <- MACS %>%
  dplyr::filter(Age %in% "Old") %>%
  dplyr::pull(CellName)
# Filtering for activation
Active_pca <- MACS %>%
  dplyr::filter(Stimulus %in% c("Active","Activated")) %>%
  dplyr::pull(CellName)
Unstim_pca <- MACS %>%
  dplyr::filter(Stimulus %in% "Unstimulated") %>%
  dplyr::pull(CellName)
activity_pca <- vector(length = ncol(PCA_Data))
activity_pca [which(grepl(paste(Active_pca, collapse = "|"), colnames(PCA_Data)))] <- "Active"
activity_pca [which(grepl(paste(Unstim_pca, collapse = "|"), colnames(PCA_Data)))] <- "Unstimulated"
age <- vector(length = ncol(PCA_Data))
age [which(grepl(paste(Young_B6, collapse = "|"), colnames(PCA_Data)))] <- "Young"
age [which(grepl(paste(Old_B6, collapse = "|"), colnames(PCA_Data)))] <- "Old"
# Calculating PCA
pca <- prcomp(t(log10(PCA_Data + 1)))
pca.df <- data.frame(pca.1 = pca$x[,1], pca.2 = pca$x[,2], age = age, activity = activity_pca)
pca.df.young <- subset(pca.df, age %in% "Young")
pca.df.old <- subset(pca.df, age %in% "Old")

# tSNE
# Extracting only old B6 mice.
Old <- metadata %>%
  dplyr::filter(Age %in% "Old") %>%
  dplyr::filter(Strain %in% "Mus musculus domesticus")

TSNE_Data <- normalized_counts %>%
  dplyr::select(Old$CellName)
# Filtering each cell type.
Cell1 <- Old %>%
  dplyr::filter(Celltype %in% "MACS-purified Naive") %>%
  dplyr::pull(CellName)
Cell2 <- Old %>%
  dplyr::filter(Celltype %in% "FACS-purified Naive") %>%
  dplyr::pull(CellName)
Cell3 <- Old %>%
  dplyr::filter(Celltype %in% "FACS-purified Effector Memory") %>%
  dplyr::pull(CellName)
# Filtering activity.
Active_tsne <- Old %>%
  dplyr::filter(Stimulus %in% c("Active","Activated")) %>%
  dplyr::pull(CellName)
Unstim_tsne <- Old %>%
  dplyr::filter(Stimulus %in% "Unstimulated") %>%
  dplyr::pull(CellName)
# Differentiating by cell type and naive vs active.
cell_type <- vector(length = ncol(TSNE_Data))
cell_type[which(grepl(paste(Cell1, collapse = "|"), colnames(TSNE_Data)))] <- "MACS-purified Naive" 
cell_type[which(grepl(paste(Cell2, collapse = "|"), colnames(TSNE_Data)))] <- "FACS-purified Naive" 
cell_type[which(grepl(paste(Cell3, collapse = "|"), colnames(TSNE_Data)))] <- "FACS-purified Effector Memory"

activity_tsne <- vector(length = ncol(TSNE_Data))
activity_tsne [which(grepl(paste(Active_tsne, collapse = "|"), colnames(TSNE_Data)))] <- "Active"
activity_tsne [which(grepl(paste(Unstim_tsne, collapse = "|"), colnames(TSNE_Data)))] <- "Unstimulated"
# Calculating tSNE.
tsne <- Rtsne(t(log10(TSNE_Data + 1)), perplexity = 10)
tsne.df <- data.frame(tsne.1 = tsne$Y[,1], tsne.2 = tsne$Y[,2], cell_type = cell_type, activity = activity_tsne)
tsne.df.FACSeff <- subset(tsne.df, cell_type %in% "FACS-purified Effector Memory")
tsne.df.FACSn <- subset(tsne.df, cell_type %in% "FACS-purified Naive")
tsne.df.MACSn <- subset(tsne.df, cell_type %in% "MACS-purified Naive")

# Correlation plot
MACS_corr <- metadata %>%
  dplyr::filter(Celltype %in% "MACS-purified Naive") %>%
  dplyr::filter(Strain %in% "Mus musculus domesticus")
Corr_counts <- normalized_counts %>%
  dplyr::select(MACS$CellName)
correlationMatrix <- cor(Corr_counts)
ann_col <- MACS %>%
  dplyr::select(CellName, Stimulus, Individuals) %>%
  column_to_rownames(var = "CellName")







# Define UI 

ui <- shinyUI(fluidPage(
  titlePanel("Hazan Final Project"),
  sidebarLayout(
    sidebarPanel(
      selectInput("plot","Plot to Display",choices = c("PCA", "tSNE", "Correlation plot")),
      conditionalPanel(
        condition = "input.plot == 'PCA'", 
        radioButtons("age","Age of mice", choices = c("Young", "Old", "All"))),
      conditionalPanel(
        condition = "input.plot == 'tSNE'", 
        radioButtons("cell","Cell Type", choices = c("FACS-purified effector memory", "FACS-purified naive", 
                                                     "MACS-purified naive", "All"))),
      conditionalPanel(
        condition = "input.plot == 'Correlation plot'", 
        numericInput("corr_row","Row divisions", value = 3)),
      conditionalPanel(
        condition = "input.plot == 'Correlation plot'", 
        numericInput("corr_col","Column divisions", value = 3)),
      uiOutput("filter_degree")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      uiOutput("plot")
    )
  )))

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  output$plot <- renderUI({
    if(input$plot=="PCA"){
      if(input$age=="Young"){
        
        output$plot1<-renderPlot({
          ggplot(data = pca.df.young, aes(pca.1, pca.2)) + 
            geom_point(size = 4, mapping = aes(fill = activity), shape = 22) +
            scale_fill_manual(values = c("tomato2", "grey40")) +
            theme_minimal() + ylab("PCA 1") + xlab("PCA 2") +
            labs(title = "PCA for cells from young mice")
        })
        plotOutput("plot1")
      }
      
      else if(input$age=="Old"){
        output$plot2<-renderPlot({
          ggplot(data = pca.df.old, aes(pca.1, pca.2)) + 
            geom_point(size = 4, mapping = aes(fill = activity), shape = 24) +
            scale_fill_manual(values = c("tomato2", "grey40")) +
            theme_minimal() + ylab("PCA 1") + xlab("PCA 2") +
            labs(title = "PCA for cells from old mice")
        })
        plotOutput("plot2")
      }
      
      else if(input$age=="All"){
        output$plot3<-renderPlot({
          ggplot(data = pca.df, aes(pca.1, pca.2)) + 
            geom_point(size = 4, mapping = aes(fill = activity, shape = age)) +
            scale_fill_manual(values = c("tomato2", "grey40")) +
            scale_shape_manual(values = c(22, 24)) +
            theme_minimal() + ylab("PCA 1") + xlab("PCA 2") +
            labs(title = "PCA for cells from young and old mice") +
            guides(fill=guide_legend(override.aes=list(shape=21)))
        })
        plotOutput("plot3")
      }
    }
    
    else if(input$plot=="tSNE"){
      if(input$cell=="FACS-purified effector memory"){
        
        output$plot4<-renderPlot({
          ggplot(data = tsne.df.FACSeff, aes(tsne.1, tsne.2)) + 
            geom_point(size = 4, mapping = aes(fill = activity), shape = 22) +
            scale_fill_manual(values = c("coral4", "tan3")) +
            theme_minimal() + ylab("tSNE 1") + xlab("tSNE 2") +
            labs(title = "tSNE for FACS-purified effector memory cells in old mice")
        })
        plotOutput("plot4")
      }
      
      else if(input$cell=="FACS-purified naive"){
        output$plot5<-renderPlot({
          ggplot(data = tsne.df.FACSn, aes(tsne.1, tsne.2)) + 
            geom_point(size = 4, mapping = aes(fill = activity), shape = 21) +
            scale_fill_manual(values = c("coral4", "tan3")) +
            theme_minimal() + ylab("tSNE 1") + xlab("tSNE 2") +
            labs(title = "tSNE for FACS-purified naive T4 cells in old mice")
        })
        plotOutput("plot5")
      }
      
      else if(input$cell=="MACS-purified naive"){
        output$plot6<-renderPlot({
          ggplot(data = tsne.df.MACSn, aes(tsne.1, tsne.2)) + 
            geom_point(size = 4, mapping = aes(fill = activity), shape = 24) +
            scale_fill_manual(values = c("coral4", "tan3")) +
            theme_minimal() + ylab("tSNE 1") + xlab("tSNE 2") +
            labs(title = "tSNE for MACS-purified naive T4 cells in old mice")
        })
        plotOutput("plot6")
      }
      else if(input$cell=="All"){
        output$plot7<-renderPlot({
          ggplot(data = tsne.df, aes(tsne.1, tsne.2)) + 
            geom_point(size = 4, mapping = aes(fill = activity, shape = cell_type)) +
            scale_fill_manual(values = c("coral4", "tan3")) +
            scale_shape_manual(values = c(22,21,24)) +
            theme_minimal() + ylab("tSNE 1") + xlab("tSNE 2") +
            labs(title = "tSNE for all cell types in old mice") +
            guides(fill=guide_legend(override.aes=list(shape=21)))
        })
        plotOutput("plot7")
      }
    }
    
    else if(input$plot=="Correlation plot"){
      output$plot8<-renderPlot({
        pheatmap(correlationMatrix,  
                 annotation_col = ann_col, 
                 cutree_cols = input$corr_col,
                 cutree_rows = input$corr_row,
                 show_rownames = FALSE,
                 show_colnames = FALSE,
                 main = "Heatmap showing gene correlation between groups")
      })
      plotOutput("plot8")
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
