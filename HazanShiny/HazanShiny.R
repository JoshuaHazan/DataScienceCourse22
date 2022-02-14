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
library(compGenomRData)
library(pheatmap)
library(stats)
library(ggplot2)
library(ggfortify)
library(stats)
library(corrplot)

#Required files
counts_file <- system.file("extdata/rna-seq/SRP029880.raw_counts.tsv",
                           package = "compGenomRData")
coldata_file <- system.file("extdata/rna-seq/SRP029880.colData.tsv",
                            package = "compGenomRData")
counts <- as.matrix(read.table(counts_file, header = T, sep = '\t'))
colData <- read.table(coldata_file, header = T, sep = '\t',
                      stringsAsFactors = TRUE)
cpm <- apply(subset(counts, select = c(-width)), 2,
             function(x) x/sum(as.numeric(x)) * 10^6)

geneLengths <- as.vector(subset(counts, select = c(width)))

  # compute rpkm
rpkm <- apply(X = subset(counts, select = c(-width)),
              MARGIN = 2,
              FUN = function(x) {
                10^9 * x / geneLengths / sum(as.numeric(x))
              })

  # Computing TPM
  #find gene length normalized values
rpk <- apply( subset(counts, select = c(-width)), 2,
              function(x) x/(geneLengths/1000))

  #normalize by the sample size using rpk values
tpm <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)
  #Selecting genes
V <- apply(tpm, 1, var)
selectedGenes <- names(V[order(V, decreasing = T)][1:100])

#transpose the matrix
M <- t(tpm[selectedGenes,])
# transform the counts to log2 scale
M <- log2(M + 1)
#compute PCA
pcaResults <- prcomp(M)
#Correlation matrix
correlationMatrix <- cor(tpm)

## Starting the UI
ui <- pageWithSidebar(
  headerPanel('Hazan Exercise 12 Shiny'),
  sidebarPanel(
    uiOutput("filter_degree")
    
  ),
  mainPanel(
    uiOutput('plot')
    
  )
)
#server.r
server <- function(input, output, session) {
  output$filter_degree<-renderUI({
    selectInput("plot","Plot to Display",choices = c("Heatmap", "PCA", "Correlation plot"))})
  
  
  output$plot <- renderUI({
    if(input$plot=="Heatmap"){
      
      output$plot1<-renderPlot({
        pheatmap(tpm[selectedGenes,], scale = 'row',
                 show_rownames = FALSE,
                 annotation_col = colData)
      })
      plotOutput("plot1")
    }
      
    else if(input$plot=="PCA"){
      output$plot2<-renderPlot({
        x <- c(1:100)
        autoplot(pcaResults, data = colData, colour = 'group') + theme_bw()
      })
      plotOutput("plot2")
    }
    
    else if(input$plot=="Correlation plot"){
      output$plot3<-renderPlot({
        
        corrplot(correlationMatrix, order = 'hclust',
                 addrect = 2, addCoef.col = 'white',
                 number.cex = 0.7)
      })
      plotOutput("plot3")
    }
  })
}

shinyApp(ui = ui, server = server)