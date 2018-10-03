# Load packages
y <- c("shiny", "plyr", "dplyr", "reshape", "shinythemes", "ggplot2", "readr", "shinyFiles",
       "shinyjs", "taRifx", "shinydashboard", "here", "plotrix", "data.table", "gridExtra",
       "rmarkdown", "ggthemes", "DT", "knitr", "here", "Cairo", "grid")
for(i in 1:length(y)){is_installed <- function(mypkg){is.element(mypkg, installed.packages()[,1])}
if(!is_installed(y[i])){install.packages(y[i], repos="http://lib.stat.cmu.edu/R/CRAN")}
library(y[i], character.only=TRUE, quietly=TRUE,verbose=FALSE)}

header <- dashboardHeader(title = "BRS Analysis")
ui <- dashboardPage(skin = "purple",
                    title="Baroreflex Sensitivity Analysis",
                    header,                    
                    dashboardSidebar(width = 300,
                                     uiOutput("fileId"),
                                     actionButton("saveData", "Save Analysis")
                    ),      
                    dashboardBody(
                      fluidRow(
                        column(12,
                               h3(strong("Sum of MSNA")))
                      ),
                      fluidRow(
                        column(9,
                               plotOutput('p1',
                                  click = "plot1_click")),
                        column(3,
                               DT::dataTableOutput("dfTable"))
                      ),
                      fluidRow(
                        tags$head(
                          tags$style(HTML('#resetP1{background-color:#6a51a3;
                                                    color:#f7fbff}'))
                        ),
                        column(2,
                               actionButton(inputId = "resetP1", label = "Reset Plot")
                        )
                      ),
                      fluidRow(
                        column(9,
                               h3(strong("Weighted Burst Probability")))
                      ),
                      fluidRow(
                        column(9,
                               plotOutput('p2',
                                          click = "plot2_click")),
                        column(3,
                               DT::dataTableOutput('dfTableWt'))
                        
                      ),
                      fluidRow(
                        tags$head(
                          tags$style(HTML('#resetP2{background-color:#6a51a3;
                                                    color:#f7fbff}'))
                        ),
                        column(2,
                               actionButton(inputId = "resetP2", label = "Reset Plot")
                               )
                      ),
                      fluidRow(
                        column(9,
                               h3(strong("Regression Summary")))
                      ),
                      fluidRow(
                        column(4,
                               h4(strong("Sum of MSNA"))
                               ),
                        column(4,
                               h4(strong("Weighted Probability"))
                               ),
                        column(4,
                               h4(strong("Unweighted Probability"))
                        )
                      ),
                      fluidRow(
                        column(4,
                               tableOutput('m1')),
                        column(4,
                               tableOutput('m2')),
                        column(4,
                               tableOutput('m3'))
                      )
                    )
)



# Define server function
server <- function(input, output, session) {
  session$onSessionEnded(stopApp)

  plotPDF <- function(plotList, ncol, filename, width, height){
    # Load packages
    cairo_pdf(file = filename, width = width, height = height)
    gpVec <- c("gp1", "gp2", "gp3", "gp4", "gp5", "gp6", "gp7", "gp8")
    p1vec <- c("p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8")
    plotNum <- length(plotList)
    names(plotList) <- p1vec[1:plotNum]
    args <- ""
    colArgs <- ""
    for (i in 1:plotNum){
      eval(parse(text = paste0(gpVec[i], "<-ggplot_gtable(ggplot_build(plotList[[i]]))")))
      args <- paste(args, ", ", gpVec[i], "$widths[2:3]", sep = "")
      colArgs <- paste(colArgs, ", ", gpVec[i], sep = "")
    }
    substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}
    args <- substrRight(args, nchar(args)-2)
    colArgs <- substrRight(colArgs, nchar(colArgs)-2)
    maxWidth = eval(parse(text = paste("unit.pmax(", args, ")", sep = "")))
    gpList <- list()
    for (i in 1:plotNum){
      eval(parse(text = paste(gpVec[i], "$widths[2:3]", "<- maxWidth", sep = "")))
      eval(parse(text = paste(gpVec[i], "$heights", "<- gp1$heights", sep = "")))
      gpList[[i]] <- eval(parse(text = gpVec[i]))
    }
    eval(parse(text = paste("grid.arrange(", colArgs, ", ncol = ", ncol, ")", sep = "")))
    dev.off()
  }  
  
  fileIndex <- reactive({
    input$updateFileList
    fileID <- list.files(path =  "./rawData")
    fileID
  })
  
  output$fileId <- renderUI({
    fileIndex <- fileIndex() 
    selectInput(inputId = "fileId", label = strong("File"),choices = fileIndex)
  })
  

cleanData <- reactive({
  df <- read.csv(paste0("./rawData/", input$fileId))
  cols <-  c("start", "duration", "bp_cmt_text", "bp_mean", "SBP", "DBP",
             "sna_cmt_text", "sna_cmt_no", "sna_Min","sna_maxMin")
  colnames(df) <- cols
  df <- df[-(1:2),]
  numCols <- cols[-c(3,7)]
  for(i in numCols){
    df[,i] <- as.numeric(as.character(df[,i]))
  }
  df
})

dfDBPBursts <- reactive({
  df <- cleanData()
  dfDBPBursts <- filter(df, sna_cmt_text == "BURST" & bp_cmt_text != "CAL" | sna_cmt_text == "burst") %>% .[!duplicated(.$sna_cmt_no), ]
  dfDBPBursts$normalized <- dfDBPBursts$sna_maxMin/max(dfDBPBursts$sna_maxMin)*100
  dfDBPBursts$bins <- cut(dfDBPBursts$DBP,seq(40,120, by = 2))
  dfDBPBursts
})  

dfRaw <- reactive({
  df <- cleanData()
  dfRaw <- df[complete.cases(df),]
  dfRaw <- filter(dfRaw, bp_cmt_text != "CAL" & bp_cmt_text != "cal")
  dfRaw$bins <- cut(dfRaw$DBP,seq(40,120, by = 2))
  dfRaw$DBPbinEnd <- as.numeric(substring(unlist(lapply(strsplit(as.character(dfRaw$bins),
                                                                 split = ',', fixed = TRUE), '[', 2)), 1,2))
  dfRaw
})
 
dfRawCount <- reactive({
  dfRaw <- dfRaw()
  dfRawCount <- dfRaw %>% count(bins)
  dfRawCount$cumsum <- cumsum(dfRawCount$n)
  dfRawCount$binFreq <- round((dfRawCount$n / sum(dfRawCount$n))*100, 0) # why round?
  dfRawCount
})  

dfBurstsCount <- reactive({
  dfDBPBursts <- dfDBPBursts()
  dfBurstsCount <- dfDBPBursts %>% count(bins)
  dfBurstsCount$cumsum <- cumsum(dfBurstsCount$n)
  dfBurstsCount$binFreq <- round((dfBurstsCount$n / sum(dfBurstsCount$n))*100, 0) # why round?
  dfBurstsCount
})  

dfCount <- reactive({
  dfRawCount <- dfRawCount()
  dfBurstsCount <- dfBurstsCount()
  dfDBPBursts <- dfDBPBursts()
  dfCount <- merge(dfRawCount, dfBurstsCount, by = "bins", suffixes = c(".all", ".bursts"))
  dfCount$burstProb <- round(dfCount$n.bursts/dfCount$n.all*100, 2)
  dfCount$DBPbinEnd <- as.numeric(substring(unlist(lapply(strsplit(as.character(dfCount$bins),
                                                                   split = ',', fixed = TRUE), '[', 2)), 1,2))
  sumDF <-  summarise_all(group_by(dfDBPBursts[,-c(3,7)], bins), funs(sum(., na.rm = TRUE)))
  dfCount <- merge(dfCount, sumDF[,c(1,10)], by = "bins")
  dfCount <- mutate(dfCount, amp_burst = normalized/n.bursts) %>%
    mutate(amp_incidence = amp_burst*burstProb)
  
  dfCount
})  

dfRawWt <- reactive({
  dfRaw <- dfRaw(); dfCount<- dfCount()
  merge(dfRaw, dfCount[,8:9], by = "DBPbinEnd", all.x = TRUE)
})  

rxVals <- reactiveValues(p1keepRows = NULL, p2keepRows = NULL)

observeEvent(input$fileId,{
  rxVals$p1keepRows <-  data.frame(matrix(1,nrow = isolate(nrow(dfCount())),
                                          ncol = isolate(ncol(dfCount())), byrow = FALSE))
  colnames(rxVals$p1keepRows) <- colnames(dfCount())
  rxVals$p2keepRows <-  data.frame(matrix(1,nrow = isolate(nrow(dfRawWt())),
                                          ncol = isolate(ncol(dfRawWt())), byrow = FALSE))
  colnames(rxVals$p2keepRows) <- colnames(dfRawWt())
})

dfCountClean <- reactive({
  dfCount <- dfCount()
  dfCount[which(rxVals$p1keepRows != 1, arr.ind = TRUE)] <- NA
  dfCount <- dfCount[which(complete.cases(dfCount)),]
  dfCount
})

dfRawWtClean <- reactive({
  dfRawWt <- dfRawWt()
  dfRawWt[which(rxVals$p2keepRows != 1, arr.ind = TRUE)] <- NA
  dfRawWt
})

dtDfTable <- reactive({
  dfCount <- dfCountClean()
  dfCount2 <- dfCount[,c(9,12)]
  dfCount2[,2] <- round(dfCount2[,2], 1)
  colnames(dfCount2) <- c("Diastolic BP", "Amp*Incidence")
  dfCount2
})
output$dfTable <- DT::renderDataTable({
  datatable(dtDfTable(), options = list(dom = 't'), rownames = FALSE,
            selection = 'none')
})

dtDfTableWt <- reactive({
  dfWt <- dfRawWtClean()
  dfCount <- dfCount()
  dfCount <- dfCount[which(dfCount$DBPbinEnd %in% unique(dfWt$DBPbinEnd)),]
  dfCount2 <- dfCount[,c(9,8)]
  colnames(dfCount2) <- c("Diastolic BP", "Burst Probability")
  dfCount2[,2] <- round(dfCount2[,2], 1)
  dfCount2
})
  
 output$dfTableWt <- DT::renderDataTable({ 
  datatable(dtDfTableWt(), options = list(dom = 't'), rownames = FALSE,
            selection = 'none')
})

m1 <- reactive({
  dfCount <- dfCountClean()
  lm(amp_incidence ~ DBPbinEnd, dfCount)
})

p1Fig <- reactive({
  m <- m1()
  dfCount <- dfCountClean()
  eq <- substitute(italic(y) == b %.% italic(x) + a*","~~italic(r)^2~"="~r2,
                   list(a = round(coef(m)[1], digits = 1),
                        b = round(coef(m)[2], digits = 3),
                        r2 = format(summary(m)$r.squared, digits = 3)))  
  p1 <- ggplot(dfCount, aes(x= DBPbinEnd, y = amp_incidence)) +
    geom_point(size = 5, color = '#984ea3') +
    stat_smooth(method = "lm", se = FALSE) +
    labs(x = "Diastolic BP", y = "Total MSNA") +
    annotate(geom = 'text', x = -Inf, y = Inf, hjust = -0.2, vjust = 2,
             label = as.character(as.expression(eq)), parse = TRUE, size = 8)+
    theme_economist(base_size = 20) 
  p1
  })
output$p1 <- renderPlot({p1Fig()})
  
m2 <- reactive({
  dfRaw <- dfRawWtClean()
  lm(burstProb ~ DBPbinEnd, dfRaw)
})

p2Fig <- reactive({
  m2 <- m2()
  dfRawWt <- dfRawWtClean()
    eq2 <- substitute(italic(y) == b %.% italic(x) + a*","~~italic(r)^2~"="~r2,
                    list(a = round(coef(m2)[1], digits = 1),
                         b = round(coef(m2)[2], digits = 3),
                         r2 = format(summary(m2)$r.squared, digits = 3)))
  
  p2 <- ggplot(dfRawWt, aes(x = DBPbinEnd, y = burstProb)) +
    geom_point(size = 5, color = '#984ea3') +
    stat_smooth(method = "lm", se = FALSE) +
    labs(x = "Diastolic BP", y = "Burst Probability") +
    annotate(geom = 'text', x = -Inf, y = Inf, hjust = -0.2, vjust = 2,
             label = as.character(as.expression(eq2)), parse = TRUE, size = 8)+
    theme_economist(base_size = 20) 
  p2
})
output$p2 <- renderPlot({p2Fig()})

m3 <- reactive({
  dfCount <- dfCountClean()
  lm(burstProb  ~ DBPbinEnd, dfCount)
})

m1Df <- reactive({
  m1 <- m1()
  sum <- summary(m1)
  df <- data.frame(matrix(nrow = 4, ncol = 2))
  df[,1] <- c("m", "b", "p", "r2")
  df[1:2,2] <- round(sum$coefficients[c(2,1),1],2)
  df[3,2] <- round(sum$coefficients[1,4],2)
  df[4, 2] <-  sum$r.squared
  df 
})
output$m1 <- renderTable(colnames = FALSE, m1Df())

m2Df <- reactive({
  m1 <- m2()
  sum <- summary(m1)
  df <- data.frame(matrix(nrow = 4, ncol = 2))
  df[,1] <- c("m", "b", "p", "r2")
  df[1:2,2] <- round(sum$coefficients[c(2,1),1],2)
  df[3,2] <- round(sum$coefficients[1,4],2)
  df[4, 2] <-  sum$r.squared
  df
})
output$m2 <- renderTable(colnames = FALSE, m2Df())

m3Df <- reactive({
  m1 <- m3()
  sum <- summary(m1)
  df <- data.frame(matrix(nrow = 4, ncol = 2))
  df[,1] <- c("m", "b", "p", "r2")
  df[1:2,2] <- round(sum$coefficients[c(2,1),1],2)
  df[3,2] <- round(sum$coefficients[1,4],2)
  df[4, 2] <-  sum$r.squared
  df
})
output$m3 <- renderTable(colnames = FALSE, m3Df())

observeEvent(input$plot1_click, {
  res <- nearPoints(dfCount(), input$plot1_click, allRows = TRUE)
  rxVals$p1keepRows[which(res$selected_==TRUE),1] <- 0
})

observeEvent(input$plot2_click, {
  res <- nearPoints(dfRawWt(), input$plot2_click, allRows = TRUE)
  rxVals$p2keepRows[which(res$selected_==TRUE),1] <- 0
})

observeEvent(input$resetP1, {
  rxVals$p1keepRows <-  data.frame(matrix(1,nrow = isolate(nrow(dfCount())),
                                          ncol = isolate(ncol(dfCount())), byrow = FALSE))
  colnames(rxVals$p1keepRows) <- colnames(dfCount())
})

observeEvent(input$resetP2, {
  rxVals$p2keepRows <-  data.frame(matrix(1,nrow = isolate(nrow(dfRawWt())),
                                          ncol = isolate(ncol(dfRawWt())), byrow = FALSE))
  colnames(rxVals$p2keepRows) <- colnames(dfRawWt())
})

outputDf <- reactive({
  m1Df <- m1Df(); m2Df <- m2Df(); m3Df <- m3Df()
  dff <- merge(m1Df, m2Df, by = "X1")
  dff <- merge(m3Df, dff, by = "X1")
  ddf <- data.frame(t(dff))[-1,]
  colnames(ddf) <- c("intercept", "mean", "p-value", "r2")
  for(i in 1:ncol(ddf)){ddf[,i] <- as.numeric(as.character(ddf[,i]))}
  ddf$fileId <- substr(input$fileId, 0, nchar(input$fileId)-4)
  ddf$model <- c("sum", "prob_wt", "prob")
  ddf
})

observeEvent(input$saveData, {
  if(!dir.exists("./output")){
    dir.create("./output")
  }
  outputDf <- outputDf()
  write.csv(outputDf, file = paste0("./output/",
              substr(input$fileId, 0, nchar(input$fileId)-4), "-lmOutput.csv"),
            row.names = FALSE)
  p1 <- p1Fig()
  p1<- p1 +ggtitle(paste0("Sum of MSNA (", substr(input$fileId, 0, nchar(input$fileId)-4), ")"))
  p2 <- p2Fig()
  p2<- p2 +ggtitle("Weighted Burst Probability")
  plotPDF(list(p1,p2), 1, paste0("./output/",
                substr(input$fileId, 0, nchar(input$fileId)-4), ".pdf"), 12, 12)
  })

}

# Create Shiny object
shinyApp(ui = ui, server = server)
