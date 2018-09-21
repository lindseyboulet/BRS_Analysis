# Load packages
y <- c("shiny", "plyr", "dplyr", "reshape", "shinythemes", "ggplot2", "readr", "shinyFiles",
       "shinyjs", "taRifx", "shinydashboard", "here", "plotrix", "data.table", "gridExtra",
       "rmarkdown", "ggthemes", "DT")
for(i in 1:length(y)){is_installed <- function(mypkg){is.element(mypkg, installed.packages()[,1])}
if(!is_installed(y[i])){install.packages(y[i], repos="http://lib.stat.cmu.edu/R/CRAN")}
library(y[i], character.only=TRUE, quietly=TRUE,verbose=FALSE)}

header <- dashboardHeader(title = "BRS Analysis")
ui <- dashboardPage(skin = "purple",
                    title="Baroreflex Sensitivity Analysis",
                    header,                    
                    dashboardSidebar(width = 150,
                                     uiOutput("fileId")
                    ),      
                    dashboardBody(
                      fluidRow(
                        column(12,
                               h3(strong("Sum of MSNA")))
                      ),
                      fluidRow(
                        column(9,
                               plotOutput('p1')),
                        column(3,
                               DT::dataTableOutput("dfTable"))
                      ),
                      fluidRow(
                        column(9,
                               h3(strong("Weighted Burst Probability")))
                      ),
                      fluidRow(
                        column(9,
                               plotOutput('p2')),
                        column(3,
                               DT::dataTableOutput('dfTableWt'))
                        
                      ),
                      fluidRow(
                        column(2,
                               actionButton(inputId = "resetP2", label = "Reset Plot")
                               )
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
  # df <- read.csv("rawData.csv")
  
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
             "sna_cmt_text", "sna_cmt_no", "sna_maxMin")
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
  processList <- list()
  j <- 1
  for(i in 1:nrow(df)){
    if(i!=nrow(df)){
      if(df$sna_cmt_no[i] != df$sna_cmt_no[i+1]){
        if(!is.na(df$DBP[i+1])){
          if(!is.na(df$bp_cmt_text[i+1] != "CMT")){
            processList[[j]] <- df[i+1, c(6,8,9)]
            j <- j+1
          }
        }
      }
    }
  }
  dfDBPBursts <- ldply(processList, data.frame)
  dfDBPBursts$normalized <- dfDBPBursts$sna_maxMin/max(dfDBPBursts$sna_maxMin)*100
  dfDBPBursts$bins <- cut(dfDBPBursts$DBP,seq(40,120, by = 2))
  dfDBPBursts
})  

dfRaw <- reactive({
  df <- cleanData()
  dfRaw <- df[complete.cases(df),]
  dfRaw$bins <- cut(dfRaw$DBP,seq(40,120, by = 2))
  dfRaw$DBPbinEnd <- as.numeric(substring(unlist(lapply(strsplit(as.character(dfRaw$bins),
                                                                 split = ',', fixed = TRUE), '[[', 2)), 1,2))
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
                                                                   split = ',', fixed = TRUE), '[[', 2)), 1,2))
  sumDF <-  summarise_all(group_by(dfDBPBursts, bins), funs(sum(., na.rm = TRUE)))
  dfCount <- merge(dfCount, sumDF[,c(1,5)], by = "bins")
  dfCount <- mutate(dfCount, amp_burst = normalized/n.bursts) %>%
    mutate(amp_incidence = amp_burst*burstProb)
  
  dfCount
})  


output$dfTable <- DT::renderDataTable({
  dfCount <- dfCount()
  dfCount2 <- dfCount[,c(9,12)]
  colnames(dfCount2) <- c("Diastolic BP", "Amp*Incidence")
  datatable(dfCount2, options = list(dom = 't'), rownames = FALSE,
            selection = 'none')
})

output$dfTableWt <- DT::renderDataTable({
  dfCount <- dfCount()
  dfCount2 <- dfCount[,c(9,8)]
  colnames(dfCount2) <- c("Diastolic BP", "Burst Probability")
  datatable(dfCount2, options = list(dom = 't'), rownames = FALSE)
})

m1 <- reactive({
  dfCount <- dfCount()
  lm(amp_incidence ~ DBPbinEnd, dfCount)
})

output$p1 <- renderPlot({
  m <- m1()
  dfCount <- dfCount()
  eq <- substitute(italic(y) == b %.% italic(x) + a*","~~italic(r)^2~"="~r2,
                   list(a = round(coef(m)[1], digits = 1),
                        b = round(coef(m)[2], digits = 3),
                        r2 = format(summary(m)$r.squared, digits = 3)))  
  p1 <- ggplot(dfCount, aes(x= DBPbinEnd, y = amp_incidence)) +
    geom_point(size = 3) +
    stat_smooth(method = "lm", se = FALSE) +
    labs(x = "Diastolic BP", y = "Total MSNA") +
    annotate(geom = 'text', x = -Inf, y = Inf, hjust = -0.2, vjust = 2,
             label = as.character(as.expression(eq)), parse = TRUE, size = 8)+
    theme_economist(base_size = 20) 
  p1
})
  
dfRawWt <- reactive({
  dfRaw <- dfRaw(); dfCount<- dfCount()
  merge(dfRaw, dfCount[,8:9], by = "DBPbinEnd", all.x = TRUE)
})  


rxVals <- reactiveValues(keepRows = NULL)

observeEvent(input$fileId,{
                 rxVals$keepRows <-  data.frame(matrix(1,nrow = isolate(nrow(dfRawWt())),
                                                         ncol = isolate(ncol(dfRawWt())), byrow = FALSE))
                 colnames(rxVals$keepRows) <- colnames(dfRawWt())
               })

observeEvent(input$dfTableWt_cell_clicked, {
  res <- input$dfTableWt_cell_clicked
  df <- dfRawWt()
  rxVals$keepRows[which(df$DBPbinEnd == res$value),] <- 0
})

dfRawWtClean <- reactive({
  dfRawWt <- dfRawWt()
  dfRawWt[which(rxVals$keepRows != 1, arr.ind = TRUE)] <- NA
  dfRawWt
})

proxy = dataTableProxy('dfTableWt')

observeEvent(input$resetP2, {
  rxVals$keepRows <-  data.frame(matrix(1,nrow = isolate(nrow(dfRawWt())),
                                        ncol = isolate(ncol(dfRawWt())), byrow = FALSE))
  colnames(rxVals$keepRows) <- colnames(dfRawWt())
  proxy %>% selectRows(NULL)
})


m2 <- reactive({
  dfRaw <- dfRawWtClean()
  lm(burstProb ~ DBPbinEnd, dfRaw)
})


output$p2 <- renderPlot({
  m2 <- m2()
  dfRawWt <- dfRawWtClean()
    eq2 <- substitute(italic(y) == b %.% italic(x) + a*","~~italic(r)^2~"="~r2,
                    list(a = round(coef(m2)[1], digits = 1),
                         b = round(coef(m2)[2], digits = 3),
                         r2 = format(summary(m2)$r.squared, digits = 3)))
  
  p2 <- ggplot(dfRawWt, aes(x = DBPbinEnd, y = burstProb)) +
    geom_point(size = 3) +
    stat_smooth(method = "lm", se = FALSE) +
    labs(x = "Diastolic BP", y = "Burst Probability") +
    annotate(geom = 'text', x = -Inf, y = Inf, hjust = -0.2, vjust = 2,
             label = as.character(as.expression(eq2)), parse = TRUE, size = 8)+
    theme_economist(base_size = 20) 
  p2
})

m3 <- reactive({
  dfCount <- dfCount()
  lm(burstProb  ~ DBPbinEnd, dfCount)
})

output$m1 <- renderTable(colnames = FALSE, {
  m1 <- m1()
  sum <- summary(m1)
  df <- data.frame(matrix(nrow = 4, ncol = 2))
  df[,1] <- c("m", "b", "p", "r2")
  df[1:2,2] <- round(sum$coefficients[c(2,1),1],2)
  df[3,2] <- round(sum$coefficients[1,4],2)
  df[4, 2] <-  sum$r.squared
  df
})
output$m2 <- renderTable(colnames = FALSE, {
  m1 <- m2()
  sum <- summary(m1)
  df <- data.frame(matrix(nrow = 4, ncol = 2))
  df[,1] <- c("m", "b", "p", "r2")
  df[1:2,2] <- round(sum$coefficients[c(2,1),1],2)
  df[3,2] <- round(sum$coefficients[1,4],2)
  df[4, 2] <-  sum$r.squared
  df
})
output$m3 <- renderTable(colnames = FALSE, {
  m1 <- m3()
  sum <- summary(m1)
  df <- data.frame(matrix(nrow = 4, ncol = 2))
  df[,1] <- c("m", "b", "p", "r2")
  df[1:2,2] <- round(sum$coefficients[c(2,1),1],2)
  df[3,2] <- round(sum$coefficients[1,4],2)
  df[4, 2] <-  sum$r.squared
  df
})
}

# Create Shiny object
shinyApp(ui = ui, server = server)
